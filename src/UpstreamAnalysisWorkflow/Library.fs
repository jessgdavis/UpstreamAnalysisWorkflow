namespace Bio.FSharp

open Bio
open Bio.IO
open Bio.IO.GenBank
open Bio.Web.Blast
open Bio.FSharp.Experimental

open System
open System.IO
open System.Threading
open System.Threading.Tasks
open System.Collections.Generic
open System.Net

// initalize the type
type GenBankSequence = GenBankProviderDev< @"../../../../lib/NC_012686-Chlamydia trachomatis.gbk" >

type SingleBlastHit = {FeatureSequence: GenBankSequence; HitSequence: GenBankSequence; Hit: Hit; HitEval: float}
type BlastOutputRecord = {QueryGene: GenBankSequence; BlastResults: seq<SingleBlastHit>}
type ExtraBlastParams = seq<KeyValuePair<string, string>>

    module UpstreamAnalysis =

        // the genbank input should be for a single gene
        let GetBlastResultForGene (gene : GenBankSequence) blastdb blastp (extraParams : ExtraBlastParams) = 

            async {
                let bp = new BlastRequestParameters([gene.ISequence], extraParams)
                bp.Database <- blastdb
                bp.Program <- blastp

                // set up web handler
                let ncbiWebHandler = new NcbiBlastWebHandler(
                    EndPoint = @"https://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
                    LogOutput = fun str -> Console.WriteLine(str)
                )

                use task = ncbiWebHandler.ExecuteAsync(bp, CancellationToken.None)
                // wait for response to come back before passing along
                Async.AwaitIAsyncResult task |> ignore
                let resStream = 
                    task.Result
                    |> (fun brs ->
                            use sr = new System.IO.StreamReader(brs)
                            sr.ReadToEnd()
                        )
                task.Result.Close()
                task.Result.Dispose()

                let fileName = Directory.GetCurrentDirectory() + gene.Accession + ".txt"

                File.CreateText(fileName) |>
                    (fun (rs:string) (f:StreamWriter) ->
                        f.AutoFlush <- true
                        f.Write(rs)
                        f.Close()
                        f.Dispose()
                    ) resStream

                let bres = 
                    new System.IO.FileStream(fileName, FileMode.Open, FileAccess.Read) |>
                        (fun fs ->
                            let bParser = new BlastXmlParser()
                            // should only be one result set in stream
                            bParser.ParseOne(fs)
                        )

                let res = 
                    bres.Records.[0].Hits 
                    |> Seq.map (fun ht -> 
                                    let seq = GenBankSequence(new Sequence(Alphabets.Protein, ht.Hsps.[0].QuerySequence) :> ISequence)
                                    // Associate the BlastHit and QueryGene with the ISequence so data is not lost
                                    seq.ISequence.ID <- ht.Accession // so the sequence to which the alignment belongs can be easily identified
                                    {FeatureSequence=null; HitSequence=seq; Hit=ht; HitEval=ht.Hsps.[0].EValue}
                                )
            // return a tuple, which has ISeq of query gene assoiated with seq of iseqs of hits
            return {QueryGene=gene; BlastResults=res}
        }

        // run a batch of blast requests in parallel
        let GetBlastResultsForAllGenes blastdb blastp (extraParams: ExtraBlastParams) (queryGenes : seq<GenBankSequence>) =
            queryGenes
             |> Seq.map (fun gene ->
                    GetBlastResultForGene gene blastdb blastp extraParams 
                )
            |> Async.Parallel
            |> Async.RunSynchronously

        // method to delete the temp files created for parsing xml
        let Cleanup filePath querySeqs = 
            Seq.mapi(fun num iseq ->
                let fn = filePath + num.ToString() + ".txt"
                System.IO.File.Delete(fn)
            ) querySeqs

        // okay so now we want to get full genomes for all our blast hits
        // get list of accession numbers for what we need to retrieve, is for all hits
        let UniqueAccessions (allBlastResults:BlastOutputRecord[]) = 
            allBlastResults
            |> Seq.map (fun brSeqs -> brSeqs.BlastResults)
            |> Seq.concat
            |> Seq.map(fun res -> 
                res.Hit.Accession
            )
            |> Seq.distinct // only want to download any given genome once
            |> Seq.toList

        // download the genome corresponding to an accession
        let GetWholeSequenceForAlignment alignmentAccession = 
            async {
                Console.WriteLine("Requesting genome for " + alignmentAccession)
                let request = WebRequest.Create("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+alignmentAccession+"&rettype=gb&retmode=text")
                request.Method <- WebRequestMethods.Http.Get
                let response = request.GetResponse()
                use str = response.GetResponseStream()
                System.Console.WriteLine("Have response for " + alignmentAccession)
                let pars = new GenBankParser()
                // returning a type provider for the downloaded genome
                return GenBankSequence(pars.ParseOne(str))
            }

        // do all requests in parallel to download genomes based on list of accessions
        let AllGenomesFromAcessions accessionList = 
            accessionList 
            |> Seq.map (fun ac -> GetWholeSequenceForAlignment ac)
            |> Async.Parallel
            |> Async.RunSynchronously

        // extract a feature sequence from a genome file
        let FindFeatureSequence (genome : GenBankSequence) (feature:FeatureItem) = 
            let geneStart = int64(feature.Location.LocationStart - 1)
            let geneEnd = int64(feature.Location.LocationEnd)
            let geneLength = geneEnd - geneStart
            let seq = genome.ISequence.GetSubSequence(geneStart, geneLength)
            //// Associate the Gene object with the ISeq so that information is not lost
            seq.Metadata.Add("Feature", feature)
            seq

        // get the full sequence of the feature for each hit, returns updated list of blastResultRecords
        let FullSeqForHits (genomes:seq<GenBankSequence>) (allBlastResults:BlastOutputRecord[]) = 
            allBlastResults
            |> Seq.map (fun orGenWBr ->
                    let updatedBlastResultSequences = 
                        orGenWBr.BlastResults
                        // map over BlastResultSequences to find full seqs for each
                        |> Seq.map (fun sbh -> 
                            // find the downloaded genome that corresponds to this record
                            let genomeForRecord = 
                                genomes |> 
                                Seq.tryFind (fun (genome:GenBankSequence) ->
                                            // is a match if the ID (accession) of the Blast result matches the accession of the downloaded record
                                            genome.Accession = sbh.Hit.Accession
                                         )
                            // try to find the blast result as a feature in that genome
                            match genomeForRecord with 
                            | None -> {FeatureSequence=null; HitSequence=null; Hit=null; HitEval=0.0}
                            | Some(genomeForRecord) ->
                                let gfrmd = genomeForRecord.GenBankMetadata
                                let hitStart = sbh.Hit.Hsps.[0].HitStart
                                let hitEnd = sbh.Hit.Hsps.[0].HitEnd
                                let featureSeq = 
                                    gfrmd.Features.All.Exists(fun f ->
                                            int64(f.Location.LocationStart) <= hitStart && int64(f.Location.LocationEnd) >= hitEnd && not (f.Key.Equals("source"))
                                        ) // identifies whether any features matching this predicate exist
                                    |> (fun exists -> 
                                            match exists with
                                            | true -> 
                                                // if it exists, grab the feature
                                                gfrmd.Features.All.Find(fun f ->
                                                        int64(f.Location.LocationStart) <= hitStart && int64(f.Location.LocationEnd) >= hitEnd && not (f.Key.Equals("source"))
                                                    ) // we have a feature
                                            | false ->
                                                // if not return nothing
                                                null
                                        )
                                    |> (fun findResult ->
                                                    match findResult with 
                                                    | null -> new Bio.Sequence(Alphabets.DNA, "") :> ISequence // if the hit didn't match any regions
                                                    | _ -> FindFeatureSequence genomeForRecord findResult
                                                ) // end mchs
                                    |> GenBankSequence 
                                // we want to update the original blast result seq so that the Sequence component is now the sequence of the feature we found
                                // new record which has old info + newly found FeatureSequence
                                {FeatureSequence=featureSeq; HitSequence=sbh.HitSequence; Hit=sbh.Hit; HitEval=sbh.HitEval}
                        )
                        // we may as well filter out any dodgy hits which had no sequences
                        |> Seq.filter(fun s -> 
                                        match s.FeatureSequence.ToString() with 
                                        | "" | null -> false
                                        | _ -> true
                                    )
                    // return updated overall record
                    {QueryGene=orGenWBr.QueryGene; BlastResults=updatedBlastResultSequences}
                ) 

        // run second blast with this filtered list of sequences, using the full feature sequence as the parmeter
        let GetBlast2BestHitsForAllGenes blastdb blastp genome (allGenes:seq<BlastOutputRecord>) =
            allGenes
            |> Seq.map (fun singleInitialBlastOutputRecord -> 
                    let ``original gene`` = singleInitialBlastOutputRecord.QueryGene
                    let blast2extraparams = [new KeyValuePair<string, string>("ENTREZ_QUERY", ``original gene``.ISequence.ID + "[Accession]")]

                    let isBiDir = 
                        singleInitialBlastOutputRecord.BlastResults
                        |> Seq.map(fun ibrs -> ibrs.FeatureSequence)
                        |> GetBlastResultsForAllGenes blastdb blastp blast2extraparams
                        |> FullSeqForHits genome
                        |> Seq.map(fun blastOrForHit -> 
                                let bestHitInGenome = 
                                    blastOrForHit.BlastResults
                                    |> Seq.sortBy (fun brs2 -> brs2.HitEval)
                                    |> Seq.toList
                                    |> (fun x -> x.[0])
                    
                                if (bestHitInGenome.FeatureSequence.ToString() = ``original gene``.ISequence.ToString()) then
                                    // associate the now known location data with the query gene for later use
                                    let fmd = bestHitInGenome.FeatureSequence.ISequence.Metadata.["Feature"] :?> Bio.IO.GenBank.FeatureItem
                                    ``original gene``.ISequence.Metadata.Add("Feature", fmd)
                                    true
                                else
                                    false
                            )

                    let bidirbesthits = 
                        Seq.map2 (fun isBiDirel blastrsel -> 
                                if isBiDirel then 
                                    blastrsel
                                else 
                                    {FeatureSequence=null; HitSequence=null; Hit=null; HitEval=0.0}
                            ) isBiDir singleInitialBlastOutputRecord.BlastResults
                        |> Seq.filter (fun x ->
                                match x.FeatureSequence with 
                                | null -> false
                                | _ -> true
                            )

                    {QueryGene=``original gene``; BlastResults=bidirbesthits}
                )

        /// <summary> Method to grab the upstream region (1000bp) of a given genome feature, ensuring correct behaviour for features on the compliment sequence </summary>
        /// <param name="genomeSeq"> Bio.ISequence of the genome which contains the given gene </param>
        /// <param name="gene"> Bio.Genbank.FeatureItem to find upstream region of </param>
        /// <returns> Bio.ISequence representing the upstream region of the feature </returns>
        let GetUpstreamRegion (genomeSeq:GenBankSequence) (gene:FeatureItem) = 
            // Helper method for getting upstream region of a gene -> grabs 1000bp upstream of specified gene start location. Ff there are not 1000bp before the feature, just grab from start of genome
            let getUpstreamSeq1000bp (genomeSeq:Bio.ISequence) upstreamStartPoint = 
                let numBPToGet = int64(1000)
                let genomeStart = int64(0)
                if (upstreamStartPoint <= genomeStart) then
                    genomeSeq.GetSubSequence(genomeStart, numBPToGet) else
                        genomeSeq.GetSubSequence(upstreamStartPoint, numBPToGet)
            // if the feature is on the other strand make sure to take upstream of reversed complemented sequence
            let upstreamRegion = 
                match gene.Location.Operator.ToString() with 
                    |"Complement" -> int64(gene.Location.LocationEnd)
                    | _ -> int64(gene.Location.LocationStart - 2)
                |> getUpstreamSeq1000bp genomeSeq
            match gene.Location.Operator.ToString() with 
                    |"Complement" -> upstreamRegion.GetReverseComplementedSequence()
                    | _ -> upstreamRegion
