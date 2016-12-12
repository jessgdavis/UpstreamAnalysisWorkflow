# Upstream Analysis Workflow 
This workflow aims to automate the majority of the work required to prepare a set of genes for upstream analysis.

## Instructions 


## Use Cases ###
### Use Case 1 - Find orthologs 
As a user, I have a transcription factor (e.g. Fur) and I know it regulates 80 genes in E.coli. I also have a list of 10 target genomes that I want to study. I want to find, of the 80 genes in E.coli, how many orthologs exists in the 10 target genomes using the reciprocal best hit method.

### Use Case 2 - Find upstream regions 
As a user, once the reciprocal best hits are found, I want to be able to examine both the upstream intergenic region, and the 1000bp upstream region, of a particular hit.

### Use Case 3 - Queries on results
As a user, I want to be able to perform queries on a result set, where a result set has the information as represented in the figure below.
The queries I want to be able to perform are:

* get the list of orthologs (the reciprocal hits) for query gene X

* get the list of orthologs for a specified query gene that all belong to the same genome G

Once I have this information, I want to be able to:

* get the upstream regions of the results in this list

* view the sequences of each of the results in this list

* filter out results which do not meet a specified evalue threshold

![workflowdata.png](https://bitbucket.org/repo/gaaxA5/images/3889694038-workflowdata.png)

## Comments and Known Issues
At this stage, the upstream workflow is very much in the experimental stage, and still needs rigorous testing, as well as the addition of some further 'nice to have' features. 
Some known issues which are currently being addressed are:

* In some cases, if a query gene or the feature associated with a BLAST hit has a particularly long sequence, HTTP 414 errors are still being thrown, which means that sometimes not all the orthologs can be verified. This is an issue that needs further investigation. For example, if the sequence is sufficiently long, would it be okay to trim the query to within the URI length limitations, as the chances of random matches at these kind of lengths is extremely unlikely?

* BLAST requests take a long time, and since the number of BLAST2 requests is equal to the number of query genes multiplied by the number of hits, there can be quite a number of BLAST2 requests being sent. This also poses a problem as these large numbers of BLAST requests are in violation of the BLAST terms of use. As such, it might be worth looking in to porting some kind of local BLAST to .NET, so the alignments for BLAST2 requests can be run locally. This makes sense given that BLAST2 requests are run against a single genome which is already stored in memory.

* Given that there can be a large number of hits on a query gene (though this can be limited, it is likely that users will want to choose how many hits they consider), there are often a large number of efetch requests that need to be sent to get the genomes for each of these hits. This has already been optimised as much as possible by ensuing that the accession numbers sent off are a distinct list, but it can still be quite large. This could lead to violations of the efetch terms of use.

* In some cases, requests to retreive genomes using efetch were timing out due to a lack of a timely response from the server. This appears to have been mostly resolved by increasing the timeout of the request locally to wait indefinitely, and modifying the asynchronous parallel algorithms. However, this needs more rigorous testing before anything can be said with confidence.

Some further features to be included are:

* Making the system more flexible, so that users can make BLASTP requests if desired. This presents some challenges around getting sequence translations, and the fact that GenBank associates data with protein sequences in different ways to the nucleotide sequences. Some infrastructure has been worked on to this effect, but it needs refinement before it is integrated with the workflow.

* An algorithm to retrieve the intergenic upstream region and not just the 1000bp upstream region needs to be developed.

# How to interact with the upstream analysis workflow 

```!fsharp

#r "../packages/NetBioCore.PCL.2.0.150722/lib/net45/Bio.Core.dll"
#r "../packages/NetBioWeb.PCL.2.0.150722/lib/Bio.WebServices.dll"

#r "System.Threading.Tasks"
#r "System.IO"
#r "System.Collections"

#load "UpstreamAnalysis.fs"
open workflowlib.workflowFunctions

// step 0 - input 
// define your working directory
let filePath = @"C:\Users\Jess\AppData\Local\Temp\Data\file20genes.fasta" // replace this with the path to your own file

// take a FASTA file with a set of genomes, or a list of Locus Tags as input (list of locus not yet implemeted)
let fp = new Bio.IO.FastA.FastAParser()
let fs = System.IO.File.OpenRead(filePath) 
let ourSeq = fp.Parse(fs) |> Seq.toList
fs.Flush()
fs.Close()
fs.Dispose()

// step 1 - define query genes and download genome for query genes
let queryGenes = ourSeq |> Seq.take 1
let ``accession number of the gene which your sequences are from`` = "U00096" // change this to the accession of the genome you are working with
let originGenome = allGenomesFromAcessions [``accession number of the gene which your sequences are from``] 

// step 2 - perform first blast with this queryseq
let database = "nr"
let program = Bio.Web.Blast.BlastProgram.Blastn
let ``maxiumum number of results`` = "2"

let ``entrez query for limiting initial blast search`` = 
    "txid1423 [ORGN] OR txid1392 [ORGN] OR txid632 [ORGN] OR txid1301 [ORGN] OR txid813 [ORGN] OR txid1313 [ORGN] OR txid1280 [ORGN]"
    |> (fun s -> s.Replace(' ', '+'))
let blast1extraparams = [new System.Collections.Generic.KeyValuePair<string, string>("MAX_NUM_SEQ", ``maxiumum number of results``)]
// you can use the entrez query to limit the search to specific organisms, at the moment, this is set to only BLAST against 
    // Bacillus subtilis (taxid:1423)
    // Bacillus anthrasis (taxid:1392)
    // Yersinia pestis (taxid:632)
    // Streptococcus (taxid:1301)
    // Chlamydia trachomastis (taxid:813)
    // streptococcus pneumonia (taxid:1313)
    // staphylococcus aureus (taxid:1280)
let blast1 = getBlastResultsForAllGenes filePath database program blast1extraparams queryGenes

// step 3 - download genomes for the results
let listOfAccessions = uniqueAccessions blast1 
let genomes = allGenomesFromAcessions listOfAccessions

// step 4 - get the full sequences for each of the hits
let fullSeqsForHits = fullSeqForHits genomes blast1 |> Seq.toList 

// step 5 - get only bidirectional best hits
let biDirectionalBestHits = getBlast2BestHitsForAllGenes filePath database program originGenome fullSeqsForHits |> Seq.toList 

// step 6 - analyse!
// each object in the biDirectionalBestHits is a "BlastOutputRecord" which has two properties
// 1. QueryGene (the gene that was sent as input to BLAST) 
// 2. BlastResultSequences (a sequence of SingleBlastHit objects which contains only the bidirectional best hits for the QueryGene) 
    // SingleBlastHit objects have 4 properties:  
    // 1. FeatureSequence (an ISequence corresponding to the full sequence matching the BLAST hit object) 
    // 2. HitSequence (an ISequence representation of the hit alignment sequence returned from BLAST, usually a subsequence of some other feature)
    // 3. Hit (Bio.Web.Blast.Hit object, contains stats about Blast result)
    //4. HitEval (The Evalue associated with this hit)

// you can get the upstream region of your input sequence like this: 
let upstreamOfInputSeq = 
    biDirectionalBestHits.[0]
    |> (fun blastoutputrecord -> 
            blastoutputrecord.QueryGene.Metadata.["Feature"] :?>  Bio.IO.GenBank.FeatureItem
        )
    |> getUpstreamRegion genomes.[0]

// you can get an upstream region for a particular hit on your input sequence like this:
let upstreamSelectedHit = 
    biDirectionalBestHits.[0]
    |> (fun blastoutputrecord -> 
            blastoutputrecord.BlastResultSequences |> Seq.toList
            |> (fun x -> x.[0].FeatureSequence.Metadata.["Feature"] :?>  Bio.IO.GenBank.FeatureItem)
        )
    |> getUpstreamRegion genomes.[0]

```
