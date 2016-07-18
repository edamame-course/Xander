#Xander

Authored by Taylor Dunivin for EDAMAME2016 based on a [previous tutorial] (https://github.com/rdpstaff/Xander_assembler) by Qiong Wang for EDAMAME2015 

[EDAMAME-2016 wiki] (https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

##Learning Objectives
* Prepare gene references for assembly
* Understand Xander assembly parameters 
* Assemble one gene from a megatenome
* Examine the abundance and diversity of genes of interest from metagenome data

---

###Citations
Wang, Q., J. A. Fish, M. Gilman, Y. Sun, C. T. Brown, J. M. Tiedje and J. R. Cole. 2015. Xander: Employing a Novel Method for Efficient Gene-Targeted Metagenomic Assembly. Microbiome. 3:32. DOI: [10.1186/s40168-015-0093-6] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6). 

Fish, J. A., B. Chai, Q. Wang, Y. Sun, C. T. Brown, J. M. Tiedje, and J. R. Cole. 2013. FunGene: the Functional Gene Pipeline and Repository. Front. Microbiol. 4: 291. DOI: [10.3389/fmicb.2013.00291] (http://www.ncbi.nlm.nih.gov/pubmed/24101916).

###Required tools
* [RDPTools] (https://github.com/rdpstaff/RDPTools)
* Python 2.7+
* Java 2.6+
* [HMMER 3.1] (http://hmmer.org) (If using HMMER 3.0 add --allcol to bin/run_xander_skel.sh )
* UCHIME 

###What is Xander?
Xander is a gene-targeted metagenome assembler. This means that it takes metagenomic data and assembles it based on your gene(s) of interest only. This can be useful if you have metagenomic data and a gene-centric study in mind. Classic approaches to metagenome assembly tend to assemble the most abundant organisms and can have more limited recovery for individual genes. This can be problematic, especially if you are interested in only a handful of specific genes. Xander avoids this problem by guiding the assembly based on your gene of interest. 

###How does Xander work?
Xander uses a profile hidden Markov models (HMMs) to guide assembly of metagenomic data. Check out the [wikipedia page] (https://en.wikipedia.org/wiki/Hidden_Markov_model) for more information on HMMs. In short, these are probabalistic graphs that predict protein sequences. Basically these models say how likely a specific sequence is to have come from a known sequence. Xander uses these to guide the assembly of metagenomic data, so you only end up with assembled contigs that resemble your gene of interest.

###What does Xander require?
1. High quality data file. Sequencing depth will greatly effect the quality of Xander results. It is important to have deep sequencing, especially for low abundance genes in order to get quality, near full length proteins. The theory of "garbage in, garbage out" also applies to Xander. Input data should already be quality filtered. 
2. Biological insight on your gene of interest. Because Xander uses existing knowledge of your gene of interest for the assembly, it is best to run well characterized genes especially proteins with crystal structures. 

The following figure is from the [Xander publication] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6) and is a good overview of how a gene-targeted metagenome assembly works. 
![Structure](https://github.com/edamame-course/Xander/blob/master/Xander_structure.png)

As you can see, Xander will require several inputs. Our metagenomic data will be the input for the de Bruijn graph assembly, but we will also need gene references to make the HMM. There are three ways to go about setting up a gene reference. 
  1. The analysis pipeline is preconfigured with the _rplB_ phylogenetic marker gene, and nitrogen cycling genes including _nirK_, _nirS_, _nifH_, _nosZ_ and _amoA_. These require no work from you. *In this tutorial we will use _rplB_ and do not need to prepare our gene reference.*
  2. Your gene of interest is on the RDP's FunGene database. This requires a bit of work and biological insight. See bottom of this page for instructions. 
  3. Your gene of interest is not preconfigured or on the FunGene database. This requires a bit of research, work, and biological insight. See bottom of this page for instructions. 
  
Overall, running Xander has three steps: _Build_, _Find_, and _Search_

1. __Build__: Make a Bloom filter (one per sequence sample) with kmers of specified length. A bloom filter helps reduce memory usage. See these papers for more information on bloom filters [1] (http://dl.acm.org/citation.cfm?doid=362686.362692) [2] (http://www.sciencedirect.com/science/article/pii/0020019094000328). 
2. __Find__: Find starting kmers for the genes using aligned reference sequences (has multithread option that allows you to find starting kmers for multiple genes at once). 
3. __Search__: Assemble contigs for your gene(s) of interest

For more detail on these steps, check out the [Xander publication] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6)!

Now that we have a basic idea of how Xander works and what we need, let's get started!

###1 Connect to Xander AMI
For this tutorial, we will use an existing AMI that contains all the necessary tools to run Xander. Search for the AMI name "RDP-Edamame-2015" or ID "ami-e973b782"

Click launch.

Select instance type: General purpose, m3.large, 7.5G

###2 Prepare gene reference

The Xander assembler, which is included in the RDPTools, can be accessed and downloaded from the [RDP staff GitHhub] (https://github.com/rdpstaff/RDPTools). Since this is already included in the AMI,  we can skip this download in the tutorial and navigate to the Xander assembler directory.

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler
ls
```

In the Xander assembler, you should see several directories including ```gene_resource```. If you navigate to this directory, you will see a list of genes that have pre-made hidden Markov models that can be used for assembly of your metagenome. We will now go through how to make your own profile hidden Markov model in case your gene of interest is not pre-processed. 

```gene_resource``` contains reference sequence files and models. The analysis pipeline is preconfigured with the _rplB_ phylogenetic marker gene, and nitrogen cycling genes including _nirK_, _nirS_, _nifH_, _nosZ_ and _amoA_.

In this tutorial, we will use the preexisting _rplB_ files and models, but we'll still go over what's in this directory. (Note: if you want to know more about making your own gene directory, see the bottom of this tutorial). 

Let's navigate to the input data for our gene of interest _rplB_. 

```
cd ~/tools/RDPTools/Xander_assembler/gene_resource/rplB/originaldata
ls
```

The ```originaldata``` directory holds the four required files for Xander: gene.seeds, gene.hmm, framebot.fa, and nucl.fa. 

You should have all of these in your directory, but we'll go over how you can get these files from the [RDP's FunGene database] (http://fungene.cme.msu.edu). You can select any gene you want for now on the home page. 
* **gene.seeds**: a small set of **full length, high quality** protein sequences in FASTA format, used to build gene.hmm, forward and reverse HMMs. These sequences are used to make FunGene databases, so if your gene of interest is in FunGene, you can download them there. If not, you would need to find high quality sequences on your own (ideally from papers with protein crystal structures) and send them to RDP to make a FunGene database. Select your gene of interest, and then select "download seeds" 
* **gene.hmm**: this is the HMM built from gene.seeds using original HMMER3. This is used to build for_enone.hmm and align contigs after assembly. You can also download this from FunGene. Select your gene of interest, and then select "download HMM"
* **framebot.fa**: a large near full length known protein set for identifying start kmers and FrameBot nearest matching. You can also download these from FunGne, but you will need to select which sequences you want. More diversity is better, more sequences means more starting points (more computational time) but less susceptible to noise than model creation. Prefer near full-length and well-annotated sequences. Select "show/hide filter options" and filter with Minimum HMM Coverage at least 80 (%). You can also filter based on number of amino acids and score (both will change depending on the gene of interest). Then click "Begin Analysis" and download the protien sequences. 
* **nucl.fa**: a large near full length known set used by UCHIME chimera check. Can be downloaded from FunGene using the same sequences used for framebot.fa. Clicl "Begin Analysis" and select "nucleotide download" before downloading.

These four files were used to make our HMMs, and they are already in our gene directory.

```
cd ..
ls
```

You should now see three files (and the ```originaldata``` directory): **ref_aligned.faa**, **for_enone.hmm**, and **rev_enone.hmm**
* **for_enone.hmm** and **rev_enone.hmm** for the forward and reverse HMMs respectively. This is used to assemble gene contigs.
* **ref_aligned.faa** file containing a set of protein reference sequences aligned with for_enone.hmm. This is used to identify starting kmers. Need to manually examine the alignment using Jalview or alignment viewing tool to spot any badly aligned sequences. If found, it is likely there are no sequences in gene.seeds close these sequences. You need to valide these problem sequences to see they are from the gene you are interested, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps.

Now that we understand our gene reference files, we can begin to set up the assembler. 

##3 Set up environment and parameters
Since you may want to run Xander multiple times, it can be useful to make a directory for each project that includes the gene name and dataset used. We will do this now. 

In our case, we will use data provided by the creators of Xander. This demo_reads file here contains a subset of reads from one of seven corn rhizosphere replicates used in the original Xander publication. This subset is enriched in reads matching rplB, nirK and nifH genes. **These paired-end reads have already been quality trimmed and merged using RDP's read assembler**.

Let's navigate to the ```Xander_assembler``` directory and make a new directory in Xander assembler for this analysis.

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler
mkdir rplB_demo
```

Now we're ready to adjust paramters for analysis!


We will need to edit one shell script to run Xander. You may want to keep the original file for future reference, so we will copy it to into ```rplB_demo``` instead of editing them directly. 

```
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/xander_setenv.sh /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
cp /home/ubuntu/tools/RDPTools/Xander_assembler/bin/run_xander_skel.sh /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
```

We need to change ```xander_setenv.sh``` to reflect our directories and gene of interest. This is also where we can adjust Xander parameters. First let's think about our parameter options. 

#####File locations
* __SEQFILE__ -- Absolute path to your sequence file(s). Can use wildcards to point to multiple files (fasta, fataq or gz format)
* __WORKDIR__ --Absolute path to where you want to include your output files. In our case, this is rplB_demo
* __REFDIR__--Absolute path to Xander_assembler directory
* __JARDIR__--Absolute path to RDPTools 
* __SAMPLE SHORTNAME__ -- a short name for your sample, prefix of final files and contig IDs (needed when pool contigs from multiple samples)

#####Input data parameters
* __K SIZE__ -- K-mer size to assemble at, must be divisible by 3 (recommend 45, maximum 63). Higher numbers yield more stringent results
* __FILTER SIZE__ -- size of the bloom filter, 2**FILTER_SIZE, 38 = 32 GB, 37 = 16 GB, 36 = 8 GB, 35 = 4 GB. Multiply by 2 if you want a count 2 bloomfilter. Increase filter size if false positive rate >1%. 
* * __MAX JVM HEAP__ -- Maximum amount of memory DBG processes can use (must be larger than FILTER_SIZE); For example, if your filter size is ```36```, you can use ```100``` for this. 
* __MIN COUNT__=1 -- Minimum kmer occurrence in SEQFILE to be included in the final bloom filter

#####Contig Merge Parameters
MIN_BITS=50 --Mimimum assembled contigs bit score. This gives the quality of contig you want. It is not recomended to go below 50. 
MIN_LENGTH=150 -- Minimum assembled protein contigs. This determines how long your final contigs will be (min length+kmer length = minimum length). In general, keep this at or above 150. You may need to reduce it for very small proteins. 

More parameter descriptions can be found in the RDP's [Xander README] (https://github.com/rdpstaff/Xander_assembler) and in greater detail in the [Xander publication] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6). 


Now that we understand the paramters a bit more, lets edit our environment. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
nano xander_setenv.sh
```

Directories must match the absolute path that we are using, and we need to name our sample (sample shortname). The sample shortname will be the beginning of every results file.  

Delete lines 14-19 and paste in the following.

```
SEQFILE=/home/ubuntu/demo/data/demo_reads.fa
WORKDIR=/home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo
REF_DIR=/home/ubuntu/tools/RDPTools/Xander_assembler
JAR_DIR=/home/ubuntu/tools/RDPTools
UCHIME=/home/ubuntu/tools/third_party_tools/uchime4.2.40_i86linux32
HMMALIGN=/usr/local/bin/hmmalign
```

Save changes to this file and exit, and remember that you can adjust other parameters as well. 

We also need to change the security of this file so that we can excecute it. 

```
chmod 755 xander_setenv.sh
```

Now we are ready to roll!

###4 Run Xander
To run Xander, we use one simple command.

```
./run_xander_skel.sh xander_setenv.sh "build find search" "rplB"
```

This should take about 15-20 minutes to run with this dataset. It will take much longer (hours) with larger datasets.  

Note: if you wanted to run multiple genes at once, you would simply list them in quotations as done below. 
```
./run_xander_skel.sh xander_setenv.sh “build find search” “rplB nirK nirS”
```

###5 Analyze results
Xander has a lot of output files, so now we will go over which ones are important for your analysis. First let's check our false positive rate.

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo/k45
cat k45_bloom_stat.txt
```

The false positive rate is at the bottom of this file. This must be smaller than 1%. If it is larger than 1%, you should increase the filter size in ```xander_setenv.sh```. 

Now that we know our false positive rate, we can examine the results. These files will be in the following directory:

```
/home/ubuntu/tools/RDPTools/Xander_assembler/rplB_demo/k45/rplB/cluster
```

Here you will find the final output files: 

* ```demo_rplB_k45_final_nucl.fasta```: Quality filtered nucleotide sequences (representative) 
* ```demo_rplB_k45_final_prot.fasta```: Quality filtered protein sequences (representative) and raw abundance (number of contigs)
* ```demo_rplB_k45_final_prot_aligned.fasta```: Aligned protein sequences and raw abundance (number of contigs)
* ```demo_rplB_k45_Taxonabund.txt```: taxonomic abundance adjusted by coverage (```coverage.txt```), grouped by lineage (phylum/class)
* ```demo_rplB_k45_Framebot.txt```: Alignment of your contig with nearest reference sequence and % amino acid identity
* ```complete.clust```: Shows how many contigs you have with different distance cutoffs. (You'll have more contigs with lower distance cutoffs)

More output file descriptions can be found in the RDP's [Xander README] (https://github.com/rdpstaff/Xander_assembler) and in greater detail in the [Xander publication] (http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6). 

-------
-------
###Preparing Gene References
  1. The analysis pipeline is preconfigured with the _rplB_ phylogenetic marker gene, and nitrogen cycling genes including _nirK_, _nirS_, _nifH_, _nosZ_ and _amoA_. 
  2. Your gene of interest is on the RDP's FunGene database. 
  3. Your gene of interest is not preconfigured or on the FunGene database. 

Below, you will find instructions for each of these scenarios. 

####Scenario 1: 
The analysis pipeline is preconfigured with the _rplB_ phylogenetic marker gene, and nitrogen cycling genes including _nirK_, _nirS_, _nifH_, _nosZ_ and _amoA_. These require no work from you. See above tutorial. 

####Scenario 2:
 Your gene of interest is on the RDP's FunGene database. This requires a bit of work and biological insight. 
For each individual gene of interest, ```mkdir genename```, navigate to this gene directory, and ```mkdir originaldata```. 

```originaldata``` will hold the four required files for Xander. This step requires biological insight!
* **gene.seeds**: a small set of **full length, high quality** protein sequences in FASTA format, used to build gene.hmm, forward and reverse HMMs. Can be downloaded from RDP's FunGene database. 
* **gene.hmm**: this is the HMM built from gene.seeds using original HMMER3. This is used to build for_enone.hmm and align contigs after assembly. Can be downloaded from RDP's FunGene database. 
* **framebot.fa**: a large near full length known protein set for identifying start kmers and FrameBot nearest matching. More diversity is better, more sequences means more starting points (more computational time) but less susceptible to noise than model creation. Prefer near full-length and well-annotated sequences. Filter with Minimum HMM Coverage at least 80 (%).
* **nucl.fa**: a large near full length known set used by UCHIME chimera check.

In your web browser, download files from the [RDP's FunGene database] (http://fungene.cme.msu.edu). Here you will find a list of gene families based on searches of the NCBI non-redundant protein database using "training sequences" or **gene.seeds**.

To download the gene.seeds and gene.hmm files, click on the links at the top left of a page. Note this may not work with Safari.
(add image of fungene with links boxed in red)

To obtain the framebot.fa and nucl.fa files, click on the show/hide filter options link at the top right of the page. Here you can limit the the score, minimum size aa, and minimum HMM coverage. 
* The score will vary by gene. You may choose to leave this area blank. You can also manually look though the sequences to see if there is a large drop off in score at a specific number and use this number as a cutoff. 
* The minimum size aa will depend on the actual protein size. You'll want to find a balance between near-full length sequences and diversity. 
* The minimum HMM coverage must be greater than 80%. What percent cutoff you choose will require biological insight. Lower % HMM coverage will increase diversity but may lower the quality of your search. Once you have set your parameters, click filter. Then you will need to select sequences. It is not recommended to blindly select all sequences. Instead, manually go through and make sure no oddballs are included. For reference, the Xander paper used over 700 near full-length sequences for their .fa files. Once your sequences are selected, click begin analysis at the top right of the page. This will take you to a page where you can download your framebot.fa (protein) and nucl.fa (nucleotide) sequences. Then change the names to framebot.fa and nucl.fa respectively. Move all of these files to ```originaldata```. 

From these files, we need to make three new files
* **for_enone.hmm** and **rev_enone.hmm** for the forward and reverse HMMs respectively. This is used to assemble gene contigs.
* **ref_aligned.faa** file containing a set of protein reference sequences aligned with for_enone.hmm. This is used to identify starting kmers. Need to manually examine the alignment using Jalview or alignment viewing tool to spot any badly aligned sequences. If found, it is likely there are no sequences in gene.seeds close these sequences. You need to valide these problem sequences to see they are from the gene you are interested, then either remove them or add some representative sequences to gene.seeds and repeat the prepartion steps.

This can be accomplished using the ```prepare_gene_ref.sh``` script, but first, we need to edit this file. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/bin
nano prepare_gene_ref.sh
```

We need to change this file to reflect our gene of interest and our paths. The paths in this initial script are generalized from the RDP, but we want them to be specific to our AMI. 
* line 4: change gene to your gene of interest
* line 9: change the jar directory to ```/home/ubuntu/tools/RDPTools```
* line 10: change the reference directory to ```/home/ubuntu/tools/RDPTools/Xander_assembler```
* line 13: change the hmmer xanderpatch location to ```hmmer_xanderpatch=/home/ubuntu/tools/third_party_tools/hmmer-3.0-xanderpatch```

Now save the changes you've made to the shell script. 

Excellent! You now have a script ready to prepare the gene reference for Xander. Let's run it. All you need to do is excecute this script and specify your gene of interest. 

```
./prepare_gene_ref.sh rplB
```

Check and see that all of your output files are in your gene directory. 

```
cd /home/ubuntu/tools/RDPTools/Xander_assembler/gene_resource/rplB
```

You should now see three new files: **ref_aligned.faa**, **for_enone.hmm**, and **rev_enone.hmm**

Remember that when running Xander on your own gene of interest, you will need to manually check ```ref_aligned.faa``` for any poorly aligned sequences. This can be done using Jalview or another alignment viewing tool. 

####Scenario 3: 
Your gene of interest is not preconfigured or on the FunGene database. This requires a bit of research, work, and biological insight. 

You need to look through existing literature/ protein databases to find high quality sequences on your own (ideally from papers with protein crystal structures) and send them to RDP to make a FunGene database.

Once you have a FunGene database, you can proceed with scenario 2.
