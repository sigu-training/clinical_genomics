---
layout: tutorial_hands_on

title: Advanced data analysis of genomic data
zenodo_link: https://doi.org/10.5281/zenodo.3531577
questions:
- What is the minimal depth of coverage of target regions in a typical clinical genomics setting?
- Which are the most common computational strategies to identify copy-number alterations in NGS experiments?

objectives:
- To calculate sequencing statistics
- To manipulate genomic regions files
- To identify genomic mosaic variants
- To use webtools for genomic variants annotation
- To analyze CNV and Regions of Homozygosity (ROH)

time_estimation: 2h

contributors:
- abrusell
- aciolfi
- gmauro
- m-giuseppe
- puva
- tommasopippucci

---


# Introduction
{:.no_toc}

This section of the tutorial will cover in more details many different aspects of data analysis and interpretation for clinical genomics. It will delve into quality control procedures, and into specific strategies to be used to identify different variations of individuals' genetic architecture, as mosaic variants, copy number variants, and regions of homozygosity.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Requirements

This tutorial requires that the [basic course](https://sigu-training.github.io/clinical_genomics/)
has been completed successfully or that the user is familiar with the [Galaxy](https://galaxyproject.org/) platform.
In particular, we'll use the European Galaxy server running at [https://usegalaxy.eu](https://usegalaxy.eu).

# Datasets 

Input datasets are available:
 - at [Zenodo](https://zenodo.org/record/3531578), an open-access repository developed under the European OpenAIRE program and operated by CERN
 - as *Shared Data Libraries* in [Galaxy](https://usegalaxy.eu/library/list): *[Galaxy courses / Sigu](https://usegalaxy.eu/library/list#folders/F3d08bb711e4e3b26)*
   

---

# Quality control


## Computation of per-base coverage depth at specific genomic intervals

Commercial next-generation sequencing platforms usually provide users with analysis programs, which include tools for the identification of low coverage regions (for instance, target regions that have a coverage depth lower than 20x).

The present tutorial is aimed to show how to perform a custom coverage analysis of NGS experiments by using tools that are available in Galaxy.

Starting material:
* Alignment (*bam*) file(s) on which you want to perform the coverage evaluation.
* A reference *bed* file listing the genomic regions for which you want to obtain coverage data. If you performed a targeted sequencing experiment by using commercial kits (either custom or from the catalogue), you should already have obtained a *bed* file listing the target regions: it should be the file you want to use.

>    > ### {% icon comment %} BED format specifications
>    > **BED files** are tab-delimited files with one line for each genomic region.
>    > The first 3 fields (columns) are required: 
>    > 1. chromosome
>    > 2. the starting position
>    > 2. the ending position
>    >
>    > Further columns may be present but are optional.
>    > Additional details may be found here: [UCSC BED format specifications](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
>    >
>    > Please be aware that in BED files the description of genomic regions follows the “0-start, half-open” coordinate system. Further details may be found here: [the "0-based, half-open" UCSC Genome Browser Coordinate Counting Systems](http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/). Practically speaking, it means that the starting position is diminished by 1 bp with respect to the actual position. For instance, the region "chr2 50149071 50149399" in a *bed* file corresponds to the genomic interval from position 50149072 (included) to position 50149399 (included) on chromosome 2 (when the first base of the chromosome is defined as position 1). If you want to indicate a single base position in a *bed* file, it should be written like this: "chr2 50149071 50149072" .
>    {: .comment}

 
> ### {% icon hands_on %} Hands-on: Compute per base coverage depth with BEDtools
>    Before starting, you have to upload the files you need for the analysis,
>    following standard Galaxy procedure. 
>
>    You may use files provided as examples with this tutorial and called
>    `Panel_alignment.bam` and `Panel_Target_regions.bed`. 
>
>    Please check that uploaded file datatypes (formats) are properly recognized by
>    selecting `edit attributes` (i.e. the pencil sign in correspondence of each file
>    you can find in your history) (indicated by the red arrow in figure 1)
>    and then the tab datatypes (the blue arrow in figure 1).
>    If the datatype is wrong, select the correct one from the drop-down list
>    and press the `change datatype` button (the green arrow in figure 1).
>
>    ---
>    ![Figure 1]({{site.baseurl}}/images/cov_fig1.png)
>    **Figure 1**
>
>    ---
>
>    Once ready, you can select the tool named `bedtools Compute both the depth and
>    breadth of coverage` in the `Operate on Genomic Intervals` section
>    (see the red arrow in figure 2).
>
> 1. Select the *bed* file listing the target regions as "file A" (blue arrow in figure 2)
>    and the *bam* file(s) you want to analyze as "file B" (green arrow in figure 2)
>    (they should be listed in the drop-down menu if they have the correct format).
>    You may analyze one or more *bam* files in a single run. 
> 2. If you want to analyze two or more .bam files, you can further choose if you want
>    all the results in a single file or one output file per input "file B" by selecting
>    the desired option under the `Combined or separate output files` menu.
> 3. Select "Yes" for the option `Report the depth at each position in each A feature`
>   (yellow arrow in figure 2) and check that all the other options are set to "No".
> 4. Star (Execute) the analysis.
>
>    ![Figure 2]({{site.baseurl}}/images/cov_fig2.png)
>    **Figure 2**
>
>    ---
>
>    Output file, which will be called "coverage_depth.bed" from now on, 
>    will contain all the fields of the original target_regions.bed file plus
>    two further columns:
> 1. the first value indicates a specific base within the reported genomic interval.
>    For instance, if the genomic interval described by the first 3 field is
>    "chr2 50149071 50149072" and the first new field reports the number 1000,
>    it means that the coverage value refers to nucleotide 50150071
>    (i.e.: 50149071 + 1000) on chromosome 2;
> 2. the second value indicates the depth of coverage for the defined base.
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Sort files
>    Since entries in the "coverage_depth.bed" may not be in the desired order, you can sort it by genomic positions.
>    For this purpose you may want to use the `Sort` tool in the `Text Manipulation` section (check out the blue arrow in figure 3).
>
>    You may sequentially sort on different columns:
> 1. first you can sort by chromosome by selecting column 1 (green arrow in figure 3) in "ascending order" and selecting the "flavor" `Natural/Version sort (-V)`, which allows for sorting chromosomes in their "natural" order (with alphabetical order chr10 will be right after chr1, chr2 after chr19 and so on);
> 2. after inserting a new column section (red arrow in figure 3), you can sort by column 2 in "ascending order" with `Fast numeric sort (-n)`;
> 3. after inserting a further column section, you can sort by column 3 in "ascending order" with `Fast numeric sort (-n)`;
> 4. after inserting a final column section, you can sort by column 3 in "ascending order" with `Fast numeric sort (-n)`.
>
>    ---
>    ![Figure 3]({{site.baseurl}}/images/cov_fig3.png)
>    **Figure 3**
>
>    ---
>
>    The obtained output file, which will be called "sorted_coverage_depth.bed"
>    from now on, will be sorted first by chromosome, then by starting position,
>    by ending position and by the actual position of the specific base in the
>    considered genomic interval.
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Remove duplicate rows
>    If your file contains duplicated lines you may want to remove them for further processing.
>    For this purpose you can use the `Unique lines assuming sorted input file` tool in the `Text Manipulation` section.
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Some manipulation of the *bed* file 
>    You may follow the following steps to manipulate the "sorted_coverage_depth.bed" and to obtain a *bed* file listing the exact base position to which each coverage value is referred.
>    For instance instead of having "chr2 50149071 50149399 NRXN1 1 2335" in the first row of your file, you will get "chr2 50149071 50149072 NRXN1 2335".
>
>    These steps will add further columns at the end of each line defining the base position with the “0-start, half-open” coordinate system.
> 
> 1. Select the `Compute an expression on every row` tool in the `Text Manipulation` section (indicated by the blue arrow in figure 4);
> 1. add the following expression "c2+c5-1" to obtain the sum of the values in columns 2 and 5 minus 1 (it will be used as the new starting position in the final file);
> 1. select the file "sorted_coverage_depth.bed";
> 1. select "yes" to `round Results?` (green arrow in figure 4);
> 1. execute;
> 1. you can rename the output as "temp1_coverage_depth.bed";
> 1. select the `Compute an expression on every row` tool in the `Text Manipulation` section;
> 1. add the following expression "c2+c5" to obtain the sum of the values in columns 2 and 5 (it will be used as the new ending position in the final file). You can also use the expression ; 
> 1. select the file "temp1_coverage_depth.bed";
> 1. select "yes" to `round Results?`;
> 1. execute;
> 1. you can rename the output as "temp2_coverage_depth.bed";
> 1. select the `Table Compute` tool in the `Text Manipulation` section (indicated by the blue arrow in figure 5);
> 1. select the file "temp2_coverage_depth.bed";
> 1. select the option `Drop, keep or duplicate rows and columns` from the drop-down menu `Type of table operation` (indicated by the green arrow in figure 5);
> 1. fill the field `List of columns to select` with "1,7,8,4,6" (the red arrow in figure 5);
> 1. unselect all the other options;
> 1. execute;
> 1. set the output file datatype to *bed*;
> 2. you can rename the output as "final_coverage_depth.bed".
>
>    ---
>    ![Figure 4]({{site.baseurl}}/images/cov_fig4.png)
>    **Figure 4**
>
>    ![Figure 5]({{site.baseurl}}/images/cov_fig5.png)
>    **Figure 5**
>
>    ---
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Select positions with low coverage depth
> The following procedure can be used to obtain a *bed* file listing base positions with a coverage depth lower than a certain threshold (for instance 20x).
>
> 1. Select the `Filter` tool in the `Filter and Sort` section (indicated by the blue arrow in figure 6);
> 1. select the file "final_coverage_depth.bed";
> 1. add the following expression "c5<20" to filter all positions with a coverage depth lower than (green arrow in figure 6)("c5" stands for the fifth column, in this case reporting the coverage depth);
> 1. execute.
>
>    ---
>    ![Figure 6]({{site.baseurl}}/images/cov_fig6.png)
>    **Figure 6**
> 
>    ---
>
>    The output file, which will be called "low_coverage_depth.bed" from now on, will only list all the positions with a depth lower than 20x.
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Merge low coverage regions
>    If you want to merge the positions with low coverage depth in larger genomic intervals to be used for further analyses (i.e.: Sanger sequencing of regions not covered by your NGS experiment), you may want to use the `bedtools MergeBED` tool in the `Operate on Genomic Intervals` section (see the blue arrow in figure 7). 
>
> 1. Select the file "low_coverage_depth.bed";
> 1. set the maximum distance between features (i.e.: the different positions listed in your file) allowed for features to be merged (green arrow in figure 7): if it is 0 only overlapping and/or book-ended features are merged, while if it is set to 10 (or any other different positive integer of your choice), features at the maximum distance of 10 bases will be merged; 
> 1. if you want, you may apply some operations to certain columns to get further information in your output file. For instance you may:
>    1. click on "Insert Applying operations to columns from merged intervals" (red arrow in figure 7), 
>    1. specify the column on which you want to perform the operation (in this case column 5), 
>    1. and select the operation from the drop-down list(in this case "Min", which calculates the minimum value of coverage depth among all the positions that will be merged in a single interval) (yellow arrow in figure 7);
> 4. you may add as many operations as you need. In this example we will also calculate the maximum value of coverage depth;
> 5. execute.	
>
>    ---
>    ![Figure 7]({{site.baseurl}}/images/cov_fig7.png)
>    **Figure 7**
>
>    ---
>    The output file will have the following fields (columns): chromosome, starting and ending positions of low coverage regions, the minimum and the maximum coverage depth in each region.
>
>    > ### {% icon warning %} BED files may have different columns
>    > Please be aware that the columns to use for calculations may be different
>    > compared to the example here considered, depending on the amount of columns
>    > of your *bed* files.
>    {: .warning}
>
>    > ### {% icon warning %} Preview of BED files
>    > Please be aware that the Galaxy preview of your file shows a header row that
>    > does not properly define columns in your files
>    > (it is just a standard header for the UCSC bed format).
>    {: .warning}
>
{: .hands_on}

### Final notes
The procedures listed above are to be taken as examples of the possible operations that can be performed on bed files with bedtools (you may check out their website to get further information:
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)) ad text manipulation tools available on Galaxy.

Furthermore, please be aware that the tool `bedtools Compute both the depth and breadth of coverage` does not perform any filtering based on read quality: if your are interested in that aspect you may want to rely on different tools. 

---
# Mosaic variants

After the generation of a high-quality set of mapped read pairs, it could be useful to identify different classes of DNA variants in the analyzed sample.
Users interested in germline variant calling can refer to related Galaxy's tutorials, e.g. [Exome sequencing data analysis for diagnosing a genetic disease](https://galaxyproject.github.io/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html). On the other hand, it is also possible to accurately detect mosaic variants in sequenced samples, without the need of matched controls, using **MuTect2** tool from **[GATK toolkit](https://gatk.broadinstitute.org/hc/en-us)**. 

In more details, this tool executes different operations:

- Determine haplotypes by local assembly of the genomic regions in which the samples being analyzed show substantial evidence of variation relative to the reference
- Evaluate the evidence for haplotypes and variant alleles
- Assigning per-sample genotypes

Users could use files provided as examples with this tutorial and called `Panel_alignment.bam` and `Panel_Target_regions.bed`, and run **Mutect2** restricting the search space on target regions with "-L" option to reduce computational burden.

The first step of this procedure is needed to create an internal database of controls (i.e. **Panel Of Normals** - PoN) to reduce the bias for somatic calls.
It runs on a single sample at time:

 - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -L Panel_target_regions.bed -I Panel_alignment_normal1.bam -O normal_genotyped1.vcf`
 - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -L Panel_target_regions.bed -I Panel_alignment_normal2.bam -O normal_genotyped2.vcf`

Then users can take advantage of GATK's *CreateSomaticPanelOfNormals* tool to generate the PoN with the following commands:

- `gatk GenomicsDBImport -L Panel_target_regions.bed -R HSapiensReference_genome_hg19.fasta --genomicsdb-workspace-path PoN_db -V normal_genotyped1.vcf -V normal_genotyped2.vcf`
- `gatk CreateSomaticPanelOfNormals -R HSapiensReference_genome_hg19.fasta -V gendb://PoN_db -O panel_of_normals.vcf`

   > ### {% icon comment %} Note
   > The --genomicsdb-workspace-path must point to a non-existent or empty directory.
   {: .comment}

Then, to effectively call somatic mutations, users can use variants contained in the **PoN** and/or other public repositories  (e.g. by means of the option *--germline-resource*, using a VCF file containing frequencies of germline variants in the general population) to exclude germline variation. Finally, to properly classify somatic variants, *FilterMutectCalls filtering* could be applied to produce the final subset annotated VCF file, as described by the following commands:

 - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -I Panel_alignment.bam --germline-resource af-only-gnomad.vcf --panel-of-normals panel_of_normals.vcf -O somatic_genotyped_unfiltered.vcf`
 
 - `gatk FilterMutectCalls -R HSapiensReference_genome_hg19.fasta -V somatic_genotyped_unfiltered.vcf -O somatic_genotyped_filtered.vcf`

The VCF file obtained with this analysis can then be annotated by means of any annotation tools, as described below.

---
# Annotation and filtering with wANNOVAR

The web tool [wANNOVAR](http://wannovar.wglab.org/index.php) allows for rapid annotation of your variants and for some basic filtering steps to find disease genes.
It is based on its command line counterpart [ANNOVAR](http://annovar.openbioinformatics.org/), but it is more user-friendly since it does not require any programming skills.

The output consist in tabular text files that can be easily manipulated with Excel or other spreadsheet programs.

The annotation is performed against some of the most commonly used databases: RefSeq, UCSC Known, ENSEMBL/Genecode, dbSNP, ClinVar, 1000genomes, ExAC, ESP6500, gnomAD (minor allele frequencies in different populations are reported) and various precalculated prediction scores for any possible single nucleotide variant in coding regions (see [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)). 

The gene-based annotation results in a single row for each input variant: only the most deleterious consequence is reported (i.e.: if a certain variant may result to be missense for one transcript and nonsense for a second transcript, only the latter consequence will be reported).

Unfortunately, unlike the command line version, wANNOVAR does neither allow for the use of custom annotation databases, nor for the selection of different pubblicly available databases.

To annotate your file with wANNOVAR you need to provide your email address, to be notified when the annotation is complete, and just upload your input file (or paste a series of variant in the designated field). Results are usually ready within minutes.
![Figure 1]({{site.baseurl}}/images/wann_fig1.png)


You will get both *csv* files or *txt* files (which can be saved as they are by clicking on the specific link, right-clicking in any point of the page and selecting `Save as`). They can be both opened with Excel or other spreadsheet programs (concerning *csv* files, please ensure that comma are set as default list separator/delimiter in your version of Excel).

You will also get both `exome summary results` (only conding variants are included) and `genome summary results` (all variants included).

You can also provide a list of Disease or Phenotype Terms that the program can use for filtering your results (only for single sample `vcf` files).
![Figure 2]({{site.baseurl}}/images/wann_fig2.png)

Finally, there are some Parameter Settings that can be modified:
 - Result duration: it can now only be set to "1 day", since your files will be automatically removed after 24 hours.
 - Reference genome: you can choose between hg19 (GRCh37) and hg38 (GRCh38).
 - Input Format: you can upload not only *vcf* files, but also other kinds of variant files.
 - Gene Definition: the database you want to use for gene function annotation. Three oprtions are available: RefSeq, UCSC Known, ENSEMBL/Gencode
 - Individual analysis: the "Individual analysis" option allows you to perform further filtering steps (based on the Disease/Phenotype terms or the Disease Model) on a single sample (if you upload a multisample *vcf* only the first sample will be considered); the "All Annotaions" option will annotate all variants in your multisample *vcf*, maintaining the original columns of your *vcf* as the last columns of your output file.
 - Disease Model: this options allows for some basic filtering of your variants based on the expected mechanism of inheritance. In mainly consider frequencies, sample genotypes and consequences at level of genic function. FIltering step are summarized among the results and they are only performed on a single sample: the program does not perform any multisample evaluation (i.e.: variant segregation in a trio) and cannot classify any variant as de novo even if you provide a multisample *vcf* with parental genotypes.

## Clinical databases for further manual variant annotation

Once you have obtained a file with the annotation of your variants, you might find useful to annotate also the involved genes, in order to know, for instance, the list of diseases that may be associated with them.

Some databases that can the jb are the following:
 1. The gene2phenotype dataset ([G2P](https://www.ebi.ac.uk/gene2phenotype/disclaimer)) integrates data on genes, variants and phenotypes for example relating to developmental disorders. In the "Download" section you will find both databases of genes related to cancer and to developmental disorders. Those files report for each gene listed: the OMIM number, the associated disease name, the disease OMIM number, the disease category (you can fing more details in the "Terminology" section), whether the disease is caused by biallelic or monoallelic variants ("allelic requirement"), the expected category of variant to be causative of the disease, and a few other details.
 1. In the "Download" section of the [OMIM database](https://www.omim.org/downloads/), if you register for research and educational use, you may obtain different lists of OMIM genes and their associated phenotypes.
 1. Among the files available for download from the [gnomAD database](https://gnomad.broadinstitute.org/downloads#constraint), you may get per-gene constraint scores (for further details, please check the paper by the Exome Aggregation Consortium on "Nature. 2016 Aug 18; 536(7616):285–291."). Those score may indicate if a specific gene is expected to be intolerant to loss-of-function variants (pLI) (haploinsufficiency), or if it is predicted to be associate to recessive diseases.

In the end, you can add these annotations to your wANNOVAR files by using the `VLOOKUP` function in Excel.



# CNV detection from targeted sequencing data

- *Copy Number Variants* (CNVs) are imbalances in the copy number of the genomic material that result in either DNA **deletions** (copy loss) or **duplications**
- *CNVs* can cause/predispose to human diseases by altering gene structure/dosage or by **position effect**
- *CNV* size ranges from 50 bp to megabases
- Classical methods to identify *CNVs* use array-based technologies (SNP/CGH)
- Computational approaches have been developed to identify *CNVs* in targeted sequencing data from **hybrid capture** experiments

## Computational approaches

There are four main methods for *CNV* identification from short-read +NGS data (see figure below):
 - **Read Count** (RC)
 - **Read Pair** (RP)
 - **Split Read** (SR)
 - **De Novo Assembly** (AS)
  
  
 - *RP* and *SR* require continuous coverage of the *CNV* region or reads encompassing *CNV* breakpoints, as in whole genome sequencing. The sparse nature and small size of exonic targets hamper the application of these methods to targeted sequencing. 
 - *RC* is the best suited method for *CNV* detection from whole exomes or gene panels where:
  - deletions appear as exonic targets devoid of reads
  - duplications appear as exonic targets characterized by excess of coverage

---

![CNV detection methods]({{site.baseurl}}/images/methods_identification_cnv.png)

**Figure 1**. Methods for detection of CNVs in short read NGS data (adapted from [Tattini et al., 2015](https://doi.org/10.3389/fbioe.2015.00092))

---

## RC method and data normalization

In targeted sequencing, a method to study DNA copy number variation by *RC* (as implemented in *EXCAVATOR* tool, [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))) is to consider the **exon mean read count** (EMRC):

*EMRC* = *RC*<sub>e</sub>/*L*<sub>e</sub>

where *RC*<sub>e</sub> is the number of reads aligned to a target genomic region e and *L*<sub>e</sub> is the size of that same genomic region in base pairs ([Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))).
Three major bias sources are known to affect *EMRC* dramatically in targeted sequencing data:
- local **GC content** percentage
- genomic **mappability**
- target region **size**

These biases contribute to non uniform read depth across target regions and, together with the sparse target nature, challenge the applicability of *RC* methods to targeted data. As shown in Figure 2 (left panel), in single-sample data the *EMRC* distributions of genomic regions characterized by different copy numbers largely overlap, revealing poor *CNV* prediction capability. 

*EMRC ratio* between two samples can be used as a normalization procedure. The effect of *EMRC ratio*-based normalization is clear in Figure 2 (right panel) as a markedly improved correspondece between the predicted and the real copy number states.

---

![Data Normalization_1]({{site.baseurl}}/images/normalization_EMRC.png)

**Figure 2**. Effect of *EMRC ratio* on DNA copy number prediction (adapted from [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))

---
    
Similarly FPKM, a normalized measure of read depth implemented in ExomeDepth tool ([Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526)), is affected by extensive exon–exon variability but comparison between pairs of exome datasets demonstrates high correlation level of the normalized read count data across samples making it possibile to use one or more combined exomes as reference set to base the CNV inference on ([Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526)).

---
    
![Data Normalization_2]({{site.baseurl}}/images/normalization_FPKM.png)

**Figure 3**. Correlation of the normalized read count data between samples (from [Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526))

---    
    
## *CNV* detection accuracy

Strong correlation is observed between Affymetrix array-SNP and exome-derived data for *CNVs* >1 Mb (Figure 4, right panel), while including CNVs of any size dramatically decreases the correlation level (Figure 4, left panel). This can be explained by the different distribution of exons and SNP probes throughout the genome. Candidate *CNV* regions as identified by *EXCAVATOR* contain a comparable number of exons and SNP probes when they are > 1 Mb (*R*=0.8), while regions <100 Kb do not (*R*=-0.02). Accordingly, [Krumm et al., 2012](https://genome.cshlp.org/content/22/8/1525.full) report 94% precision in detecting *CNVs* across three or more exons.

---    
    
![Correlation_Exome_Affymetrix]({{site.baseurl}}/images/corr_exome_affymetrix.png)

**Figure 4**. Correlation between array-SNP and exome-derived *CNVs* including all (left panel) or > 1 Mb *CNVs* (adapted from [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))

---    
    
To increase chances to detect *CNVs* encompassing few exons or in non-coding regions from exome data, [D'Aurizio et al., 2016](https://academic.oup.com/nar/article/44/20/e154/2607979) have proposed an extension to *EXCAVATOR* to exploit off-target read count (*EXCAVATOR2*). Similar approaches taking advantage of off-target reads have been described by [Bellos et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4147927/) or [Talevich et al., 2016](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873).

## Tools for *CNV* detection from gene panels or exome data

A number of tools or pipelines implementing modified versions of previously published tools have been reported to detect single-exon *CNVs* in clinical gene panels:
- [CoNVaDING](https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22969)
- [DeCON](https://wellcomeopenresearch.org/articles/1-20/v1)
- [ExomeDepth v1.1.6](https://www.nature.com/articles/ejhg201742)
- [Atlas-CNV](https://www.nature.com/articles/s41436-019-0475-4)
- [DeviCNV](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6192323/)

In addition to the above-mentioned methods, many tools have been developed that detect *CNVs* from exomes. Evaluation of these tools is not straightforward as there is lack of a gold standard for comparison. As a consequence, there is no consensus level on pipelines as high as for single nucleotide variants. [Zare et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452530/pdf/12859_2017_Article_1705.pdf) have reviewed this topic with focus on cancer, reporting poor overlap of different tools and added challenges for somatic variant calling.

 
---    
  
> ### {% icon hands_on %} Hands-on: CNV calling with ExomeDepth
>    ExomeDepth uses pairs of case and control exomes to identify copy number imbalances in case(s) compared to control(s).
>    Here we provide a small data set including a pair of case/control BAM files and an exome target BED file restricted to a specific polymorphic region on chromosome 2q.
>
> - CNV_case.bam.
> - CNV_control.bam.
> - CNV_TruSeq_Chr2.bam
>
>    Your aim is to identify a large polymorphic deletion in case.
{: .hands_on}

# Regions of Homozygosity

- *Runs of Homozygosity* (ROHs) are sizeable stretches of consecutive homozygous markers encompassing genomic regions where the two haplotypes are identical
- Haplotypes in *ROHs* can underlie identity either by **state** (IBS) or by **descent** (IBD). IBD occurs when two haplotypes originate from a common ancestor (Figure 1), a condition whcih we refer to as **autozygosity**
- *Autozygosity* in an individual is the hallmark of high levels of genomic inbreeding (g*F*), often occurring in the offspring from consanguineous unions where parents are related as second-degree cousins or closer
- The most well-known medical impact of parental consanguinity is the increased risk of rare autosomal recessive diseases in the progeny, with excess risk inversely proportional to the disease-allele frequency ([Bittles, 2001](https://onlinelibrary.wiley.com/doi/full/10.1034/j.1399-0004.2001.600201.x?sid=nlm%3Apubmed))

---

![homozygosity_mapping]({{site.baseurl}}/images/homozygosity_mapping.png)
**Figure 1**. Schematic of haplotype flow along a consanguineous pedigrees to form autozygous *ROHs* in the progeny (from [McQuillan et al., 2008](https://www.sciencedirect.com/science/article/pii/S000292970800445X?via%3Dihub))

---

## Homozygosity mapping

- We refer to *homozygosity mapping* ([Lander and Botstein, 1987](https://science.sciencemag.org/content/236/4808/1567.long)) as an approach that infers *autozygosity* from the detection of long ROHs to map genes associated with recessive diseases
- A classical strategy to identify *ROHs* is based on SNP-array data and analysis with [PLINK](http://zzz.bwh.harvard.edu/plink/)
- The combination of *homozygosity mapping* with exome sequencing has boosted genomic research and diagnosis of recessive diseases in recent years ([Alazami et al., 2015](https://www.sciencedirect.com/science/article/pii/S2211124714010444?via%3Dihub); [Monies et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5502059/))
- *ROHs* may be of clinical interest also since they can unmask **uniparental isodisomy**
- Computational approaches have been developed to identify *ROHs* directly in exome data, allowing **simultaneous detection of variants and surrounding ROHs** from the same datasets.


## Computational approaches

As for CNVs, the sparse nature and small size of exonic targets are a challenge for the identification of continuous strethces of homozygous markers throughout the genome. This may affect the performance of tools originally tailored to SNPs when applied to exome data. 

In Figure 2 a synthetic overview is given of algorithms, input data and output file types in bioinformatic tools for exome-based *ROH* detection. *ROHs* are divided in three size classes reflecting different mechanisms that have shaped them (according to [Pemberton et al., 2012](https://www.sciencedirect.com/science/article/pii/S0002929712003230?via%3Dihub)). 

Usually, long ROH (approximately >1.5 Mb) arise as a result of close parental consanguinity, but also short-medium *ROHs* can be of medical interest in populations as involved in disease susceptibility, natural selection and founder effects ([Ceballos et al., 2018](https://www.nature.com/articles/nrg.2017.109)).

---

![roh_exome_methods]({{site.baseurl}}/images/roh_exome_methods.png)
**Figure 2**. Summary of the tools for *ROH* detection from exome data (from [Pippucci et al., 2014](https://www.karger.com/Article/Pdf/362412))

---

### H<sup>3</sup>M<sup>2</sup>

*H<sup>3</sup>M<sup>2</sup>* ([Magi et al., 2014](https://academic.oup.com/bioinformatics/article/30/20/2852/2422169)) is based on an heterogeneous [*Hidden Markov Model*](https://en.wikipedia.org/wiki/Hidden_Markov_model) (HMM) that incorporates inter-marker distances to detect *ROHs* from exome data.

*H<sup>3</sup>M<sup>2</sup>* calculates B-allele frequencies of a set of polymorphic sites throughout the exome as the ratio between allele *B* counts and the total read count at site *i*:

*BAF*<sub>i</sub> = *N<sub>b</sub>/N*

*H<sup>3</sup>M<sup>2</sup>* retrieves *BAF*<sub>i</sub> **directly from BAM files** to predict the heterozygous/homozygous genotype state at each polymorphic position *i*:

- when *BAF*<sub>i</sub> \~ 0 the predicted genotype is **homozygous reference**
- when *BAF*<sub>i</sub> \~ 0.5 the predicted genotype is **heterozygous**
- when *BAF*<sub>i</sub> \~ 1 the predicted genotype is **homozygous alternative**

*H<sup>3</sup>M<sup>2</sup>* models *BAF* data by means of the *HMM* algorithm to discriminate between regions of homozygosity and non-homozygosity according to the *BAF* distribution along the genome (Figure 3).

---
  
![h3m2_baf]({{site.baseurl}}/images/h3m2_baf.png)
**Figure 3**. ***BAF* data distribution** Panels a, b and c show the distributions of BAF values against the genotype calls generated by the HapMap consortium on SNP-array data ( a ), the genotype calls made by SAMtools ( b ) and the genotype calls made by GATK ( c ). For each genotype caller, the distribution of BAF values is reported for homozygous reference calls (HMr), heterozygous calls (HT) and homozygous alternative calls (HMa). R is the Pearson correlation coefficient. Panels d–g show the distribution of BAF values in all the regions of the genome ( d ), in heterozygous regions ( e ), in homozygous regions ( f ) and in the X chromosome of male individuals ( g ). For each panel, the main plot reports the zoomed histogram, the left subplot shows the BAF values against genomic positions, whereas the right subplot shows the entire histogram of BAF values. (from [Magi et al., 2014](https://academic.oup.com/bioinformatics/article/30/20/2852/2422169))


> ### {% icon trophy %} **Congratulations!**
> You successfully completed the advanced tutorial.
{: .comment}


# Contributors
{:.no_toc}

 * [Tommaso Pippucci](https://www.aosp.bo.it/content/curriculum?E=154659) - Sant’Orsola-Malpighi University Hospital, Bologna, Italy
 * [Alessandro Bruselles](https://www) - Istituto Superiore di Sanità, Rome, Italy
 * [Andrea Ciolfi](http://www.ospedalebambinogesu.it) - Ospedale Pediatrico Bambino Gesù, IRCCS, Rome, Italy
 * [Gianmauro Cuccuru](https://gmauro.github.io) - Albert Ludwigs University, Freiburg, Germany
 * [Giuseppe Marangi](http://www) - Institute of Genomic Medicine, Fondazione Policlinico Universitario A. Gemelli IRCCS, Università Cattolica del Sacro Cuore, Roma, Italy
 * [Paolo Uva](https://www.researchgate.net/profile/Paolo_Uva) - IRCCS G. Gaslini, Genoa, Italy

