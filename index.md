---
layout: tutorial_hands_on

title: Data analysis and interpretation for clinical genomics
zenodo_link: https://doi.org/10.5281/zenodo.3531577
questions:
- What are the specific challenges for the interpretation of sequencing data in the clinical setting?
- How can you annotate variants in a clinically-oriented perspective?
objectives:
- Perform in-depth quality control of sequencing data at multiple levels (fastq, bam, vcf)
- Call, classify and annotate variants with information extracted from public databases for clinical interpretation
- Analyze CNV and Regions of Homozygosity (ROH)

time_estimation: 6h
key_points:
- Focus on clinical interpretation of variants.
- Provides real uses cases.
- Use public tools and free annotation servers.
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

In years 2018-2019, on behalf of the Italian Society of Human Genetics ([SIGU](https://www.sigu.net/)) an itinerant [Galaxy](https://usegalaxy.eu/)-based “hands-on-computer” training activity entitled “Data analysis and interpretation for clinical genomics” was held four times on invitation from different Italian institutions (Università Cattolica del Sacro Cuore in Rome, University of Genova, SIGU 2018 annual scientific meeting in Catania, University of Bari) and was offered to about 30 participants each time among clinical doctors, biologists, laboratory technicians and bioinformaticians. Topics covered by the course were NGS data quality check, detection of variants, copy number alterations and runs of homozygosity, annotation and filtering and clinical interpretation of sequencing results.

Realizing the constant need for training on NGS analysis and interpretation of sequencing data in the clinical setting, we designed an on-line [Galaxy](https://usegalaxy.eu/)-based training resource articulated in presentations and practical assignments by which students will learn how to approach NGS data quality at the level of fastq, bam and VCF files and clinically-oriented examination of variants emerging from sequencing experiments.

This training course is not to be intended as a tutorial on NGS pipelines and variant calling. This on-line training activity is indeed focused on data analysis for clinical interpretation. If you are looking for training on variant calling, visit this **Galaxy** tutorial on [Exome sequencing data analysis for diagnosing a genetic disease](https://galaxyproject.github.io/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html).


>    > ### {% icon comment %} SIGU
>    > The *Italian Society of Human Genetics* (**SIGU**) was established on November 14,
>    > 1997, when the pre-existing Italian Association of Medical Genetics and the Italian
>    > Association of Medical Cytogenetics joined.
>    > SIGU is one of the 27 member societies of FEGS (Federation of European Genetic
>    > Societies).
>    > Animated by a predominant scientific spirit, SIGU wants to be reference for
>    > all health-care issues involving human genetics in all its applications.
>    > Its specific missions are to develop quality criteria for medical genetic
>    > laboratories, to promote writing of guidelines in the field of human genetics
>    > and public awareness of the role and limitations of genetic diagnostic techniques. 
>    > SIGU coordinates activities of several working groups: Clinical Genetics,
>    > Cytogenetics, Prenatal Diagnosis, Neurogenetics, Fingerprinting, Oncological 
>    > Genetics, Immunogenetics, Genetic Counseling, Quality Control, Medical Genetics 
>    > Services, Bioethics. More than 1000 medical geneticists and biologists are active
>    > members of the society.
>    {: .comment}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Requirements

This tutorial is based on the [Galaxy](https://galaxyproject.org/) platform,
therefore a basic knowledge of Galaxy is required to get most out of the course.
In particular, we'll use European Galaxy server running at [https://usegalaxy.eu](https://usegalaxy.eu).

Registration is **free**, and you get access to **250GB** of disk space for your analysis.

1. Open your browser. We recommend Chrome or Firefox (please don't use Internet Explorer or Safari).
1. Go to [https://usegalaxy.eu](https://usegalaxy.eu)
   - If you have previously registered on this server just log in:
     - On the top menu select: User -> Login
     - Enter your user/password
     - Click Submit
   - If you haven’t registered on this server, you’ll need to do now.
     - On the top menu select: User -> Register
     - Enter your email, choose a password, repeat it and add a one word name (lower case)
     - Click Submit

To familiarize with the Galaxy interface (e.g. working with histories, importing dataset),
we suggest to follow the [Galaxy 101](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-101/tutorial.html) tutorial.


# Datasets 

Input datasets used in this course are available:
 - at [Zenodo](https://zenodo.org/record/3531578), an open-access repository developed under the European OpenAIRE program and operated by CERN
 - as *Shared Data Libraries* in [Galaxy](https://usegalaxy.eu/library/list): *[Galaxy courses / Sigu](https://usegalaxy.eu/library/list#folders/F3d08bb711e4e3b26)*
   
## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a meaningful name (e.g. Clinical genomics)
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import files from [Zenodo](https://zenodo.org/record/3531578) or Shared Data Library:
>
>    ```
>    https://zenodo.org/record/3531578/files/CNV_case.bam
>    https://zenodo.org/record/3531578/files/CNV_control.bam
>    https://zenodo.org/record/3531578/files/CNV_TruSeq_Chr2.bed
>    https://zenodo.org/record/3531578/files/HighQuality_Reads.fastq.gz
>    https://zenodo.org/record/3531578/files/LowQuality_Reads.fastq.gz
>    https://zenodo.org/record/3531578/files/Panel_alignment.bam
>    https://zenodo.org/record/3531578/files/Panel_target_regions.bed
>    https://zenodo.org/record/3531578/files/Sample1.all_exons.hg19.vcf
>    ```
>    {% include snippets/import_via_link.md %} 
>
>    The same files may be available on the Galaxy server
>    through a *Shared Data Libraries* in 
>    *[Galaxy courses/Sigu](https://usegalaxy.eu/library/list#folders/F3d08bb711e4e3b26)*.
>    You may prefer to import the data directly from there.
>
>    {% include snippets/import_from_data_library.md %}
>
>    > ### {% icon comment %} Note
>    > All the files are based on `hg19` reference genome which is
>    > available with pre-built indexes for widely used tools such as 
>    > *bwa-mem* **and** *samtools* by selecting `hg19` version as an option under
>    > *"(Using) reference genome"*).
>    {: .comment}
>
> 3. In case you import datasets from Zenodo, check that all datasets in your history
>    have their datatypes assigned correctly, and fix it when necessary.
>    For example, to assign BED datatype do the following:
>     
>    {% include snippets/change_datatype.md datatype="bed" %}
>
> 4. Rename the datasets
>
>    For datasets uploaded via a link, Galaxy will use the link
>    as the dataset name. In this case you may rename datasets.
>
>    {% include snippets/rename_dataset.md %}
>
>
{: .hands_on}


# Next Generation Sequencing

*Next (or Second) Generation Sequencing* (NGS/SGS) is an umbrella-term covering a number of approaches to DNA sequencing that have been developed after the first, widespread and for long time most commonly used Sanger sequencing.

*NGS* is also known as *Massive Parallel Sequencing* (MPS), a term that makes explicit the paradigm shared by all these technologies, that is to sequence in parallel a massive library of spatially separated and clonally amplified DNA templates. 

For a comprehensive review of the different *NGS* technologies see [Goodwin et al., 2016](https://www.nature.com/articles/nrg.2016.49), which also includes an introduction to the third generation methods allowing sequencing of long single-molecule reads.

## NGS in the clinic

In the span of less than a decade, NGS approaches have pervaded clinical laboratories revolutionizing genomic diagnostics and increasing yield and timeliness of genetic tests.

In the context of disorders with a recognized strong genetic contribution such as neurogenetic diseases, *NGS* has been firmly established as the strategy of choice to rapidly and efficiently diagnose diseases with a Mendelian basis. A general diagnostic workflow for these disorders currently embraces different *NGS*-based diagnostic options as illustrated in **Figure 1**.

---

![ngs_neuro_diagnostics]({{site.baseurl}}/images/ngs_neuro_diagnostics.png)

 **Figure 1**. General workflow for genetic diagnosis of neurological diseases. (\*If considering high-yield single-gene testing of more than 1–3 genes by another sequencing method, note that next-generation sequencing is often most cost-effective. †Genetic counselling is required before and after all genetic testing; other considerations include the potential for secondary findings in genomic testing, testing parents if inheritance is sporadic or recessive, and specialty referral.) From [Rexach et al., 2019](https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(19)30033-X/fulltext)

---

Currently, most common *NGS* strategies in clinical laboratories are the so-called **targeted sequencing** methods that, as opposed to **genome sequencing** covering the whole genomic sequence, focus on a pre-defined set of regions of interest (the **targets**). The targets can be selected by **hybrid capture** or **amplicon sequencing**, and the target-enriched libraries is then sequenced. The most popular target designs are:

 - **gene panels** where the coding exons of only a clinically-relevant group of genes are targeted
 - **exome sequencing** where virtually all the protein-coding exons in a genome are simultaneously sequenced

## Basics of NGS bioinformatic analysis

Apart from the different width of the target space in exome and gene panels, these two approaches usually share the same experimental procedure for NGS library preparation. After clonal amplification, the fragmented and adapter-ligated DNA templates are sequenced from both ends (**paired-end** sequencing) of the *insert* to produce short reads in opposite (**forward** and **reverse**) orientation. 

Bioinformatic analysis of NGS data usually follows a general three-step workflow to variant detection. Each of these three steps is marked by its "milestone" file type containing sequence data in different formats and metadata describing sequence-related information collected during the analysis step that leads to generation of that file.

| NGS workflow step | File content | File format | File Size (individual exome) |
| --- | --- | --- | --- |
| Sample to reads | Unaligned reads and qualities | **fastQ** | gigabytes |
| Reads to alignments | Aligned reads and metadata | **BAM** | gigabytes |
| Alignments to variants | Genotyped variants and metadata | **VCF** | megabytes |


Here are the different formats explained: 

- **[fastQ](https://en.wikipedia.org/wiki/FASTQ_format)** (sequence with quality): the *de facto* standard for storing the output of high-throughput sequencing machines
  - Usually not inspected during data analysis
- **[BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)** (binary sequence alignment/map): the most widely used TAB-delimited file format to store alignments onto a reference sequence
  - Aligned reads
- **[VCF](http://samtools.github.io/hts-specs/VCFv4.3.pdf)** (variant call format): the standard TAB-delimited format for genotype information associated with each reported genomic position where a variant call has been recorded

Another useful file format is [BED](https://www.ensembl.org/info/website/upload/bed.html), to list genomic regions of interest such as the exome or panel targets.

The steps of the ***reads-to-variants*** workflow can be connected through a bioinformatic **pipeline** ([Leipzig et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429012/)), consisting in read alignments, post-alignment BAM processing and variant calling.


## Alignment

As generated by the sequencing machines, *paired-end* reads are written to two *fastQ* files in which *forward* and *reverse* reads are stored separately together with their qualities. *FastQ* files are taken as input files by tools (the **aligners**) that align the reads onto a reference genome. One of the most used aligners is [BWA](http://bio-bwa.sourceforge.net/) among the many that have been developed (**Figure 2**).

---

![aligners_timeline]({{site.baseurl}}/images/aligners_timeline.png)

**Figure 2**. Aligners timeline 2001-2012 (from [Fonseca et al., 2012](https://academic.oup.com/bioinformatics/article/28/24/3169/245777))

---
 
During the bioinformatic process, *paired-end* reads from the two separate *fastQ* files are re-connected in the alignment, where it is expected that they will:
- map to their correct location in the genome
- be as distant as the insert size of the fragment they come from
- be in opposite orientations
a combination which we refer to as **proper pairing**. All these data about *paired-end* reads mapping are stored in the BAM file and can be used to various purposes, from alignment quality assessment to structural variant detection.

In **Figure 3**, the [Integrative Genomic Viewer (IGV)](http://software.broadinstitute.org/software/igv/) screenshot of an exome alignment data over two adjacent *ASXL1* exons is shown. Pink and violet bars are *forward* and *reverse* reads, respectively. The thin grey link between them indicates that they are *paired-end* reads. The stack of reads is concentrated where exons are as expected in an exome, and the number of read bases covering a given genomic location *e* (depicted as a hill-shaped profile at the top of the figure) defines the **depth of coverage (DoC)** over that location:

*DoC*<sub>e</sub>=*number of read bases over e/genomic length of e*

---

![coverage]({{site.baseurl}}/images/coverage.png)
**Figure 3**. Exome data visualization by [*IGV*](http://software.broadinstitute.org/software/igv/)

---

## Post-alignment BAM processing

Regarding post-alignment *pipelines*, the most famous for germline SNP and InDel calling is probably that developed as part of the GATK toolkit (**Figure 4**).

---

![gatk_germline_pipeline]({{site.baseurl}}/images/gatk_germline_snps_indels.png)

**Figure 4**. After-alignment pipeline for germline SNP and InDel variant calling according to [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145)

---

    
According to GATK best practices, in order to be ready for variant calling the BAM file should undergo the following processing:
- [marking duplicate reads](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php) to flag (or discard) reads that are mere optical or PCR-mediated duplications of the actual reads
- [recalibrating base quality scores](https://software.broadinstitute.org/gatk/documentation/article?id=11081) to correct known biases of the native base-related qualities
While GATK BAM processing is beyond doubt important to improve data quality, it has to be noticed that it is not needed to obtain variant calls and that non GATK-based pipelines may not use it or may use different quality reparametrization schemes. Duplicate flagging or removal is not recommended in *amplicon* sequencing experiments.

## Variant calling

The process of variant detection and genotyping is performed by *variant callers*. These tools use probabilistic approaches to collect evidence that non-reference read bases accumulating over a given locus support the presence of a variant, usually differing in algorithms, filtering strategies, recommendations ([Sandmann et al., 2017](https://doi.org/10.1038/srep43169)). To be confident that a variant is a true event, its supporting evidence should be significantly stronger than chance; e.g. the C>T on the left of the screenshot in **Figure 5** is supported by all its position-overlapping reads, claiming for a variant. In contrast, the C>A change on the right of the screenshot is seen only once over many reads, challenging its interpretation as a real variant. In fact, DNA variants that occur in germ cells (i.e., **germline/constitutional variants** that can be passed on to offspring) are either diploid/biallelic, so expected alternative allele frequency is 50% for a heterozygous change. On the other hand, if only a smaller subset of aligned reads indicates variation, that could result from technology bias or be a **mosaicism**, i.e. an individual which harbour two or more populations of genetically distinct cells as a result of postzygotic mutation. Postzygotic *de novo* mutations may result in **somatic mosaicism** (potentially causing a less severe and/or variable phenotype compared with the equivalent constitutive mutation) **and/or germline mosaicism** (hence enabling transmission of a pathogenic variant from an unaffected parent to his affected offspring) ([Biesecker et al., 2013](https://doi.org/10.1038/nrg3424)). To identify mosaicism, a probabilistic approach should consider deviation of the proband variant allele fraction (VAF, defined as the number of alternative reads divided by the total read depth) from a binomial distribution centred around 0.5.

---

![igv_screenshot_variant]({{site.baseurl}}/images/igv_variant.png)

**Figure 5**. Variant visualization by *IGV*

---

The GATK variant calling pipeline first produces a **genomic VCF ([gVCF](https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf))**, whose main difference with *VCF* is that it records all sites in a target whether there is a variant or not, while *VCF* contains only information for variant sites, preparing multiple samples for **joint genotyping** and creation of a **multi-sample** *VCF* whose variants can undergo quality **filtering** in order to obtain the final set of quality-curated variants ready to be annotated.

In downstream analyses, annotations can be added to *VCF* files themselves or information in *VCF* files can be either annotated in TAB- or comma- deimited files to be visually inspected for *clinical* variant searching or used as input to **prioritization** programs.


# Quality control

In-depth quality control (QC) of data generated during an NGS experiment is crucial for an accurate interpretation of the results. For example an accurate QC could help in identifying poor quality experiments, sequence contamination or genomic regions with low sequence coverage, and all these factors have a large impact on the downstream processing. 

Most of the programs used during an NSG workflow do not include steps for quality control, therefore artifacts needs to be detected using ad-hod developed tools for QC. The table  summarizes the main tools available at [https://usegalaxy.eu](https://usegalaxy.eu) for quality checking, at each step of the analysis.

| NGS workflow step | File format | Tools for quality control |
| --- | --- | --- |
| Sample to reads | fastQ | **FastQC** |
| Reads to alignments |  BAM | General statistics: **bam.io.bio**, **samtools**; target coverage: **Picard CollectHSMetrics**; per base coverage depth: **bedtools** |
| Alignments to variants | VCF | **vcf.io.bio** |

---

## Quality control of FASTQ files

Before starting the analysis workflow, you should identify possible issues that could affect alignment and variant calling. This first step of quality control is based on the raw sequence data (fastQ) generated by the sequencer. Common issues with sequence quality can be easily addressed by further processing your original sequences to trim or remove low-quality reads. In presence of severe artefacts you should consider to repeat the experiment instead of starting the downstream analysis that will generate poor quality results, according to the rule 'garbage in, garbage out'.

Here we'll use the **FastQC** software for a standard quality check, using the two FASTQ files *HighQuality_Reads.fastq* and *LowQuality_Reads.fastq*.

*FastQC* is relatively easy to use. The output of **FastQC** consists of multiple modules analysing a specific aspect of the quality of the data. A detailed help can be found in the [help manual](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

The names of the modules are preceded by an icon that reflects the quality of the data, and indicates whether the results of the module are:
 - normal (green)
 - slightly abnormal (orange)
 - very unusual (red)

> ### {% icon comment %} Note on FastQC interpretation
> These evaluations must be taken in the context of what you are expecting from your
> dataset. For **FastQC** a *normal* sample includes random sequences with high diversity. 
> If your experiment generates biased libraries (e.g. low complexity
> libraries) you should interpret the report with attention. In general, you should
> concentrate on the icons different from green and try to understand the reasons for
> this behaviour.
{: .comment}


> ### {% icon hands_on %} Hands-on: Compute sequence quality with FastQC
> 1. Run **FastQC** {% icon tool %} on your fastq datasets *HighQuality_Reads.fastq* and
> *LowQuality_Reads.fastq*. You can select both datasets with the  **Multiple datasets**
> option.
>
>    {% include snippets/select_multiple_datasets.md %}
>
>    For each input file you will get two datasets, one with the raw QC statistics and
>    another with an HTML report with figures.
>
> 1. Using the **MultiQC** {% icon tool %} software, you can aggregate multiple raw
> **FastQC** output in one unique report. This helps in comparing multiple samples at
> the same time in order to quickly identify low quality samples that will be displayed
> as outliers.  
>     - *"Which tool was used generate logs?"*: `FastQC`
>       - In *"FastQC output"*
>          - *"Type of FastQC output?"*: `Raw data`
>          - {% icon param-files %} *"FastQC output"*: the two *RawData*
>            outputs of **FastQC** {% icon tool %}
>
> 1. Inspect **MultiQC** report
>
> For a detailed explanation of the different analysis modules of **FastQC** you may refer
> to the [Quality control](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html) tutorial.
> 
>    > ### {% icon question %} Questions
>    >
>    > 1. Based on the **MultiQC** report, check which modules highlight differences in sequence quality betwen the two datasets
>    {: .question}
{: .hands_on}

## Quality control of BAM files

### Fast quality check with bam.iobio.io

BAM files are binary files containing information on the sequences aligned onto a reference genome. Exploring BAM files you can address several questions, e.g.:

- which is the amount of duplicated sequences? For non PCR-free protocols, it should be < 15%. Duplicated sequences are not used in downstream analysis to identify variants, therefore should be kept at a minimum to avoid waste of reagents.
- which is the fraction of unmapped reads? It should be < 2%. If higher, you should ask why so many reads are not properly mapped onto the reference genome. One possible reason could be sample contamination.

In Galaxy, BAM files can be explored using the *bam.iobio.io* web app. Leveraging on random subsampling of reads, **bam.iobio.io** quickly draws several quality control reports (**Figure 6**). 

On top of each plot, clicking on the question mark you can open a window with a detailed explanation of the expected output. The number of reads sampled is shown at the top-right of the page, and can be increased by clicking on the arrow.

![bam_iobio]({{site.baseurl}}/images/bamiobio.png)

**Figure 6**. BAM quality control using *bam.iobio.io*

> ### {% icon hands_on %} Hands-on: Compute BAM quality with bam.iobio.io
> 1. Run **bam.iobio.io** {% icon tool %} on a BAM dataset. To start **bam.iobio.io**
>    click on the link `diplay at bam.iobio.io` in the dataset section. Please note
>    that the link will be visible only for datasets with the appropriate database 
>    field set to `hg19`
> 
>    {% include snippets/change_dbkey.md dbkey="hg19" %}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which is the amount of duplicated sequences?
>    > 1. And the fraction of aligned reads?
>    {: .question}
{: .hands_on}

### Compute statistics with Picard CollectHsMetrics

**Collect Hybrid Selection (HS) Metrics** {% icon tool %} tool computes a set of metrics that are specific for sequence datasets generated through hybrid-selection. Hybrid selection is the commonly used protocol to capture specific sequences for targeted experiments such as exome sequencing.

In order to run this tool you need a file with the aligned sequences in BAM format, and files with the intervals corresponding to bait and target regions. These files can be generally obtained from the website of the kit manufacturer. 

> ### {% icon comment %} Note
> Please note that interval files are generally available as BED files,
> and must be converted in Picard interval_list format using Picard's
> **BedToInterval** {% icon tool %} before running CollectHsMetrics - see
> the hands on below for details.
{: .comment}

Metrics generated by CollectHsMetrics are grouped into three classes:
* Basic sequencing metrics: genome size, the number of reads, the number of aligned reads, the number of unique reads, etc.
* Metrics for evaluating the performance of the wet-lab protocol: number of bases mapping on/off/near baits, number of bases mapping on target, etc.
* Metrics for coverage estimation: mean bait coverage, mean and median target coverage, the percentage of bases covered at different coverage (e.g. 1X, 2X, 10X, 20X, ...), the percentage of filtered bases, etc.

For a detailed description of the output see [Picard's CollectHsMetrics](http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics)

In the next tutorial we will compute hybrid-selection metrics for BAM files containing 
aligned sequences from an exome sequencing experiment.


> ### {% icon hands_on %} Hands-on: Compute BAM statistics with Picard CollectHsMetrics
>
> 1. Before computing the statistics, we first need to convert the BED files
>    with bait and target regions, in Picard interval_list format.
>    Remember that you can select multiple datasets
>    with the **Multiple datasets** option.
>    
>    {% include snippets/select_multiple_datasets.md %}
>
>    Run **BedToIntervalList** {% icon tool %} to convert BED files.
>     - *"Load picard dictionary file from?"*: `Local cache`
>       - In *"Use dictionary from the list"*: `Human (Homo sapiens): hg19`
>     - *"Select coordinate dataset or dataset collection?"*: your BED file to be converted
>
> 1. Run **CollectHsMetrics** {% icon tool %} using as input the BAM files and the intervals
>    in Picard interval_list format, corresponding to the bait and target regions,
>    generated in the previous step.
>
> 1. Use **MultiQC** {% icon tool %} to aggregate CollectHsMetrics output in one unique
>    report to facilitate the comparison across multiple samples.  
>     - *"Which tool was used generate logs?"*: `Picard`
>       - In *"Picard output"*
>          - *"Type of Picard output?"*: `HS Metrics`
>          - {% icon param-files %} *"Picard output"*: the output of ***CollectHsMetrics***
>            {% icon tool %}
>
> 1. Inspect **MultiQC** {% icon tool %} report
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which is the average target coverage?
>    > 1. And the fraction of bases covered at least 10X?
>    {: .question}
{: .hands_on}


### BAM statistics with samtools

Samtools is a widely used suite of programs for manipulating alignments in the SAM/BAM/CRAM format.
Among the different tools, we will focus on **samtools flagstat** and **samtools stats** for computing BAM statistics.
Both programs take as input a file with the aligned sequences, and generate an output in text format that can be
visualized with **MultiQC**.

Here is an example of simple statistics obtained with **samtools flagstat**:

```
492739 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
476344 + 0 mapped (96.67% : N/A)
492739 + 0 paired in sequencing
243667 + 0 read1
249072 + 0 read2
462914 + 0 properly paired (93.95% : N/A)
469984 + 0 with itself and mate mapped
6360 + 0 singletons (1.29% : N/A)
5061 + 0 with mate mapped to a different chr
2034 + 0 with mate mapped to a different chr (mapQ>=5)
```

For more detailed statistics, use **samtools stats**:
```
SN	raw total sequences:	492739																																										
SN	filtered sequences:	0																																										
SN	sequences:	492739																																										
SN	is sorted:	1																																										
SN	1st fragments:	243667																																										
SN	last fragments:	249072																																										
SN	reads mapped:	476344																																										
SN	reads mapped and paired:	469984	# paired-end technology bit set + both mates mapped																																									
SN	reads unmapped:	16395																																										
SN	reads properly paired:	462914	# proper-pair bit set																																									
SN	reads paired:	492739	# paired-end technology bit set																																									
SN	reads duplicated:	0	# PCR or optical duplicate bit set																																									
SN	reads MQ0:	111721	# mapped and MQ=0																																									
SN	reads QC failed:	0																																										
SN	non-primary alignments:	0																																										
SN	total length:	49766639	# ignores clipping																																									
SN	total first fragment length:	24610367	# ignores clipping																																									
SN	total last fragment length:	25156272	# ignores clipping																																									
SN	bases mapped:	48110744	# ignores clipping																																									
SN	bases mapped (cigar):	45368693	# more accurate																																									
SN	bases trimmed:	0																																										
SN	bases duplicated:	0																																										
SN	mismatches:	254006	# from NM fields																																									
SN	error rate:	5.598707e-03	# mismatches / bases mapped (cigar)																																									
SN	average length:	101																																										
SN	average first fragment length:	101																																										
SN	average last fragment length:	101																																										
SN	maximum length:	101																																										
SN	maximum first fragment length:	101																																										
SN	maximum last fragment length:	101																																										
SN	average quality:	31.5																																										
SN	insert size average:	247.9																																										
SN	insert size standard deviation:	75.8																																										
SN	inward oriented pairs:	231769																																										
SN	outward oriented pairs:	452																																										
SN	pairs with other orientation:	240																																										
SN	pairs on different chromosomes:	2530																																										
SN	percentage of properly paired reads (%):	93.9
```

> ### {% icon hands_on %} Hands-on: Compute alignment statistics with samtools
> 1. Run **samtools flagstat** {% icon tool %} on your BAM dataset *CNV_case.bam*.
>
> 1. Using the **MultiQC** {% icon tool %} software, you can aggregate and visualize results
> obtained with **samtools flagstat** to identify low quality samples that will be displayed
> as outliers.  
>     - *"Which tool was used generate logs?"*: `Samtools`
>       - In *"Samtools output"*
>          - *"Type of Samtools output?"*: `flagstat`
>          - {% icon param-files %} *"Samtools flagstat output"*: the output 
>            of **Samtools flagstat** {% icon tool %}
>
> 1. Inspect **MultiQC** report
> 
>    > ### {% icon question %} Questions
>    >
>    > 1. Which is the percentage of mapped reads? And
>    >    the percentage of properly mapped reads?
>    {: .question}
{: .hands_on}


### Computation of per-base coverage depth at specific genomic intervals

Commercial next-generation sequencing platform usually provide users with analysis programs that include tools for the identification of low coverage regions (for instance, target regions that have a coverage depth lower than 20x).

The present tutorial is aimed to show how to perform a custom coverage analysis
of NGS experiment(s) by using tools that are available in Galaxy.

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

#### Final notes
The procedures listed above are to be taken as examples of the possible operations that can be performed on bed files with bedtools (you may check out their website to get further information:
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)) ad text manipulation tools available on Galaxy.

Furthermore, please be aware that the tool `bedtools Compute both the depth and breadth of coverage` does not perform any filtering based on read quality: if your are interested in that aspect you may want to rely on different tools. 

# Variant calling and classification

After the generation of a high-quality set of mapped read pairs, we can proceed to call different classes of DNA variants.
Users interested in germline variant calling can refer to related Galaxy's tutorials, e.g. [Exome sequencing data analysis for diagnosing a genetic disease](https://galaxyproject.github.io/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html).
To accurately detect mosaic variants in sequencing data without matched controls we will use **MuTect2** tool from **[GATK toolkit](https://gatk.broadinstitute.org/hc/en-us)**. 

In more details, this tool executes different operations:

- Determine haplotypes by local assembly of the genomic regions in which the samples being analyzed show substantial evidence of variation relative to the reference;
- Evaluate the evidence for haplotypes and variant alleles;
- Assigning per-sample genotypes.

> ### {% icon hands_on %} Hands-on: Mosaic Variant calling.
> 
> You may use files provided as examples with this tutorial and called
>    `Panel_alignment.bam` and `Panel_Target_regions.bed`. 
> 
> Run **Mutect2** {% icon tool %} restricting the search space on target regions with "-L" option to reduce computational burden.
>    The first step is needed to create an internal database of controls (i.e. **Panel Of Normals** - PoN) to reduce bias for somatic calls. It runs on a single sample at time:
>
>  - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -L Panel_target_regions.bed -I Panel_alignment_normal1.bam -O normal_genotyped1.vcf`
>  - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -L Panel_target_regions.bed -I Panel_alignment_normal2.bam -O normal_genotyped2.vcf`
>
>  Then use GATK's *CreateSomaticPanelOfNormals* tool to generate the PoN:
>
>  - `gatk GenomicsDBImport -L Panel_target_regions.bed -R HSapiensReference_genome_hg19.fasta --genomicsdb-workspace-path PoN_db -V normal_genotyped1.vcf -V normal_genotyped2.vcf`
>  - `gatk CreateSomaticPanelOfNormals -R HSapiensReference_genome_hg19.fasta -V gendb://PoN_db -O panel_of_normals.vcf`
>
>    > ### {% icon comment %} Note
>    > The --genomicsdb-workspace-path must point to a non-existent or empty directory.
>    {: .comment}
>
> 
> Then, to effectively call somatic mutations, we can use variants contained in the **PoN** and/or other public repositories  (e.g. by means of the option *--germline-resource*, using a VCF file containing frequencies of germline variants in the general population) to exclude germline variation. Finally, to properly classify somatic variants, we apply *FilterMutectCalls filtering*, which produces the final subset annotated VCF file. To this aim, we can run the following commands:
>
> - `gatk Mutect2 -R HSapiensReference_genome_hg19.fasta -I Panel_alignment.bam --germline-resource af-only-gnomad.vcf --panel-of-normals panel_of_normals.vcf -O somatic_genotyped_unfiltered.vcf`
> 
> - `gatk FilterMutectCalls -R HSapiensReference_genome_hg19.fasta -V somatic_genotyped_unfiltered.vcf -O somatic_genotyped_filtered.vcf`
>
 

# Variant annotation

Once called, variants (SNPs and InDels) need to be annotated.

We want to know for example if a variant is located in a gene, if it’s in the coding portion of that gene, if it causes an aminoacid substitution, if that substitution is deleterious for the encoded protein function.

**Variant annotation** is the process of attaching technical and/or biological information from multiple available source (e.g. **public databases**) to variants

---
## Gene model

The choice of gene model is essential for variant downstream variant annotation: it describe genomic positions of genes and each exon-intron exact locations

Different gene models can give different annotations:
![BRCA1]({{site.baseurl}}/images/brca1_var.jpg)

**Figure 1.** Variant indicated by the red dashed line can be annotated as *intronic* or *exonic* (on one of the UCSC transcript variants), depending on the adopted gene model:

- [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
- [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)
- [Gencode](https://www.gencodegenes.org/human/)

---
## Sequence variant nomenclature

Variant nomenclature should be described univocally:

- [HGVS](https://varnomen.hgvs.org/)
- [HGMC](https://www.genenames.org/)

## Variant class

- [Sequence Ontology](http://www.sequenceontology.org/)
- ![Sequence Ontology]({{site.baseurl}}/images/seqOnt.png)

## Population sequencing db

- [gnomAD](https://gnomad.broadinstitute.org/)
- [ExAC](http://exac.broadinstitute.org/)
- [1000 Genomes Project](https://www.internationalgenome.org/)
- [NHLBI-ESP 6500 exomes](https://evs.gs.washington.edu/EVS/)
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)

## Variant-disease/gene-disease db

- [Human Gene Mutation Database](http://www.hgmd.cf.ac.uk/ac/index.php)
- [clinically relevant variants known in a gene](https://www.ncbi.nlm.nih.gov/clinvar/)
- [Simple ClinVar](http://simple-clinvar.broadinstitute.org/)
- [Online Mendelian Inheritance in Man](https://www.omim.org/)

- [**dbNSFP**](https://sites.google.com/site/jpopgen/dbNSFP), a big database of curated annotations and precomputed functional predictions for all potential non-synonymous and splice-site single-nucleotide variants in the human genome

## Annotation Software and tools

- [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
- [Single Nucleotide Polymorphism Effect](http://snpeff.sourceforge.net/)
- [KGGSeq](http://grass.cgs.hku.hk/limx/kggseq/)

- Web Interface:
   - [wAnnovar](http://wannovar.wglab.org)
   - [VEP](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP)

---
## Annotation and filtering with wANNOVAR

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

# Variant prioritization

Once annotated, variants need to be filtered and prioritized

 - WES: Tens of thousands
 - WGS: Millions

**No universal filters, they depend on the experimental features**

## Variant impact

First of all you usually want to filter variants by consequence on the encoded protein, keeping those which have an higher impact on protein:
 
- Missense
- Nonsense
- Splice sites
- Frameshift indels
- Inframe indels

## Variant frequency
 
- Common variants are unlikely associated with a clinical condition
- A rare variant will probably have a higher functional effect on the protein
- Frequency cut-off have to be customized on each different case
- Typical cut-offs: 1% - 0.1%
- Allele frequencies may differ a lot between different populations

## Variant effect prediction Tools

- Tools that predict consequences of amino acid substitutions on protein function
- They give a score and/or a prediction in terms of "Tolerated", "Deleterious" (SIFT) or "Probably Damaging", "Possibly Damaging", "Benign" (Polyphen2)
  
  - fitCons
  - GERP++
  - SIFT
  - PolyPhen2
  - CADD
  - DANN
  - Condel
  - fathmm
  - MutationTaster
  - MutationAssessor
  - REVEL

## ACMG/AMP 2015 guidelines

The *American College of Medical Genetics* and the *Association for Molecular Pathology* published guidelines for the interpretation of sequence variants in May of 2015 [(Richards S. et al, 2015)](https://www.nature.com/articles/gim201530). This report describes updated standards and guidelines for classifying sequence variants by using criteria informed by expert opinion and experience

- 28 evaluation criteria for the clinical interpretation of variants. Criteria falls into 3 sets:
  - pathogenic/likely pathogenic (P/LP)
  - benign/likely benign (B/LB)
  - variant of unknown significance (VUS)

- [**Intervar**](http://wintervar.wglab.org): software for automatical interpretation of the 28 criteria. Two major steps:
 - automatical interpretation by 28 criteria
 - manual adjustment to re-interpret the clinical significance

## Prioritization

Phenotype-based prioritization tools are methods working by comparing the phenotypes of a patient with gene-phenotype known associations.

- [Phenolyzer](http://phenolyzer.wglab.org/): Phenotype Based Gene Analyzer, a tool to prioritize genes based on user-specific disease/phenotype terms


> ### {% icon hands_on %} Hands-on: Variant prioritization
>
>    Using knowledge gained on [Genomic databases and variant annotation]({{site.url}}{{site.baseurl}}/lectures/04.db.html) section, try to ***annotate***, ***filter*** and ***prioritize*** an example exome variant data, using two disease terms
>
> - Download file Sample1.all_exons.hg19.vcf from the *dataset*, as learned in the [Home]({{site.url}}{{site.baseurl}}/lectures/01.home.html) section
>
> - Go to [wANNOVAR](http://wannovar.wglab.org/)
>
> - Use the vcf file as **input file** and *hearing loss* and *deafness autosomal recessive*, as disease terms to prioritize results
>
> - Choose *rare recessive Mendelian disease* as Disease Model in the ***Parameter Settings*** section
>
> - Provide an istitutional email address and submit the Job
>
> - ... *wait for results* ...
>
>    In the results page you can navigate and download results.
>    Click ***Show*** in the *Network Visualization* section to see
>    Phenolyzer prioritization results
{: .hands_on}

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


# Contributors
{:.no_toc}

 * [Tommaso Pippucci](https://www.aosp.bo.it/content/curriculum?E=154659) - Sant’Orsola-Malpighi University Hospital, Bologna, Italy
 * [Alessandro Bruselles](https://www) - Istituto Superiore di Sanità, Rome, Italy
 * [Andrea Ciolfi](http://www.ospedalebambinogesu.it) - Ospedale Pediatrico Bambino Gesù, IRCCS, Rome, Italy
 * [Gianmauro Cuccuru](https://gmauro.github.io) - Albert Ludwigs University, Freiburg, Germany
 * [Giuseppe Marangi](http://www) - Institute of Genomic Medicine, Fondazione Policlinico Universitario A. Gemelli IRCCS, Università Cattolica del Sacro Cuore, Roma, Italy
 * [Paolo Uva](http://www.crs4.it/peopledetails/183/paolo-uva) - Centro di Ricerca, Sviluppo e Studi Superiori in Sardegna (CRS4), Pula, Cagliari, Italy

