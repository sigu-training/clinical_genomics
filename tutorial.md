---
layout: tutorial_hands_on

title: "NGS analysis using Galaxy"
zenodo_link: "https://doi.org/10.5281/zenodo.3054169"
questions:
  - "How do you explore NGS data without command line tools?"
  - "How do you identify genetic variants in samples based on exome sequencing
    data?"
  - "How do you, among the set of detected variants, identify candidate
    causative variants for a given phenotype/disease?"
objectives:
  - "Learn how to use Galaxy: load datasets, start an analysis"
  - "Use variant annotation and inheritance pattern to identify candidate
     causative variants and to prioritize them"
time_estimation: "2h"
key_points:
  - "Galaxy is an on line free platform for NGS analysis" 
  - "Exome sequencing is an efficient way to identify disease-relevant genetic
    variants."
  - "Variant annotation and being able to exploit genotype information across
    family members is key to identifying candidate disease variants"
  - "SnpEff and GEMINI are powerful tools for variant annotation and filtering"
contributors:
  - puva
  - gmauro
---


# Introduction
{:.no_toc}

Exome sequencing is a method that enables the selective sequencing of the
exonic regions of a genome - that is the transcribed parts of the genome present
in mature mRNA, including protein-coding sequences, but also untranslated
regions (UTRs).

In humans, there are about 180,000 exons with a combined length of ~ 30 million
base pairs (30 Mb). Thus, the exome represents only 1% of the human genome, but
has been estimated to harbor up to 85% of all disease-causing variants ([Choi
et al., 2009](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2768590/)).

Exome sequencing, thus, offers an affordable alternative to whole-genome
sequencing in the diagnosis of genetic disease, while still covering far more
potential disease-causing variant sites than genotyping arrays. This is of
special relevance in the case of rare genetic diseases, for which the causative
variants may occur at too low a frequency in the human population to be
included on genotyping arrays.

Of note, a recent study focusing on the area of clinical pediatric neurology
indicates that the costs of exome sequencing may actually not be higher even
today than the costs of conventional genetic testing ([Vissers et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5589982/)).

> ### {% icon details %} Exome sequencing *vs* whole-genome sequencing
>
> In principle, the steps illustrated in this tutorial are suitable also for
> the analysis of whole-genome sequencing (WGS) data. At comparable mean
> coverage, however, WGS datasets will be much larger than exome sequencing
> ones and their analysis will take correspondingly more time.
>
> The obvious benefit of WGS compared to exome-sequencing, of course, is that
> it will allow variant detection in even more regions of the genome.
> As a less apparent advantage, the more complete information of WGS data can
> make it easier to detect *copy number variation (CNV)* and
> *structural variants* such as translocations and inversions (although such
> detection will require more sophisticated analysis steps, which are not
> covered by this tutorial).
>
> Very generally, one could argue that exome-sequencing captures most of the
> information that can be analyzed with standard bioinformatical tools today at
> reasonable costs. WGS, on the other hand, captures as much information as
> today's sequencing technology can provide, and it may be possible to
> reanalyze such data with more powerful bioinformatical software in the
> future to exploit aspects of the information that were not amenable to
> analysis at the time of data acquisition.
{: .details}

The identification of causative variants underlying any particular genetic
disease is, as we will see in this tutorial, not just dependent on the
successful detection of variants in the genome of the patient, but also on
variant comparison between the patient and selected relatives. Most often
family trio data, consisting of the genome sequences of the patient and their
parents, is used for this purpose. With multisample data like this it becomes
possible to search for variants following any kind of Mendelian inheritance
scheme compatible with the observed inheritance pattern of the disease,
e.g. recessive, *de-novo*.

> ### {% icon details %} Related tutorials
>
> This tutorial focuses on the practical aspects of analyzing real-world
> patient data, and is a reduced version of the tutorial
> [Exome sequencing data analysis for diagnosing a genetic disease](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html). 
> If you are more interested in the theoretical aspects of
> variant calling, you may want to have a look at the related tutorial on
> [Calling variants in diploid systems](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/dip/tutorial.html).
>
> The tutorial on [Somatic variant calling](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/somatic-variants/tutorial.html)
> follows an analysis workflow that is rather similar to the one here, but
> tries to identify tumor variants by comparing a tumor sample to healthy
> tissue from the same patient. Together, the two tutorials are intended to get
> you started with genomics medicine using Galaxy.
>
{: .details}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Log in and data preparation

## Log in into Galaxy

Several Galaxy servers are available worldwide, with differ in the set of tools available. 
In this tutorial we'll use the European Galaxy server [usegalaxy.eu](https://usegalaxy.eu).


> ### {% icon hands_on %} Hands-on: Log in or register
> 1. Open your favorite browser (Chrome/Chromium, Safari or Firefox, but not Internet Explorer/Edge!)
> 2. Browse to [usegalaxy.eu](https://usegalaxy.eu)
> 3. Choose *Login or Register* from the navigation bar at the top of the page
> 4. If you have previously registered an account with this particular instance of Galaxy (user accounts are *not* shared between public servers!), proceed by logging in with your registered *public name*, or email address, and your password.  
If you need to create a new account, click on *Register here* instead.
> 5. After you logged in, click on the following URL: [https://usegalaxy.eu/join-training/obicourse](https://usegalaxy.eu/join-training/obicourse)  
For the training time, your user will be added to a dedicated training group and put into a private queue, which should be a bit faster than the regular queue. 
>>
{: .hands_on}


## Data Preparation

In this tutorial, we are going to analyze exome sequencing data from a family
trio, in which the boy child is affected by the disease
[osteopetrosis](https://ghr.nlm.nih.gov/condition/osteopetrosis/), while both
parents, who happen to be consanguineous, are unaffected. Our goal is to
identify the genetic variation that is responsible for the disease.

At the end of this tutorial you'll learn how to:

- perform a qualty control of the original sequenced reads in `fastq`
  format
- work with the final files of a typical exome analysis workflow: `bam` (mapped reads),
  and `vcf` (list of variants)

Therefore we assume that the complete analysis, from raw `fastq` to `bam` and `vcf`
has been already performed. This shortened analysis may suit your time frame
better, and is probably closer to how you would analyze your own
data.

The following hands-on section will guide you through obtaining the right
data for the tutorial.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a meaningful name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Obtain the raw `fastq` sequencing data
>
>
>    The whole set of raw sequence data for the family trio is available at
>    [Zenodo](https://zenodo.org/record/3054169).
>    To reduce the time of the analysis, we'll only use the raw sequence reads
>    of the proband
>    ```
>    https://zenodo.org/record/3243160/files/proband_R1.fq.gz
>    https://zenodo.org/record/3243160/files/proband_R2.fq.gz
>    ```
>
>    {% include snippets/import_via_link.md format="fastqsanger.gz" %}
>
>
> 3. Check that the newly created datasets in your history have their
>    datatypes assigned correctly to `fastqsanger.gz`, and fix any missing or
>    wrong datatype assignment
>
>    {% include snippets/change_datatype.md datatype="fastqsanger.gz" %}
>
> 4. Obtain the mapped sequencing data in `bam` format 
>
>    As in the previous step, the whole dataset is avalable from
>    [Zenodo](https://zenodo.org/record/3054169). Here we will only use the proband's data:
>    ```
>    https://zenodo.org/record/3243160/files/mapped_reads_proband.bam
>    ```
>
>    {% include snippets/import_via_link.md format="bam" %}
>
>
> 5. Check that the newly created dataset in your history has its 
>    datatype assigned correctly to `bam`, and fix any missing or wrong
>    datatype assignment
>
>    {% include snippets/change_datatype.md datatype="bam" %}
>
> 6. Obtain the variants in vcf format
>
>    Now you'll get the list of variants for the whole trio - father,
>    mother and proband:
>    ```
>    https://zenodo.org/record/4030079/files/Post-processed_FreeBayes_calls.vcf
>    ```
>
>    {% include snippets/import_via_link.md format="vcf" %}
>
> 7. Check that the newly created datasets in your history has its 
>    datatype assigned correctly to `vcf`, and fix any missing or wrong
>    datatype assignment
>
>    {% include snippets/change_datatype.md datatype="vcf" %}
>
> 8. Specify the genome version that was used for mapping and variant calling
>
>    Change the database/build (dbkey) for your `bam` and `vcf` datasets
>    to `hg19`.
>
>    {% include snippets/change_dbkey.md dbkey="Human Feb. 2009 (GRCh37/hg19) (hg19)" %}
>
>    > ### {% icon details %} Why specify genome versions
>    > When you are starting with sequencing data that has already been mapped
>    > to a particular genome version (human hg19 in this case), it is good
>    > practice to attach this information as metadata to the datasets.
>    >
>    > Doing so helps prevent accidental use of a different version of the
>    > reference genome in later steps like variant calling, which would
>    > typically lead to nonsensical results because of base position changes
>    > between the different versions.
>    {: .details}
>
>
> 9. Rename the datasets
>
>    For datasets that you upload via a link, Galaxy will pick the link
>    address as the dataset name, which you will likely want to shorten to
>    just the file names.
>
>    {% include snippets/rename_dataset.md %}
>
{: .hands_on}

> ### {% icon question %} Break
> {% icon trophy %} **Congratulations!** You successfully loaded all the datasets we'll use in this tutorial.
>
> Please stop here and move to the main room for discussion.
{: .question}

# Quality control

This step serves the purpose of identifying possible issues with the raw
sequenced reads input data before embarking on any "real" analysis steps.

Some of the typical problems with NGS data can be mitigated by preprocessing
affected sequencing reads before trying to map them to the reference genome.
Detecting some other, more severe problems early on may at least save you a lot
of time spent on analyzing low-quality data that is not worth the effort.

Here, we will perform a standard quality check on our input data and only point
out a few interesting aspects about that data. For a more thorough explanation
of NGS data quality control, you may want to have a look at the dedicated
tutorial on [Quality control](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html).

> ### {% icon hands_on %} Hands-on: Quality control of the input datasets
> 1. Run **FastQC** {% icon tool %} on each of your two `fastq` datasets
>       - {% icon param-files %} *"Short read data from your current history"*: all 2 FASTQ datasets selected with **Multiple datasets**
>
>    {% include snippets/select_multiple_datasets.md %}
>
>    When you start this job, four new datasets (one with the calculated raw
>    data, another one with an html report of the findings for each input
>    dataset) will get added to your history.
>
> 2. Use **MultiQC** {% icon tool %} to aggregate the raw **FastQC** data of all input datasets into one report
>      - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-files %} *"FastQC output"*: all two *RawData*
>             outputs of **FastQC** {% icon tool %})
>
> 3. Inspect the *Webpage* output produced by the tool
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Based on the report, do you think preprocessing of the reads
>    >    (trimming and/or filtering) will be necessary before mapping?
>    > 2. Why do all samples show a non-normal GC content distribution, and
>    >    should you be worried?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. Sequence quality is quite good overall. If anything you might
>    > >    consider trimming the 3' ends of reads (base qualities decline
>    > >    slightly towards the 3' ends) or to filter out the small fraction
>    > >    of reads with a mean base quality < 5.
>    > >    Feel free to run, *e.g.*, **Trimmomatic** {% icon tool %} on the
>    > >    fastq datasets if you want to, but don't expect this to have a big
>    > >    effect on the analysis given the high overall quality of the data
>    > >    of all samples.
>    > > 2. In whole-genome sequencing, a non-normal distribution of the GC
>    > >    content of the reads from a sample is typically considered to hint
>    > >    at possible contamination.
>    > >    Here, however, we are dealing with sequencing data from captured
>    > >    exomes, *i.e*, the reads are not representing random sequences from
>    > >    a genome, but rather a biased selection.
>    > >
>    > >    A bimodal GC content distribution, like for the samples at hand, is
>    > >    a characteristic feature of many exome capture methods and has also
>    > >    been observed with an Illumina Nextera Rapid Capture exome kit
>    > >    before (compare Fig. 3D in
>    > >    [Shigemizu et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4522667/)
>    > >    ).
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}


> ### {% icon question %} Break
> {% icon trophy %} **Congratulations!** You successfully performed the quality 
> control of your raw sequence reads.
>
> Please stop here and move to the main room for discussion.
{: .question}


# Variant annotation

A list of variants detected in a set of samples is a start, but to discover
biologically or clinically relevant information in it is almost impossible
without some additional tools and data. In particular, the
variants in the list need to be:

- **annotated** with respect to their potential relevance for the biological
  / clinical phenotype that is studied

  Even with exome sequencing, only a fraction of the detected variants will
  have a clear impact on the function of a protein (many variants will
  introduce silent mutations, or reside in intronic regions still covered by
  the exome-enriched sequencing data). Of these, many will have been observed
  before in healthy individuals arguing against them playing an important role
  in an adverse phenotype.

- **filtered** based on the inheritance pattern expected for a causative
  variant

  A multisample VCF file records the most likely genotypes of all samples at
  every variant site. Knowing which individuals (samples) are affected by a
  phenotype we can exclude variants with inheritance patterns that are
  incompatible with the observed inheritance of the phenotype.

## Get data

Our workhorse for annotating and reporting variants and the genes affected by
them will be the **GEMINI** framework. GEMINI comes bundled with a wealth of
annotation data for human variants from many different sources. These can be
used to annotate any list of human variants conveniently, without the need for
separate downloads and conversion between different annotation data formats.

> ### {% icon details %} Custom annotations and GEMINI
> Beyond its bundled annotation data, GEMINI also provides (limited) support
> for using custom annotations.
> The [Somatic variant calling tutorial](../somatic-variants/tutorial.html)
> demonstrates the use of **GEMINI annotate** {% icon tool %} for this purpose.
{: .details}

The only additional annotation tool we need, for the purpose of this
tutorial, is the tool **SnpEff**, which can annotate variants with their
calculated effects on known genomic features. Because SnpEff is a generic tool
that can be used on variants found in the genome of any organism we need to
provide it with a so-called **SnpEff genome file** that holds the annotated
features (genes, transcripts, translated regions, *etc.*) for our genome of
interest.

While annotated variants are all we need to *prioritize* them as described
above, *filtering based on inheritance patterns* requires a way to inform
GEMINI about the relationship between our samples and their observed
phenotypes. This is done through a so-called **pedigree file** in PED format,
which is rather simple to generate manually.

> ### {% icon hands_on %} Hands-on: Obtain GEMINI pedigree files
> 1. Create a PED-formatted pedigree dataset describing our single-family sample trio:
>
>    ```
>    #family_id    name     paternal_id    maternal_id    sex    phenotype
>    FAM           father   0              0              1      1
>    FAM           mother   0              0              2      1
>    FAM           proband  father         mother         1      2
>    ```
>
>    and set its datatype to `tabular`.
>
>    {% include snippets/create_new_file.md format="tabular" %}
>
>    > ### {% icon warning %} Remember those sample names
>    >
>    > Names in the above pedigree file should be exactly `father`, `mother`,
>    > `proband` as those names were used as sample names at the read mapping step,
>    > and have been propagated to the VCF dataset of variants.
>    > It is **important** that
>    > you use matching sample names in the pedigree and in the VCF dataset, or
>    > GEMINI will not be able to connect the information in them.
>    >
>    {: .warning}
>
>    > ### {% icon details %} More on PED files
>    >
>    > The PED format is explained in the help section of **GEMINI load**
>    > {% icon tool %}.
>    >
>    > Take a moment and try to understand the information that is encoded in
>    > the PED dataset we are using here.
>    {: .details}
{: .hands_on}

## Variant annotation with functional genomic effects

We need to start annotating our variants with SnpEff simply because Gemini
knows how to parse SnpEff-annotated VCFs, while GEMINI output cannot be used
with SnpEff.

> ### {% icon hands_on %} Hands-on: Adding annotations with SnpEff
> 1. **SnpEff eff** {% icon tool %} (not the one for the SARS-CoV-2 pipeline)
>    - {% icon param-file %} *"Sequence changes (SNPs, MNPs, InDels)"*: the
>      uploaded **VCF** file
>    - *"Input format"*: `VCF`
>    - *"Output format"*: `VCF (only if input is VCF)`
>    - *"Genome source"*: `Locally installed snpEff database`
>       - *"Genome"*: `Homo sapiens: hg19` (or a similarly named option)
>
>    - *"Produce Summary Stats"*: `Yes`
>
{: .hands_on}

Running the above job will produce two datasets. One is a *Summary Stats* HTML
report, which contains some interesting general metrics such as a distribution
of variants across gene features. The other one is the main annotation result -
a VCF like the input, but with annotations of variant effects added to the INFO
column.

> ### {% icon hands_on %} Optional hands-on: Inspect the Summary Stats output produced by SnpEff
>
> 1. Display the dataset:
>    - Click the {% icon galaxy-eye %} icon next to the HTML dataset generated
>      by SnpEff to display its contents.
>
>    > ### {% icon question %} Question
>    >
>    > One section in the report is **Number of effects by type and region**.
>    > Given that you are analyzing exome data, what is the most surprising
>    > aspect in this section? Do you have an idea how to explain it?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > According to the report, intronic variants make up 50% of all
>    > > variants detected!
>    > >
>    > > Exome capture kits are designed to capture exons plus a bit of
>    > > surrounding sequence to ensure proper coverage of the exon ends
>    > > including splice junction sites.
>    > >
>    > > Thus, even though intronic sequences are underrepresented in exome
>    > > sequencing data, not all of them are eliminated. Since mutations
>    > > in most intron bases are neutral, they can accumulate at higher
>    > > frequency than most mutations in exons and, thus, still represent a
>    > > relevant fraction of all detected variants.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

> ### {% icon question %} Break
> {% icon trophy %} **Congratulations!** You successfully annotated your variants.
>
> Please stop here and move to the main room for discussion.
{: .question}


# Variant filtering

## Generating a GEMINI database of variants for further annotation and efficient variant queries

Next, we are going to use the SnpEff-annotated VCF as the basis for more
exhaustive annotation with GEMINI. Unlike SnpEff, GEMINI does not just add
annotations to a list of variants in VCF format. Instead the framework
extracts the variants from the VCF input and stores them, together with newly
added annotations, in an SQL database. It then lets you formulate queries for
retrieving and reporting subsets of variants. The combined variant
extraction/annotation/storage step is performed by the **GEMINI load** tool. In
addition, that same tool can be used to incorporate sample pedigree info into
the database.

> ### {% icon hands_on %} Hands-on: Creating a GEMINI database from a variants dataset
>
> 1. **GEMINI load** {% icon tool %} with
>    - {% icon param-file %} *"VCF dataset to be loaded in the GEMINI database"*:
>      the output of **SnpEff eff** {% icon tool %}
>    - *"The variants in this input are"*: `annotated with snpEff`
>    - *"This input comes with genotype calls for its samples"*: `Yes`
>
>      Sample genotypes were called by Freebayes for us.
>
>    - *"Choose a gemini annotation source"*: select the latest available annotations snapshot (most likely, there will be only one)
>    - *"Sample and family information in PED format"*: the pedigree file
>      prepared above
>    - *"Load the following optional content into the database"*
>      - {% icon param-check %} *"GERP scores"*
>      - {% icon param-check %} *"CADD scores"*
>      - {% icon param-check %} *"Gene tables"*
>      - {% icon param-check %} *"Sample genotypes"*
>      - {% icon param-check %} *"variant INFO field"*
>
>      Leave **unchecked** the following:
>      - *"Genotype likelihoods (sample PLs)"*
>
>         Variants have been generated by Freebayes which does not generate these values
>
>      - *"only variants that passed all filters"*
>
>        This setting is irrelevant for our input because Freebayes did not
>        apply any variant filters.
{: .hands_on}

Running this job generates a GEMINI-specific database dataset, which can only
be processed with other GEMINI tools. The benefit, however, is that we now have
variants, rich annotations and pedigree info stored in a format that enables
flexible and highly efficient queries, which will greatly simplify our actual
task to identify the variant responsible for the child's disease!

> ### {% icon details %} The GEMINI suite of tools
>
> For a beginner, the sheer number of GEMINI tools may be a bit daunting and
> give the impression that this framework adds a lot of complexity.
> In practice, however, you will likely only need a very limited number of
> these tools for any given analysis. To help you get an overview, here is a
> list of the most general-purpose tools and their function:
>
> - Every analysis will always start with **GEMINI load** {% icon tool %} to
>   build the GEMINI database from a variant list in VCF format and the
>   default annotations bundled with GEMINI.
> - **GEMINI amend** {% icon tool %} and **GEMINI annotate** {% icon tool %}
>   let you add additional information to a GEMINI database after it has been
>   created. Use **GEMINI amend** to overwrite the pedigree info in a database
>   with the one found in a given PED input, and **GEMINI annotate** to add
>   custom annotations based on a VCF or BED dataset.
> - **GEMINI database info** {% icon tool %} lets you peek into the content of
>   a GEMINI database to explore, *e.g*, what annotations are stored in it.
> - **GEMINI query** {% icon tool %} and **GEMINI gene-wise** {% icon tool %}
>   are the two generic and most flexible tools for querying, filtering and
>   reporting variants in a database. There is hardly anything you cannot do
>   with these tools, but you need to know or learn at least some SQL-like
>   syntax to use them.
> - **GEMINI inheritance pattern** {% icon tool %} simplifies variant queries
>   based on their inheritance pattern. While you could use the *query* and
>   *gene-wise* tools for the same purpose, this tool allows you to perform a
>   range of standard queries without typing any SQL. **This tool is what we
>   are going to use here to identify the causative variant.**
> - The remaining tools serve more specialized purposes, which are beyond the
>   scope of this tutorial.
>
> The [Somatic variant calling tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/somatic-variants/tutorial.html)
> demonstrates the use of the **GEMINI annotate** and **GEMINI query** tools,
> and tries to introduce some essential bits of GEMINI's SQL-like syntax.
>
> For a thorough explanation of all tools and functionality you should consult
> the [GEMINI documentation](https://gemini.readthedocs.io).
{: .details}

## Candidate variant detection

Let us now try to identify variants that have the potential to explain the
boy child's osteopetrosis phenotype. Remember that the parents are
consanguineous, but both of them do not suffer from the disease. This
information allows us to make some assumptions about the inheritance pattern of the causative variant.

> ### {% icon question %} Question
>
> Which inheritance patterns are in line with the phenotypic observations
> for the family trio?
>
> Hint: GEMINI easily lets you search for variants fitting any of the following
> inheritance patterns:
>    - Autosomal recessive
>    - Autosomal dominant
>    - X-linked recessive
>    - X-linked dominant
>    - Autosomal de-novo
>    - X-linked de-novo
>    - Compound heterozygous
>    - Loss of heterozygosity (LOH) events
>
> Think about which of these might apply to the causative variant.
>
> > ### {% icon solution %} Solution
> > - Since both parents are unaffected the variant cannot be dominant and
> >   inherited.
> > - A de-novo acquisition of a dominant (or an X-linked recessive) mutation
> >   is, of course, possible.
> > - A recessive variant is a possibility, and a more likely one given
> >   the parents' consanguinity.
> > - For both the de-novo and the inherited recessive case, the variant could
> >   reside on an autosome or on the X chromosome. Given that we provided you
> >   with only the subset of sequencing reads mapping to chr8, an X-linked
> >   variant would not be exactly instructive though ;)
> > - A compound heterozygous combination of variant alleles affecting the
> >   same gene is possible, but less likely given the consanguinity of the
> >   parents (as this would require two deleterious variant alleles in the
> >   gene circulating in the same family).
> > - A loss of heterozygosity (LOH) turning a heterozygous recessive variant
> >   into a homozygous one could be caused by uniparental disomy or by an LOH
> >   event early in embryonic development, but both these possibilities have
> >   an exceedingly small probability.
> >
> > Based on these considerations it makes sense to start looking for
> > inherited autosomal recessive variants first. Then, if there is no
> > convincing candidate mutation among them, you could extend the search to
> > de-novo variants, compund heterozygous variant pairs and LOH events -
> > probably in that order.
> {: .solution}
{: .question}

Since our GEMINI database holds the variant and genotype calls for the
family trio and the relationship between the family members, we can make use
of **GEMINI inheritance pattern** {% icon tool %} to report all variants
fitting any specific inheritance model with ease.

Below is how you can perform the query for inherited autosomal recessive
variants. Feel free to run analogous queries for other types of variants that
you think could plausibly be causative for the child's disease.

> ### {% icon hands_on %} Hands-on: Finding and reporting plausible causative variants
> 1. **GEMINI inheritance pattern** {% icon tool %}
>    - *"GEMINI database"*: the GEMINI database of annotated variants; output
>      of **GEMINI load** {% icon tool %}
>    - *"Your assumption about the inheritance pattern of the phenotype of interest"*:
>      `Autosomal recessive`
>      - {% icon param-repeat %} *"Additional constraints on variants"*
>        - *"Additional constraints expressed in SQL syntax"*:
>          `impact_severity != 'LOW'`
>
>          This is a simple way to prioritize variants based on their
>          functional genomic impact. Variants with *low impact severity* would
>          be those with no obvious impact on protein function (*i.e.*, silent
>          mutations and variants outside coding regions)
>      - *"Include hits with less convincing inheritance patterns"*: `No`
>
>        This option is only meaningful with larger family trees to account
>        for errors in phenotype assessment.
>      - *"Report candidates shared by unaffected samples"*: `No`
>
>        This option is only meaningful with larger family trees to account
>        for alleles with partial phenotypic penetrance.
>    - *"Family-wise criteria for variant selection"*: keep default settings
>
>      This section is not useful when you have data from just one family.
>    - In *"Output - included information"*
>      - *"Set of columns to include in the variant report table"*:
>        `Custom (report user-specified columns)`
>        - *"Choose columns to include in the report"*:
>          - {% icon param-check %} *"alternative allele frequency (max_aaf_all)"*
>        - *"Additional columns (comma-separated)"*:
>          `chrom, start, ref, alt, impact, gene, clinvar_sig,
>          clinvar_disease_name, clinvar_gene_phenotype, rs_ids`
>
>    > ### {% icon comment %} Variant- *vs* gene-centric ClinVar information
>    > Most annotation data offered by GEMINI is variant-centric, *i.e.*,
>    > annotations reported for a variant are specific to this exact variant.
>    > A few annotation sources, however, also provide gene-centric
>    > information, which applies to the gene affected by a variant, not
>    > necessarily the variant itself.
>    >
>    > In the case of [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), the
>    > annotation fields/columns `clinvar_sig` and `clinvar_disease_name` refer
>    > to the particular variant, but `clinvar_gene_phenotype` provides
>    > information about the affected gene.
>    >
>    > Including the gene phenotype in the report can be crucial because a
>    > gene may be well known to be disease-relevant, while a particular
>    > variant may not have been clinically observed or been reported before.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Question
>
> From the GEMINI reports you generated, can you identify the most likely
> candidate variant responsible for the child's disease?
{: .question}

> ### {% icon details %} More GEMINI usage examples
>
> While only demonstrating command line use of GEMINI, the following tutorial
> slides may give you additional ideas for variant queries and filters:
>
> - [Introduction to GEMINI](https://s3.amazonaws.com/gemini-tutorials/Intro-To-Gemini.pdf)
> - [Identifying *de novo* mutations with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-DeNovo-Tutorial.pdf)
> - [Identifying recessive gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf)
> - [Identifying dominant gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Dominant-Tutorial.pdf)
{: .details}


## Displaying data in UCSC genome browser

A good way to proceed with candidate variants is to look at their coverage. 
This can be done by using genome browsers to display the aligned reads in the position of the candidate variant.
Aligned reads are stored in `bam` files, and Galaxy can display `bam` launching a genome browser such as IGV on 
your local machine, and it can connect to online genome browsers as well.
An example of such an online genome browser is the UCSC Genome Browser.

> ### {% icon hands_on %} Hands-on: UCSC genome browser
>
> 1. First, check that the **database** of your `bam` dataset is `hg19`. If not, click on the {% icon galaxy-pencil %} pencil icon and modify the **Database/Build:** field to `Human Feb. 2009 (GRCh37/hg19) (hg19)`.
>
>    {% include snippets/change_dbkey.md dbkey="hg19" %}
>>
> 2. To **visualize the data in UCSC genome browser**, click on `display at UCSC main` option visible when you expand the history item.
>
>    ![`display at UCSC main` link]({{site.baseurl}}/images/101_displayucsc.png)
>
>    This will upload the data to UCSC as custom track. To see your data look at the `User Track` near the top.
>    You can enter the coordinates of one of your variants at the top to jump to that location.
>
>    ![`User Track` shown in the UCSC genome browser]({{site.baseurl}}/images/101_21.png)
{: .hands_on}

UCSC provides a large number of tracks that can help you get a sense of your genomic area, it contains common SNPs, repeats, genes, and much more (scroll down to find all possible tracks).


# Conclusion
{:.no_toc}

It was not hard to find the most likely causative mutation for the child's
disease (you did find it, right?).

Even though it will not always provide as strong support for just one specific
causative variant, analysis of whole-exome sequencing data of family trios (or
other related samples) can often narrow down the search for the cause of a
genetic disease to just a very small, manageable set of candidate variants, the
relevance of which can then be addressed through standard methods.

> ### {% icon question %} Break
> {% icon trophy %} **Congratulations!** You successfully completed the tutorial.
>
> Please move to the main room for discussion.
{: .question}


# Credits
{:.no_toc}

This tutorial is a short version of the Galaxy tutorial [Exome sequencing data analysis for diagnosing a genetic disease](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html).

