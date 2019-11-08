---
layout: page
title: Bioinformatic analysis in NGS
summary: "Basics of file formats and bioinformatic workflows for NGS data analysis"
sidebar_link: true
---

# Next Generation Sequencing


- *Next (or Second) Generation Sequencing* (NGS/SGS) is an umbrella-term covering a number of approaches to DNA sequencing that have been developed after the first, widespread and for long time most commonly used Sanger sequencing.
- *NGS* is also known as *Massive Parallel Sequencing* (MPS), a term that makes explicit the paradigm shared by all these technologies, that is to sequence in parallel a massive library of spatially separated and clonally amplified DNA templates. 
    - For a comprehensive review of the different *NGS* technologies see [Goodwin et al., 2016](https://www.nature.com/articles/nrg.2016.49), which also includes an introduction to the third generation methods allowing sequencing of long single-molecule reads.
    
## NGS in the clinic

In the span of less than a decade, NGS approaches have pervaded clinical laboratories revolutionizing genomic diagnostics and increasing yield and timeliness of genetic tests.

In the context of disorders with a recognized strong genetic contribution such as neurogenetic diseases, *NGS* has been firmly established as the strategy of choice to rapidly and efficiently diagnose diseases with a Mendelian basis. A general diagnostic workflow for these disorders currently embraces different *NGS*-based diagnostic options as illustrated in Figure 1.

[![ngs_neuro_diagnostics]({{site.url}}{{site.baseurl}}/images/ngs_neuro_diagnostic.jpg)]
**Figure 1**. General workflow for genetic diagnosis of neurological diseases. (\*If considering high-yield single-gene testing of more than 1–3 genes by another sequencing method, note that next-generation sequencing is often most cost-effective. †Genetic counselling is required before and after all genetic testing; other considerations include the potential for secondary findings in genomic testing, testing parents if inheritance is sporadic or recessive, and specialty referral. From [Rexach et al., 2019](https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(19)30033-X/fulltext))

Currently, the most common *NGS* strategies in clinical laboratories are:

- **gene panels** where the coding exons of a clinically-relevant group of genes is interrogated 
- **exome sequencing** where virtually all the protein-coding exons in a genome are simultaneously sequenced

Usually, the 


## CNV detection from targeted sequencing data


### Computational approaches


- There are four main methods for *CNV* identification from short-read +NGS data (see figure below):
  - **Read Count** (RC)
  - **Read Pair** (RP)
  - **Split Read** (SR)
  - **De Novo Assembly** (AS)
  
  
- *RP* and *SR* require continuous coverage of the *CNV* region or reads encompassing *CNV* breakpoints, as in whole genome sequencing. The sparse nature and small size of exonic targets hamper the application of these methods to targeted sequencing. 
- *RC* is the best suited method for *CNV* detection from whole exomes or gene panels where:
  - deletions appear as exonic targets devoid of reads
  - duplications appear as exonic targets characterized by excess of coverage


[![CNV detection methods]({{site.url}}{{site.baseurl}}/images/methods_identification_cnv.png)]
**Figure 1**. Methods for detection of CNVs in short read NGS data (adapted from [Tattini et al., 2015](https://doi.org/10.3389/fbioe.2015.00092))


### RC method and data normalization


In targeted sequencing, a method to study DNA copy number variation by *RC* (as implemented in *EXCAVATOR* tool, [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))) is to consider the **exon mean read count** (EMRC):

*EMRC = RCe/Le*

where *RCe* is the number of reads aligned to a target genomic region e and *Le* is the size of that same genomic region in base pairs ([Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))).
Three major bias sources are known to affect *EMRC* dramatically in targeted sequencing data:
- local **GC content** percentage
- genomic **mappability**
- target region **size**

These biases contribute to non uniform read depth across target regions and, together with the sparse target nature, challenge the applicability of *RC* methods to targeted data. As shown in Figure 2 (left panel), in single-sample data the *EMRC* distributions of genomic regions characterized by different copy numbers largely overlap, revealing poor *CNV* prediction capability. 

*EMRC ratio* between two samples can be used as a normalization procedure. The effect of *EMRC ratio*-based normalization is clear in Figure 2 (right panel) as a markedly improved correspondece between the predicted and the real copy number states.

[![Data Normalization_1]({{site.url}}{{site.baseurl}}/images/normalization_EMRC.png)]
**Figure 2**. Effect of *EMRC ratio* on DNA copy number prediction (adapted from [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))

Similarly FPKM, a normalized measure of read depth implemented in ExomeDepth tool ([Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526)), is affected by extensive exon–exon variability but comparison between pairs of exome datasets demonstrates high correlation level of the normalized read count data across samples making it possibile to use one or more combined exomes as reference set to base the CNV inference on ([Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526)).

[![Data Normalization_2]({{site.url}}{{site.baseurl}}/images/normalization_FPKM.png)]
**Figure 3**. Correlation of the normalized read count data between samples (from [Plagnol et al., 2012](https://doi.org/10.1093/bioinformatics/bts526))

### *CNV* detection accuracy

Strong correlation is observed between Affymetrix array-SNP and exome-derived data for *CNVs* >1 Mb (Figure 4, right panel), while including CNVs of any size dramatically decreases the correlation level (Figure 4, left panel). This can be explained by the different distribution of exons and SNP probes throughout the genome. Candidate *CNV* regions as identified by *EXCAVATOR* contain a comparable number of exons and SNP probes when they are > 1 Mb (*R*=0.8), while regions <100 Kb do not (*R*=-0.02). Accordingly, [Krumm et al., 2012](https://genome.cshlp.org/content/22/8/1525.full) report 94% precision in detecting *CNVs* across three or more exons.

[![Correlation_Exome_Affymetrix]({{site.url}}{{site.baseurl}}/images/corr_exome_affymetrix.png)]
**Figure 4**. Correlation between array-SNP and exome-derived *CNVs* including all (left panel) or > 1 Mb *CNVs* (adapted from [Magi et al., 2013](http://genomebiology.com/2013/14/10/R120))

To increase chances to detect *CNVs* encompassing few exons or in non-coding regions from exome data, [D'Aurizio et al., 2016](https://academic.oup.com/nar/article/44/20/e154/2607979) have proposed an extension to *EXCAVATOR* to exploit off-target read count (*EXCAVATOR2*). Similar approaches taking advantage of off-target reads have been described by [Bellos et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4147927/) or [Talevich et al., 2016](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873).

### Tools for *CNV* detection from gene panels or exome data

A number of tools or pipelines implementing modified versions of previously published tools have been reported to detect single-exon *CNVs* in clinical gene panels:
- [CoNVaDING](https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22969)
- [DeCON](https://wellcomeopenresearch.org/articles/1-20/v1)
- [ExomeDepth v1.1.6](https://www.nature.com/articles/ejhg201742)
- [Atlas-CNV](https://www.nature.com/articles/s41436-019-0475-4)
- [DeviCNV](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6192323/)

In addition to the above-mentioned methods, many tools have been developed that detect *CNVs* from exomes. Evaluation of these tools is not straightforward as there is lack of a gold standard for comparison. As a consequence, there is no consensus level on pipelines as high as for single nucleotide variants. [Zare et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5452530/pdf/12859_2017_Article_1705.pdf) have reviewed this topic with focus on cancer, reporting poor overlap of different tools and added challenges for somatic variant calling.

---