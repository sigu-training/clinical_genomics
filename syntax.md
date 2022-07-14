---
layout: tutorial_hands_on

title: Computational genomics platform

---


# Section title

## Section subtitle

### Section sub sub title (absent in left menu)
This is an introduction. Examples are extracted from
[Galaxy exome-seq tutorial](https://galaxyproject.github.io/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html).

> ### {% icon details %} Details (e.g. Algorithms for normalization)
> Details about normalization.
>
> Thus, normalization based on the flag has two consequences:
>
> - first
> - second
>
>
> > ### {% icon hands_on %} Hands-on within details...wow!!
> >
> > 1. **Filter SAM or BAM, output SAM or BAM** {% icon tool %}:
> >   - {% icon param-files %} *"SAM or BAM file to filter"*: use as input the output of
> >     the previous step
> >   - *"Filter on bitwise flag"*: `yes`
> >     - *"Only output alignments with all of these flag bits set"*:
> >       - {% icon param-check %} *"Read is mapped in a proper pair"*
> >     - *"Skip alignments with any of these flag bits set"*:
> >       - {% icon param-check %} *"The read is unmapped"*
> >
> {: .hands_on}
>
{: .details}

> ### {% icon comment %} This is a note/comment
> The *Italian Society of Human Genetics* (**SIGU**) was established on November 14,
> 1997, when the pre-existing Italian Association of Medical Genetics and the Italian
> Association of Medical Cytogenetics joined.
> SIGU is one of the 27 member societies of FEGS (Federation of European Genetic
> Societies).
{: .comment}

> ### {% icon warning %} This is a warning
> Please be aware that ...
{: .warning}


### This is a list

1. Item 1
1. Item 2
1. Item 3
   - Sub item
     - Sub sub item

> ### {% icon hands_on %} Hands-on: TITLE
>
>    Introduction to hands-on
>
> 1. First step
> 1. Second step
>
>    > ### {% icon comment %} This is a note
>    > All the files are based on `hg19` reference genome
>    {: .comment}
>
>    > ### {% icon warning %} This is a warning
>    > Please be aware that the columns to use for calculations may be different
>    > compared to the example here considered, depending on the amount of columns
>    > of your *bed* files.
>    {: .warning}
>
>    The following is a snippet (from snippets folder):
>    {% include snippets/create_new_history.md %}
>
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Based on the report, do you think preprocessing of the reads
>    >    (trimming and/or filtering) will be necessary before mapping?
>    > 1. Which is your candidate gene?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. Solution for Sequence quality is quite good overall. If anything you might
>    > >    consider trimming the 3' ends of reads (base qualities decline
>    > >    slightly towards the 3' ends) or to filter out the small fraction
>    > >    of reads with a mean base quality < 5.
>    > > 1. The gene was ...
>    > >
>    > {: .solution}
>    {: .question}
>
>
{: .hands_on}

