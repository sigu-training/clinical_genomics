## On-line training tutorial on [“Data analysis and interpretation for clinical genomics” using Galaxy](https://sigu-training.github.io/clinical_genomics)

In years 2018-2019, on behalf of the Italian Society of Human Genetics ([SIGU](https://www.sigu.net/)) an itinerant [Galaxy](https://usegalaxy.eu/)-based “hands-on-computer” training activity entitled “Data analysis and interpretation for clinical genomics” was held four times on invitation from different Italian institutions (Università Cattolica del Sacro Cuore in Rome, University of Genova, SIGU 2018 annual scientific meeting in Catania, University of Bari) and was offered to about 30 participants each time among clinical doctors, biologists, laboratory technicians and bioinformaticians. Topics covered by the course were NGS data quality check, detection of variants, copy number alterations and runs of homozygosity, annotation and filtering and clinical interpretation of sequencing results.

Realizing the constant need for training on NGS analysis and interpretation of sequencing data in the clinical setting, we designed an on-line [Galaxy](https://usegalaxy.eu/)-based training resource articulated in presentations and practical assignments by which students will learn how to approach NGS data quality at the level of fastq, bam and VCF files and clinically-oriented examination of variants emerging from sequencing experiments.

### Updating GitHub Pages - for Admins only
Automatic update of GitHub pages through Jekyll is disabled. To update pages you have to build the website locally and push the updated pages:

 * Clone this repository with `git clone https://github.com/sigu-training/clinical_genomics.git`
 * (If not done yet) Set up the conda environment - see this [tutorial](https://galaxyproject.github.io/training-material/topics/contributing/tutorials/running-jekyll/tutorial.html) for details:
   * Move to the cloned repository with `cd clinical_genomics`
   * Install conda, if not already installed: `make install conda`
   * Create the **galaxy_training_material** conda environment: `make create-env`
   * Install Jekyll and related modules into the conda environment: `make install`
 * Build the website locally with `make build` (static) or `make serve` (interactive). This will generate the new HTML Pages in the *docs* folder.
 * Push the updated files. 

When updating text, please follow the [formatting instructions](https://sigu-training.github.io/clinical_genomics/syntax.html).

### Contributors

 * [Tommaso Pippucci](https://www.aosp.bo.it/content/curriculum?E=154659) - Sant’Orsola-Malpighi University Hospital, Bologna, Italy
 * [Alessandro Bruselles](https://www) - Istituto Superiore di Sanità, Rome, Italy
 * [Andrea Ciolfi](http://www.ospedalebambinogesu.it) - Ospedale Pediatrico Bambino Gesù, IRCCS, Rome, Italy
 * [Gianmauro Cuccuru](http://www) - Albert Ludwigs University, Freiburg, Germany
 * [Giuseppe Marangi](http://www) - Institute of Genomic Medicine, Fondazione Policlinico Universitario A. Gemelli IRCCS, Università Cattolica del Sacro Cuore, Roma, Italy
 * [Paolo Uva](https://www.researchgate.net/profile/Paolo_Uva) - IRCCS G. Gaslini, Genoa, Italy

### Italian Society of Human Genetics

The *Italian Society of Human Genetics* (**SIGU**) was established on November 14, 1997, when the pre-existing Italian Association of Medical Genetics and the Italian Association of Medical Cytogenetics joined. SIGU is one of the 27 member societies of FEGS (Federation of European Genetic Societies).
Animated by a predominant scientific spirit, SIGU wants to be reference for all health-care issues involving human genetics in all its applications.
Its specific missions are to develop quality criteria for medical genetic laboratories, to promote writing of guidelines in the field of human genetics and public awareness of the role and limitations of genetic diagnostic techniques. 
SIGU coordinates activities of several working groups: Clinical Genetics, Cytogenetics, Prenatal Diagnosis, Neurogenetics, Fingerprinting, Oncological Genetics, Immunogenetics, Genetic Counseling, Quality Control, Medical Genetics Services, Bioethics. More than 1000 medical geneticists and biologists are active members of the society.
