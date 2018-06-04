# sRNACharP

A  Pipeline for small RNA Characterization (sRNACharP) written in the [Nextflow DSL](http://nextflow.io).

The pipeline obtains for each sRNA provided:
* the free energy of its secondary structure using centroidFold,
* the distance to the closest promoter predicted by bprom, 
* the distance to the closest terminator predicted by transterm,
* the distance to the closest ORFs listed in the genome annotation.


## 1. Requisites

## 1.1 Nextflow

Install [[nextflow]](http://nextflow.io) with the following command:
```
curl -fsSL get.nextflow.io | bash
```
Nextflow can be installed on any POSIX system (UNIX-like system such as Linux, OS X) with Java 7 or 8.

## 1.2 Other software

sRNACharP  requires several pieces of software. There are two ways to install all of the software packages required: natively, or by using [Docker](https://www.docker.com/) (recommended).

## *Docker
This is certainly the easier (and more reproducible) method. In order to run the pipeline with Docker, you need to install [Docker](https://www.docker.com/) 18.03 (or higher). See the included [Dockerfile](Dockerfile) for the configuration details of the Docker image we have built. Note that this Dockerfile is only included for information only as we have already generated the Docker image. To pull the docker image:

```
docker pull penacastillolab/srnacharp@sha256:c2a07f176d7cfe8cea3530bd76da05b30b182cdfe4d4b878f7d90e81f2d6a5f3
```

## *Natively

The pipeline can also be used without Docker by installing the following software components on your system:

* [Bedtools](http://bedtools.readthedocs.io/en/latest/index.html) version 2.26 (or higher)
* [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) version 2.4 (or higher)
* [CentroidFold](https://github.com/satoken/centroid-rna-package) (requires Boost libraries)
* [TransTermHP](http://transterm.cbcb.umd.edu/index.php) version 2.09
* [Bprom](http://www.softberry.com/berry.phtml?topic=fdp.htm&no_menu=on)
* [R](https://www.r-project.org/) version 3.4 (or higher)

## 2. Pipeline usage

Launching the pipeline with the `--help ` parameter shows the help message:

```
nextflow run sRNACharP.nf  --help
```

```
N E X T F L O W  ~  version 0.27.0
Launching `sRNACharP.nf` [adoring_volta] - revision: c85fc3d135

sRNACharP: sRNA Characterization Pipeline
----------------------------------------------------

Usage: 
    sRNACharP [options]

Options:
--org           ORGANISM         the organism name [required]
--dir           DIRECTORY        the directory containing the input files [required]
--genome        GENOME_FILE      the genome file in FASTA format [required]
--annot         ANNOTATION_FILE  the annotation file in BED format (including only protein coding genes)[required]
--sRNAs         SEQUENCE_FILE    the sRNAs file in BED format [required]
--transtermFile TRANSTERM_FILE   the TransTerm predictions file [optional].

```

## 3. Running the pipeline
sRNACharP is configured to run using the [Docker](https://www.docker.com/) container engine by default (see [Nextflow config file](nextflow.config)). If you have installed the required software natively, you also need to modify the [Nextflow config file](nextflow.config).

A [sample data set](test_data) is provided. To run sRNACharP with the test data, replace ADD_PATH below with the correct path in your system and run the following command on your terminal:

```
nextflow run sRNACharP.nf  --org="R_capsulatus" --dir="/ADD_PATH/sRNACharP/test_data" --genome="R_capsulatus_genome.fasta"  --sRNAs="R_capsulatus_sRNAs.bed" --genomeAnnotation="R_capsulatus_genomeAnnotation_proteincoding.bed"

```

On your terminal you should see something like this:

```
N E X T F L O W  ~  version 0.27.0
Launching `sRNACharP.nf` [crazy_heyrovsky] - revision: a560734176
[warm up] executor > local
[36/f67e03] Submitted process > createCRDFile
[80/28804e] Submitted process > getFASTAsRNAs
[72/6591a1] Submitted process > getUpstreamSequences (1)
[84/f7b852] Submitted process > getFreeEnergySS
[cb/bd7b45] Submitted process > runTransterm
[5d/6d813e] Submitted process > getPromoterSites (1)
[0b/eafc90] Submitted process > getPromoterSites (2)
[be/757ee4] Submitted process > getPromoterSites (3)
[b9/db95c2] Submitted process > getPromoterSites (4)
[b8/ad8cab] Submitted process > getPromoterSites (5)
[e6/5f8f87] Submitted process > parseTranstermResults (1)
[b5/53ee8a] Submitted process > getDistances (1)
[bb/b6a9ab] Submitted process > createAttributeTable (1)

    Pipeline execution summary
    ---------------------------
    Completed at: Fri May 04 11:29:32 CEST 2018
    Duration    : 17.3s
    Success     : true
    workDir     : /home/lpena/sRNACharP/work
    exit status : 0
```

## 4. Pipeline results

Analyses results are saved into the working directory.

Output files are the following (replace ORGANISM with the value of the --org parameter):

* `ORGANISM_FeatureTable.tsv` - table containing the characteristics per sRNA. Columns are: sequence identifier (sRNA ID), free energy of predicted secondary structure, distance to the -10 position of the closest predicted promoter, distance to the closest predicter terminator, distance to the
closest upstream ORF, a flag indicating whether the sRNA is on the same strand as the closest upstream ORF, distance to the closest downstream ORF, and a flag indicating whether the sRNA is on the same strand as the closest downstream ORF.
Here are some lines of a feature table generated:

```
STnc1010      -23.3   -3      0       -5      1       16      1
Stnc1020      -21.6   -15     0       -6      1       168     1
STnc470       -16.3   -55     1000    -103    0       34      1
```

* `ORGANISMRNAsBpromRes.txt`- compilation of Bprom results. Columns are: sequence identifier (sRNA ID), TSS Position, LDF, -10 Position, -10 Sequence, -10 Score, -35 Position, -35 Sequence, -35 Score.
Here are some lines of a predicted promoters table generated:

```
STnc1010        162     5.16    147     TGGAATAAT       58      130     CTGTCA  18
STnc1160        160     2.76    145     TTGAATTAT       37      124     TTGTCG  47
STnc1390        150     3.93    135     TGCCATAAT       62      115     TTGATT  53
```

* `ORGANISMRNAsTranstermRes.gtf` - compilation of TransTermHP results. Columns are: Replicon, Program, Terminator ID, Start, End, Score, Strand, and an empty field. 
Here are some lines of a predicted terminators table generated:

```
NC_016810.1     TransTermHP     TERM2   13530   13572   83      +       .
NC_016810.1     TransTermHP     TERM3   13538   13564   100     +       .
NC_016810.1     TransTermHP     TERM4   14760   14779   100     +       .

```

* Two more files are generated in the working directory: `ORGANISMRNAsSS.txt` and `ORGANISMGenomelength.txt` containing the free energy of the predicted secondary structure and the number of nucleotides in the genome, respectively. 

## 5. Citing
If you use the pipeline, please cite:

Eppenhof EJ, Peña-Castillo L. (2018) [Prioritizing bona fide bacterial small RNAs with machine learning classifiers. PeerJ Preprints 6:e26974v1](https://doi.org/10.7287/peerj.preprints.26974v1)

