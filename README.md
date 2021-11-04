# sRNACharP

A  Pipeline for small RNA Characterization (sRNACharP) written in the [Nextflow DSL](http://nextflow.io).

The pipeline obtains for each sRNA provided:

* the free energy of its secondary structure using centroidFold,
* the distance to the closest promoter predicted by Promotech,
* the distance to the closest terminator predicted by transterm,
* the distance to the closest ORFs listed in the genome annotation, and
* a vector with the reverse tetranucleotide (tetranucleotide + its reverse complement) frequencies.

## 1. Requisites

## 1.1 Nextflow

Install [[nextflow]](http://nextflow.io) with the following command:

```
curl -fsSL get.nextflow.io | bash
```

Nextflow can be installed on any POSIX system (UNIX-like system such as Linux, OS X) with Java 7 or 8.

## 1.2 Other software

sRNACharP  requires several pieces of software. There are two ways to get sRNACharP to run: 1) install all of the software packages required natively, or 2) use sRNACharP [Docker](https://www.docker.com/) container (recommended).

## *Docker

This is certainly the easier (and more reproducible) method. In order to run the pipeline with Docker, you need to install [Docker](https://www.docker.com/) 18.03 (or higher). See the included [Dockerfile](Dockerfile) for the configuration details of the Docker image we have built. Note that this Dockerfile is included for information only as we have already generated the Docker image. To pull the docker image:

```
docker pull penacastillolab/srnacharpv2
```

## *Natively

The pipeline can also be used without Docker by installing the following software components on your system:

* [Bedtools](http://bedtools.readthedocs.io/en/latest/index.html) version 2.27 (or higher)
* [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) version 2.4.18 (or higher)
* [CentroidFold](https://github.com/satoken/centroid-rna-package) (requires Boost libraries)
* [TransTermHP](http://transterm.cbcb.umd.edu/index.php) version 2.09
* [R](https://www.r-project.org/) version 3.4 (or higher)

## 2. Pipeline usage

After downloading the file [sRNACharP.nf](sRNACharP.nf), one can launch the pipeline with the `--help ` parameter (make sure you run this command on the directory where sRNACharP.nf is located). This should show the help message:

```
nextflow run sRNACharP.nf  --help
```

```
N E X T F L O W  ~  version 20.10.0
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
--promoterFile  PREDICTED_PROMOTER_FILE the Promoter predictions file [required]

```

The organism name is any alphanumeric string without whitespaces.

The annotation file is a BED file with seven columns: Sequence, Start, End, Name, Score, Strand, and Type. For instance,

```
NC_008060.1	277	1393	Bcen_0001	·	-	gene
NC_008060.1	1407	2622	Bcen_0002	·	-	gene
NC_008060.1	2665	2992	Bcen_0003	·	-	gene
```

The sRNAs BED file requires six columns: Sequence, Start, End, Name, Score, Strand. A seventh column with the Type may also be provided. The name in the fourth column must uniquely identify each sRNA. For instance,

```
NC_008060.1	208288	208699	s208288	412	+
NC_008060.1	604653	605021	s604653	369	-
NC_008060.1	2044980	2045161	s2044980	182	+
NC_008062.1	409323	409390	s409323	68	+
```

NOTE: Make sure that the sequence identifier is the same in both BED files and in the genome FASTA file (in the example above the sequence identifier is NC_008060.1). The sequence identifier should not contain any whitespace.

The promoter file provides the location of promoter sequences in the genome. You can use our programs [Promotech](https://github.com/BioinformaticsLabAtMUN/PromoTech) or [Promotech-CNN](https://github.com/BioinformaticsLabAtMUN/Promotech-CNN) to obtain promoter predictions. However, prior knowledge or othe software can be used to obtain promoter locations. The promoter file contains six columns: Sequence, start of promoter, end of promoter, probability of being a promoter, strand and promoter sequence.  For instance, 

```
NC_008060.1      3647120 3647159 0.63685 -       CGTTTTCCGTCCGCAGGCTCCGGTTGGCACGATGGATGCG
NC_008060.1      3651194 3651233 0.65932 +       CCGCTGGAATTCAACATCGCCTGACCCTAGATTCAAGGCG
NC_008060.1      3652212 3652251 0.68253 +       CCCTGTTCCCCGCCCCCGTGCCGGAGATAACCTTTCCGCA
NC_008060.1      3656315 3656354 0.67378 +       CCGCCGCCGCATCGCGCAGGTTGGTGCCATAATCGTCCCG
```

## 3. Running the pipeline

sRNACharP is configured to run using the [Docker](https://www.docker.com/) container engine by default (see [Nextflow config file](nextflow.config)). If you are using Docker, download the [Nextflow config file](nextflow.config) and save it on the same directory as the sRNACharP.nf file. If you have installed the required software natively, you will need to modify the [Nextflow config file](nextflow.config).

A [sample data set](test_data) is provided. To run sRNACharP with the test data, replace ADD_PATH below with the correct path in your system and run the following command on your terminal (make sure you are executing nextflow on the directory where the sRNACharP.nf file and the Nextflow config file are located):

```
nextflow run sRNACharP.nf  --org="R_capsulatus" --dir="/ADD_PATH/sRNACharP/test_data" --genome="R_capsulatus_genome.fasta"  --sRNAs="R_capsulatus_sRNAs.bed" --genomeAnnotation="R_capsulatus_genomeAnnotation_proteincoding.bed" --promoterFile="R_capsulatus_genome_promoter_predictions.csv"

```

On your terminal you should see something like this:

```
N E X T F L O W  ~  version 20.10.0
Launching `sRNACharP.nf` [astonishing_yalow] - revision: d038e223d1
executor >  local (9)
[99/a861f0] process > getFASTAsRNAs                        [100%] 1 of 1 ✔
[c2/e039ce] process > reorderAndSortPromoterPredictions    [100%] 1 of 1 ✔
[b4/eafc8f] process > getFreeEnergySS                      [100%] 1 of 1 ✔
[d6/2d2435] process > createCRDFile                        [100%] 1 of 1 ✔
[56/51aa8c] process > runTransterm                         [100%] 1 of 1 ✔
[cb/72eb37] process > parseTranstermResults                [100%] 1 of 1 ✔
[2b/c2bdf0] process > getDistances (1)                     [100%] 1 of 1 ✔
[7d/9244ec] process > createAttributeTable (1)             [100%] 1 of 1 ✔
[ec/9eec60] process > createTetranucleotideRC_features (1) [100%] 1 of 1 ✔

    Pipeline execution summary
executor >  local (9)
[99/a861f0] process > getFASTAsRNAs                        [100%] 1 of 1 ✔
[c2/e039ce] process > reorderAndSortPromoterPredictions    [100%] 1 of 1 ✔
[b4/eafc8f] process > getFreeEnergySS                      [100%] 1 of 1 ✔
[d6/2d2435] process > createCRDFile                        [100%] 1 of 1 ✔
[56/51aa8c] process > runTransterm                         [100%] 1 of 1 ✔
[cb/72eb37] process > parseTranstermResults                [100%] 1 of 1 ✔
[2b/c2bdf0] process > getDistances (1)                     [100%] 1 of 1 ✔
[7d/9244ec] process > createAttributeTable (1)             [100%] 1 of 1 ✔
[ec/9eec60] process > createTetranucleotideRC_features (1) [100%] 1 of 1 ✔
```

## 4. Pipeline results

Analyses results are saved into the working directory.

Output files are the following (replace ORGANISM with the value of the --org parameter):

* `ORGANISM_FeatureTableNew.tsv` - table containing the characteristics per sRNA. Columns are: sequence identifier (sRNA ID), free energy of predicted secondary structure, distance to the closest predicted promoter, distance to the closest predicted terminator, distance to the closest upstream ORF, a flag indicating whether the sRNA is on the same strand as the closest upstream ORF, distance to the closest downstream ORF, a flag indicating whether the sRNA is on the same strand as the closest downstream ORF, and Tetranucleotide_RC frequencies. Here are some lines of a feature table generated:

```
	SS	PromoterDistance	DistTerm	Distance	sameStrand	DownDistance	sameDownStrand	AAAA/TTTT	AAAC/GTTT	AAAG/CTTT	AAAT/ATTT	AACA/TGTT	AACC/GGTT	AACG/CGTT	AACT/AGTT	AAGA/TCTT	AAGC/GCTT	AAGG/CCTT	AAGT/ACTT	AATA/TATT	AATC/GATT	AATG/CATT	AATT/AATT	ACAA/TTGT	ACAC/GTGT	ACAG/CTGT	ACAT/ATGT	ACCA/TGGT	ACCC/GGGT	ACCG/CGGT	ACCT/AGGT	ACGA/TCGT	ACGC/GCGT	ACGG/CCGT	ACGT/ACGT	ACTA/TAGT	ACTC/GAGT	ACTG/CAGT	AGAA/TTCT	AGAC/GTCT	AGAG/CTCT	AGAT/ATCT	AGCA/TGCT	AGCC/GGCT	AGCG/CGCT	AGCT/AGCT	AGGA/TCCT	AGGC/GCCT	AGGG/CCCT	AGTA/TACT	AGTC/GACT	AGTG/CACT	ATAA/TTAT	ATAC/GTAT	ATAG/CTAT	ATAT/ATAT	ATCA/TGAT	ATCC/GGAT	ATCG/CGAT	ATGA/TCAT	ATGC/GCAT	ATGG/CCAT	ATTA/TAAT	ATTC/GAAT	ATTG/CAAT	CAAA/TTTG	CAAC/GTTG	CAAG/CTTG	CACA/TGTG	CACC/GGTG	CACG/CGTG	CAGA/TCTG	CAGC/GCTG	CAGG/CCTG	CATA/TATG	CATC/GATG	CATG/CATG	CCAA/TTGG	CCAC/GTGG	CCAG/CTGG	CCCA/TGGG	CCCC/GGGG	CCCG/CGGG	CCGA/TCGG	CCGC/GCGG	CCGG/CCGG	CCTA/TAGG	CCTC/GAGG	CGAA/TTCG	CGAC/GTCG	CGAG/CTCG	CGCA/TGCG	CGCC/GGCG	CGCG/CGCG	CGGA/TCCG	CGGC/GCCG	CGTA/TACG	CGTC/GACG	CTAA/TTAG	CTAC/GTAG	CTAG/CTAG	CTCA/TGAG	CTCC/GGAG	CTGA/TCAG	CTGC/GCAG	CTTA/TAAG	CTTC/GAAG	GAAA/TTTC	GAAC/GTTC	GACA/TGTC	GACC/GGTC	GAGA/TCTC	GAGC/GCTC	GATA/TATC	GATC/GATC	GCAA/TTGC	GCAC/GTGC	GCCA/TGGC	GCCC/GGGC	GCGA/TCGC	GCGC/GCGC	GCTA/TAGC	GGAA/TTCC	GGAC/GTCC	GGCA/TGCC	GGCC/GGCC	GGGA/TCCC	GGTA/TACC	GTAA/TTAC	GTAC/GTAC	GTCA/TGAC	GTGA/TCAC	GTTA/TAAC	TAAA/TTTA	TACA/TGTA	TAGA/TCTA	TATA/TATA	TCAA/TTGA	TCCA/TGGA	TCGA/TCGA	TGAA/TTCA	TGCA/TGCA	TTAA/TTAA
sRNA00508	-33.2	-1000	1000	-119	1	315	1	0.011494253	0	0.011494253	0.011494253	0	0.011494253	0	0.011494253	0	0.011494253	0.011494253	0	0.011494253	0.011494253	0	0	0	0.022988506	0	0	0.022988506	0.011494253	0.011494253	0	0	0.011494253	0.022988506	0	0	0	0.011494253	0	0	0	0	0	0.022988506	0	0	0	0.011494253	0.011494253	0	0	0	0	0.011494253	0	0	0.011494253	0.011494253	0	0	0.011494253	0.011494253	0	0.011494253	0	0.011494253	0	0.011494253	0.011494253	0.034482759	0.011494253	0	0.011494253	0	0	0.011494253	0	0.022988506	0.022988506	0	0	0.022988506	0.034482759	0	0.034482759	0.022988506	0.011494253	0	0	0	0	0.011494253	0.011494253	0.022988506	0.011494253	0.045977011	0.022988506	0	0	0.011494253	0	0	0	0.011494253	0.011494253	0	0	0.011494253	0.011494253	0	0	0	0	0	0	0	0	0.034482759	0.011494253	0	0.011494253	0	0.022988506	0	0.022988506	0.022988506	0.022988506	0	0.011494253	0	0	0.011494253	0.011494253	0	0.011494253	0	0	0	0	0	0.011494253	0	0
sRNA00526	-36.2	-1000	1000	-64	1	475	1	0	0	0.012345679	0	0.012345679	0	0.012345679	0.024691358	0	0	0.012345679	0.012345679	0	0.012345679	0	0	0	0.012345679	0.012345679	0	0	0.012345679	0.012345679	0	0	0	0.037037037	0	0	0	0.012345679	0	0.012345679	0	0	0	0.012345679	0	0	0	0.024691358	0.012345679	0	0	0	0	0	0	0	0	0.012345679	0	0	0	0	0	0	0.012345679	0	0.012345679	0.012345679	0.012345679	0	0.024691358	0.012345679	0.012345679	0.012345679	0	0	0	0	0.012345679	0.012345679	0.012345679	0.024691358	0.037037037	0.012345679	0.024691358	0.024691358	0	0.012345679	0	0.012345679	0	0.024691358	0	0	0.049382716	0.024691358	0	0	0	0	0	0	0.012345679	0	0.024691358	0	0	0.012345679	0.049382716	0	0.024691358	0	0	0	0	0.012345679	0	0.012345679	0.037037037	0	0	0	0.037037037	0	0.012345679	0.037037037	0.012345679	0	0	0	0	0	0	0	0	0	0	0.024691358	0	0	0.024691358	0	0
sRNA00627	-25.8	-1000	1000	-263	0	59	1	0	0	0.013157895	0	0	0.026315789	0	0	0	0.013157895	0.013157895	0	0	0	0	0	0	0	0.013157895	0	0.013157895	0.013157895	0	0.013157895	0	0.026315789	0	0	0	0.013157895	0	0.013157895	0	0	0.013157895	0.026315789	0.013157895	0	0	0	0	0.039473684	0	0.013157895	0	0	0	0	0	0.013157895	0.013157895	0	0.013157895	0.026315789	0	0	0	0	0.013157895	0	0	0	0.013157895	0	0.013157895	0.026315789	0.026315789	0	0.013157895	0.013157895	0	0	0.026315789	0	0.013157895	0.026315789	0	0.013157895	0.026315789	0.013157895	0	0	0.013157895	0	0	0.026315789	0	0.039473684	0	0.013157895	0.013157895	0	0.013157895	0	0.013157895	0.013157895	0.026315789	0	0	0.013157895	0	0.026315789	0.013157895	0	0.013157895	0	0	0.013157895	0	0.013157895	0.013157895	0.039473684	0.013157895	0	0	0	0.026315789	0.013157895	0	0.013157895	0	0	0	0	0	0	0	0	0	0	0.013157895	0	0	0.026315789	0	0

```

* `ORGANISMRNA_PromoterDistances.txt`- sRNA distances to predicted promoters. Columns are: sRNA Sequence, Start, End, sRNA Name, Score, Strand, PromoterSequence, Start, End, Promoter Name, Score, Strand, distance
  Here are some lines of promoter distances table generated:

```
Chromosome      40220   40387   sRNA00747       15.2512228906   -       Chromosome      40202   40241   CCTCTGGTGATCGGGGAGGCATGTGCTATCCTCCCCGACA        0.61461 -       0
Chromosome      260228  260433  sRNA00798       149.01043287    -       Chromosome      264063  264102  GCGGCGGGGGCGGCGCTGCTTTCCTGGTATGATGCCCAGG        0.75134 -       -3631
Chromosome      354469  354612  sRNA00123       12.6197304188   +       Chromosome      351554  351593  CCTGACGCGACCTAGATCGTCGATTGCTATGATTGTTTCT        0.60641 +       -2877
Chromosome      397555  397640  sRNA00822       9.92096705552   -       Chromosome      401074  401113  CAGGGCGAAACCGGCCCGTTCTGGGGGCACAATGCCATCG        0.61885 -       -3435
```

* `ORGANISMRNAsTranstermRes.gtf` - compilation of TransTermHP results. Columns are: Replicon, Program, Terminator ID, Start, End, Score, Strand, and an empty field.
  Here are some lines of a predicted terminators table generated:

```
NC_016810.1     TransTermHP     TERM2   13530   13572   83      +       .
NC_016810.1     TransTermHP     TERM3   13538   13564   100     +       .
NC_016810.1     TransTermHP     TERM4   14760   14779   100     +       .

```

* Three more files are generated in the working directory: `ORGANISMRNAsSS.txt`, `ORGANISMGenomelength.txt`, and `ORGANISMRNA_sortedPromoterPredictions.bed` containing the free energy of the predicted secondary structure, the number of nucleotides in the genome, and sorted promoter predictions respectively.

## 5. Citing

If you use the pipeline, please cite:

Eppenhof EJ, Peña-Castillo L. (2019) [Prioritizing bona fide bacterial small RNAs with machine learning classifiers. PeerJ 7:e6304](https://doi.org/10.7717/peerj.6304)
