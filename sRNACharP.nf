#!/usr/bin/env nextflow

/*
* This sRNA characterization pipeline obtains for each sRNA:
* the free energy of each sRNA secondary structure using centroidFold,
* the distance to the closest promoter predicted by bprom, 
* the distance to the closest terminator predicted by transterm,
* the distance to the closest ORFs listed in the genome annotation.
*
* The pipeline required input is:
* 1) a FASTA file with the genome (genome),
* 2) a BED file with the sRNAs genomic coordinates (sRNAs), and 
* 3) a BED file with the genome annotation (only protein coding genes) (genomeAnnotation).
* 4) a CSV file with predicted promoters
* To get a BED file from a GFF file one can use:
* awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";'
*
* Optional input includes:
* 1) A transterm predictions file (transtermFile)
* It can be downloaded from http://transterm.cbcb.umd.edu/cgi-bin/transterm/predictions.pl. 
* If not provided, the pipeline will run TransTerm
*
* Other required parameters:
* org = organism name
* dir = full path to the directory with the input datasets: genome, sRNAs, genomeAnnotation, and transtermFile.
*
* You can redistribute this pipeline and/or modify it under the terms of the 
* GNU General Public License as published by the Free Software Foundation, 
* either version 3 of the License, or (at your option) any later version
* (see <http://www.gnu.org/licenses/>)
*
* This Nextflow pipeline is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*
* Author: Lourdes Peña-Castillo
* Version 1.0.0. 2018
*/

params.transtermFile = null
params.help = false

//print usage
if (params.help) {
  log.info ''
  log.info 'sRNACharP: sRNA Characterization Pipeline'
  log.info '----------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info "    ${workflow.projectDir.baseName} [options]"
  log.info ''
  log.info 'Options:'
  log.info '--org           ORGANISM         the organism name [required]'
  log.info '--dir           DIRECTORY        the directory containing the input files [required]'
  log.info '--genome        GENOME_FILE      the genome file in FASTA format [required]'
  log.info '--annot         ANNOTATION_FILE  the annotation file in BED format (including only protein coding genes)[required]'
  log.info '--sRNAs         SEQUENCE_FILE    the sRNAs file in BED format [required]'
  log.info '--transtermFile TRANSTERM_FILE   the TransTerm predictions file [optional].'
  log.info '--promoterFile  PREDICTED_PROMOTER_FILE the Promoter predictions file [required]' 
  exit 1
}

/*
* Verify that the files provided exist
*/
genomeFile = file("${params.dir}/${params.genome}")
if( !genomeFile.exists() ) {
  exit 1, "The specified input genome file does not exist: ${params.genome}"
}

bedFile = file("${params.dir}/${params.sRNAs}")
if( !bedFile.exists() ) {
  exit 1, "The specified input BED file does not exist: ${params.sRNAs}"
}

annotationFile = file("${params.dir}/${params.genomeAnnotation}")
if( !annotationFile.exists() ) {
  exit 1, "The specified input annotation file does not exist: ${params.genomeAnnotation}"
}

exptermFile = file("${TERM_DATA}")

promoterFile = file("${params.dir}/${params.promoterFile}")
if( !promoterFile.exists() ) {
  exit 1, "The specified input promoter prediction file does not exist: ${params.promoterFile}"
}
/*
* Set global variables
*/
name = "name"    
runTransterm = true
if (params.transtermFile != null ) {  //TransTerm predictions are provided
  transtermFile = file("${params.dir}/${params.transtermFile}")
  if( !transtermFile.exists() ) {
    exit 1, "The transterm results file does not exist: ${params.transtermFile}"
  }
  runTransterm = false
} //TransTerm will be run thus verify that the expterm.dat file exists



/*************************************************************************************
* GENERATE OTHER REQUIRED INPUT FILES
*************************************************************************************/
/*
* Create a file with the length of the genome. 
* This file contains one line with the genome ID and the genome length
*/
lengthFile = file("${params.org}GenomeLength.txt")

Channel
   .from(genomeFile)
   .splitFasta( record: [id: true, seqString: true ])
   .collectFile(name: lengthFile) { record -> record.id + "\t" + record.seqString.length() + "\n"}
   .set{lengthGenome}

/*
* Get a FASTA file with the sRNAs sequences
*/
process getFASTAsRNAs{
   input:
      file genomeFile
      file bedFile

   output:
      file "${params.org}_sRNAs.fasta" into sRNAsFASTA
      file "sequences.txt" into sequences

   script:
   """
   fastaFromBed -fi $genomeFile -bed $bedFile -fo ${params.org}_sRNAs.fasta -s -$name
   fastaFromBed -fi $genomeFile -bed $bedFile -fo sequences.txt -s -$name -tab
   """
}

// sequences
//   .collectFile(name: file("${params.org}_sequences.txt"))


/*
*Reorder and sort Promoter Prediction File
*/
process reorderAndSortPromoterPredictions{
   input:
	file promoterFile
   output:
	file "${params.org}_sortedPromoterPredictions.bed" into sortedPromoterPredictions 
   script:
	"""
	awk -v OFS="\t" 'NR>1{print \$1,\$2,\$3,\$6,\$4,\$5}' ${promoterFile} > ${params.org}_promoterPredictions.bed
	sortBed -i ${params.org}_promoterPredictions.bed > ${params.org}_sortedPromoterPredictions.bed

	"""
}

sortedPromoterPredictions
  .collectFile(name: file("${params.org}_sortedPromoterPredictions.bed"))


/*************************************************************************************
* GET FREE ENERGY OF PREDICTED SECONDARY STRUCTURE
*************************************************************************************/
/*
* Run centroidFold to get the free energy of the predicted secondary structure (SS)
* then process the output to obtain a tab-delimited file with the sRNA ID and the free Energy of the SS
*/
process getFreeEnergySS{
   input:
      file  "${params.org}_sRNAs.fasta" from sRNAsFASTA

   output:
      file "${params.org}_sRNAs_freeEnergy.txt" into freeEnergy4sRNA

   script:
   """
    centroid_fold -e "CONTRAfold" -g 4 "${params.org}_sRNAs.fasta" > ${params.org}_CentroidFold_out.txt
    grep "e=" ${params.org}_CentroidFold_out.txt | perl -p -e 's/.*e=//' | perl -p -e 's/\\)\$//' > freeEnergy.txt
    grep ">" ${params.org}_CentroidFold_out.txt > sequencesNames.txt
    paste sequencesNames.txt freeEnergy.txt > ${params.org}_sRNAs_freeEnergy.txt
    perl -pi -e 's/>//' ${params.org}_sRNAs_freeEnergy.txt
   """
}

freeEnergy4sRNA
  .collectFile(name: file("${params.org}sRNAsSS.txt"))
  .set{sRNAsEnergy}

/*************************************************************************************
* GET GENOMIC COORDINATES OF PREDICTED RHO-INDEPENDENT TERMINATORS
*************************************************************************************/
/* 
* If a transterm file has been provided then
* generate a GTF file with the results
* otherwise
* execute a transterm job giving as input the FASTA file with the genome,
* the expterm.dat file distributed with transterm, and a file with the genomic
* coordinates of the protein-coding genes
* then generate a GTF file with the results
*/

process createCRDFile{
 	input:
		file annotationFile
		
   output:
   	file "transtermAnnotation.crd" into CRDFile
   
   script:	
   """
   #transterm start coordinates are 1-based
   awk -F '\\t' '\$6 == "+" {print \$4,\$2+1,\$3,\$1} \$6 == "-" {print \$4,\$3,\$2+1,\$1}' $annotationFile > transtermAnnotation.crd
   """
}

process runTransterm{
 	input:
 	  file genomeFile
 	  file "transtermAnnotation.crd" from CRDFile
 
   output:
     file "transtermRes.txt" into predictedTerminators

   script:
   if (runTransterm) {
    """
      transterm -p $exptermFile $genomeFile transtermAnnotation.crd > transtermRes.txt
    """
   } else {
   """
      ln -s "${params.dir}/${params.transtermFile}" transtermRes.txt
   """
   }
}

process parseTranstermResults{
  input:
     file "transtermRes.txt" from predictedTerminators
     
  output:
     file "transtermRes.gtf" into transtermGTF

   script:
    lengthFile = file("${params.org}GenomeLength.txt")
	(id) = (lengthFile.text=~ /^(\S+)/)[0]
   """  
     grep TERM transtermRes.txt | perl -p -e 's/ +/\\t/g' | perl -p -e 's/\\|.*\$//' | perl -p -e 's/TERM\\t/TERM/' \
     | cut -f2,3,5- \
     | awk -F "\\t" '\$3 > \$2 {print "${id}", "TransTermHP", \$1, \$2, \$3, \$6, \$4, "."} \$2 > \$3 {print "${id}", "TransTermHP", \$1, \$3, \$2, \$6, \$4, "." }'  \
     | perl -p -e 's/ +/\\t/g' |  sortBed > transtermRes.gtf
   """
}


transtermGTF
  .collectFile(name: file("${params.org}sRNAsTranstermRes.gtf"))
  .set{sRNAsTerminators}


/*************************************************************************************
* GET DISTANCE TO ORFs, TERMINATORS AND PREDICTED PROMOTERS
*************************************************************************************/
process getDistances{
  input:
     file bedFile
     file annotationFile
     file predTerms from sRNAsTerminators
     file sortedPromoterFile from sortedPromoterPredictions

  output:
          file "sRNAsSorted.bed" into sortedBed
  	  file "ClosestDownstreamTranscript.txt" into downORFs
  	  file "ClosestUpstreamTranscript.txt" into upORFs
  	  file "ClosestDownstreamTerm.txt" into downTerminator
          file "ClosestPromoterDistance.txt" into closestPromoterDistances 
  script: 
   """
      sortBed -i $bedFile > sRNAsSorted.bed
      bedtools closest -a sRNAsSorted.bed -b ${annotationFile} -D "ref" -iu -k 2 > ClosestDownstreamTranscript.txt
      bedtools closest -a sRNAsSorted.bed -b ${annotationFile} -D "ref" -id -k 2 > ClosestUpstreamTranscript.txt
      bedtools closest -a sRNAsSorted.bed -b ${predTerms}  -D "a" -iu  > ClosestDownstreamTerm.txt
      cut -f 1-6 sRNAsSorted.bed > firstSixColumnsFromsRNA.bed
      bedtools closest -D a -s -id -a firstSixColumnsFromsRNA.bed -b ${sortedPromoterFile} > ClosestPromoterDistance.txt
    """
}

closestPromoterDistances
  .collectFile(name: file("${params.org}_closestPromoterDistances.txt"))
  .set{promoterDistances}


/*************************************************************************************
* CREATE FINAL DATASET WITH R SCRIPT
* Reguires output from getDistances,sRNAsEnergy
*************************************************************************************/

process createAttributeTable{
  input:
         file energySS from sRNAsEnergy
         file sRNAs from sortedBed
  	 file dORFs from downORFs
  	 file uORFs from upORFs
  	 file terminators from downTerminator
         file promoters from promoterDistances

  output:
  	 file "FeatureTable.tsv" into attributesTable 

  script:
    """
#!/usr/bin/env Rscript
    
selectClosestORF <- function(m, up = TRUE){
	# There is only one close ORF
	if (nrow(m) < 2) {
        	if (m[,"Distance"] == -1 && m[,7] == ".") { ## there is no close Feature Replicon == . and BedTools has returned -1 as distance, then make distance 1000
            		tmp <-m
            		tmp[,"Distance"] <- 1000
            		return(tmp)
        	} else {
            		return(m)
        	}
	}
	#order by distance
	o <- order(abs(m[,"Distance"]), decreasing = FALSE) 
	m <- m[o,]
	    # To be an overlap upstream the start of the sRNA has to be within the ORF start and end
	if ( (( up && m[1,"ORFStart"] <= m[1,"Start"] && m[1,"ORFEnd"] >= m[1,"Start"] ) || 
		# To be an overlap downstream the end of the sRNA has to be within the ORF start and end
	      (!up && m[1,"ORFStart"] <= m[1,"End"] && m[1,"ORFEnd"] >= m[1,"End"] ) || 
	      # There is an overlap with a smaller ORF totally enclosed by the sRNA
	      ( m[1,"ORFStart"] >= m[1,"Start"]  && m[1,"ORFEnd"] <= m[1,"End"] ))  ){
			#there is an overlap with first ORF
			m[1,] #Return first ORF
	} else if (abs(m[1,"Distance"]) > 0) { 
		#there is not overlap and the distance returned by bedtools is greater than zero, then
		#return the closest ORF
			m[1,]
	} else {
		#there is not overlap in the right end but bedtools is returning an overlap, then return the second
		#closest ORF
			m[2,]
	}
}

#sRNAs
sRNAs <- read.table( "$sRNAs", header = FALSE, sep = "\\t", stringsAsFactors = FALSE)
#the first 6 columns are required
nCols <- ncol(sRNAs)
if (nCols < 6 || nCols > 7){
stop("sRNAs BED file must contain 6 or 7 columns: Replicon,Start, End, ID, Score, Strand, Type. The last column Type is optional.")
}

#Take only the first columns
sRNAs<- sRNAs[,1:6]
colnames(sRNAs) <- c("Replicon","Start", "End", "ID", "Score", "Strand")
row.names(sRNAs ) <- sRNAs[,"ID"]


#SS
sRNA_E <- read.table( "$energySS", header = FALSE, sep = "\\t", stringsAsFactors = FALSE)
colnames(sRNA_E) <- c("ID", "Energy")
row.names(sRNA_E ) <- gsub("::.*\$", "", sRNA_E[,"ID"], perl = TRUE)

#Distance to Promoters
promotersRaw <- read.table("$promoters", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
colnames(promotersRaw) <- c("Sequence","Start", "End", "ID", "Score", "Strand", "PromoterSequence",  "ORFStart", "ORFEnd", "PromoterID", "PromoterScore", "PromoterStrand", "Distance")
closestPromoterDistance <- by(promotersRaw, promotersRaw[,"ID"], selectClosestORF, TRUE, simplify = FALSE)
closestPromoterDistance <- do.call("rbind", closestPromoterDistance)
row.names(closestPromoterDistance) <- closestPromoterDistance[,"ID"]
closestPromoterDistance[["Distance"]] <- ifelse((abs(closestPromoterDistance[["Distance"]])>=1000 | closestPromoterDistance[["Distance"]] == -1),-1000, closestPromoterDistance[["Distance"]])

#Distance to ORFs
upstreamRaw <- read.table("$uORFs", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
downstreamRaw <- read.table("$dORFs", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)

#A BED file has to be used to provide the genome annotation. Otherwise the number of column number is not correct!
if (nCols == 7 && ncol(upstreamRaw) == 15){ #obtained when using a bed file with the genome annotation
downstreamRaw <- downstreamRaw[,-7]
upstreamRaw <- upstreamRaw[,-7]
colnames(downstreamRaw) <- colnames(upstreamRaw) <- c("Replicon","Start", "End", "ID", "Score", "Strand", "RepliconORF",  "ORFStart", "ORFEnd", "ORFDescription", "ORFScore", "ORFStrand", "ORFType", "Distance")
} else if ( nCols == 6 && ncol(upstreamRaw) == 14 ) {
colnames(downstreamRaw) <- colnames(upstreamRaw) <- c("Replicon","Start", "End", "ID", "Score", "Strand", "RepliconORF",  "ORFStart", "ORFEnd", "ORFDescription", "ORFScore", "ORFStrand", "ORFType", "Distance")
} else{ # unexpected number of columns
stop("Unexpected number of columns in neighbor ORF files")
}

upstreamC <- by(upstreamRaw, upstreamRaw[,"ID"], selectClosestORF, TRUE, simplify = FALSE)
upstreamC <- do.call("rbind", upstreamC)
row.names(upstreamC) <- upstreamC[,"ID"]

downstreamC <- by(downstreamRaw, downstreamRaw[,"ID"], selectClosestORF, FALSE, simplify = FALSE)
downstreamC <- do.call("rbind", downstreamC)
row.names(downstreamC) <- downstreamC[,"ID"]

#Terminator
sRNA_closestTerm <- read.table("$terminators", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)

if (nCols == 7 ){
sRNA_closestTerm <- sRNA_closestTerm[,-7]
}

if (ncol(sRNA_closestTerm) != 15) {
stop("Unexpected number of columns in closest Term file")
}

colnames(sRNA_closestTerm) <-  c("Replicon","Start", "End", "ID", "Score", "Strand","RepliconTerm", "TermSource","TermID", "ORFStart", "ORFEnd", "TermScore", "TermStrand", "TermExtra", "Distance")

sRNA_closestTerm <- by(sRNA_closestTerm, sRNA_closestTerm[,"ID"], selectClosestORF, FALSE, simplify = FALSE)
sRNA_closestTerm <- do.call("rbind", sRNA_closestTerm)

#Cap distance to terminator to 1000. Erik experiments show that having very large distances to terminator decrease the performance of the classifiers
sRNA_closestTerm[["Distance"]] <- ifelse(sRNA_closestTerm[["Distance"]]>1000, 1000, sRNA_closestTerm[["Distance"]])

#Create dataset
Data <- cbind(sRNAs[,"Strand"], sRNA_E[row.names(sRNAs), "Energy"], closestPromoterDistance[row.names(sRNAs),"Distance"],sRNA_closestTerm[row.names(sRNAs),"Distance"], 
upstreamC[row.names(sRNAs), c("Distance", "ORFStrand")], downstreamC[row.names(sRNAs), c("Distance", "ORFStrand")])

colnames(Data) <- c("Strand", "SS", "PromoterDistance","DistTerm", "Distance", "ORFStrand","DownDistance",  "DownORFStrand")

Data[["sameStrand"]] <- ifelse(Data[["Strand"]] == Data[["ORFStrand"]], 1, 0)
Data[["sameDownStrand"]] <- ifelse(Data[["Strand"]] == Data[["DownORFStrand"]], 1, 0)
Data[["DistTerm"]] <- ifelse(is.na(Data[["DistTerm"]]), 1000, Data[["DistTerm"]])

DataF <- Data[,c("SS", "PromoterDistance","DistTerm", "Distance", "sameStrand", "DownDistance", "sameDownStrand")]
write.table(DataF, file = "FeatureTable.tsv", sep = "\\t", row.names = TRUE, col.names = TRUE)

    """
}


process createTetranucleotideRC_features{
  input:
    file  "sequences.txt" from sequences
    file "FeatureTable.tsv" from attributesTable 

  output:
  file 'featureTableNew.tsv' into featureTableNew
    
  '''
  #!/usr/bin/env python3
  import os
  import itertools
  import pandas as pd
  from skbio import Sequence
  from skbio import DNA

  featureTable = pd.read_csv('Featuretable.tsv', "\t")
  Sequences = pd.read_csv('sequences.txt', header=None, sep="\t")
  Sequences.iloc[:,0] = Sequences.iloc[:,0].str.split("(",expand=True).iloc[:,0]

  # Number of sequences
  rowsCount = Sequences.shape[0]

  # Number of Nucleotide
  NucleotideNum = 4

  # Appending NucleotidesColumn(new features) and initialize with zeros
  iter = itertools.product('ACGT', repeat=NucleotideNum)
  iterJoin = []

  for i in iter:
    colLable = "".join(i)
    iterJoin.append(colLable)
    colValues_zeros = [0]*rowsCount
    featureTable[colLable] = colValues_zeros


  # Filling NucleotidesColumn with their frequency for each sequence
  for idIndex in range(rowsCount):
    id = Sequences.iloc[idIndex,0]
    seq = Sequences.iloc[idIndex,1]
    s = Sequence(seq)
    freqs = s.kmer_frequencies(NucleotideNum, relative=True, overlap=True)
    for nucleotide in freqs:
      if nucleotide in iterJoin :
        featureTable.loc[id , nucleotide] = freqs[nucleotide]


  ###### Creating Reverse Complement Features ######

  # Extracting Tetranucleotides' name
  featuresName = [ x for x in featureTable.keys()[7:]]

  # Calculating new features based on tetraNucleotides & their reverse complement
  while len(featuresName)>0 : 
    tetraNucleotideName = featuresName[0]

    seq = DNA(tetraNucleotideName)
    tetraNucleotide_ReverseComp = str(seq.reverse_complement())

    if(tetraNucleotide_ReverseComp == tetraNucleotideName):
        featureTable[f'{tetraNucleotideName}/{tetraNucleotide_ReverseComp}'] = featureTable[tetraNucleotideName]
    else:
        featureTable[f'{tetraNucleotideName}/{tetraNucleotide_ReverseComp}'] = featureTable[tetraNucleotideName] + featureTable[tetraNucleotide_ReverseComp]
        featureTable.drop(tetraNucleotide_ReverseComp, axis=1, inplace=True)
        featuresName.remove(tetraNucleotide_ReverseComp)

    featureTable.drop(tetraNucleotideName, axis=1, inplace=True)
    featuresName.remove(tetraNucleotideName)   

  featureTable.to_csv("featureTableNew.tsv", sep="\t")
  '''
}

featureTableNew
  .collectFile(name: file("${params.org}_featureTableNew.tsv"))

  

/*************************************************************************************
* END of Workflow
*************************************************************************************/
workflow.onComplete { 
   println( 
   """
    Pipeline execution summary
    ---------------------------
    Run as      : ${workflow.commandLine}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """)
}

workflow.onError = {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
