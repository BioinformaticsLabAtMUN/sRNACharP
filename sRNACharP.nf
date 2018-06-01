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

   script:
   """
   fastaFromBed -fi $genomeFile -bed $bedFile -fo ${params.org}_sRNAs.fasta -s -$name
   """
}

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
* GET GENOMIC COORDINATES OF PREDICTED PROMOTERS
*************************************************************************************/
/*
* Get the 150 nt upstream sequences of the sRNAs into a FASTA file
*/
process getUpstreamSequences{
   input:
      file "${params.org}GenomeLength.txt" from lengthGenome
      file bedFile
      file genomeFile
    
   output:
      file "${params.org}_sRNAs_150ntUp.fasta" into sRNAsUpstreamFASTA
      val id into genomeID 
      
   script:
	 /*Get the genome ID*/
	lengthFile = file("${params.org}GenomeLength.txt")
	(id) = (lengthFile.text=~ /^(\S+)/)[0]
   """
    slopBed -i $bedFile -g ${params.org}GenomeLength.txt -l 150 -r 0 -s > ${params.org}_sRNAs_150ntUp.gtf
    fastaFromBed -fi $genomeFile -bed ${params.org}_sRNAs_150ntUp.gtf -fo ${params.org}_sRNAs_150ntUp.fasta -s -$name
   """
}


/*
* Split the upstream sequences FASTA file into single-sequence FASTA files 
* then execute a bprom job for each sequence
*/
sRNAsUpstreamFASTA
  .splitFasta(by:1, file:true)
  .set{singleSequences}


/* 
* Execute a bprom job to predict promoter sites in each of the upstream sequences
* and obtain a tabular representation of the promoter site associated
* to each sRNA
* Directory to access Bprom data has to be set in the nextflow.config file
*/
process getPromoterSites{
   input:
      file "seq.fa" from singleSequences //automatically does it for all the files
      
   output:
      file "bpromRes.tab" into predictedPromoters
      
   script:
   '''
    #create a multiline fasta file with maximum 70 characteres per line
    awk 'BEGIN{RS=">";FS="\\n"}NR>1{seq="";for (i=2;i<=NF;i++) seq=seq""$i;a[$1]=seq;b[$1]=length(seq)}END{for (i in a) {k=sprintf("%d", (b[i]/70)+1); printf ">%s\\n",i;for (j=1;j<=int(k);j++) printf "%s\\n", substr(a[i],1+(j-1)*70,70)}}' seq.fa> seq2.fa
    bprom seq2.fa bpromRes.txt
	 grep -E '>|LDF|-10|-35' bpromRes.txt \
    | awk -F "\\n" '!/^>/ {printf "\\t%s", $0; n= "\\n"} /^>/ {printf "\\n%s", $0; n = ""} END {printf "\\t%s", n}' \
    | grep Promoter | perl -p -e 's/ +/\\t/g' | perl -p -e 's/\\t+/\\t/g' | cut -f1,4,6,11,12,14,19,20,22 \
    | perl -p -e 's/^>//'  > bpromRes.tab
   '''
}

predictedPromoters
  .collectFile(name: file("${params.org}sRNAsBpromRes.txt"))
  .set{sRNAsPromoters}


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
   awk -F '\\t' '\$6 == "+" {print \$4,\$2,\$3,\$1} \$6 == "-" {print \$4,\$3,\$2,\$1}' $annotationFile > transtermAnnotation.crd
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
     val ID from genomeID
     
  output:
     file "transtermRes.gtf" into transtermGTF

   script:
   """  
     grep TERM transtermRes.txt | perl -p -e 's/ +/\\t/g' | perl -p -e 's/\\|.*\$//' | perl -p -e 's/TERM\\t/TERM/' \
     | cut -f2,3,5- \
     | awk -F "\\t" '\$3 > \$2 {print "${ID}", "TransTermHP", \$1, \$2, \$3, \$6, \$4, "."} \$2 > \$3 {print "${ID}", "TransTermHP", \$1, \$3, \$2, \$6, \$4, "." }'  \
     | perl -p -e 's/ +/\\t/g' |  sortBed > transtermRes.gtf
   """
}


transtermGTF
  .collectFile(name: file("${params.org}sRNAsTranstermRes.gtf"))
  .set{sRNAsTerminators}


/*************************************************************************************
* GET DISTANCE TO ORFs AND TERMINATORS
*************************************************************************************/
process getDistances{
  input:
     file bedFile
     file annotationFile
     file predTerms from sRNAsTerminators

  output:
     file "sRNAsSorted.bed" into sortedBed
  	  file "ClosestDownstreamTranscript.txt" into downORFs
  	  file "ClosestUpstreamTranscript.txt" into upORFs
  	  file "ClosestDownstreamTerm.txt" into downTerminator
  
  script:
    """
      sortBed -i $bedFile > sRNAsSorted.bed
      bedtools closest -a sRNAsSorted.bed -b ${annotationFile} -D "ref" -iu -k 2 > ClosestDownstreamTranscript.txt
      bedtools closest -a sRNAsSorted.bed -b ${annotationFile} -D "ref" -id -k 2 > ClosestUpstreamTranscript.txt
      bedtools closest -a sRNAsSorted.bed -b ${predTerms}  -D "a" -iu  > ClosestDownstreamTerm.txt
    """
}


/*************************************************************************************
* CREATE FINAL DATASET WITH R SCRIPT
* Reguires output from getDistances, sRNAsPromoters, sRNAsEnergy
*************************************************************************************/

process createAttributeTable{
  input:
     file promoters from sRNAsPromoters
     file energySS from sRNAsEnergy
     file sRNAs from sortedBed
  	  file dORFs from downORFs
  	  file uORFs from upORFs
  	  file terminators from downTerminator

  output:
  	  file "${params.org}_FeatureTable.tsv" into attributesTable 

  script:
    """
#!/usr/bin/env Rscript
    
selectClosestORF <- function(m, up = TRUE){
	# if there is only one close ORF, return it
	if (nrow(m) < 2) {
		return(m) 
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
colnames(sRNAs) <- c("Replicon","Start", "End", "ID", "Score", "Strand", "Type")
row.names(sRNAs ) <- sRNAs[,"ID"]

#SS
sRNA_E <- read.table( "$energySS", header = FALSE, sep = "\\t", stringsAsFactors = FALSE)
colnames(sRNA_E) <- c("ID", "Energy")
row.names(sRNA_E ) <- gsub("::.*\$", "", sRNA_E[,"ID"], perl = TRUE)

#Promoter sites
sRNA_bprom <- read.table("$promoters", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
colnames(sRNA_bprom) <- c("ID", "TSSPos", "LDF", "Pos10", "Seq10", "Score10","Pos35","Seq35","Score35")
row.names(sRNA_bprom) <- gsub("::.*\$", "", sRNA_bprom[,"ID"], perl = TRUE)

#Get distance to promoters
sRNA_bprom[["Pos10wrtsRNAStart"]] <- sRNA_bprom[["Pos10"]] - 150
sRNA_bprom[["Pos35wrtsRNAStart"]] <- sRNA_bprom[["Pos35"]] - 150

#Distance to ORFs
upstreamRaw <- read.table("$uORFs", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
downstreamRaw <- read.table("$dORFs", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)

#A BED file has to be used to provide the genome annotation. Otherwise the number of columns is not 15!
if (ncol(upstreamRaw) == 15){ #obtained when using a bed file with the genome annotation
	colnames(downstreamRaw) <- colnames(upstreamRaw) <- c("Replicon","Start", "End", "ID", "Score", "Strand", "Type", "RepliconORF",  "ORFStart", "ORFEnd", "ORFDescription", "ORFScore", "ORFStrand", "ORFType", "Distance")
} else{ # unexpected number of columns
	stop("Unexpected number of columns in neighbors files")	
}

upstreamC <- by(upstreamRaw, upstreamRaw[,"ID"], selectClosestORF, TRUE, simplify = FALSE)
upstreamC <- do.call("rbind", upstreamC)
row.names(upstreamC) <- upstreamC[,"ID"]

downstreamC <- by(downstreamRaw, downstreamRaw[,"ID"], selectClosestORF, FALSE, simplify = FALSE)
downstreamC <- do.call("rbind", downstreamC)
row.names(downstreamC) <- downstreamC[,"ID"]

#Terminator
sRNA_closestTerm <- read.table("$terminators", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
colnames(sRNA_closestTerm) <-  c("Replicon","Start", "End", "ID", "Score", "Strand", "Type", "RepliconTerm", "TermSource","TermID", "ORFStart", "ORFEnd", "TermScore", "TermStrand", "TermExtra", "Distance")

sRNA_closestTerm <- by(sRNA_closestTerm, sRNA_closestTerm[,"ID"], selectClosestORF, FALSE, simplify = FALSE)
sRNA_closestTerm <- do.call("rbind", sRNA_closestTerm)


#Cap distance to terminator to 1000. Erik experiments show that having very large distances to terminator decrease the performance of the classifiers
sRNA_closestTerm[["Distance"]] <- ifelse(sRNA_closestTerm[["Distance"]]>1000, 1000, sRNA_closestTerm[["Distance"]])


#Create dataset
Data <- cbind(sRNAs[,"Strand"], sRNA_E[row.names(sRNAs), "Energy"], sRNA_bprom[row.names(sRNAs), "Pos10wrtsRNAStart"], sRNA_closestTerm[row.names(sRNAs),"Distance"], 
upstreamC[row.names(sRNAs), c("Distance", "ORFStrand")], downstreamC[row.names(sRNAs), c("Distance", "ORFStrand")])

colnames(Data) <- c("Strand", "SS", "Pos10wrtsRNAStart","DistTerm", "Distance", "ORFStrand","DownDistance",  "DownORFStrand")

Data[["sameStrand"]] <- ifelse(Data[["Strand"]] == Data[["ORFStrand"]], 1, 0)
Data[["sameDownStrand"]] <- ifelse(Data[["Strand"]] == Data[["DownORFStrand"]], 1, 0)
Data[["Pos10wrtsRNAStart"]] <- ifelse(is.na(Data[["Pos10wrtsRNAStart"]]), -1000, Data[["Pos10wrtsRNAStart"]])
Data[["DistTerm"]] <- ifelse(is.na(Data[["DistTerm"]]), 1000, Data[["DistTerm"]])

DataF <- Data[,c("SS", "Pos10wrtsRNAStart","DistTerm", "Distance", "sameStrand", "DownDistance", "sameDownStrand")]
write.table(DataF, file = "${params.org}_FeatureTable.tsv", sep = "\\t", row.names = TRUE, col.names = TRUE)

    """
}


attributesTable
  .collectFile(name: file("${params.org}_FeatureTable.tsv"))


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
