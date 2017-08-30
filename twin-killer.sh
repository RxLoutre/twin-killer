#!/bin/bash
#A script designed to reduce duplication level within genome assembly based on Busco result
#Giving a Busco result (full table file from busco result directory), we are looking for contigs
#that can be removed because they contain redondant informations.
#When gene are "duplicated" in Busco results, they appear twice or more on several contigs
#This script analyze if the smallest contigs (they can be several) contains only redondant informations
#or also contaisn uniq informations. To see wheather it carries uniq informations or not, we are using ORF finder
#to discover all ORF contained in one contig
#If all ORF of one contig can be retrieved on the rest of the assembly (the genome file - the concerned contig)
#then the contig is considered as redundant and then eliminated from the assembly

#Usage:
#./manual_duplication_cleaning.sh busco_result.tsv my_assembly.fai my/output/dir my_assembly.fasta

#/!\ Fasta file must be indexed and masked

#b : Busco full table tsv output
#o : Output directory where to write results
#f : Fasta file to work with. Must be indexed and masked, where these files can be found in the same directory
#t : Tag of the analysism while be usefull to name output files
#m : Minimum number of ORF to be found in shortest contig to be kept
#h : Print available options and their meanings
#r : Redo a short analysis, skipping the generation of the resume file, allow to play quickly with -m parameter


#Example : ./twin-killer.sh -b '/home/loutre/buscov2/80-3clades-busco-analysis.tsv' -o '/media/loutre/SUZUKII/assembly/2-duplication_removal/2017/test' -f '/media/loutre/SUZUKII/assembly/merged/3-suzukii-polished-80-merged-renamed.fasta' -t twin-killer-test -m 10 -r


SKIP=0

while getopts b:o:f:t:m:hr OPTION
do
	case "${OPTION}"
	in
	b)	FULL_TABLE_FILE=${OPTARG};;
	o)	OUTPUT_DIR=${OPTARG};;
	f)	FASTA_FILE=${OPTARG};;
	t)	TAG=$OPTARG;;
	m)	TEMP=${OPTARG}
		NB_MIN_ORF=$((TEMP+0))
		;;
	h)	echo "
		******Twin killer V1.0******
		Usage :
		./twin-killer.sh -b busco-full-table.tsv -o path/to/my/output/dir -f my/assembly.fasta -t myFancyTag -m myThreshold (-r -h)
		Options :
		-b : Busco full table tsv output
		-o : Output directory where to write results
		-f : Fasta file to work with. Must be indexed and masked, where these files can be found in the same directory
		-t : Tag of the analysis while be usefull to name output files
		-m : Minimum number of ORF to be found in shortest contig to be kept
		-h : Print available options and their meanings
		-r : Redo a short analysis, skipping the generation of the resume file, allow to play quickly with -m parameter

		"
		exit;;
	r)	SKIP=1;;
	esac
done

FASTA_FILE_MASKED=$FASTA_FILE".masked"
INDEX_FILE=$FASTA_FILE".fai"

#Check existence of index and masked file for given fasta
if [ ! -f $FASTA_FILE_MASKED ]; then
	echo "Error : the given fasta file is not masked."
	echo "You should run repeat masker before using twin-killer."
	echo "Please check a .masked file with the same name and path as your fasta file exist."
	echo "Given fasta file : $FASTA_FILE"
	exit
fi
if [ ! -f $INDEX_FILE ]; then
	echo "Warning : the given fasta file is not indexed."
	echo "Creating index..."
	samtools faidx $FASTA_FILE
fi

#Name of the main output
RESUME=$OUTPUT_DIR/"resume.txt"
TO_DEL=$OUTPUT_DIR/"contigs_to_delete.txt"
CLEANED_FASTA=$OUTPUT_DIR/$TAG".fasta"
ALT_FASTA=$OUTPUT_DIR/$TAG"_alt.fasta"	
if [ ! -f $TO_DEL ]; then
	touch $TO_DEL
else
	> $TO_DEL
fi


#If the -r option is not present, run a full analysis
if [ $SKIP = 0 ]
then

	DUPLICATED_FILE=$OUTPUT_DIR/"duplicated.txt"
	CONFLICT_FILE=$OUTPUT_DIR/"CONFLICT_to_blast.txt"

	#1) Keep only the lines containing "Duplicated" in Busco full table
	cat $FULL_TABLE_FILE | awk '$2 == "Duplicated" {print $0}' > $DUPLICATED_FILE

	#2) Call python script that is able to determine which contig conflict to compare
	#The output file is a list of contig names conflict
	#Example :
	#Longest_contig Shortest_contig
	#tig00000001 tig000000002
	#tig00000250 tig000025501
	python src/busco_duplicated.py $DUPLICATED_FILE $INDEX_FILE $CONFLICT_FILE
	sort $CONFLICT_FILE | uniq > $CONFLICT_FILE'_sorted'

	ORF_DIR=$OUTPUT_DIR/"ORF"
	ORF_UNIQ=$ORF_DIR/"uniq"
	BLAST_DIR=$OUTPUT_DIR/"Blast"
	CONTIG_DIR=$OUTPUT_DIR/"contigs_fasta"
	
	

	# CREATING OUTPUT DIRS AND FILES AND CHECKING IF THEY ALREADY EXIST
	if [[ ! -d $CONTIG_DIR ]]; then	
		mkdir $CONTIG_DIR
	fi
	if [[ ! -d $ORF_DIR ]]; then	
		mkdir $ORF_DIR
	fi
	if [[ ! -d $BLAST_DIR ]]; then	
		mkdir $BLAST_DIR
	fi
	if [[ ! -d $ORF_UNIQ ]]; then	
		mkdir $ORF_UNIQ
	fi
	if [ ! -f $RESUME ]; then
		touch $RESUME	
	else
		> $RESUME
	fi
	echo "#Short	Long	NbORF	NbUniqORF	Deleted" > $RESUME
	
	clear
	echo "*********************ANALYSIS START*************************"
	#Read the conflict file previously created
	while read LINE; do
		echo "Analyzing conflict $LINE";
		LONG=$(echo $LINE | awk '{print $1;}')
		SHORT=$(echo $LINE | awk '{print $2;}')
		echo "     Long : " $LONG
		echo "     Short : "$SHORT
		ORF_FASTA=$ORF_DIR/$SHORT"_orf.fasta"
		ORF_FASTA_FILTERED=$ORF_DIR/$SHORT"_orf_filtered.fasta"
		BLAST_FILE=$BLAST_DIR/$SHORT"_orf.txt"
		SHORT_ONLY=$ORF_UNIQ/$SHORT"_exclusives.txt"
		#3)Extract fasta sequence for the 2nd contig (aka the smallest one)
		SEQ_SHORT=$(samtools faidx $FASTA_FILE_MASKED $SHORT)	
		#If the file of this contig doesn't exist
		if [ ! -f $CONTIG_DIR/$SHORT.txt ]; then
			#create the file
			echo $SEQ_SHORT > $CONTIG_DIR/$SHORT.txt
			sed s/' '/'\n'/g $CONTIG_DIR/$SHORT.txt > $CONTIG_DIR/$SHORT".fasta"
			rm $CONTIG_DIR/$SHORT.txt
		fi
		#4 Analyze ORF for each contig file	
		echo "      Launching ORFfinder...."
		#Outfmt 1 is for fasta nucleotide sequence
		src/ORFfinder -in $CONTIG_DIR/$SHORT".fasta" -outfmt 1 > $ORF_FASTA
		#5 Filter ORF so it allows only ORF > 200 bases
		echo "      Filtering ORFs...."
		perl src/removesmalls.pl 200 $ORF_FASTA > $ORF_FASTA_FILTERED
		#6 Create an index for ORF file with all ORF names
		echo "      Indexing ORFs...."
		grep '^>' $ORF_FASTA_FILTERED | awk '{print substr($1,2)}' > $ORF_FASTA_FILTERED".idx"
		#7 Retrieve analyzed contig from assembly file before blast
		python src/Remove_seqs_from_fasta_using_name.py $FASTA_FILE_MASKED $SHORT > $BLAST_DIR"/to_blast.fasta"
		#8 Blast ORF found in short against all assembly	
		echo "      Blasting ORF against assembly...."
		blastn -query $ORF_FASTA_FILTERED -subject $BLAST_DIR"/to_blast.fasta" -out $BLAST_FILE -outfmt 6
		#9 Filter blast results
		python src/filter_blast_hit.py $BLAST_FILE $BLAST_FILE"_filtered"
		#10 Check if every ORF is present in filtered blast file
		echo "      Checking ORF redundancy...."
		NB_ORF=$(cat $ORF_FASTA_FILTERED".idx" | wc -l)
		ORF_MISSING=0
		while read ORF_ID; do
			COUNT=$(grep -w $ORF_ID $BLAST_FILE"_filtered" | wc -l)
			if [ $COUNT = 0 ]
			then
				ORF_MISSING=$((ORF_MISSING + 1))
				echo $ORF_ID >> $SHORT_ONLY
			fi
		done < $ORF_FASTA_FILTERED".idx"
			
		if [ "$ORF_MISSING" -gt "$NB_MIN_ORF" ]
		then
			echo "      $ORF_MISSING ORF of contig $SHORT are missing in the assembly."
			echo "      Check $SHORT_ONLY file for more details"
			echo "      Contig $SHORT Kept in assembly."
			
		
		else
			echo "      $ORF_MISSING ORF of $SHORT are missing in the assembly."
			echo "      $SHORT is added to redundant contigs."
			echo $SHORT >> $TO_DEL
			echo "$SHORT	$LONG	$NB_ORF	$ORF_MISSING	1" >> $RESUME
		fi
	

	done < $CONFLICT_FILE'_sorted'


	echo "*********************ANALYSIS DONE*************************"

#If the -r option is present, run a short analysis using previously created in the same output dir as specified resume.txt file
else
	if [ ! -f $RESUME ]; then
		echo "Error : $RESUME file do not exist."
		echo "You should have run twin-killer at least once in the specified directory."
		exit
	else
		#Read line by line the resume file
		while read LINE; do
			case "$LINE" in \#*) continue ;; esac
			NB_UNIQ_ORF=$(echo $LINE | cut -f 4 -d" ")
			CONTIG=$(echo $LINE | cut -f 1 -d" ")
			if [ "$NB_UNIQ_ORF" -gt "$NB_MIN_ORF" ]
			then
				echo "      $NB_UNIQ_ORF ORF of contig $CONTIG are missing in the assembly."
				echo "      Contig $CONTIG Kept in assembly."
			else
				echo "      $NB_UNIQ_ORF ORF of contig $CONTIG are missing in the assembly."
				echo "      $CONTIG is added to redundant contigs."
				echo $CONTIG >> $TO_DEL
			fi
			
		done < $RESUME	
	fi

fi

# 11 Delete all redundant contigs from non masked assembly
echo "Deleting bad contigs from assembly...."
python src/Remove_seqs_from_list.py $FASTA_FILE $TO_DEL > $CLEANED_FASTA
# 12 Isolated all contigs that have been deletd in the assembly
echo "Gathering deleted contigs in alternate file...."
python src/Keep_seqs_from_list.py $FASTA_FILE $TO_DEL > $ALT_FASTA
# 13 Printing a nice resume at the end
echo "Deleted contigs : "
cat $TO_DEL
echo "Numbers of deleted contigs : "
cat $TO_DEL | wc -l
echo "Cleaned file : "$CLEANED_FASTA
echo "Alternative file : "$ALT_FASTA
