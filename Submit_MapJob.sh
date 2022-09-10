
#!/bin/bash

#check for required command line argument

if [[ -z "$1" ]]
then
		echo "
		################
		ERROR: You must include the path to your accession file in the command line call.
		eg. sh submit.sh Path/To/Your/AccessionFile.txt
		##################
		"

		exit
fi

#iterates through list of accessions and passes to mapping script

accession="/scratch/zlewis/JGI_Sept2022/AccessionLists/BatchOneSraAccessionsOnly.txt"
fastqPath="/scratch/zlewis/JGI_Sept2022/BatchOne/downloadSRA/FastqFiles"
outdir="/scratch/zlewis/JGI_Sept2022/BatchOne/ParTestt"

mkdir ${outdir}
mkdir ${outdir}/logs

#make output file folders
trimmed="${outdir}/TrimmedFastQs/${line}"
mkdir "${outdir}/TrimmedFastQs"
mkdir $trimmed
mkdir "${outdir}/BamFiles"
mkdir "${outdir}/counts"

while read -r line

	do
	sleep 10
	echo "$line mapping job submitted"
	sbatch --export=ALL,line="${line}" parallel.sh & done <"$1"
