#!/bin/bash
#The following script uses ReadSimulator.py to simulate Illumina short reads (PE 150bp) with the specified error sources and ratios, then generates the corresponding genone assemblies all.
#Finally quast is run to collect the genome assembly quality metrics.
#bbmap and MutationSimulator should be placed in the current directory.

ids="ls *fna | cut -f 1 -d_" #assumes all reference genome sequences to be located in the current directory. extension should be .fna.
covs="25 50 75 100 125 150"
err="0 0.01 0.025 0.05"
odups="0.01 0.05 0.15 0.30"
pdups="0.01 0.05 0.15 0.30"
rep="1 2 3 4 5"

sed -i 's/^>/>GENOMESEQ_/' *fna

for r in $rep
do
	for n in $ids;
	do
		for i in $covs
		do
			for l in $err
			do
				echo "readSimulator.py --input ${n}*fna --simulator wgsim  --outdir ${n}_cov${i}_err${l}_rep${r} --iterations 10 --readlen 150 --depth $i --opts '-e $l -r 0 -R 0 -X 0 -h -S 5'"
			done
		done
	done
done > readsim.sh

sh readsim.sh

for l in $ids
do
	for i in $( ls ${l}* | grep -v "gz\|html\|zip\|fna" )
	do
		ii=ids_${i}; echo "mv $i $ii"
	done
done | sed 's/://g' | parallel -j48

mkdir shs_reps #to run created scripts in parallel

for r in $rep
do
	for n in $ids;
	do
		for i in $covs
		do
			for l in $err
			do
				echo "mkdir ${n}_cov${i}_err${l}_rep${r}

				zcat ids_${n}_cov${i}_err${l}_rep${r}/*1.fastq.gz | awk 'sub(\"GENOMESEQ\",cnt+1, \$0){cnt++}1' | gzip > ${n}_cov${i}_err${l}_rep${r}/${n}_1.fastq.gz

				zcat ids_${n}_cov${i}_err${l}_rep${r}/*2.fastq.gz | awk 'sub(\"GENOMESEQ\",cnt+1, \$0){cnt++}1' | gzip > ${n}_cov${i}_err${l}_rep${r}/${n}_2.fastq.gz" > shs_reps/${n}_cov"${i}"_err"${l}"_rep"${r}"_make_ids.sh
			done
		done
	done
done

parallel -j48 bash :::: <(ls shs_reps/*make_ids.sh)

for r in $rep
do
	for n in $ids
	do
		for i in $covs
		do
			for l in $err
			do
				for o in $odups
				do

					echo "zgrep "^@" ${n}_cov${i}_err${l}_rep${r}/${n}_*1.fastq.gz | sed 's/\/1$//' > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids" > shs_reps/extract_ids_${n}_cov"${i}"_err"${l}"_odup"${o}"_pdup0_rep"${r}".sh

				done
			done
		done
	done
done

parallel -j48 bash :::: <(ls shs_reps/extract_ids_*.sh)

for r in $rep
do
	for n in $ids
	do
		for i in $covs
		do
			for l in $err
			do
				for o in $odups
				do

					echo "no_reads=\`cat ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids | wc -l\`

					echo \$no_reads

					no_sub=\`awk -v x="\$no_reads" -v y=$o 'BEGIN {print x*y}' | awk '{ printf(\"%d\", \$1 + 0.5) }'\`
					echo \$no_sub


					cat ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids | shuf | head -n \$no_sub > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.randids

					zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.randids ${n}_cov${i}_err${l}_rep${r}/*1.fastq.gz --no-group-separator | sed 's/^@/@odup/' | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads1.gz
					zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.randids ${n}_cov${i}_err${l}_rep${r}/*2.fastq.gz --no-group-separator | sed 's/^@/@odup/' | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads2.gz

					mkdir ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}

					zcat ${n}_cov${i}_err${l}_rep${r}/${n}*1.fastq.gz ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads1.gz | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/${n}_1.fastq.gz

					zcat ${n}_cov${i}_err${l}_rep${r}/${n}*2.fastq.gz ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads2.gz | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/${n}_2.fastq.gz" > ./shs_reps/merge_${n}_cov"${i}"_err"${l}"_odup"${o}"_pdup0_rep"${r}".sh

				done
			done
		done
	done
done

parallel -j48 bash :::: <(ls shs_reps/merge_*.sh)

pdups="0.01 0.05 0.15 0.30"

for r in $rep
do
for n in $ids
do
for i in $covs
do
for l in $err
do
for o in $odups
do
for p in $pdups
do

echo "no_reads=\`cat ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids | wc -l\`

echo \$no_reads

no_sub=\`awk -v x="\$no_reads" -v y=$p 'BEGIN {print x*y}' | awk '{ printf(\"%d\", \$1 + 0.5) }'\`
echo \$no_sub

cat ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids | shuf | head -n \$no_sub > ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}.prandids

zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}.prandids ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/*1.fastq.gz --no-group-separator | sed 's/^@/@pdup/' | seqtk seq -a > temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1.fa
zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}.prandids ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/*2.fastq.gz --no-group-separator | sed 's/^@/@pdup/' | seqtk seq -a > temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2.fa


./Mutation-Simulator/mutation-simulator.py temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1.fa args -sn 0.01
./Mutation-Simulator/mutation-simulator.py temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2.fa args -sn 0.01

./bbmap/reformat.sh in=temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1_ms.fa out=${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1_ms.fq qfake=40
./bbmap/reformat.sh in=temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2_ms.fa out=${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2_ms.fq qfake=40

gzip ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1_ms.fq
gzip ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2_ms.fq

mkdir ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}

zcat ${n}_cov${i}_err${l}_rep${r}/*1.fastq.gz ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads1.gz ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads1_ms.fq.gz | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}/${n}_1.fastq.gz
zcat ${n}_cov${i}_err${l}_rep${r}/*2.fastq.gz ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}_dupreads2.gz ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}_pdupreads2_ms.fq.gz | gzip > ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}/${n}_2.fastq.gz

rm temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}*" > shs_reps/mutate_${n}_cov"${i}"_err"${l}"_odup"${o}"_pdup"${p}"_rep"${r}".sh


done
done
done
done
done
done

parallel -j48 bash :::: <(ls shs_reps/mutate_cbot*.sh)

#odup=0 and pdup=$pdup
for r in $rep
do
	for n in $ids
	do
		for i in $covs
		do
			for l in $err
			do
			for p in $pdups
			do

				echo "no_reads=\`cat ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}.ids | wc -l\`

				echo \$no_reads

				no_sub=\`awk -v x="\$no_reads" -v y=$p 'BEGIN {print x*y}' | awk '{ printf(\"%d\", \$1 + 0.5) }'\`
				echo \$no_sub

				cat ${n}_cov${i}_err${l}_odup0.01_pdup0_rep${r}.ids | shuf | head -n \$no_sub > ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}.prandids

				zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}.prandids ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/*1.fastq.gz --no-group-separator | sed 's/^@/@pdup/' | seqtk seq -a > temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1.fa
				zgrep -A 3 -F -f ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}.prandids ${n}_cov${i}_err${l}_odup${o}_pdup0_rep${r}/*2.fastq.gz --no-group-separator | sed 's/^@/@pdup/' | seqtk seq -a > temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2.fa


				./Mutation-Simulator/mutation-simulator.py temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1.fa args -sn 0.01
				./Mutation-Simulator/mutation-simulator.py temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2.fa args -sn 0.01

				./bbmap/reformat.sh in=temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1_ms.fa out=${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1_ms.fq qfake=40
				./bbmap/reformat.sh in=temp_${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2_ms.fa out=${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2_ms.fq qfake=40

				gzip ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1_ms.fq
				gzip ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2_ms.fq

				mkdir ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}

				zcat ${n}_cov${i}_err${l}_rep${r}/*1.fastq.gz ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads1_ms.fq.gz | gzip > ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}/${n}_1.fastq.gz
				zcat ${n}_cov${i}_err${l}_rep${r}/*2.fastq.gz ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}_pdupreads2_ms.fq.gz | gzip > ${n}_cov${i}_err${l}_odup0_pdup${p}_rep${r}/${n}_2.fastq.gz

				rm temp_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}*" > shs_reps/mutate_${n}_cov"${i}"_err"${l}"_odup0_pdup"${p}"_rep"${r}".sh

			done
		done
	done
done

parallel -j48 bash :::: <(ls shs_reps/mutate_llac*odup0_*.sh)

renames=$(find ./ -name "*err*" -type d | grep -v "odup\|ids") #move assemblies with zero error

for i in $renames
do
	ii=$(echo "$i" | sed 's/_rep/_odup0_pdup0_rep/')
	echo "mv ${i} ${ii}"
done

err="0 0.01 0.025 0.05" #reassign to assemble zero error reads
odups="0 0.01 0.05 0.15 0.30"
pdups="0 0.01 0.05 0.15 0.30"
rep="1 2 3 4 5"

for r in $rep
do
	for n in $ids
	do
		for i in $covs
		do
			for l in $err
			do
				for o in $odups
				do
					for p in $pdups
					do
						echo "{ time spades.py -1 ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}/${n}_1.fastq.gz -2 ${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}/${n}_2.fastq.gz -o spades_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r} -t 16 --only-assembler ; } 2> time_spades_${n}_cov${i}_err${l}_odup${o}_pdup${p}_rep${r}.txt" > shs_reps/run_spades_${n}_cov"${i}"_err"${l}"_odup"${o}"_pdup"${p}"_"${r}".sh
					done
				done
			done
		done
	done
done

parallel -j3 bash :::: <(ls shs_reps/run_spades_*.sh)

for n in $ids
do
	for i in $rep
	do
		find ./ -name "contigs.fasta" | grep "$n" | grep rep"${i}" | xargs quast.py -o "${n}"_rep"${i}"_quast --glimmer -t 16 --min-contig 200 -r "${n}"*.fna
	done
done

for n i $ids
do

	find ./${n}_rep*_quast -name "transposed_report.tsv" | xargs cat | grep -v Assembly > temp
	cat <(grep Assembly ${n}_rep1_quast/transposed_report.tsv) temp > ${n}_reps_quast_transposed_report.tsv

done
