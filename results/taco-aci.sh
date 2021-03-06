#!/bin/bash

#id="sra-hisat-ref-A02"
id="gtex2377-SC2-taco02"
maxid=2376
gtfdir=/storage/home/mxs2589/shared/projects/aletsch-test/results/gtex2377-SC2

cores="20"
threads="40"

dir=`pwd`
cur=$dir/$id
mkdir -p $cur

ref=/storage/home/mxs2589/shared/data/gencode/GRCh38/gencode.v33.annotation.gtf
#ref=/gpfs/group/mxs2589/default/shared/data/ensembl/release-97/GRCh38/Homo_sapiens.GRCh38.97.gtf

pbsfile=$cur/taco.pbs
gtflist=$cur/gtflist

rm -rf $gtflist
for k in `seq 0 $maxid`
do
	if [ -s $gtfdir/$k.gtf ]; then
		echo "$gtfdir/$k.gtf" >> $gtflist
	fi
done

echo "#!/bin/bash" > $pbsfile
echo "#PBS -l nodes=1:ppn=$cores" >> $pbsfile
echo "#PBS -l mem=120gb" >> $pbsfile
echo "#PBS -l walltime=200:00:00" >> $pbsfile
echo "#PBS -A mxs2589_b_g_sc_default" >> $pbsfile
echo "module load python/2.7.14-anaconda5.0.1" >> $pbsfile
echo "$dir/taco.sh $threads $gtflist $cur $ref cov" >> $pbsfile

cd $cur
qsub $pbsfile
cd -
