mkdir 01_cutadapt
mkdir 02_bowtie2
mkdir 03_sort_uni_bam
mkdir 04_rmdup_bam
mkdir 05_bw

for i in `ls *.2.fq.gz`
do 
base=$(basename $i ".2.fq.gz")
/media/xionglab/data1/02_Software/Anaconda/bin/cutadapt -u 46 -o ./01_cutadapt/${base}.cut.1.fq.gz ${base}.1.fq.gz
/media/xionglab/data1/02_Software/Anaconda/bin/cutadapt -u 52 -o ./01_cutadapt/${base}.cut.2.fq.gz ${base}.2.fq.gz
/media/xionglab/data1/02_Software/Anaconda/bin/cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 30 --max-n 0.1 --trim-n -o ./01_cutadapt/${base}.trimmed.1.fq.gz -p ./01_cutadapt/${base}.trimmed.2.fq.gz ./01_cutadapt/${base}.cut.1.fq.gz ./01_cutadapt/${base}.cut.2.fq.gz
/media/xionglab/data1/02_Software/Anaconda/bin/bowtie2 -x  /media/xionglab/data1/01_Database/02_bowtie2/hg19/hg19 -1 ./01_cutadapt/${base}.trimmed.1.fq.gz -2 ./01_cutadapt/${base}.trimmed.2.fq.gz --very-sensitive-local --no-unal -X 3000 -S ./02_bowtie2/${base}.hg19.sam -p 28 2>./02_bowtie2/${base}.align.log
/media/xionglab/data1/02_Software/samtools-1.11/samtools view -hbS -q 20 ./02_bowtie2/${base}.hg19.sam | samtools sort -T -> ./03_sort_uni_bam/${base}.hg19.unique.bam
java -jar /media/xionglab/data1/02_Software/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=./03_sort_uni_bam/${base}.hg19.unique.bam o=./04_rmdup_bam/${base}.hg19.rmdup.bam M=./04_rmdup_bam/${base}.hg19.txt
/media/xionglab/data1/02_Software/samtools-1.11/samtools index ./04_rmdup_bam/${base}.hg19.rmdup.bam
done

cd ./03_sort_uni_bam
for x in ./*.bam
do 
base=$(basename $x ".bam")
samtools view -c  ./${base}.bam  >./${base}.totolreads.log
done

cd ../04_rmdup_bam
for x in ./*.bam
do 
base=$(basename $x ".bam")
samtools view -c  ./${base}.bam  >./${base}.totolreads.log
done

cd ..

conda activate py39 
for i in ./04_rmdup_bam/*.hg19.rmdup.bam
do
base=$(basename $i ".hg19.rmdup.bam")
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i   -o   ./05_bw/${base}.10M.bw -p 13
#bamCoverage --scaleFactor  $res -b  $i  -e 300  --smoothLength 500 -o  ./05_bw/${base}.ext500.smo200.bw
done
