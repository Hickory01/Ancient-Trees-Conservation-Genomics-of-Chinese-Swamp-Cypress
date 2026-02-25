#parameter#
id=$1
min=$2
max=$3

#input#
refseq="/home/ZhangWP/water_pine/Gpen_ref_index/Gpen_refgenome_100k.fa"
#bampath="/home/ZhangWP/water_pine/mappingout/fasta_bam/bam_psmc"
bampath="/home/ZhangWP/water_pine/demography/0_psmc/psmc54_bam"
sampath="/home/ZhangWP/water_pine/consensus/samtools_output_Gpen"

#output#
#psmcpath="/home/ZhangWP/water_pine/demography/0_psmc/Gpen_psmc2023"
psmcpath="/home/ZhangWP/water_pine/demography/0_psmc/Gpen_psmc2025"
mkdir -p $psmcpath

#workdir#
scripts="/home/ZhangWP/water_pine/consensus/scripts"
psmc="/home/ZhangWP/software/psmc-master"

#gvcf
bcftools mpileup -C 50 -q 20 -Q 20 -f $refseq $bampath/$id.bam | bcftools call --ploidy 2 -c -V indels -o $sampath/$id.samtool.gvcf 

#VCF for PSMC
vcfutils.pl vcf2fq -d 10 -D 100 $sampath/$id.samtool.gvcf | gzip > $psmcpath/$id.fq.gz
#vcf2fq: VCF->fastq

#Convert fq to psmcfa
$psmc/utils/fq2psmcfa -q20 $psmcpath/$id.fq.gz >$psmcpath/$id.psmcfa

#cd $psmcpath
$psmc/psmc -N30 -t15 -r5 -p "4+25*2+4+6" -o $id.psmc $id.psmcfa

