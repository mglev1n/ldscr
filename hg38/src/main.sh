#!/bin/bash

fx () {

  POP=$1

  ## bedfile for LiftOver

  zcat ../../inst/extdata/${POP}/UKBB.${POP}.l2.ldscore.gz.norsid \
    | sed '1d' \
    | awk -v OFS="\t" '{split($2,arr,":");
      print arr[1],arr[2]-1,arr[2],arr[3],arr[4],$2
    }' > ../tmp/UKBB.${POP}.l2.ldscore.hg19.bed

  bedtools getfasta -tab -fi ~/publicdata/reference/Homo_sapiens_assembly19.fasta \
    -bed ../tmp/UKBB.${POP}.l2.ldscore.hg19.bed > ../tmp/UKBB.${POP}.l2.ldscore.hg19.ref

  ## Confirm hg19 All ref/alt match

  paste ../tmp/UKBB.${POP}.l2.ldscore.hg19.bed ../tmp/UKBB.${POP}.l2.ldscore.hg19.ref \
    | awk '$4!=$8'

  ## Make a VCF for LiftOver

  cat ../tmp/UKBB.${POP}.l2.ldscore.hg19.bed \
    | awk -v OFS="\t" '{print "chr"$1,$3,$6,$4,$5,".",".","."}'  \
    | sed '1s/^/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n/' \
    | sed '1s/^/##fileformat=VCFv4.2\n/' \
    | gzip -c > ../tmp/UKBB.${POP}.l2.ldscore.hg19.vcf.gz

  java \
    -Djava.io.tmpdir=/broad/hptmp/skoyama \
    -jar ~/repo/picard/picard.jar LiftoverVcf \
    I=../tmp/UKBB.${POP}.l2.ldscore.hg19.vcf.gz \
    O=../tmp/UKBB.${POP}.l2.ldscore.hg38.vcf.gz \
    CHAIN=~/repo/liftOver/hg19ToHg38.over.chain.gz \
    REJECT=../tmp/UKBB.${POP}.l2.ldscore.hg19.rejected.vcf.gz \
    RECOVER_SWAPPED_REF_ALT=true \
    R=~/publicdata/reference/Homo_sapiens_assembly38.fasta

  zcat ../tmp/UKBB.${POP}.l2.ldscore.hg38.vcf.gz \
    | grep -v "#" \
    | awk -v OFS="\t" '{print $3,$1,$1":"$2":"$4":"$5,$2}' \
    | sed 's/chr//' > ../tmp/UKBB.${POP}.l2.ldscore.hg19tohg38.linker

  join -t $'\t' \
    <(cat ../tmp/UKBB.${POP}.l2.ldscore.hg19tohg38.linker | sort -k 1,1) \
    <(zcat ../../inst/extdata/${POP}/UKBB.${POP}.l2.ldscore.gz.norsid \
    | sed '1d' | awk -v OFS="\t" '{print $2,$4}' | sort -k 1,1) \
    | cut -f 2- \
    | sort -V -S 12G -T . \
    | sed '1s/^/CHR\tSNP\tBP\tL2\n/' \
    | gzip -c > ../../inst/extdata/${POP}/UKBB.${POP}.l2.ldscore.gz.norsid_hg38

  rm ../tmp/UKBB.${POP}.l2.ldscore.{hg19,hg38}.vcf.gz
  rm ../tmp/UKBB.${POP}.l2.ldscore.hg38.vcf.gz.tbi
  rm ../tmp/UKBB.${POP}.l2.ldscore.hg19.{ref,bed}
  rm ../tmp/UKBB.${POP}.l2.ldscore.hg19.rejected.vcf.gz
  rm ../tmp/UKBB.${POP}.l2.ldscore.hg19tohg38.linker

}

fx AFR
fx AMR
fx EUR
fx EAS
fx CSA
fx MID

