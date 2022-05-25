#!/usr/bin/env bash
# Author: Junhao Chen
# Date: 2022-05-05
# Title: Stripe-seq pipeline for wheat, running on SLUHPC
# Version: 0.1.0

#########config
# The input file address
sample_dir=/data/projects/junhao.chen/02wheat/01rawData/all_fastq
# The suffix of the samples
suffix=.R[12].fastq.gz
# work dir
wkdir=${PWD}
ref_dir=/data/projects/junhao.chen/02wheat/run_all_fastq_by_pipeline/00ref
ref_name=wheat
max_intronlen=42000
# Threads of fastqc
fastqc_threads=38
fastp_threads=16
hisat2_threads=40
samtools_threads=40
##########Dont edit it if you dont know the fuction it is below this line!
sample_list=sample_list 

# constom config:
module load cuda10.2/toolkit/10.2.89
######### generate the config file
if [ ! -f ${sample_list} ]
then
  ls ${sample_dir} | sed "s/${suffix}//" | sort | uniq > ${sample_list} \
  && echo "${sample_list} generated!" \
  && sample_number=$(cat ${sample_list}| wc -l ) \
  && echo -e "Total sample number is: " ${sample_number}
else
  if [ -s ${sample_list} ] 
  then
    echo "The sample_list file exist!"\
    && sample_number=$(cat ${sample_list}| wc -l ) \
    && echo -e "Total sample number is: " ${sample_number}
  else
    echo "The sample_list file exist but it seems an empty file! Please check it!"
  fi
fi

# generate the qc report
if [ ! -f rawData_fastqc_ok ]
then
  if [ ! -d ${wkdir}/01QC ]
  then
    mkdir -p ${wkdir}/01QC
    if [ ! -d ${wkdir}/01QC/rawData ]
    then
      mkdir -p ${wkdir}/01QC/rawData 
    else
      echo "${wkdir}/01QC/rawData alreay exist!"
    fi 
  else
    echo "${wkdir}/01QC alreay exist!"
    if [ ! -d ${wkdir}/01QC/rawData ]
    then
      mkdir -p ${wkdir}/01QC/rawData 
    else
      echo "${wkdir}/01QC/rawData alreay exist!"
    fi 
  fi \
  && echo "=====Start fastqc!=====" \
  && fastqc  -t ${fastqc_threads} -o ${wkdir}/01QC/rawData ${sample_dir}/* && touch rawData_fastqc_ok
else
  echo "fastqc: The rawData_fastqc_ok already exist! skip fastqc step!"
fi

# multiqc
if [ ! -f rawData_multiqc_ok ]
then
  if [ -f rawData_fastqc_ok ]
  then
    if [ -d ${wkdir}/01QC/rawData ]
    then
      echo "=====Start rawData multiqc!=====" \
      && multiqc -f ${wkdir}/01QC/rawData -o ${wkdir}/01QC/rawData \
      && touch rawData_multiqc_ok
      # generate the qc list
        if [ -f ${wkdir}/01QC/rawData/multiqc_data/multiqc_fastqc.txt ]
        then
          cat ${wkdir}/01QC/rawData/multiqc_data/multiqc_fastqc.txt | cut -f 1,5 | sed '1d' | tr '\t' ',' | sort -V | uniq > ${PWD}/rawData.qc.txt \
          && echo "The rawData.qc.txt already generate!" 
        else
          echo "${wkdir}/01QC/rawData/multiqc_data/multiqc_fastqc.txt file may not exist! Please check if the multiqc run correctly"
        fi
    else
      echo "multiqc: ${wkdir}/01QC/rawData may not exist! "
    fi
  else
    echo "multiqc: the file rawData_fastqc_ok does not exist! Please check the fastqc step!"
  fi
else
  echo "multiqc: The rawData_multiqc_ok already exist! skip multiqc step!"
fi

# fastp
if [ ! -f fastp_ok ]
then
  if [ -f rawData_multiqc_ok ]
  then
    if [ ! -d ${wkdir}/02fastp ]
    then
      mkdir -p ${wkdir}/02fastp
    else
      echo "fastp: ${wkdir}/02fastp already exist!"
    fi \
    && echo "=====Start run fastp!=====" \
    && for i in `cat sample_list`
       do
         if [ -f "${wkdir}/02fastp/${i}.R1.dedup.fastq.gz" -a -f "${wkdir}/02fastp/${i}.R2.dedup.fastq.gz" ]
         then
           echo "The ${wkdir}/02fastp/${i}.R1.dedup.fastq.gz file already exist! Skip the ${i}!"
        else
           echo -n "${i} start at:" \
         && date \
         && fastp \
            -i ${sample_dir}/${i}.R1.fastq.gz -o ${wkdir}/02fastp/${i}.R1.dedup.fastq.gz \
            -I ${sample_dir}/${i}.R2.fastq.gz -O ${wkdir}/02fastp/${i}.R2.dedup.fastq.gz \
            --umi --umi_loc read1 \
            --umi_len 15 \
            --trim_poly_g \
            --json ${wkdir}/02fastp/${i}.json \
            --html ${wkdir}/02fastp/${i}.html \
            --thread ${fastp_threads} \
            --report_title "${i}" \
            --dedup --dup_calc_accuracy 5 \
          && echo "${i} fastp done!"
        fi 
       done \
       && touch fastp_ok
  else
    echo "fastp: the file rawData_multiqc_ok does not exist! Please check the multiqc step!"
  fi
else
  echo "fastp: The fastp_ok file already exist! skip fastp step!"
fi

# generate the duplicate list and 
if [ -f fastp_ok ]
then
  if [ ! -f fastp_stat_ok ]
  then
    rm -rf duplicate.rate.txt \
    && rm -rf cleanData.qc.txt \
    && for i in `cat sample_list`
      do 
        echo -n "${i}" >> duplicate.rate.txt \
        && cat ${PWD}/02fastp/${i}.json | grep -w "rate" | cut -d ':' -f 2 >> duplicate.rate.txt \
        && echo -n "${i}" >> cleanData.qc.txt \
        && cat ${PWD}/02fastp/${i}.json | grep -w "read1_after_filtering" -A 1 | tail -1 | cut -d ':' -f 2 | tr -d ','  >> cleanData.qc.txt 
 
      done \
        && cat duplicate.rate.txt | sort -V | sed 'p'  > duplicate.rate.txt.tmp \
        && mv duplicate.rate.txt.tmp duplicate.rate.txt \
        && cat cleanData.qc.txt | sort -V | sed 'p'  > cleanData.qc.txt.tmp \
        && mv cleanData.qc.txt.tmp cleanData.qc.txt \
        && touch fastp_stat_ok
    else
      echo "fastp_stat: The fastp_stat_ok file already exist! skip fastp step!"
    fi
else
 echo "It seems fastp_ok file does not exist. Please check fastp stap!"
fi

# fastp_qc 
# 感觉不需要这个，可以直接从fastp里获得数据。
if [ ! -f fastp_qc_ok ]
then
  if [ ! -d ${PWD}/01QC/cleanData ]
  then
    mkdir -p ${PWD}/01QC/cleanData \
    && fastqc -t ${fastqc_threads} -o ${PWD}/01QC/cleanData \
    && multiqc -f ${PWD}/01QC/cleanData -o ${PWD}/01QC/cleanData \
    && cat ${wkdir}/01QC/rawData/multiqc_data/multiqc_fastqc.txt | cut -f 1,5 > ${PWD}/rawData.qc.txt \
    && echo "The rawData.qc.txt already generate!"
  else
    echo "The cleanData folder does not exist!"
  fi \
  && touch fastp_qc_ok
else
  echo "fastp_qc: The fastp_qc_ok file already exist! skip fastp_qc step!"
fi



# ribodetector
# run ribodetector.sh
if [ ! -f ribodetector_ok ]
then
sbatch ribodetector.sh
else
  echo "ribodetector: The ribodetector_ok file already exist! skip ribodetector step!"
fi

# hisat2
if [ ! -d ${PWD}/04hisat2 ]
then
  mkdir -p ${PWD}/04hisat2 \
  && echo "The 04hisat2 folder created!"
else
  echo "The 04hisat2 folder already exist!"
fi

if [ ! -f hisat2_ok ]
then
  echo "Start run hisat2"
  for i in `cat ${sample_list}`
  do
    echo -n "start mapping ${i} with hisat2 at:" \
    && date \
    && hisat2 -p ${hisat2_threads} \
    --new-summary --summary-file ${PWD}/04hisat2/${i}.log \
    -x ${ref_dir}/${ref_name} \
    --no-softclip --max-intronlen ${max_intronlen} \
     -1 ${PWD}/03ribodetector/ribo/${i}.1.ribo.fq.gz -2 ${PWD}/03ribodetector/ribo/${i}.2.ribo.fq.gz | samtools sort -@ ${samtools_threads} -o ${PWD}/04hisat2/${i}.sort.bam - \
    && echo "finished ${i} at:" \
    && date
  done \
  && touch hisat2_ok \
  && echo "Hisat2 mapping finished!"
else
  echo "hisat2: The hisat2_ok file already exist! skip hisat2 step!"
fi

# bam-level dedup
#if [ ! -d ${PWD}/05bamLevelDedup ]
#then
#  mkdir ${PWD}/05bamLevelDedup \
#  && echo "The folder ${PWD}/05bamLevelDedup created!"
#else
#  echo "The folder ${PWD}/05bamLevelDedup created! Do the next step..."
#fi
#
#if [ ! -f picard_ok ]
#then
#  echo "Start run the picard bam-level dedup" \
#  && for i in `cat ${PWD}/sample_list`
#  do 
#    picard MarkDuplicates \
#    I=${PWD}/04hisat2/${i}.sort.bam  O=${PWD}/05bamLevelDedup/${i}.sort.dedup.bam  REMOVE_DUPLICATES=true CREATE_INDEX=ture M=${PWD}/05bamLevelDedup/${i}_marked_dup_metrics.txt 
#  done
#else
#  echo "Picard: The picard_ok file already exist! skip picard step!"
#fi

