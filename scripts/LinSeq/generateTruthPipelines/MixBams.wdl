workflow MixBams {
  File tumor_sample
  File normal_sample
  Int spike_pct_as_int
  String cloud = "true"
  String sample_name
  
  File ref_fasta
  File ref_fasta_index
  File wgs_coverage_interval_list

  Float spike_pct = spike_pct_as_int/100.0
  call SpikeInBam {
    input:
      spikein_bam = tumor_sample,
      normal_bam = normal_sample,
      spike_pct = spike_pct,
      cloud = cloud
  }
  
  call CollectWgsMetrics {
  	input:
      input_bam = SpikeInBam.spiked_bam,
      input_bam_index = SpikeInBam.spiked_bam_index,
  	  metrics_filename = sample_name + ".txt",
  	  wgs_coverage_interval_list = wgs_coverage_interval_list,
  	  ref_fasta = ref_fasta,
  	  ref_fasta_index = ref_fasta_index
  }

  output {
    File mixed_bam = SpikeInBam.spiked_bam
    File mixed_bai = SpikeInBam.spiked_bam_index
    File wgs_metrics = CollectWgsMetrics.metrics
    Float mixed_af = SpikeInBam.allele_fraction
  }
}

task SpikeInBam {
   File spikein_bam
   File normal_bam
   Float spike_pct
   Float desired_cov = 1515607380
   String cloud # TODO: move it to the workflow level

   command <<<
      set -e

      samtools view -c ${spikein_bam} > spike.count
      read SPIKE_COUNT < spike.count
      echo "Number of reads in spikein bam: $SPIKE_COUNT"

      samtools view -c ${normal_bam} > normal.count
      read NORMAL_COUNT < normal.count
      echo "Number of reads in normal bam: $NORMAL_COUNT"

	  awk -v cov=${desired_cov} -v pct=${spike_pct} 'BEGIN {print (cov * pct)}' > TTCOV
      read TOTAL_TUMOR_COV < TTCOV
      
      awk -v cov=${desired_cov} -v pct=${spike_pct} 'BEGIN {print (cov * (1 - pct))}' > TNCOV
      read TOTAL_NORMAL_COV < TNCOV

      awk -v tcov=$TOTAL_TUMOR_COV -v sc=$SPIKE_COUNT 'BEGIN {print tcov / sc }' > SDS
      read SPIKE_DOWNSAMPLE_TMP < SDS
      
      awk -v a=$SPIKE_DOWNSAMPLE_TMP -v b=1 'BEGIN{print (a>b)?b:a}' > min_s
      read SPIKE_DOWNSAMPLE < min_s

      echo "Downsampling spike in bam by " $SPIKE_DOWNSAMPLE
      samtools view -s $SPIKE_DOWNSAMPLE ${spikein_bam} -o spike_downsampled.bam

      awk -v nc=$NORMAL_COUNT -v tcov=$TOTAL_NORMAL_COV 'BEGIN {print tcov / nc}' > NDS
      read NORMAL_DOWNSAMPLE_TMP < NDS
      
      awk -v a=$NORMAL_DOWNSAMPLE_TMP -v b=1 'BEGIN{print (a>b)?b:a}' > min_n
      read NORMAL_DOWNSAMPLE < min_n

      echo "Downsampling normal bam by " $NORMAL_DOWNSAMPLE
      samtools view ${normal_bam} -s $NORMAL_DOWNSAMPLE -o normal_downsampled.bam

	  awk -v tc=$SPIKE_COUNT -v nc=$NORMAL_COUNT -v tds=$SPIKE_DOWNSAMPLE -v nds=$NORMAL_DOWNSAMPLE 'BEGIN{print (tds * tc) / ((tds * tc) + (nds * nc))}' > af.txt
      read AF < af.txt
      echo "Allele fraction should be around " $AF
        
      echo "Merging downsampled bam into final spiked bam"
      samtools merge -f contaminated.bam normal_downsampled.bam spike_downsampled.bam

      # prem: /seq/software/picard/current/bin/picard-private.jar
      # cloud: /usr/gitc/picard.jar

      picard_jar="/usr/gitc/picard.jar"
      if [ ${cloud} != "true" ]; then
        picard_jar="/seq/software/picard/current/bin/picard-private.jar"
      fi

      java -jar $picard_jar AddOrReplaceReadGroups \
        I=contaminated.bam \
        O=rg_fixed.contaminated.bam \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=mixture \
        CREATE_INDEX=true

   >>>

   runtime {
      memory: "5 GB"
      disks: "local-disk " + 500 + " HDD"
      docker: "mshand/genomesinthecloud:2.3.3-1523386344"
   }
   output {
      File spiked_bam = "rg_fixed.contaminated.bam"
      File spiked_bam_index = "rg_fixed.contaminated.bai"
      Float allele_fraction = read_float("af.txt")
   }
}

task CollectWgsMetrics {
  File input_bam
  File input_bam_index
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int read_length = 250

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=${wgs_coverage_interval_list} \
      OUTPUT=${metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=${read_length}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}
