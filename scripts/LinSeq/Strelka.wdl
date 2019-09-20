workflow RunAllPossibleTumorNormalPairs {
  Array[File] bams
  Array[File] bam_indices
  File ref_fasta = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File ref_fasta_index = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
  File ref_dict = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File calling_bed = "gs://broad-dsp-spec-ops/scratch/jsoto/strelka/wgs_calling_regions.hg38.bed.gz"
  File calling_bed_index = "gs://broad-dsp-spec-ops/scratch/jsoto/strelka/wgs_calling_regions.hg38.bed.gz.tbi"
  Array[Pair[File, File]] tumor_normal_pairs = cross(bams, bams)
  scatter(pair in tumor_normal_pairs) {
    if(pair.left != pair.right) {
      String basename = basename(pair.left, ".bam") + "_" + basename(pair.right, ".bam")
      call Strelka {
        input:
          tumor_bam = pair.left,
          normal_bam = pair.right,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          calling_bed = calling_bed,
          bam_indices = bam_indices,
          basename = basename,
          calling_bed_index = calling_bed_index
      }
    }
  }
}

task Strelka {
  File tumor_bam
  File normal_bam
  Array[File] bam_indices
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File calling_bed
  String basename
  File calling_bed_index
  Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB"))  * 2) + 10

  command {
    set -e

    /strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam ${normal_bam} \
        --tumorBam ${tumor_bam} \
        --ref ${ref_fasta} \
        --runDir workspace \
        --callRegions ${calling_bed}

    workspace/runWorkflow.py -m local -j 4

    cp workspace/results/variants/somatic.snvs.vcf.gz ${basename}.snvs.vcf.gz
    cp workspace/results/variants/somatic.snvs.vcf.gz.tbi ${basename}.snvs.vcf.gz.tbi
    cp workspace/results/variants/somatic.indels.vcf.gz ${basename}.indels.vcf.gz
    cp workspace/results/variants/somatic.indels.vcf.gz.tbi ${basename}.indels.vcf.gz.tbi
  }

  runtime {
    docker:"jsotobroad/strelka:latest"
    cpu: 4
    memory: "11 GB"
    preemptible: 0
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File snp_vcf = "${basename}.snvs.vcf.gz"
    File snp_vcf_index = "${basename}.snvs.vcf.gz.tbi"
    File indel_vcf = "${basename}.indels.vcf.gz"
    File indel_vcf_index = "${basename}.indels.vcf.gz.tbi"
  }
}


