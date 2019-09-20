version 1.0

workflow BamToGvcf {

  input {
    File calling_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Array[File] input_bams
    Array[File] input_bams
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp_vcf
    File dbsnp_vcf_index

  }

  String base_file_name = "blah"
  String final_gvcf_base_name = "lin_seq.vcf.gz"

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate GVCF by interval
    call HaplotypeCaller_GATK4_VCF_many_bams as HaplotypeCaller {
      input:
        contamination = 0.0007,
        input_bams = input_bams,
        interval_list = ScatterIntervalList.out[index],
        vcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = true,
        preemptible_tries = 0
     }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_vcf_name = final_gvcf_base_name + ".g.vcf.gz",
      preemptible_tries = 0
  }

  Float gvcf_size = size(MergeVCFs.output_vcf, "GB")

  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }
}


# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
    Int break_bands_at_multiples_of
  }

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    memory: "2 GB"
  }
}

task HaplotypeCaller_GATK4_VCF_many_bams {
  input {
    Array[String] input_bams
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Int preemptible_tries
    Int hc_scatter
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  parameter_meta {
    input_bams: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    gatk --java-options "-Xms11000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{sep=' -I ' input_bams} \
      -L ~{interval_list} \
      -O ~{vcf_basename}.vcf.gz \
      -contamination ~{default=0 contamination} \
      -bamout bamout.bam \
      -new-qual
    ls
  >>>
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "13 GB"
    cpu: "2"
    disks: "local-disk 100 HDD"
  }
  output {
    File output_vcf = "~{vcf_basename}.vcf.gz"
    File output_vcf_index = "~{vcf_basename}.vcf.gz.tbi"
    File bamout = "bamout.bam"
    File bamout_index = "bamout.bai"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcfs, "GB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}
