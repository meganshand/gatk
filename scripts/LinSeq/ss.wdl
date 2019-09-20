version 1.0

import "pipelines/dna_seq/germline/reprocessing/wgs/WholeGenomeReprocessing.wdl" as reprocess
import "structs/dna_seq/germline/GermlineStructs.wdl"


workflow reprocess_alot {

input {
  Array[File] bams
  String unmapped_bam_suffix

  File wgs_coverage_interval_list

  PapiSettings papi_settings
  GermlineSingleSampleReferences references
}
scatter (bam in bams) {
  String basename = basename(bam, ".bam")
  call reprocess.WholeGenomeReprocessing as WholeGenomeReprocessing {
    input:
      input_bam = bam,
      sample_name = basename,
      base_file_name = basename,
      final_gvcf_base_name = basename,
      unmapped_bam_suffix = unmapped_bam_suffix,
      references = references,
      papi_settings = papi_settings,
      wgs_coverage_interval_list = wgs_coverage_interval_list
  }
 }
}
