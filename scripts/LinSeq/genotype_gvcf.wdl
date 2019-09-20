version 1.0

import "tasks/GermlineVariantDiscovery.wdl" as Calling


workflow reprocess_alot_gvcfs {

input {
  Array[File] gvcfs
  File ref_dict
  File ref_fasta
  File ref_fasta_index
}

scatter (gvcf in gvcfs) {
  String basename = basename(gvcf, ".g.vcf.gz")
  call Calling.HaplotypeCaller_Genotype_GVCFs {
    input:
    input_gvcf = gvcf,
    vcf_basename = basename,
    ref_dict = ref_dict,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index
  }
 }
}
