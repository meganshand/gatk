import subprocess
import glob

# function that runs each command
def run_command(command):
    print(command)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    if process.returncode != 0:
        print("Something went wrong.")
        exit(process.returncode)

commands = []

# list of samples you want to iterate over
samples = ["HT6W34", "HT6W45", "HT6W47", "HT6W48", "HT6W49", "HT6W54", "HT6W63", "HT115_42112", "HT115_42221", "HT115_5211", "HT115_7221"]

# path to the joint called vcf
jc_vcf = "~/Downloads/joint_called_vcf/hc_vcf_no_genotypegvcf/lin_seq_samples.vcf.gz"
# path to output the filtered joint called vcf
filtered_jc_vcf_path = "~/Downloads/joint_called_vcf/hc_vcf_no_genotypegvcf/joint_called_with_cnv_filters.vcf"
# path to cnv seg file with copy ratio
cnv_seg_copy_ratio_path = "good_branch_cnvs.seg"

# directory where you want the individual hc vcfs to go when selected out of the joint called vcf
hc_vcf_dir = "~/Downloads/hc_vcfs/"

# directory where all strelka vcfs have been downloaded
# all strelka vcfs for a given sample are expected to be found in the structure {strelka_base_dir}/{sample}/
strelka_base_dir = "~/Downloads/strelka_vcfs/"

# path to gatk4 jar with filters in both LinSeqFilter and CombineGVCFs
gatk_4_with_filters_path = "../../build/libs/gatk.jar"

# any gatk3 jar, needs the combinevariants tool
gatk_3_jar = "~/Downloads/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"


# filter joint called vcf by ad/gq/no calls
command = "java -jar " + gatk_4_with_filters_path + " LinSeqFilter -V " + jc_vcf + " -O " + filtered_jc_vcf_path + " --cnv " + cnv_seg_copy_ratio_path
#run_command(command)

for sample in samples:

     path = strelka_base_dir + sample + "/"

     # combine all strelka indel and snv files for sample with filters - strelka
     snv_vcfs = glob.glob(path + "*.snvs.vcf.gz")
     indel_vcfs = glob.glob(path + "*.indels.vcf.gz")
     java_command = "java -jar " + gatk_4_with_filters_path + " CombineGVCFs -R ~/Downloads/Homo_sapiens_assembly38.fasta"
     priority = " -priority "
     inputs = ""
     for file in snv_vcfs:
         if "_HT115." in file or "_HT6W37." in file:
             continue
         ROD = file.split("/")[-1].split(".")[0]
         priority += ROD + ","
         inputs += " -V " + file
     priority = priority[:-1]
     output = " -O %s.snvs.combined.vcf" % (path + sample)
     snv_combine_command = java_command + inputs
     snv_combine_command += output
     run_command(snv_combine_command)

     inputs = ""
     for file in indel_vcfs:
         if "_HT115." in file or "_HT6W37" in file:
             continue
         ROD = file.split("/")[-1].split(".")[0]
         inputs += " -V " + file
     output = " -O %s.indels.combined.vcf" % (path + sample)
     indel_combine_command = java_command + inputs
     indel_combine_command += output
     run_command(indel_combine_command)


     # combine indel and snv strelka vcfs - strelka
     combined_vcfs = glob.glob(path + "*.combined.vcf")
     java_command = "java -jar " + gatk_3_jar + " -T CombineVariants -R ~/Downloads/Homo_sapiens_assembly38.fasta -genotypeMergeOptions PRIORITIZE"
     inputs = ""
     priority = " -priority snvs,indels"
     for file in combined_vcfs:
         if "indels" in file:
             inputs += " -V:indels " + file
         elif "snvs" in file:
             inputs += " -V:snvs " + file
         else:
             continue

     output = " -o %s.all.combined.vcf" % (path + sample)
     run_command(java_command + inputs + output + priority)

     # removed filtered calls from combined vcf - strelka
     input_vcf = " %s.all.combined.vcf" % (path + sample)
     output_vcf = "%s.all.filtered.combined.vcf" % (path + sample)
     filter_combined_strelka_command = "java -jar " + gatk_4_with_filters_path + " SelectVariants -V " + input_vcf + " --exclude-filtered true -O " + output_vcf
     run_command(filter_combined_strelka_command)

    # select variants for each sample from joint called vcf - haplotypecaller
    output_vcf = hc_vcf_dir + sample + ".vcf.gz"
    command = "java -jar " + gatk_4_with_filters_path + " SelectVariants -V " + filtered_jc_vcf_path + " --exclude-filtered false --exclude-non-variants true -O " + output_vcf + " -sn " + sample
    run_command(command)


# combine all combined strelka vcfs across samples - strelka
 inputs = ""
 priority = " -priority "
 for sample in samples:
     path = "~/Downloads/strelka_vcfs/%s/" % sample
     vcf = path + sample + ".all.filtered.combined.vcf"
     inputs += " -V:" + sample + " " + vcf
     priority += sample + ","
 command = "java -jar " + gatk_3_jar + " -T CombineVariants -R ~/Downloads/Homo_sapiens_assembly38.fasta -genotypeMergeOptions PRIORITIZE " \
           + inputs + priority[:-1] + " -o ~/Downloads/strelka_vcfs/everything.vcf"
 run_command(command)

# combine variants across all hc vcfs - haplotypecaller
inputs = ""
priority = " -priority "
for sample in samples:
    vcf = "~/Downloads/hc_vcfs/%s.vcf.gz" % sample
    inputs += " -V:" + sample + " " + vcf
    priority += sample + ","
command = "java -jar " + gatk_3_jar + " -T CombineVariants -R ~/Downloads/Homo_sapiens_assembly38.fasta -genotypeMergeOptions PRIORITIZE " \
          + inputs + priority[:-1] + " -o ~/Downloads/hc_vcfs/everything_cnv_filter.vcf"
run_command(command)

