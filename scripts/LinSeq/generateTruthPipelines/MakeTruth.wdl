workflow MakeTruth {
    Array[String] sample_names = ["HT6W34", "HT6W45", "HT6W47", "HT6W48", "HT6W49", "HT6W54", "HT6W63", "HT115_42112", "HT115_42221", "HT115_5211", "HT115_7221"]
    File annotated_combined_seg_file = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/annotated_combined.seg"
    File joint_called_vcf = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/lin_seq_samples.vcf.gz"
    File joint_called_vcf_idx = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/lin_seq_samples.vcf.gz.tbi"
    File gatk4_jar_with_special_filters = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/gatk4_special_filters.jar"
    File gatk3_jar = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/gatk3.jar"
    File picard_jar = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/picard.jar"
    File lin_seq_samples = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/lin_seq_samples"
    File header = "gs://fc-5d8d0763-22af-45d3-ad83-5a28c7ba4d8d/inputs_for_generating_truth/header"
    File full_genome_intervals
    File ref_fasta
    File ref_fasta_fai
    File ref_dict

    File normal_bam
    File tumor_bam
    String tumor_sample_name = basename(tumor_bam, ".bam")
    String normal_sample_name = basename(normal_bam, ".bam")

    call getCNVIntervalList {
        input:
            annotated_combined_seg_file = annotated_combined_seg_file,
            tumor_sample_name = tumor_sample_name,
            ref_dict = ref_dict
    }

    call filterHCCalls {
        input:
            joint_called_vcf = joint_called_vcf,
            joint_called_vcf_idx = joint_called_vcf_idx,
            good_branch_cnvs_seg = getCNVIntervalList.good_branch_cnvs,
            gatk4_jar_with_special_filters = gatk4_jar_with_special_filters,
    }

    scatter(sample in sample_names) {
        call splitVcf {
            input:
                gatk4_jar_with_special_filters = gatk4_jar_with_special_filters,
                filtered_vcf = filterHCCalls.filtered_vcf,
                filtered_vcf_idx = filterHCCalls.filtered_vcf_idx,
                sample_name = sample
        }
    }

    call recombineVcfs {
        input:
            sample_vcfs = splitVcf.sample_vcf,
            sample_vcf_idxs = splitVcf.sample_vcf_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            gatk3_jar = gatk3_jar
    }

    call filterToTree {
        input:
            everything_vcf = recombineVcfs.everything_vcf,
            everything_vcf_idx = recombineVcfs.everything_vcf_idx,
            lin_seq_samples = lin_seq_samples,
            header = header,
            tumor_sample = tumor_sample_name,
            normal_sample = normal_sample_name
    }

    call getHiConfRegion {
        input:
        	picard_jar = picard_jar,
            tumor_leaf_vcf = filterToTree.tumor_leaf_variants,
            tumor_leaf_vcf_idx = filterToTree.tumor_leaf_variants_idx,
            normal_leaf_vcf = filterToTree.normal_leaf_variants,
            normal_leaf_vcf_idx = filterToTree.normal_leaf_variants_idx,
            full_genome_intervals = full_genome_intervals,
            bad_branch_cnvs = getCNVIntervalList.bad_branch_intervals,
            tumor_sample_name = tumor_sample_name,
            normal_sample_name = normal_sample_name
    }

    output {
        File truth_vcf = filterToTree.tumor_good_variants
        File truth_vcf_idx = filterToTree.tumor_good_variants_idx
        File high_conf_region = getHiConfRegion.intervals
        Float total_bases_in_hi_conf = getHiConfRegion.total_bases_in_hi_conf
    }

}

task getHiConfRegion {
	 File picard_jar
     File tumor_leaf_vcf
     File tumor_leaf_vcf_idx
     File normal_leaf_vcf
     File normal_leaf_vcf_idx
     File full_genome_intervals
     File bad_branch_cnvs
     String tumor_sample_name
     String normal_sample_name

     command {
        java -jar ${picard_jar} IntervalListTools ACTION=SUBTRACT I=${full_genome_intervals} \
        	SI=${tumor_leaf_vcf} \
            SI=${normal_leaf_vcf} \
            SI=${bad_branch_cnvs} \
            COUNT_OUTPUT=bases.count \
            OUTPUT_VALUE=BASES \
            O=${tumor_sample_name}_and_${normal_sample_name}.interval_list
     }
     output {
        File intervals = "${tumor_sample_name}_and_${normal_sample_name}.interval_list"
        Float total_bases_in_hi_conf = read_float("bases.count")
     }
     runtime{
        docker: "mshand/genomesinthecloud:2.3.3-1523386344"
     }

}

task filterToTree {
    File everything_vcf
    File everything_vcf_idx
    File lin_seq_samples
    File header
    String tumor_sample #in HT6W47 format
    String normal_sample #in HT6W47 format

    command <<<
python3 <<CODE
######## categorize variants into leaf, valid branch, invalid branch, root variants
class Node:

    def __init__(self, parent):
        self.parent = parent
        self.values = set()


def create_node(parent, values):
    node = Node(parent)
    node.values.update(values)
    return node


top_left = create_node(None, ["T54", "T47"])

lower_left_left = create_node(top_left, ["T54"])
lower_left_right = create_node(top_left, ["T47"])


top_right = create_node(None, ["T34", "T44", "T49", "T63", "T48", "T45", "T38", "T56", "T57"])

mid_right_left = create_node(top_right, ["T34", "T44", "T49", "T63"])
lower_right_left_left = create_node(mid_right_left, ["T34", "T44"])
lower_right_left_right = create_node(mid_right_left, ["T49", "T63"])


mid_right_right = create_node(top_right, ["T48", "T45", "T38", "T56", "T57"])
mid_low_right_right_left = create_node(mid_right_right, ["T48", "T45", "T38"])
lower_right_right_left_left = create_node(mid_low_right_right_left, ["T48", "T45"])
lower_right_right_left_right = create_node(mid_low_right_right_left, ["T38"])
mid_low_right_right_right = create_node(mid_right_right, ["T56", "T57"])
lower_right_right_right_left = create_node(mid_low_right_right_right, ["T56"])
lower_right_right_right_right = create_node(mid_low_right_right_right, ["T57"])

leafs = [lower_left_left, lower_left_right, lower_right_left_left, lower_right_left_right, lower_right_right_left_left, lower_right_right_left_right, lower_right_right_right_left, lower_right_right_right_right]

sample_to_leaf = {}
with open("${lin_seq_samples}") as f:
    for line in f:
        s = line.split("\t")
        sample = s[-1].split("/")[-1].split(".")[0]
        sample_to_leaf[sample] = s[0]

leaf_variants = 0
good_branch_variants = 0
bad_branch_variants = 0
germline_variants = 0

printnow = False

leaf_variants_vcf = open("leaf_variants_cnv_aware.vcf", "w")
germline_variants_vcf = open("germline_variants_cnv_aware.vcf", "w")
good_branch_variants_vcf = open("good_branch_variants_cnv_aware.vcf", "w")
bad_branch_variants_vcf = open("bad_branch_variants_cnv_aware.vcf", "w")
all_non_germline_variants_vcf = open("all_non_germline_variants_cnv_aware.vcf", "w")

bad_groupings = {}
good_groupings = {}

with open("${everything_vcf}") as f:
    for line in f:
        if "#" in line:
            continue
        vcf_columns = line.split("\t")
        info = vcf_columns[7]
        sample_set = None
        annotations = info.split(";")
        for an in annotations:
            if "set=" in an:
                sample_set = an.split("=")[1].split("-")
        if len(sample_set) == 1 and (sample_set[0] == "Intersection" or sample_set[0] == "FilteredInAll"):
            germline_variants += 1
            germline_variants_vcf.write(line)
            continue
        if len(sample_set) == 1:
            leaf_variants += 1
            leaf_variants_vcf.write(line)
            all_non_germline_variants_vcf.write(line)
            continue

        leaf_from_sample = sample_to_leaf[sample_set[0]]

        variant_leafs = set()
        for sample in sample_set:
            variant_leafs.add(sample_to_leaf[sample])

        string = "".join(sorted(variant_leafs))

        leaf_to_start = None
        for leaf in leafs:
            if leaf_from_sample in leaf.values:
                leaf_to_start = leaf

        current_node = leaf_to_start
        while current_node is not None:
            if current_node.values == variant_leafs:
                if vcf_columns[0] == "chr20" and printnow:
                    print(str(vcf_columns))
                good_branch_variants += 1
                good_branch_variants_vcf.write(line)
                all_non_germline_variants_vcf.write(line)

                if string in good_groupings:
                    good_groupings[string] = good_groupings[string] + 1
                else:
                    good_groupings[string] = 1
                break
            if current_node.values < variant_leafs:
                current_node = current_node.parent
                continue

            if string in bad_groupings:
                bad_groupings[string] = bad_groupings[string] + 1
            else:
                bad_groupings[string] = 1

            bad_branch_variants += 1
            bad_branch_variants_vcf.write(line)
            all_non_germline_variants_vcf.write(line)

            break
        if current_node is None:
            if string in bad_groupings:
                bad_groupings[string] = bad_groupings[string] + 1
            else:
                bad_groupings[string] = 1

            bad_branch_variants += 1
            bad_branch_variants_vcf.write(line)
            all_non_germline_variants_vcf.write(line)


leaf_variants_vcf.close()
germline_variants_vcf.close()
good_branch_variants_vcf.close()
bad_branch_variants_vcf.close()
all_non_germline_variants_vcf.close()

total = leaf_variants + germline_variants + good_branch_variants + bad_branch_variants
print(f"Total Variants: {total}")
float_total = total * 1.0
print(f"Leaf Variants: {leaf_variants}, {leaf_variants / float_total}")
print(f"Germline Variants: {germline_variants}, {germline_variants /float_total}")
print(f"Good Branch Variants: {good_branch_variants}, {good_branch_variants / float_total}")
print(f"Bad Branch Variants: {bad_branch_variants}, {bad_branch_variants / float_total}")

for key, value in sorted(bad_groupings.items(), key=lambda p: (p[1], p[0])):
    print("%s: %s" % (key, value))
CODE

        echo "Python code done"
        grep ${tumor_sample} leaf_variants_cnv_aware.vcf > ${tumor_sample}_leaf_variants.vcf.tmp
        cat ${header} ${tumor_sample}_leaf_variants.vcf.tmp > ${tumor_sample}_leaf_variants.vcf
        
        grep ${normal_sample} leaf_variants_cnv_aware.vcf > ${normal_sample}_leaf_variants.vcf.tmp
        cat ${header} ${normal_sample}_leaf_variants.vcf.tmp > ${normal_sample}_leaf_variants.vcf

        grep ${tumor_sample} good_branch_variants_cnv_aware.vcf | grep -v ${normal_sample} > ${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf.tmp
        cat ${header} ${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf.tmp > ${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf

        gatk IndexFeatureFile -F ${tumor_sample}_leaf_variants.vcf
        gatk IndexFeatureFile -F ${normal_sample}_leaf_variants.vcf
        gatk IndexFeatureFile -F ${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf
    >>>
    output {
        File leaf_variants = "leaf_variants_cnv_aware.vcf"
        File good_branch_variants = "germline_variants_cnv_aware.vcf"
        File bad_branch_variants = "bad_branch_variants_cnv_aware.vcf"
        File all_non_germline_variants = "all_non_germline_variants_cnv_aware.vcf"
        File tumor_leaf_variants = "${tumor_sample}_leaf_variants.vcf"
        File tumor_leaf_variants_idx = "${tumor_sample}_leaf_variants.vcf.idx"
        File normal_leaf_variants = "${normal_sample}_leaf_variants.vcf"
        File normal_leaf_variants_idx = "${normal_sample}_leaf_variants.vcf.idx"
        File tumor_good_variants = "${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf"
        File tumor_good_variants_idx = "${tumor_sample}_no_${normal_sample}_good_branch_variants.vcf.idx"
    }
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.1.2.0"
    }
}

task filterHCCalls {
    File joint_called_vcf
    File joint_called_vcf_idx
    File good_branch_cnvs_seg
    File gatk4_jar_with_special_filters

    command <<<
        java -jar ${gatk4_jar_with_special_filters} ExampleVariantWalker -V ${joint_called_vcf} -O with_cnv_filters.vcf --cnv ${good_branch_cnvs_seg}
    >>>
    output {
        File filtered_vcf = "with_cnv_filters.vcf"
        File filtered_vcf_idx = "with_cnv_filters.vcf.idx"
    }
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.1.2.0"
    }
}

task splitVcf {
    File gatk4_jar_with_special_filters
    File filtered_vcf
    File filtered_vcf_idx
    String sample_name

    command <<<
        java -jar ${gatk4_jar_with_special_filters} SelectVariants -V ${filtered_vcf} --exclude-filtered false --exclude-non-variants true -O ${sample_name}.vcf -sn ${sample_name}
    >>>
    output {
        File sample_vcf = "${sample_name}.vcf"
        File sample_vcf_idx = "${sample_name}.vcf.idx"
    }
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.1.2.0"
    }
}

task recombineVcfs {
    Array[File] sample_vcfs
    Array[File] sample_vcf_idxs
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File gatk3_jar

    command <<<
    ln -s ${sep=' .\nln -s ' sample_vcfs} .

    ln -s ${sep=' .\nln -s ' sample_vcf_idxs} .

python3 <<CODE
import subprocess

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

inputs = ""
priority = " -priority "
for sample in samples:
    vcf = "%s.vcf" % sample
    inputs += " -V:" + sample + " " + vcf
    priority += sample + ","
command = "java -Xmx7G -jar ${gatk3_jar} -T CombineVariants -R ${ref_fasta} -genotypeMergeOptions PRIORITIZE " + inputs + priority[:-1] + " -o everything_cnv_filter.vcf --disable_auto_index_creation_and_locking_when_reading_rods"
run_command(command)
CODE
    >>>
    output {
        File everything_vcf = "everything_cnv_filter.vcf"
        File everything_vcf_idx = "everything_cnv_filter.vcf.idx"
    }
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.1.2.0"
        disks: "local-disk 100 HDD"
        memory: "8 GB"
    }
}

task getCNVIntervalList {
    File annotated_combined_seg_file
    String tumor_sample_name #in HT6W47 format
    File ref_dict

    command <<<
        Rscript --no-save -<<'RCODE'
            library(plyr)

            combined = read.table("${annotated_combined_seg_file}", header=T, comment.char = "@", stringsAsFactors = F)

            bad_branches = subset(combined,
                                  ((FollowsTreeAmp == "BadBranch" | FollowsTreeDel == "BadBranch") & (Call_${tumor_sample_name} == "+" | Call_${tumor_sample_name}== "-")) |
                                    (FollowsTreeAmp == "Leaf" & Call_${tumor_sample_name}=="+") | (FollowsTreeDel == "Leaf" & Call_${tumor_sample_name}== "-"))
            write.table(data.frame(CHROM=bad_branches[,'CONTIG'], START=bad_branches[,'START'], END=bad_branches[,'END'], STRAND="+", name="."),
                        "bad_branch.interval_list.tmp", col.names = F, row.names = F, quote = F, sep="\t")

            #Only need amplifications because deletions will only cause HomVars which already won't be filtered.
            good_branches = subset(combined,
                                   ((FollowsTreeAmp == "GoodBranch" | FollowsTreeAmp == "Germline") & Call_${tumor_sample_name} == "+"))
            #round to nearest .5. Segment_Mean is in log_2
            good_branches[,'COPY_RATIO'] = round((2^(good_branches[,'Segment_Mean_${tumor_sample_name}']))*2)/2
            write.table(data.frame(CONTIG=good_branches[,'CONTIG'], START=good_branches[,'START'], END=good_branches[,'END'], COPY_RATIO=good_branches[,'COPY_RATIO']),
                        "good_branch_cnvs.seg", row.names = F, quote = F, sep="\t")

        RCODE

        cat ${ref_dict} bad_branch.interval_list.tmp > bad_branch.interval_list
    >>>
    output {
        File bad_branch_intervals = "bad_branch.interval_list"
        File good_branch_cnvs = "good_branch_cnvs.seg"
    }
    runtime {
        docker: "skwalker/rscripting:rwithdatatable"
        disks: "local-disk 50 HDD"
    }

}