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
with open("lin_seq_samples") as f:
    for line in f:
        s = line.split("\t")
        sample = s[-1].split("/")[-1].split(".")[0]
        sample_to_leaf[sample] = s[0]

leaf_variants = 0
good_branch_variants = 0
bad_branch_variants = 0
germline_variants = 0

printnow = False

leaf_variants_vcf = open("/~/Downloads/strelka_vcfs/leaf_variants.vcf", "w")
germline_variants_vcf = open("~/Downloads/strelka_vcfs/germline_variants.vcf", "w")
good_branch_variants_vcf = open("~/Downloads/strelka_vcfs/good_branch_variants.vcf", "w")
bad_branch_variants_vcf = open("~/Downloads/strelka_vcfs/bad_branch_variants.vcf", "w")

leaf_variants_vcf = open("~/Downloads/hc_vcfs/leaf_variants_cnv_aware.vcf", "w")
germline_variants_vcf = open("~/Downloads/hc_vcfs/germline_variants_cnv_aware.vcf", "w")
good_branch_variants_vcf = open("~/Downloads/hc_vcfs/good_branch_variants_cnv_aware.vcf", "w")
bad_branch_variants_vcf = open("~/Downloads/hc_vcfs/bad_branch_variants_cnv_aware.vcf", "w")
all_non_germline_variants_vcf = open("~/Downloads/hc_vcfs/all_non_germline_variants_cnv_aware.vcf", "w")

bad_groupings = {}
good_groupings = {}

#with open("~/Downloads/strelka_vcfs/everything.vcf") as f:
with open("~/Downloads/hc_vcfs/everything_cnv_filter.vcf") as f:
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

            # if vcf_columns[0] == "chr20" and printnow:
            #     print(str(vcf_columns))

            if string in bad_groupings:
                bad_groupings[string] = bad_groupings[string] + 1
            else:
                bad_groupings[string] = 1

            bad_branch_variants += 1
            bad_branch_variants_vcf.write(line)
            all_non_germline_variants_vcf.write(line)

            break
        if current_node is None:
            # if vcf_columns[0] == "chr20" and printnow:
            #     print(str(vcf_columns))

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

# for key, value in sorted(good_groupings.items(), key=lambda p: (p[1], p[0])):
#     print("%s: %s" % (key, value))

