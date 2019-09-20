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

updated_seg = open("annotated_combined.seg", "w")

bad_groupings = {}
good_groupings = {}

sample_columns_in_seg_file = {}

def getCall(sample_set, leaf_variants, good_branch_variants, bad_branch_variants):
    if len(sample_set) == 1:
        leaf_variants += 1
        return "Leaf"
    if len(sample_set) == 0:
        return "NA"

    leaf_from_sample = sample_set[0]

    variant_leafs = set()
    for sample in sample_set:
        variant_leafs.add(sample)

    string = "".join(sorted(variant_leafs))

    leaf_to_start = None
    for leaf in leafs:
        if leaf_from_sample in leaf.values:
            leaf_to_start = leaf

    current_node = leaf_to_start
    while current_node is not None:
        if current_node.values == variant_leafs:
            if string in good_groupings:
                good_groupings[string] = good_groupings[string] + 1
            else:
                good_groupings[string] = 1
            good_branch_variants += 1
            return "GoodBranch"
        if current_node.values < variant_leafs:
            current_node = current_node.parent
            continue

        if string in bad_groupings:
            bad_groupings[string] = bad_groupings[string] + 1
        else:
            bad_groupings[string] = 1

        bad_branch_variants += 1
        return "BadBranch"

    if current_node is None:

        if string in bad_groupings:
            bad_groupings[string] = bad_groupings[string] + 1
        else:
            bad_groupings[string] = 1

        bad_branch_variants += 1
        return "BadBranch"


with open("combined.seg") as f:
    for line in f:
        if "@" in line:
            continue
        if 'CONTIG' in line:
            header = line.split("\t")
            new_line = line.rstrip() + "\tFollowsTreeAmp\tFollowsTreeDel\n"
            updated_seg.write(new_line)
            for i in range(3, 15):
                if i != 8:
                    sample_columns_in_seg_file[i] = header[i].split("Call_")[-1].rstrip()
            continue
        print(sample_columns_in_seg_file)

        vcf_columns = line.split("\t")
        sample_set_amp = []
        sample_set_del = []
        for i in range(3,15):
            if i != 8:
                if vcf_columns[i].rstrip() == '+':
                    sample_set_amp.append(sample_to_leaf[sample_columns_in_seg_file[i]])
                if vcf_columns[i].rstrip() == '-':
                    sample_set_del.append(sample_to_leaf[sample_columns_in_seg_file[i]])
        if len(sample_set_amp) == 11 and len(sample_set_del) == 0:
            germline_variants += 1
            new_line = line.rstrip() + "\tGermline\tNA\n"
            updated_seg.write(new_line)
            continue
        if len(sample_set_del) == 11 and len(sample_set_amp) == 0:
            germline_variants +=1
            new_line = line.rstrip() + "\tNA\tGermline\n"
            updated_seg.write(new_line)
            continue

        amp_call = getCall(sample_set_amp, leaf_variants, good_branch_variants, bad_branch_variants)
        del_call = getCall(sample_set_del, leaf_variants, good_branch_variants, bad_branch_variants)
        new_line = line.rstrip() + "\t" + amp_call + "\t" + del_call + "\n"
        updated_seg.write(new_line)

updated_seg.close()

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
