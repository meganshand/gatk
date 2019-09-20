library(plyr)
library(ggplot2)
seg.files = dir(".", "seg$", full.names = T)
segments = ldply(seg.files, function(x){read.table(x,header=T,stringsAsFactors = F, sep="\t")})
segments$Length = segments$End - segments$Start


combined = read.table("annotated_combined.seg", header=T, comment.char = "@", stringsAsFactors = F)
cleaned_combined = data.frame(matrix(ncol = 29, nrow = 0))
colnames(cleaned_combined) <- colnames(combined)
current_end = combined[1,]$END
current_start = combined[1,]$START

for (i in 1:(length(combined[,'CONTIG']) - 1)) {
  row = combined[i,]
  next_row = combined[i+1,]
  if (sum(row[,4:15] != next_row[,4:15]) == 0 & row$END + 1 == next_row$START & row$CONTIG == next_row$CONTIG) {
    current_end = next_row$END
  } else {
    new_row = row
    new_row$END = current_end
    new_row$START = current_start
    current_end = next_row$END
    current_start = next_row$START
    cleaned_combined = rbind(cleaned_combined, new_row)
  }
}

cleaned_combined = rbind(cleaned_combined, combined[length(combined$CONTIG),])
cleaned_combined$LENGTH = cleaned_combined$END - cleaned_combined$START
ggplot(subset(cleaned_combined, (!is.na(FollowsTreeAmp) | !is.na(FollowsTreeDel)) & LENGTH > 0), 
       aes(x=CONTIG, y=LENGTH, color=CONTIG=="chr4")) + geom_point() + facet_grid(FollowsTreeAmp ~ FollowsTreeDel) + 
  labs(title="Amplification on y axis, Deletion on x axis") + scale_y_log10()

cleaned_combined$AMP_COUNT = apply(cleaned_combined, 1, function(x){sum(x[4:15] == "+")})
cleaned_combined$DEL_COUNT = apply(cleaned_combined, 1, function(x){sum(x[4:15] == "-")})
ggplot(cleaned_combined, aes(x=AMP_COUNT, y=DEL_COUNT)) + geom_point()

bad_branches = subset(cleaned_combined, 
                      ((FollowsTreeAmp == "BadBranch" | FollowsTreeDel == "BadBranch") & (Call_HT6W47 == "+" | Call_HT6W47== "-")) |
                        (FollowsTreeAmp == "Leaf" & Call_HT6W47=="+") | (FollowsTreeDel == "Leaf" & Call_HT6W47== "-"))
write.table(data.frame(CHROM=bad_branches$CONTIG, START=bad_branches$START, END=bad_branches$END, STRAND="+", name="."), 
            "bad_branch.interval_list.tmp", col.names = F, row.names = F, quote = F, sep="\t")

#Only need amplifications because deletions will only cause HomVars which already won't be filtered.
good_branches = subset(cleaned_combined,
                       ((FollowsTreeAmp == "GoodBranch" | FollowsTreeAmp == "Germline") & Call_HT6W47 == "+"))
#round to nearest .5. Segment_Mean is in log_2
good_branches$COPY_RATIO = round((2^(good_branches$Segment_Mean_HT6W47))*2)/2
good_branches$TEST = 2^(good_branches$Segment_Mean_HT6W47)
write.table(data.frame(CONTIG=good_branches$CONTIG, START=good_branches$START, END=good_branches$END, COPY_RATIO=good_branches$COPY_RATIO),
            "good_branch_cnvs.seg", row.names = F, quote = F, sep="\t")
