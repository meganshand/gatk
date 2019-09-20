#! /bin/bash

samples=(HT115_42112 HT115_42221 HT115_5211 HT115_7221 HT6W34 HT6W37 HT6W45 HT6W47 HT6W48 HT6W49 HT6W54 HT6W63)

for (( i=0; i<${#samples[@]}; i+=2)); do
	java -jar ~/Github/gatk/build/libs/gatk.jar CombineSegmentBreakpoints \
		--segments ${samples[i]}.called.igv.seg --segments ${samples[i+1]}.called.igv.seg \
		--labels ${samples[i]} --labels ${samples[i+1]} \
		-O ${samples[i]}.${samples[i+1]}.seg \
		--columns-of-interest Call --columns-of-interest Segment_Mean \
		-R /Volumes/seq_references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta
done

for (( i=0; i<${#samples[@]}; i+=4)); do
	java -jar ~/Github/gatk/build/libs/gatk.jar CombineSegmentBreakpoints \
		--segments ${samples[i]}.${samples[i+1]}.seg --segments ${samples[i+2]}.${samples[i+3]}.seg \
		-O combined_four_${i}.seg \
		--columns-of-interest Call_${samples[i]} --columns-of-interest Call_${samples[i+1]} --columns-of-interest Call_${samples[i+2]} --columns-of-interest Call_${samples[i+3]} \
		--columns-of-interest Segment_Mean_${samples[i]} --columns-of-interest Segment_Mean_${samples[i+1]} --columns-of-interest Segment_Mean_${samples[i+2]} --columns-of-interest Segment_Mean_${samples[i+3]} \
		-R /Volumes/seq_references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta
done

java -jar ~/Github/gatk/build/libs/gatk.jar CombineSegmentBreakpoints \
		--segments combined_four_0.seg --segments combined_four_4.seg \
		-O combined_eight.seg \
		--columns-of-interest Call_${samples[0]} --columns-of-interest Call_${samples[1]} --columns-of-interest Call_${samples[2]} --columns-of-interest Call_${samples[3]} \
		--columns-of-interest Call_${samples[4]} --columns-of-interest Call_${samples[5]} --columns-of-interest Call_${samples[6]} --columns-of-interest Call_${samples[7]} \
		--columns-of-interest Segment_Mean_${samples[0]} --columns-of-interest Segment_Mean_${samples[1]} --columns-of-interest Segment_Mean_${samples[2]} --columns-of-interest Segment_Mean_${samples[3]} \
		--columns-of-interest Segment_Mean_${samples[4]} --columns-of-interest Segment_Mean_${samples[5]} --columns-of-interest Segment_Mean_${samples[6]} --columns-of-interest Segment_Mean_${samples[7]} \
		-R /Volumes/seq_references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta

java -jar ~/Github/gatk/build/libs/gatk.jar CombineSegmentBreakpoints \
		--segments combined_four_8.seg --segments combined_eight.seg \
		-O combined.seg \
		--columns-of-interest Call_${samples[0]} --columns-of-interest Call_${samples[1]} --columns-of-interest Call_${samples[2]} --columns-of-interest Call_${samples[3]} \
		--columns-of-interest Call_${samples[4]} --columns-of-interest Call_${samples[5]} --columns-of-interest Call_${samples[6]} --columns-of-interest Call_${samples[7]} \
		--columns-of-interest Call_${samples[8]} --columns-of-interest Call_${samples[9]} --columns-of-interest Call_${samples[10]} --columns-of-interest Call_${samples[11]} \
		--columns-of-interest Segment_Mean_${samples[0]} --columns-of-interest Segment_Mean_${samples[1]} --columns-of-interest Segment_Mean_${samples[2]} --columns-of-interest Segment_Mean_${samples[3]} \
		--columns-of-interest Segment_Mean_${samples[4]} --columns-of-interest Segment_Mean_${samples[5]} --columns-of-interest Segment_Mean_${samples[6]} --columns-of-interest Segment_Mean_${samples[7]} \
		--columns-of-interest Segment_Mean_${samples[8]} --columns-of-interest Segment_Mean_${samples[9]} --columns-of-interest Segment_Mean_${samples[10]} --columns-of-interest Segment_Mean_${samples[11]} \
		-R /Volumes/seq_references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta
