version 1.0

workflow Benchmark {
    input{
        File stratIntervals
        String stratLabel
        File refDict
        String gatkTag = "4.0.11.0"
        File confidenceInterval
        File evalVCF
        File evalVCFIndex
        String evalLabel
        File truthVCF
        File truthVCFIndex
        File tumor_bam
        String truth_sample_name = basename(tumor_bam, ".bam")
        File reference
        File refIndex
        Int? threadsVcfEval=2     #threads (and cpu cores) to use for vcfeval.  Defaults to 2
        String referenceVersion = "hg38"
        String output_name
    }

    call CompressVcf as CompressTruthVcf {
    	input:
        	vcf = truthVCF,
        	sample_name = truth_sample_name
    }

    call CompressVcf as CompressEvalVcf {
        input:
            vcf = evalVCF,
            sample_name = "mixture"
    }

    call ConvertIntervals as StratConvertIntervals {
        input:
            inputIntervals=stratIntervals,
            refDict=refDict,
            gatkTag=gatkTag
    }

    call ConvertIntervals as ConfidenceConvertIntervals {
        input:
            inputIntervals=confidenceInterval,
            refDict=refDict,
            gatkTag=gatkTag
    }

    call CheckForVariants as CheckForVariantsEval {
        input:
            vcf=evalVCF,
            vcfIndex=evalVCFIndex,
            confidenceIL=ConfidenceConvertIntervals.intervalList,
            stratIL=stratIntervals,
            gatkTag=gatkTag
    }

    call CheckForVariants as CheckForVariantsTruth {
        input:
            vcf=truthVCF,
            vcfIndex=truthVCFIndex,
            confidenceIL=ConfidenceConvertIntervals.intervalList,
            stratIL=stratIntervals,
            gatkTag=gatkTag
    }
    if (CheckForVariantsTruth.variantsFound && CheckForVariantsEval.variantsFound) {
        call VcfEval as StandardVcfEval {
            input:
                truthVCF=CompressTruthVcf.compressedVcf,
                truthVCFIndex=CompressTruthVcf.compressedVcfIndex,
                evalVCF=CompressEvalVcf.compressedVcf,
                evalVCFIndex=CompressEvalVcf.compressedVcfIndex,
                confidenceBed=ConfidenceConvertIntervals.bed,
                stratBed=StratConvertIntervals.bed,
                ref=reference,
                refDict=refDict,
                refIndex=refIndex,
                outputPre="vcfeval",
                threads=threadsVcfEval
        }

        call WriteXMLfile as VcfEvalWriteXMLfile {
            input:
                input_files=[StandardVcfEval.outVcf,ConfidenceConvertIntervals.bed,StratConvertIntervals.bed],
                input_names=["vcfeval","confidence_intervals",stratIntervals],
                reference_version=referenceVersion,
                file_name="vcfeval"
        }

        call CountUNKVcfEval {
            input:
                vcf=StandardVcfEval.outVcf,
                vcfIndex=StandardVcfEval.outVcfIndex,
                gatkTag=gatkTag
        }
    }

    String areVariants=if(CheckForVariantsTruth.variantsFound && CheckForVariantsEval.variantsFound) then "yes" else "no"
    call SummariseVcfEval {
        input:
            evalLabel=evalLabel,
            truthLabel="LinSeq",
            stratLabel=stratIntervals,
            summaryFile=StandardVcfEval.outSummary,
            igvSession=VcfEvalWriteXMLfile.igv_session,
            areVariants=areVariants,
            unkSNP=CountUNKVcfEval.UNK_SNP,
            unkINDEL=CountUNKVcfEval.UNK_INDEL,
            output_name=output_name
    }

    output {
        File vcfEvalSummary = SummariseVcfEval.summaryOut
    }
}

# compress the truth vcf
task CompressVcf {
    input {
        File vcf
        String sample_name
    }
    String bn = basename(vcf)

    command<<<
        gatk --java-options "-Xmx7G" SelectVariants -V ~{vcf} -sn ~{sample_name} -O ~{bn}.gz
    >>>

    runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.2.0"
      memory: "8 GB"
    }

    output {
        File compressedVcf = "~{bn}.gz"
        File compressedVcfIndex = "~{bn}.gz.tbi"
    }
}

#takes in either a .bed or .intervallist and returns both a .bed and .intervallist version of the input
task ConvertIntervals {
    input {
        File inputIntervals
        Int? preemptible
        Int? memoryMaybe
        File refDict
        String gatkTag
    }

    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2

    command <<<
        if [[ ~{inputIntervals} == *.interval_list ]]; then
            echo "0" > inputType.txt
            gatk --java-options "-Xmx~{memoryJava}G" IntervalListToBed -I ~{inputIntervals} -O intervals.bed
            touch intervals.interval_list
        elif [[ ~{inputIntervals} == *.bed || ~{inputIntervals} == *.bed.gz ]]; then
            echo "1" > inputType.txt
            gatk --java-options "-Xmx~{memoryJava}G" BedToIntervalList -I ~{inputIntervals} -O intervals.interval_list -SD ~{refDict}
            touch intervals.bed
        fi
    >>>

    runtime {
        docker: "broadinstitute/gatk:"+gatkTag
        preemptible: select_first([preemptible,0])
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: "16"
        memory: memoryRam + " GB"
    }

    output {
        File bed=if read_int("inputType.txt")==1 then inputIntervals else "intervals.bed"
        File intervalList=if read_int("inputType.txt")==0 then inputIntervals else "intervals.interval_list"
    }
}

#Check to see if there are variants in the given vcf which overlap the confidence and stratification intervals
task CheckForVariants {
    input{
        File vcf
        File vcfIndex
        File confidenceIL
        File? stratIL
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
        }
        Int memoryDefault=16
        Int memoryJava=select_first([memoryMaybe,memoryDefault])
        Int memoryRam=memoryJava+2


    command <<<
    nVariants="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V ~{vcf} -L ~{confidenceIL} ~{"-L " + stratIL} -isr INTERSECTION | tail -1)"
    if [ "$nVariants" -gt "0" ]; then echo "true" > outBool.txt; else echo "false" > outBool.txt; fi
    >>>

    runtime {
            docker: "broadinstitute/gatk:"+gatkTag
            preemptible: select_first([preemptible,0])
            disks: "local-disk 100 SSD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }

    output {
        Boolean variantsFound=read_boolean("outBool.txt")
    }
}

#Evaluate evalVCF against truthVCF using vcfeval
task VcfEval {
    input{
        File truthVCF
        File truthVCFIndex
        File evalVCF
        File evalVCFIndex
        File confidenceBed
        File? stratBed
        File ref
        File refDict
        File refIndex
        String outputPre
        Int? preemptible
        String? memUser
        Int? threads
    }
    String memDefault="16 GB"
    String mem=select_first([memUser,memDefault])

    Int cpu=select_first([threads,1])

    command <<<
    /bin/rtg-tools/rtg format -o rtg_ref ~{ref}
    /bin/rtg-tools/rtg vcfeval -b ~{truthVCF} -c ~{evalVCF} -e ~{confidenceBed} ~{"--bed-regions " + stratBed} --output-mode combine --decompose -t rtg_ref ~{"--threads "+threads} -o output_dir
    for f in output_dir/*; do mv $f ~{outputPre}_"$(basename "$f")"; done
    python3 -<<"EOF" ~{outputPre}_snp_roc.tsv.gz ~{outputPre}_non_snp_roc.tsv.gz ~{outputPre}_summary.csv
    import gzip
    import sys

    indel_sensitivity=0
    indel_precision=0
    indel_fscore=0
    indel_TP_Base=0
    indel_TP_Eval=0
    indel_FP=0
    indel_FN=0

    snp_sensitivity=0
    snp_precision=0
    snp_fscore=0
    snp_TP_Base=0
    snp_TP_Eval=0
    snp_FP=0
    snp_FN=0

    with gzip.open(sys.argv[1],"rt") as f_snp:
        for line in f_snp:
            try:
                snp_sensitivity=float(line.split()[6])
                snp_precision=float(line.split()[5])
                snp_fscore=float(line.split()[7])
                snp_TP_Eval=float(line.split()[3])
                snp_TP_Base=float(line.split()[1])
                snp_FP=float(line.split()[2])
                snp_FN=float(line.split()[4])
            except ValueError:
                continue
            except IndexError:
                continue
        f_snp.close()
    with gzip.open(sys.argv[2],"rt") as f_indel:
        for line in f_indel:
            try:
                indel_sensitivity=float(line.split()[6])
                indel_precision=float(line.split()[5])
                indel_fscore=float(line.split()[7])
                indel_TP_Eval=float(line.split()[3])
                indel_TP_Base=float(line.split()[1])
                indel_FP=float(line.split()[2])
                indel_FN=float(line.split()[4])
            except ValueError:
                continue
            except IndexError:
                continue
        f_indel.close()

    str_indel_sensitivity=str(indel_sensitivity)
    str_indel_precision=str(indel_precision)
    str_indel_fscore=str(indel_fscore)
    str_snp_sensitivity=str(snp_sensitivity)
    str_snp_precision=str(snp_precision)
    str_snp_fscore=str(snp_fscore)


    if indel_TP_Eval+indel_FP==0:
        str_indel_precision="NA"
    if indel_TP_Base+indel_FN==0:
        str_indel_sensitivity="NA"
    if str_indel_sensitivity=="NA" or str_indel_precision=="NA":
        str_indel_fscore="NA"

    if snp_TP_Eval+snp_FP==0:
        str_snp_precision="NA"
    if snp_TP_Base+snp_FN==0:
        str_snp_sensitivity="NA"
    if str_snp_sensitivity=="NA" or str_snp_precision=="NA":
        str_snp_fscore="NA"



    with open(sys.argv[3],"wt") as f_out:
        f_out.write(",".join(["Type","Precision","Recall","F1_Score","TP_Eval","TP_Base","FP","FN"])+"\n")
        f_out.write(",".join(["SNP",str_snp_precision,str_snp_sensitivity,str_snp_fscore,str(snp_TP_Eval),str(snp_TP_Base),str(snp_FP),str(snp_FN)])+"\n")
        f_out.write(",".join(["INDEL",str_indel_precision,str_indel_sensitivity,str_indel_fscore,str(indel_TP_Eval),str(indel_TP_Base),str(indel_FP),str(indel_FN)])+"\n")
        f_out.close()
    EOF
    >>>

    runtime {
        docker: "ckachulis/rtg-tools:0.1"
        preemptible: select_first([preemptible,0])
        memory: mem
        cpu: cpu
    }

    output {
        Array[File] outs=glob("${outputPre}_*")
        File outSummary="${outputPre}_summary.csv"
        File outVcf="${outputPre}_output.vcf.gz"
        File outVcfIndex="${outputPre}_output.vcf.gz.tbi"
    }
}

# creates an IGV session
# given a list of IGV compatible file paths
task WriteXMLfile {
    input {
        Array[String] input_files
        String reference_version
        String file_name

        Array[String]? input_names=""
        Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []
    }
    command {
        bash /usr/writeIGV.sh ~{reference_version} ~{sep=" " input_files} ~{sep=" " input_names_prefix}  > "~{file_name}.xml"
    }
    runtime {
        docker: "quay.io/mduran/generate-igv-session_2:v1.0"
    }
    output {
        File igv_session = "${file_name}.xml"
    }
}

#Count number of variants which were outside confidence region based on vcfeval annotated vcf
task CountUNKVcfEval {
    input {
        File? vcf=""
        File? vcfIndex=""
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
    }
    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2

    command <<<
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcf} -O selected.unk.snp.vcf.gz -select "(CALL == 'OUT')" --select-type-to-include SNP
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcf} -O selected.unk.indel.vcf.gz -select "(CALL == 'OUT')" --select-type-to-include INDEL

        UNK_SNP="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.snp.vcf.gz | tail -1)"
        UNK_INDEL="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.indel.vcf.gz | tail -1)"

        echo "$UNK_SNP" > unk_snp.txt
        echo "$UNK_INDEL" > unk_indel.txt
    >>>


    runtime {
                docker: "broadinstitute/gatk:"+gatkTag
                preemptible: select_first([preemptible,0])
                disks: "local-disk 100 SSD"
                bootDiskSizeGb: "16"
                memory: memoryRam + " GB"
    }

    output {
        Int UNK_SNP=read_int("unk_snp.txt")
        Int UNK_INDEL=read_int("unk_indel.txt")
    }
}

#Convert vcfeval output statistics to final output format
task SummariseVcfEval {
    input {
        String evalLabel
        String truthLabel
        String areVariants
        String? igvSession
        String? stratLabel
        File? summaryFile
        Int? unkSNP
        Int? unkINDEL
        Int? preemptible
        String output_name
    }

    command <<<
        Rscript -<<"EOF" ~{evalLabel} ~{truthLabel} ~{default="" summaryFile} ~{default="" stratLabel} ~{default="" igvSession} ~{default="" unkSNP} ~{default="" unkINDEL} ~{areVariants}
        args <- commandArgs(trailingOnly = TRUE)
        if (args[length(args)]=="yes") {
            table_vcfeval <- read.csv(args[3])
            if (length(args)==7) {
                table_vcfeval$Stratifier <- NA
                table_vcfeval$IGV_Session <- args[4]
                table_vcfeval$UNK[table_vcfeval$Type=="SNP"]=args[5]
                table_vcfeval$UNK[table_vcfeval$Type=="INDEL"]=args[6]
            } else {
                table_vcfeval$Stratifier <- args[4]
                table_vcfeval$IGV_Session <- args[5]
                table_vcfeval$UNK[table_vcfeval$Type=="SNP"]=args[6]
                table_vcfeval$UNK[table_vcfeval$Type=="INDEL"]=args[7]
            }
        } else {
            types <- c("INDEL","SNP")
            recall <- c(NA,NA)
            precision <- c(NA,NA)
            f1_score <- c(NA,NA)
            tp <- c(0,0)
            fp <- c(0,0)
            fn <- c(0,0)
            unk <- c(0,0)
            igv_session <- c(NA,NA)
            table_vcfeval <- data.frame("Type"=types,"Recall"=recall,"Precision"=precision,"F1_Score"=f1_score,"TP_Base"=tp,"TP_Eval"=tp,"FP"=fp,"FN"=fn,"UNK"=unk,"IGV_Session"=igv_session)
            if (length(args)==3) {
                table_vcfeval$Stratifier <- NA
            } else {
                table_vcfeval$Stratifier <- args[3]
            }
        }
        table_vcfeval$Name <- args[1]
        table_vcfeval$Truth_Set <- args[2]
        table_vcfeval$Summary_Type <- "summary"
        table_vcfeval$Comparison_Engine <-"VcfEval"
        write.csv(table_vcfeval,"~{output_name}.vcfeval.summary.csv",row.names=FALSE)
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        preemptible: select_first([preemptible,0])
    }

    output{
        File summaryOut="~{output_name}.vcfeval.summary.csv"
    }
}