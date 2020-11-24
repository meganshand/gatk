version 1.0

workflow varScanSomatic {
	input {
		File tumor_bam
		File tumor_bai
		File normal_bam
		File normal_bai
		File ref_fasta
		File ref_fasta_fai
		File ref_dict
		File varscan_jar
		String basename
		File header
	}

	Array[String] chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

	scatter(chr in chrs) {

		call pileup as normalPileup {
			input:
				bam = normal_bam,
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_dict = ref_dict,
				filename = "normal",
				chr = chr
		}

		call pileup as tumorPileup {
			input:
				bam = tumor_bam,
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_dict = ref_dict,
				filename = "tumor",
				chr = chr
		}

		call varScan {
			input:
				tumor_pileup = tumorPileup.pileup_out,
				normal_pileup = normalPileup.pileup_out,
				varscan_jar = varscan_jar,
				filename = basename
		}

		call filter as filterSnps {
			input:
				variants = varScan.snps,
				varscan_jar = varscan_jar
		}

		call filter as filterIndels {
			input:
				variants = varScan.indels,
				varscan_jar = varscan_jar
		}
	}

	call gather as gatherSnps {
		input:
			files = filterSnps.somatic,
			filename = basename + ".somatic_snps"
	}

	call gather as gatherIndels {
		input:
			files = filterIndels.somatic,
			filename = basename + ".somatic_indels"
	}

	call combineVariants {
		input:
			snps = gatherSnps.out,
			indels = gatherIndels.out,
			basename = basename,
			header = header
	}

	output {
		File out_vcf = combineVariants.vcf
	}
}

task combineVariants {
	input {
		File snps
		File indels
		String basename
		File header
	}
	command {
		set -e

R --vanilla <<CODE
	snps = read.table("~{snps}")
	vcf = data.frame(CHROM=snps[,'V1'], POS=snps[,'V2'], ID=".", REF=snps[,'V3'], ALT=snps[,'V4'], QUAL=round(-10*(log10(snps[,'V15']))), FILTER=".", INFO=".", FORMAT="GT", sample="0/1")
	write.table(vcf, "base_snps.txt", col.names = F, row.names = F, quote = F, sep="\t")

	indels = read.table("~{indels}", stringsAsFactors = F)
	indels[,'type'] = substr(indels[,'V4'], 1, 1)
	indels[,'bases'] = substring(indels[,'V4'], 2)
	indels[,'ref'] = ifelse(indels[,'type'] == "-", paste(indels[,'V3'], indels[,'bases'], sep=""), indels[,'V3'])
	indels[,'alt'] = ifelse(indels[,'type'] == "-", indels[,'V3'], paste(indels[,'V3'], indels[,'bases'], sep=""))
	vcf = data.frame(CHROM=indels[,'V1'], POS=indels[,'V2'], ID=".", REF=indels[,'ref'], ALT=indels[,'alt'], QUAL=round(-10*(log10(indels[,'V15']))), FILTER=".", INFO=".", FORMAT="GT", sample="0/1")
	write.table(vcf, "base_indels.txt", col.names = F, row.names = F, quote = F, sep="\t")

CODE

	cat ~{header} base_snps.txt > base_snps.vcf
	cat ~{header} base_indels.txt > base_indels.vcf

	gatk SortVcf -I base_snps.vcf -I base_indels.vcf -O ~{basename}.vcf

	}
	output {
		File vcf = "~{basename}.vcf"
	}
	runtime {
	    disks: "local-disk 100 HDD"
	    memory: "2 GB"
	    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0" 
  	}
}

task gather {
	input {
		Array[File] files
		String filename
	}
	command {
		cat ~{sep = " " files} > ~{filename}
	}
	output {
		File out = "~{filename}"
	}
	runtime {
		disks: "local-disk 100 HDD"
		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
	}
}

task filter {
	input {
		File variants
		File varscan_jar
	}
	String basename = basename(variants)
	command {
		set -e
		ln -s ~{variants} .
		java -jar ~{varscan_jar} processSomatic ~{basename}
		tail -n +2 ~{basename}.Somatic.hc > tmp.txt
		mv tmp.txt ~{basename}.Somatic.hc
	}
	output {
		File somatic_all_confidence = "~{basename}.Somatic"
		File germline = "~{basename}.Germline"
		File loh = "~{basename}.LOH"
        File somatic = "~{basename}.Somatic.hc"
	}
	runtime {
		disks: "local-disk 100 HDD"
		memory: "6 GB"
		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
	}
}

task pileup {
	input {
		File bam
		File ref_fasta
		File ref_fasta_fai
		File ref_dict
		String filename
		String chr
	}

	command {
		set -e
		samtools index ~{bam}
		samtools mpileup -q 1 -r ~{chr} -f ~{ref_fasta} ~{bam} > ~{filename}.pileup
	}
	output {
		File pileup_out = "~{filename}.pileup"
	}
	runtime {
		disks: "local-disk 200 HDD"
		memory: "6 GB"
		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
	}
}

task varScan {
	input {
		File tumor_pileup
		File normal_pileup
		File varscan_jar
		String filename
	}
	command {
		java -jar ~{varscan_jar} somatic ~{normal_pileup} ~{tumor_pileup} ~{filename}
	}
	output {
		File snps = "~{filename}.snp"
		File indels = "~{filename}.indel"
	}
	runtime {
		disks: "local-disk 300 HDD"
		memory: "8 GB"
		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
	}
}