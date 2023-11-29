'''
    Pipeline for RNA Seq processing of Parkinsonian Tissue Data
    Tools Used:
        - FastQC
        - Trimmomatic
        - Hisat2
        - Samtools
        - StringTie
        - DESeq2
    Above generates counts.txt.
    Differential Analysis of gene expression profiles using DeSeq.py

    To Run:
        To run the entire pipeline, run:
        $ python3 pe_rna_seq.py --runall
'''

import argparse
import time
import os

if __name__ == "__main__":

    input_dir = "se_data/"
    fastq_dir = "fastq/"
    trim_dir = "trim3/"
    hisat_dir = "hisat23/"
    stringtie_dir = "stringtie2/"
    ref_genome = "hg38/hg38"
    # deseq_dir = "deseq/"

    argparser = argparse.ArgumentParser(description="RNA Seq Pipeline")

    argparser.add_argument("--fastqc", help="Run FastQC", default=False)
    argparser.add_argument("--trimmomatic", help="Run Trimmomatic", default=False)
    argparser.add_argument("--hisat2", help="Run Hisat2", default=False)
    argparser.add_argument("--samtools", help="Run Samtools", default=False)
    argparser.add_argument("--stringtie", help="Run StringTie", default=False)
    #argparser.add_argument("--deseq2", help="Run DESeq2", default=False)
    argparser.add_argument("--runall", help="Run all the tools")
    argparser.add_argument("--SRR", nargs='+', help="List of SRRs to run the pipeline on")

    args = argparser.parse_args()
    srr_list = args.SRR
    
    for SRR in srr_list:
        if args.runall:
            args.fastqc = False
            args.trimmomatic = True
            args.hisat2 = True
            args.samtools = True
            args.stringtie = True
            #args.deseq2 = True
            
            print(f"******************************************")
            print(f"RUNNING COMPLETE ANALYSIS ON {SRR}")
            print(f"******************************************")


        if args.fastqc:
            # FastQC
            start_time = time.time()
            print(f"Running FastQC on single-end reads of SRR:{SRR}")
            os.system(f"fastqc {input_dir}{SRR}.fastq.gz -o {fastq_dir}")
            print(f"FastQC completed in {time.time() - start_time} seconds")
            print(f"******************************************")

        if args.trimmomatic:
            # Trimmomatic
            start_time = time.time()
            print(f"Running Trimmomatic on single-end reads of SRR:{SRR}")
            os.system(f"java -Xmx27g -jar trimmomatic-0.39.jar SE -threads 16 {input_dir}{SRR}.fastq.gz {trim_dir}{SRR}_trim.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 MINLEN:50 SLIDINGWINDOW:4:20")
            print(f"Trimmomatic completed in {time.time() - start_time} seconds")
            print(f"******************************************")

            # FastQC for trimmed reads
            #start_time = time.time()
            #print(f"Running FastQC on trimmed paired and unpaired reads of SRR:{SRR}")
            #os.system(f"fastqc {trim_dir}{SRR}_1P.fastq.gz {trim_dir}{SRR}_1U.fastq.gz {trim_dir}{SRR}_2P.fastq.gz {trim_dir}{SRR}_2U.fastq.gz -o {fastq_dir}")
            #print(f"FastQC completed in {time.time() - start_time} seconds")
            #print(f"******************************************")
            
        if args.hisat2:
            # Hisat2
            start_time = time.time()
            print(f"Running Hisat2 on trimmed reads of SRR:{SRR}")
            # Currently running on paired trimmed reads and not using unpaired output reads
            os.system(f"hisat2 --phred33 --dta -x {ref_genome} -U {trim_dir}{SRR}_trim.fastq.gz -S {hisat_dir}{SRR}_alignment.sam -p 12")
            print(f"Hisat2 completed in {time.time() - start_time} seconds")
            print(f"******************************************")

        if args.samtools:
            # Samtools
            start_time = time.time()
            print(f"Running Samtools on Hisat2 output of SRR:{SRR}")
            os.system(f"samtools sort -o {hisat_dir}{SRR}_alignment.bam {hisat_dir}{SRR}_alignment.sam")
            os.system(f"samtools index {hisat_dir}/{SRR}_alignment.bam")
            print(f"Samtools completed in {time.time() - start_time} seconds")
            print(f"******************************************")
            
        if args.stringtie:
            # StringTie
            new_line = f"{SRR} {stringtie_dir}{SRR}.gtf"
            with open("lst.txt", "a") as file:
                file.write(new_line + "\n")
            start_time = time.time()
            print(f"Running StringTie on Samtools output of SRR:{SRR}")
            os.system(f"stringtie -p 8 -e -G hg38.knownGene.gtf {hisat_dir}{SRR}_alignment.bam -o {stringtie_dir}{SRR}.gtf")
            print(f"StringTie completed in {time.time() - start_time} seconds")
            print(f"******************************************")
            
            
        
        #if args.deseq2:
            # DeSeq2
            #start_time = time.time()
            #print(f"Running PrepDE on StringTie output of SRR:{SRR}")
            #os.system(f"python3 prepDE.py -i {stringtie_dir} -g {deseq_dir}/gene_count_matrix.csv -t ./prepDE/{deseq_dir}transcript_count_matrix.csv")
            #print(f"PrepDE completed in {time.time() - start_time} seconds")
            #print(f"******************************************")

            #start_time = time.time()
            #print(f"Running DESeq2 on PrepDE output of SRR:{SRR}")
            #os.system(f"Rscript DESeq2.R")
            #print(f"DESeq2 completed in {time.time() - start_time} seconds")
            #print(f"******************************************")
            #print(f"COMPLETED WITH ENTIRE ANALYSIS PIPELINE")
            #print(f"******************************************")
            
        if args.runall:
            print(f"COMPLETED WITH ENTIRE ANALYSIS PIPELINE")
            print(f"******************************************")
            start_time = time.time()
            print(f"Deleting all unnecessary files of SRR:{SRR}")
            os.system(f"rm -rf {trim_dir}")
            os.system(f"mkdir {trim_dir}")
            print(f"Deleted 4 trimmed reads from {trim_dir}")
            
            os.system(f"rm -rf {hisat_dir}")
            os.system(f"mkdir {hisat_dir}")
            print(f"Deleted bam and sam alignment files from {hisat_dir}")
            
            os.system(f"cd ..")

            print(f"All deletions completed in {time.time() - start_time} seconds")
            print(f"******************************************")
        

        
	
