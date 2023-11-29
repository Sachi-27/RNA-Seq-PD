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

    input_dir = "data/"
    fastq_dir = "fastq/"
    trim_dir = "trim/"
    hisat_dir = "hisat2/"
    stringtie_dir = "stringtie/"
    ref_genome = "hg38/hg38"
    # deseq_dir = "deseq/"

    argparser = argparse.ArgumentParser(description="RNA Seq Pipeline")

    argparser.add_argument("--download", help="Download FASTQ Reads", default=False)
    argparser.add_argument("--fastqc", help="Run FastQC", default=False)
    argparser.add_argument("--trimmomatic", help="Run Trimmomatic", default=False)
    argparser.add_argument("--hisat2", help="Run Hisat2", default=False)
    argparser.add_argument("--samtools", help="Run Samtools", default=False)
    argparser.add_argument("--stringtie", help="Run StringTie", default=False)
    #argparser.add_argument("--deseq2", help="Run DESeq2", default=False)
    argparser.add_argument("--runall", help="Run all the tools")
    argparser.add_argument("--SRR", nargs='+', help="List of SRRs to run the pipeline on")
    argparser.add_argument("--i", help="To interpret SRR as <start> <end> <offset>", default=False)

    args = argparser.parse_args()
    srr_list = args.SRR
    if(args.i):
        srr_start = int(srr_list[0][3:])
        srr_end = int(srr_list[1][3:])
        offset = int(srr_list[2])
        srr_list = [f"SRR{x}" for x in range(srr_start, srr_end+1, offset)]
        
    for SRR in srr_list:
        if args.runall:
            args.download = True
            args.fastqc = False
            args.trimmomatic = True
            args.hisat2 = True
            args.samtools = True
            args.stringtie = True
            #args.deseq2 = True
            
            print(f"******************************************")
            print(f"RUNNING COMPLETE ANALYSIS ON {SRR}")
            print(f"******************************************")


        if args.download:
            start_time = time.time()
            print(f"Downloading Paired-end Reads of {SRR}")
            os.system(f"prefetch -p {SRR} -O {input_dir} && fasterq-dump -pe 24 {input_dir}{SRR}/{SRR}.sra -O {input_dir}")
            os.system(f"rm -rf {input_dir}{SRR}")
            print(f"Download completed in {time.time() - start_time} seconds")
            print(f"******************************************")  

            

        if args.fastqc:
            # FastQC
            start_time = time.time()
            print(f"Running FastQC on paired-end reads of SRR:{SRR}")
            os.system(f"fastqc {file1} {file2} -o {fastq_dir}")
            print(f"FastQC completed in {time.time() - start_time} seconds")
            print(f"******************************************")

        if args.trimmomatic:

            # Trimmomatic
            start_time = time.time()
            print(f"Running Trimmomatic on paired-end reads of SRR:{SRR}")
            os.system(f"java -Xmx27g -jar trimmomatic-0.39.jar PE -threads 16 {input_dir}{SRR}_1.fastq {input_dir}{SRR}_2.fastq -baseout {trim_dir}{SRR}.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:8:TRUE MINLEN:36 TRAILING:3 LEADING:3 SLIDINGWINDOW:4:15")
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
            print(f"Running Hisat2 on trimmed paired reads of SRR:{SRR}")
            # Currently running on paired trimmed reads and not using unpaired output reads
            os.system(f"hisat2 --phred33 --dta -x {ref_genome} -1 {trim_dir}{SRR}_1P.fastq.gz -2 {trim_dir}{SRR}_2P.fastq.gz -S {hisat_dir}{SRR}_alignment.sam -p 12")
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
            os.system(f"stringtie -p 16 -e -G hg38.knownGene.gtf {hisat_dir}{SRR}_alignment.bam -o {stringtie_dir}{SRR}.gtf")
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
            
            os.system(f"rm -rf {input_dir}")
            os.system(f"mkdir {input_dir}")
            print(f"Removed fastq reads from {input_dir}")

            print(f"All deletions completed in {time.time() - start_time} seconds")
            print(f"******************************************")
        

        
	
