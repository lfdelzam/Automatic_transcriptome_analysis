# list of samples (input files)
import os
file_list = os.listdir("./data")
samples = sorted([i.split(".")[0] for i in file_list if i.find(".fastq") !=-1]) #list of fastq files

rule all:
    input:
        expand("results/1_quality/{s}_fastqc.html", s=samples),
        expand("results/3_mapping/genome.{num}.ht2", num=["1","2","3","4","5","6","7","8"]),
        expand("results/5_differential_expression/diffExp{cut}tab", cut=[".0.01.",".0.05.","."]),
        "results/5_differential_expression/normalized.count",
        expand("results/6_visualize/fig{names}", names=["1.pdf","2.pdf","3.pdf","4.png"]),
        "results/5_differential_expression/topmorexpressed",
        "results/5_differential_expression/morexpresedseq.fasta"

rule quality:
    input: "data/{s}.fastq"
    output:"results/1_quality/{s}_fastqc.html"
    version:"v0.11.8"
    shell: "bin/FastQC/fastqc --threads 6 {input} -o results/1_quality 2>> running.log"

#-o name output directory

rule trimming:
    input: "data/{s}.fastq"
    output: "results/2_trimming/{s}.clean.fastq"
    version: "0.38"
    shell: "java -jar ./bin/trimmomatic-0.38.jar SE {input} {output} AVGQUAL:28 MINLEN:46 CROP:46 SLIDINGWINDOW:8:20 2>> running.log"

#SE stands for single end
#AVGQUAL is the average quality threshold. Sequences with a quality below this will be discarded.
#MINLEN is the minimum length threshold. Shorter sequences will be discarded (that is the length after clipping).
#CROP means that the output sequences is clipped to a length of 46 (in this case).
#SLIDINGWINDOW here means that a window size of 8
    #is used and if the average quality of this window is below 20, the sequences will be clipped from where the windows started

rule index:
    input: "data/genome.fna"
    output: expand("results/3_mapping/genome.{num}.ht2", num=["1","2","3","4","5","6","7","8"])
    version: "2.1.0"
    shell: "bin/hisat2-2.1.0/hisat2-build {input} results/3_mapping/genome 2>> running.log"
    #Hisat Creates a genome index called genome (results/3_mapping/genome)

rule mapping:
    input:"results/2_trimming/{s}.clean.fastq"
    output:"results/3_mapping/{s}.sam"
    version: "2.1.0"
    shell: "bin/hisat2-2.1.0/hisat2 -p 4 --max-intronlen 5000 -U {input} -x results/3_mapping/genome -S {output} 2>> running.log"

#-p Number of threads
#-max-intronlen Maximum intron size
#-U Input file (unpaired reads)
#-x Base name of index
#-S Output file with aligned reads in SAM format

rule readcounts:
    input:
        sam = "results/3_mapping/{s}.sam",
        gff = "data/genes.gff"
    output: "results/4_Read_counts/{s}.count"
    version: "0.11.2"
    shell:
        """htseq-count -s no -t CDS -i name -m intersection-nonempty {input.sam} {input.gff} | grep -v "^__" | tr -d "#" > {output} && 2>> running.log"""

# We have removed the hash sign from genes names (# is a special character in R)

# -i name option refers to the id of the genomic feature. In the last column of the gff file, the id is named name
# -s <yes/no/reverse>, --stranded=<yes/no/reverse> whether the data is from a strand-specific assay (default: yes)
   # For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same
   # or the opposite strand as the feature.
# -t <feature type>, --type=<feature type>
  #feature type (3rd column in GFF file) to be used, all features of other type are ignored
# -m <mode>, --mode=<mode>
  #Mode to handle reads overlapping more than one feature.
  #Possible values for <mode> are union, intersection-strict and intersection-nonempty
#Take a look at https://htseq.readthedocs.io/en/release_0.11.1/count.html to see what intersection-nonempty means.

rule DifferentialExpression:
        input: expand("results/4_Read_counts/{s}.count", s=samples)
        output: data = expand("results/5_differential_expression/diffExp{cut}tab", cut=[".0.01.",".0.05.","."]),
                count = "results/5_differential_expression/normalized.count",
                figures = expand("results/6_visualize/fig{names}", names=["1.pdf","2.pdf","3.pdf","4.png"])
        script:"scr/DifferentialExpression.R"

rule top_more_expressed:
        input:"results/5_differential_expression/diffExp.0.01.tab"
        output: "results/5_differential_expression/topmorexpressed"
        shell:"""
                cat {input} | cut -f 1,2 | head -1 > {output}
                cat {input} | cut -f 2 | sort -n -r -u | head -20 | while read line; do grep $line {input} | cut -f 1,2 ; done >> {output}
              """

rule sequences:
    input: gf="data/genes.gff",
           ge="data/genome.fna",
           more="results/5_differential_expression/topmorexpressed"
    output: "results/5_differential_expression/morexpresedseq.fasta"
    shell: "./scr/seqfromgenome.py -g {input.gf} -i {input.ge} -t {input.more} -o {output}"
#the program reads a file containing the name of the more expressed genes, and an anotation file (gff).
#The correspondig CDSs positions in the genome are extracted.
#The sequences are then extracted from the genome and a fasta file with the sequences is generated.

#optional arguments:
 # -g gff file
 # -i input file containing the genome, type .fna
 # -t input file containing the top genes more expressed
 # -o output fasta file containing the list of extracted gene sequences
