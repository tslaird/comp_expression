SRA_SAMPLES=["SRX3312086","SRX3312085","SRX4042153","SRX4042154","SRX4042155","SRX4042156","SRX4042157","SRX4042158"]
READS=[1,2]

rule all:
    input:
        expand("infofiles/{sample}.info", sample=SRA_SAMPLES),
        expand("fastq_files/{sample}_1.fastq.gz", sample = SRA_SAMPLES),
        expand("fastq_files/{sample}_2.fastq.gz", sample = SRA_SAMPLES)
        # expand("cleaned_reads/{sample}_1.clean.fastq", sample=SAMPLES.split(" ")),
        # "genome/GCF_000737145.1_ASM73714v1_genomic.gbff.gz",
        # "genome/GCF_000737145.1_ASM73714v1_genomic.gbff",
        # "genome/GCF_000737145.1_ASM73714v1_genomic_transcripts.fasta",
        # expand("cleaned_reads/{sample}_{read}.clean_fastqc.html", sample=SAMPLES.split(" "),read=READS.split(" ")),
        # directory("salmon_transcriptome_index"),
        # expand("{sample}_quant/quant.sf", sample=SAMPLES.split(" "))

rule download_infofiles:
    output:"infofiles/{sample}.info"
    run:
        url= "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession="+wildcards.sample
        urllib.request.urlretrieve(url,"infofiles/"+wildcards.sample+".info")

rule download_fastq:
    input:"infofiles/{sample}.info"
    output:
        "fastq_files/{sample}_1.fastq.gz",
        "fastq_files/{sample}_2.fastq.gz"
    run:
        info_df=pd.read_table(str(input), sep='\t')
        ftp_links=info_df['fastq_ftp'][0].split(";")
        md5_sums=info_df['fastq_md5'][0].split(";")
        for i in [0,1]:
            path_fastq = "http://"+ftp_links[i]
            name = 'fastq_files/'+path_fastq.split('/')[-1]
            md5= md5_sums[i]
            retries=20
            while(retries > 0):
                try:
                    urllib.request.urlretrieve(path_fastq,name)
                    with open(name, 'rb') as file_to_check:
                        data = file_to_check.read()
                        md5_returned = str(hashlib.md5(data).hexdigest())
                    if md5_returned==md5:
                        print("Fetched " + name)
                        break
                except:
                    print("Retrying download from " + path_fastq)
                    retries = retries - 1
                    continue


# rule uncompress_genome:
#     input: "genome/GCF_000737145.1_ASM73714v1_genomic.gbff.gz"
#     output: "genome/GCF_000737145.1_ASM73714v1_genomic.gbff"
#     shell:
#         "gunzip -k {input}"
#
# rule bbduk_trim:
#     input:
#     	read1="raw_reads/{sample}_1.fastq",
#     	read2="raw_reads/{sample}_2.fastq",
#     output:
#     	read1="cleaned_reads/{sample}_1.clean.fastq",
#     	read2="cleaned_reads/{sample}_2.clean.fastq",
#     shell: """bbduk.sh in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} ref=/anaconda3/opt/bbmap-38.22-1/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 trimq=20"""
#
# rule fastqc_clean:
#     input: "cleaned_reads/{sample}.clean.fastq"
#     output: "cleaned_reads/{sample}.clean_fastqc.html"
#     shell: """fastqc {input}"""
#
# rule generate_transcriptome:
#     input: "genome/GCF_000737145.1_ASM73714v1_genomic.gbff"
#     output: "genome/GCF_000737145.1_ASM73714v1_genomic_transcripts.fasta"
#     shell: """./src/parse_gbff_transcripts.py {input}"""
#
# rule index_transcriptome:
#     input: "genome/GCF_000737145.1_ASM73714v1_genomic_transcripts.fasta"
#     output: directory("salmon_transcriptome_index")
#     conda: "envs/salmon.yaml"
#     shell: """salmon index -t {input} -i {output} -k 31"""
#
# rule map_reads:
#     input:
#     	read1="cleaned_reads/{sample}_1.clean.fastq",
#     	read2="cleaned_reads/{sample}_2.clean.fastq",
#     params:
#         outdir="{sample}_quant"
#     output: "{sample}_quant/quant.sf"
#     conda: "envs/salmon.yaml"
#     shell: """salmon quant -i "salmon_transcriptome_index" -l A -1 {input.read1} -2 {input.read2} -o {params.outdir}"""
