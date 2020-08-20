import pandas as pd
import urllib.request
import hashlib

ESOLI_SAMPLES=["SRR7120554","SRR7120555","SRR7120556","SRR7120557","SRR7120558","SRR7120559"]
ABAU_SAMPLES=["SRR6202159","SRR6202160"]

rule all:
    input:
        expand("Abau/{abau}.info", abau= ABAU_SAMPLES),
        expand("Esoli/{esoli}.info", esoli= ESOLI_SAMPLES),
        expand("Esoli/fastq_files/{esoli}.fastq.gz", esoli= ESOLI_SAMPLES),
        expand("Abau/fastq_files/{abau}_1.fastq.gz",abau= ABAU_SAMPLES),
        expand("Abau/fastq_files/{abau}_2.fastq.gz", abau= ABAU_SAMPLES),
        expand("Esoli/fastq_files/{esoli}.fastq", esoli= ESOLI_SAMPLES),
        expand("Abau/fastq_files/{abau}_1.fastq",abau= ABAU_SAMPLES),
        expand("Abau/fastq_files/{abau}_2.fastq", abau= ABAU_SAMPLES),
        "Abau/genome/GCF_009759685.1_ASM975968v1_genomic.gbff",
        "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic.gbff",
        "Abau/genome/GCF_009759685.1_ASM975968v1_genomic_transcripts.fasta",
        "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic_transcripts.fasta",
        directory("Esoli/salmon_transcriptome_index"),
        directory("Abau/salmon_transcriptome_index"),
        expand("Abau/{abau}_quant/quant.sf", abau= ABAU_SAMPLES),
        expand("Esoli/{esoli}_quant/quant.sf", esoli=ESOLI_SAMPLES)
        #expand("cleaned_reads/{sample}_1.clean.fastq", sample=SAMPLES.split(" ")),
        # "genome/GCF_000737145.1_ASM73714v1_genomic.gbff.gz",
        # "genome/GCF_000737145.1_ASM73714v1_genomic.gbff",
        # "genome/GCF_000737145.1_ASM73714v1_genomic_transcripts.fasta",
        # expand("cleaned_reads/{sample}_{read}.clean_fastqc.html", sample=SAMPLES.split(" "),read=READS.split(" ")),
        # directory("salmon_transcriptome_index"),
        # expand("{sample}_quant/quant.sf", sample=SAMPLES.split(" "))

rule download_abau_infofiles:
    output: expand("Abau/{abau}.info", abau=ABAU_SAMPLES)
    run:
        print(wildcards.abau)
        url= "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession="+wildcards.abau
        urllib.request.urlretrieve(url,"Abau/"+wildcards.abau+".info")

rule download_esoli_infofiles:
    output: expand("Esoli/{esoli}.info", esoli=ESOLI_SAMPLES)
    run:
        print(wildcards.esoli)
        url= "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession="+wildcards.esoli
        urllib.request.urlretrieve(url,"Esoli/"+wildcards.esoli+".info")

rule download_abau_fastq:
    input:"Abau/{abau}.info"
    output:
        "Abau/fastq_files/{abau}_1.fastq.gz",
        "Abau/fastq_files/{abau}_2.fastq.gz"
    run:
        info_df=pd.read_table(str(input), sep='\t')
        ftp_links=info_df['fastq_ftp'][0].split(";")
        md5_sums=info_df['fastq_md5'][0].split(";")
        for i in [0,1]:
            path_fastq = "http://"+ftp_links[i]
            name = 'Abau/fastq_files/'+path_fastq.split('/')[-1]
            link="wget -c --tries=50 -O "+name+" "+path_fastq
            print(link)
            md5= md5_sums[i]
            retries=2
            while(retries > 0):
                try:
                    shell(link)
                    with open(name, 'rb') as file_to_check:
                        data = file_to_check.read()
                        md5_returned = str(hashlib.md5(data).hexdigest())
                        print(md5_returned+" : "+md5)
                    if md5_returned==md5:
                        print("Fetched " + name)
                        break
                except:
                    print("Retrying download from " + path_fastq)
                    retries = retries - 1
                    continue

rule download_esoli_fastq:
    input:"Esoli/{esoli}.info"
    output:
        "Esoli/fastq_files/{esoli}.fastq.gz",
    run:
        info_df=pd.read_table(str(input), sep='\t')
        ftp_links=info_df['fastq_ftp']
        md5_sums=info_df['fastq_md5']
        path_fastq = "http://"+ftp_links[0]
        name = 'Esoli/fastq_files/'+path_fastq.split('/')[-1]
        link="wget -c --tries=50 -O "+name+" "+path_fastq
        print(link)
        shell(link)
        md5= md5_sums[0]
        retries=2
        while(retries > 0):
            try:
                shell(link)
                with open(name, 'rb') as file_to_check:
                    data = file_to_check.read()
                    md5_returned = str(hashlib.md5(data).hexdigest())
                    print(md5_returned+" : "+md5)
                if md5_returned==md5:
                    print("Fetched " + name)
                    break
            except:
                print("Retrying download from " + path_fastq)
                retries = retries - 1
                continue

rule uncompress_fastq_Esoli:
    input: "Esoli/fastq_files/{esoli}.fastq.gz"
    output: "Esoli/fastq_files/{esoli}.fastq"
    shell:
        "gunzip -k {input}"

rule uncompress_fastq_Abau:
    input:
        A="Abau/fastq_files/{abau}_1.fastq.gz",
        B= "Abau/fastq_files/{abau}_2.fastq.gz"
    output:
        "Abau/fastq_files/{abau}_1.fastq",
        "Abau/fastq_files/{abau}_2.fastq"
    shell:
        """
        gunzip -k {input.A}
        gunzip -k {input.B}
        """

rule download_genomes:
    output:
        "Abau/genome/GCF_009759685.1_ASM975968v1_genomic.gbff",
        "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic.gbff",
    shell:
        '''
        wget -c --tries=50 -O Abau/genome/GCF_009759685.1_ASM975968v1_genomic.gbff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/759/685/GCF_009759685.1_ASM975968v1/GCF_009759685.1_ASM975968v1_genomic.gbff.gz
        wget -c --tries=50 -O Esoli/genome/GCF_000224675.1_ASM22467v1_genomic.gbff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/675/GCF_000224675.1_ASM22467v1/GCF_000224675.1_ASM22467v1_genomic.gbff.gz
        gunzip -k Esoli/genome/*.gz
        gunzip -k Abau/genome/*.gz
        '''
rule parse_gbff_transcripts_Abau:
    input:
        "Abau/genome/GCF_009759685.1_ASM975968v1_genomic.gbff",
    output:
        "Abau/genome/GCF_009759685.1_ASM975968v1_genomic_transcripts.fasta",
    shell:
        '''
        python src/parse_gbff_transcripts.py {input}
        '''

rule parse_gbff_transcripts_Esoli:
    input:
        "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic.gbff"
    output:
        "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic_transcripts.fasta",
    shell:
        '''
        python src/parse_gbff_transcripts.py {input}
        '''
rule index_transcriptome_Esoli:
    input: "Esoli/genome/GCF_000224675.1_ASM22467v1_genomic_transcripts.fasta",
    output: directory("Esoli/salmon_transcriptome_index")
    shell: """/home/tslaird/programs/salmon-latest_linux_x86_64/bin/salmon index -t {input} -i {output} -k 31"""

rule index_transcriptome_Abau:
    input: "Abau/genome/GCF_009759685.1_ASM975968v1_genomic_transcripts.fasta"
    output: directory("Abau/salmon_transcriptome_index")
    shell: """/home/tslaird/programs/salmon-latest_linux_x86_64/bin/salmon index -t {input} -i {output} -k 31"""

rule map_reads_Esoli:
    input:
    	read="Esoli/fastq_files/{esoli}.fastq",
        transcriptome="Esoli/salmon_transcriptome_index"
    params:
        outdir="Esoli/{esoli}_quant"
    output: "Esoli/{esoli}_quant/quant.sf"
    shell: """/home/tslaird/programs/salmon-latest_linux_x86_64/bin/salmon quant -i {input.transcriptome} -l A -r {input.read} -o {params.outdir} """

rule map_reads_Abau:
    input:
    	read1="Abau/fastq_files/{abau}_1.fastq",
    	read2="Abau/fastq_files/{abau}_2.fastq",
        transcriptome="Abau/salmon_transcriptome_index"
    params:
        outdir="Abau/{abau}_quant"
    output: "Abau/{abau}_quant/quant.sf"
    shell: """/home/tslaird/programs/salmon-latest_linux_x86_64/bin/salmon quant -i {input.transcriptome} -l A -1 {input.read1} -2 {input.read2} -o {params.outdir} """

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
