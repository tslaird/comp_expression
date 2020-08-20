#import libraries needed

library("tximport")
library(tidyr)
library(DESeq2)
library(stringr)

#import data
ESOLI_SAMPLES= paste0(c("SRR7120554","SRR7120555","SRR7120556","SRR7120557","SRR7120558","SRR7120559"),"_quant")
ABAU_SAMPLES= paste0(c("SRR6202159","SRR6202160"), "_quant")

Esoli_files<- file.path('/home/tslaird/leveau_lab/comp_expression/Esoli', ESOLI_SAMPLES, "quant.sf" )
Abau_files<- file.path('/home/tslaird/leveau_lab/comp_expression/Abau', ABAU_SAMPLES, "quant.sf" )

txi_Esoli <- tximport(Esoli_files, type = "salmon", txOut = T)
txi_Abau<- tximport(Abau_files,type="salmon", txOut = T)

#analyze data Abau
tpm_Abau<-as.data.frame(txi_Abau$counts)
colnames(tpm_Abau)<-ABAU_SAMPLES
ft<-apply(x,1, function(i) 
  fisher.test(rbind(c(as.numeric(i[2]),1e6-as.numeric(i[2])),
                    c(as.numeric(i[1]),1e6-as.numeric(i[1]) )), alternative = 'two.sided'))

tpm_Abau$Name <-rownames(tpm_Abau)
tpm_Abau<-separate(tpm_Abau,Name, into=c('assembly','locus','locus_tag','old_locus_tag','coords','product','protein_id','pseudogene'), sep = "!!", remove = TRUE,
         convert = FALSE, extra = "warn", fill = "warn")
tpm_Abau$p.value<-sapply(ft, function(i) i$p.value)
tpm_Abau$p.adj<-p.adjust(tpm_Abau$p.value,method = 'fdr')
tpm_Abau$FC<- (tpm_Abau$SRR6202159_quant /tpm_Abau$SRR6202160_quant)
tpm_Abau$log2FC<- log2(tpm_Abau$SRR6202159_quant)-log2(tpm_Abau$SRR6202160_quant)

#analyze data Esoli

sampleTable <- data.frame(condition = factor(rep(c("control", "IAA"), each = 3)))
rownames(sampleTable) <- colnames(txi_Esoli$counts)

desq2_Esoli <- DESeqDataSetFromTximport(txi_Esoli, sampleTable, ~condition)
desq2_Esoli  <- DESeq(desq2_Esoli)
desq2_Esoli_results<- results(desq2_Esoli )
desq2_Esoli_results_ordered <- desq2_Esoli_results[order(desq2_Esoli_results$pvalue),]
desq2_Esoli_results_ordered_df<-as.data.frame(desq2_Esoli_results_ordered)
desq2_Esoli_results_ordered_df$name<-rownames(desq2_Esoli_results_ordered_df)
desq2_Esoli_results_ordered_df<-desq2_Esoli_results_ordered_df %>% tidyr::separate(name, sep='!!',
                                                                                   c('assemby','accession',
                                                                                     'locus_tag','old_locus_tag',
                                                                                     'coords','protein','protein_id',"pseudogene"))

desq2_Esoli_final<-desq2_Esoli_results_ordered_df

#compare homologs
cdhit_out<-read.csv('/home/tslaird/leveau_lab/comp_expression/cdhit_data/all_proteins.fasta.cdhit90-80-70-60-50-40.clstr.id_matrix')
homologs<-cdhit_out[(!cdhit_out$GCF_000224675.1 =='') & ( !cdhit_out$GCF_009759685.1 =='' ),]


#threshold of .01 p.adj and > 1 l2fc
Abau_sig<- tpm_Abau[tpm_Abau$p.adj < 0.01 & abs(tpm_Abau$log2FC) >0.5,]
Esoli_sig<- desq2_Esoli_final[(desq2_Esoli_final$padj < 0.05) & (abs(desq2_Esoli_final$log2FoldChange) > 0.5),]

Esoli_h<-apply(homologs, 1, function(x) paste( sapply(strsplit(x[4],' '), function(y) ifelse(y %in% Esoli_sig$locus_tag, ifelse(Esoli_sig[Esoli_sig$locus_tag%in%y,]$log2FoldChange>1,"UP","DOWN"),'NO') )))
Abau_h<-apply(homologs, 1, function(x) paste( sapply(strsplit(x[5],' '), function(y) ifelse(y %in% Abau_sig$locus_tag, ifelse(Abau_sig[Abau_sig$locus_tag%in%y,]$log2FC >1,"UP","DOWN"),'NO') )))

homologs$Esoli_sig<-lapply(Esoli_h, function(x) paste(x[1:length(x)],collapse=" "))
homologs$Abau_sig<-lapply(Abau_h, function(x) paste(x[1:length(x)],collapse=" "))

t<-homologs[(sapply(homologs$Esoli_sig, function(x) str_detect(x,'UP|DOWN'))==TRUE) & (sapply(homologs$Abau_sig, function(x) str_detect(x,'UP|DOWN'))==TRUE),]
