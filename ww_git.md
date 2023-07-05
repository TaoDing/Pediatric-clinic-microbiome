### 1 16s Data

##### 16s Data quality control**

```bash
 cat *_1.fq.gz > 1.fg.gz
 cat *_2.fq.gz > 2.fq.gz
 mkdir qc
 fastqc -t 2 1.fg.gz 2.fq.gz -o qc
```

#####  qiime** 

```bash
source activate qiime2-2019.7

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.csv --input-format PairedEndFastqManifestPhred64 --output-path paired-end-demux64.qza

qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux64.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 200 --p-trunc-len-r 200 --o-representative-sequences double-end-rep-seqs.qza --o-table double-end-table.qza --o-denoising-stats double-end-stat.qza

qiime tools export --input-path double-end-table.qza --output-path exported-feature-table

cd exported-feature-table
biom convert -i feature-table.biom -o otu_table.tsv --to-tsv
```

#####  **\## taxonomy**

```bash
qiime feature-classifier classify-sklearn --i-classifier silva-132-99-515-806-nb-classifier.qza --i-reads double-end-rep-seqs.qza --o-classification taxonomy-paired-end.qza

qiime tools export --input-path taxonomy-paired-end.qza --output-path  exported-taxonomy-paired-end-table
```

##### **## relative abundance of genus level——Rstudio****

```bash
otu_genus<-aggregate(otu_taxonomy[,8:38],by=list(genus=otu_taxonomy$genus),FUN=sum)
write.table(otu_genus,"otu_genus.txt",sep='\t',quote=F,row.names = F)

otu_genus<-read.table("otu_genus.txt",row.names=1,header = T,sep='\t')
otu_genusP = t(t(otu_genus)/colSums(otu_genus))
colSums(otu_genusP)
otu_genusP1<-as.data.frame(otu_genusP)

otu_genusP1$sum <- rowSums(otu_genusP1)
otu_genusP1<- otu_genusP1[order(otu_genusP1$sum, decreasing = TRUE), ]
otu_genusP1_top10 <- otu_genusP1[1:12, -ncol(otu_genusP1)]
otu_genusP1_top10['Others', ] <- 1 - colSums(otu_genusP1_top10)
write.csv(otu_genusP1_top10, 'otu_genusP1_top10.csv', quote = FALSE)
```

##### **##alpha diversity**

```bash
library(vegan)
library(picante)

alpha <- function(x, tree = NULL, base = exp(1)) { est <- estimateR(x) Richness <- est[1, ] Chao1 <- est[2, ] ACE <- est[4, ] Shannon <- diversity(x, index = 'shannon', base = base) Simpson <- diversity(x, index = 'simpson') Pielou <- Shannon / log(Richness, base) goods_coverage <- 1 - rowSums(x == 1) / rowSums(x) result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)  if (!is.null(tree)) { PD_whole_tree <- pd(x, tree, include.root = FALSE)[1] names(PD_whole_tree) <- 'PD_whole_tree' result <- cbind(result, PD_whole_tree)} result}

otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
tree <- read.tree('tree.nwk')

alpha_all <- alpha(otu, tree, base = 2)
write.csv(alpha_all, 'otu_alpha.csv', quote = FALSE)
```

##### **## correlation analysis**

```bash
library(psych)
cor.test(otu_speciesP1T_top10[,'bacterial1'], otu_speciesP1T_top10[,'Acidobacteria'], method = 'pearson')

corr_species_spearman <- corr.test(otu_speciesP1T_top10, method = 'spearman')
corr_species_clinical_spearman <- corr.test(otu_speciesP1T_top10,clinical indicators, method = 'spearman')
corr_species_spearman$r
corr_species_spearman$p  
corr_species_p <- corr_species_spearman$p

write.csv(corr_species_p, 'corr_species_p.csv', quote = FALSE)
```

### **## 2. Bacterial whole genome sequencing data quality control**

#####  quality control**

```bash
cat *.R1.fq.gz > 1.fg.gz
cat *.R2.fq.gz > 2.fq.gz
mkdir qc
fastqc -t 2 1.fg.gz 2.fq.gz -o qc

fastp -i ${i}.R1.fq.gz -I ${i}.R2.fq.gz -o clean${i}.R1.fq.gz -O clean${i}.R2.fq.gz
```

##### **## Spades**

```bash
spades.py -k 77,87,97,117 --careful --only-assembler -1 clean${i}.R1.fq.gz -2 clean${i}.R2.fq.gz -o ${i}.assembly -t 30
```

##### **## gtdb**

```bash
source activate /public/home/wangwan/.conda/envs/gtdbtk

gtdbtk classify_wf --genome_dir ./${i}/contigs.fasta --extension fa --out_dir ./gtdb_28 --cpus 40
```

##### **## FastANI**

```bash
fastANI --ql design.txt --rl design.txt -o fastANI_output2.txt --fragLen 1000 -t 10 --matrix
```

### **## 3 SourceTracker**

```bash
library(vegan) 
otu <- read.table('otu_table.txt', header=T, sep="\t", quote = "", row.names=1,comment.char="",stringsAsFactors = FALSE)
otu_ID_Flattening = as.data.frame(t(rrarefy(t(otu), min(colSums(otu)))))
write.table (otu_ID_Flattening, file ="otu_ID_Flattening.csv",sep =",", quote =FALSE)

metadata  <- read.table('group_table2.txt',sep='\t',header=TRUE,row.names=1,check=F,comment='')
otus <- read.table("otu_ID_Flattening.txt",sep='\t',header=TRUE,row.names=1,check=F,comment='')
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]

if(length(common.sample.ids) <= 3) {
    message <- paste(sprintf('Error: there are %d sample ids in common ',length(common.sample.ids)),
                    'between the metadata file and data table')
    stop(message)
}

train.ix <- which(metadata$SourceSink=='Source')
test.ix <- which(metadata$SourceSink=='Sink')
envs <- metadata$Env

source('SourceTracker.r')
alpha1 <- alpha2 <- 0.001
st <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth=1000)

results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)
```


