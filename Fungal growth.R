# libraries
library(tidyverse)
library(data.table)
library(Amelia) 
library(seqinr) # import fasta data file
library(factoextra)
library(cluster)
library(stringr)
library(corrplot)
library(broom) # pca
library(EMT) # exact multinomial test
library(msgl) # group lasso
library(glmnet)
library(caret)
# load data set
tag <- fread('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 2 - Fungal growth\\Data\\CryptoWakeUpSampleSheetPlusTags.txt')
  
raw.gene <- fread('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 2 - Fungal growth\\Data\\CW-kallisto-abundance-foldchange-long-bygene.txt')

isoforms <- fread('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 2 - Fungal growth\\Data\\CW-kallisto1-abundance-long-bysample-byisoform.txt')

sequence <- read.fasta('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 2 - Fungal growth\\Data\\H99_allorfs_p500.fasta')

fivemers <- fread("C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 2 - Fungal growth\\Data\\H99_allorfs_promoter500nt_5mercounts.txt")

# EDA 1 part
## missing values
### count of missing values
sum(is.na(isoforms)) # there is no missing values in isoforms

sum(is.na(raw.gene))
colSums(is.na(raw.gene)) # there is some missing values in FoldChange column

missmap(raw.gene, main = 'missmap of features in data set 1')


# Data preprocessing
gene <- copy(raw.gene)
## remove all genes with missing values, zero, and _1
ind <- unique(gene[gene$FoldChange == 0 | is.na(gene$FoldChange) | str_detect(gene$Gene, '_1')]$Gene)
gene <- gene[!gene$Gene %in% ind]
## Systematic bias of replicates
gene$Replicate <- str_sub(gene$Code,-1)
gene$OtherCon <- str_sub(gene$Code,1,-2)
A <- gene[gene$Replicate == 'A']
B <- gene[gene$Replicate == 'B']
### check the order of data frame A and B
all(A$Gene == B$Gene) # It should be True
all(A$OtherCon == A$OtherCon) # It should be True
### plot of bias
y1 <- log(A$FoldChange)
y2 <- log(B$FoldChange)
ss = smooth.spline(y1, y2-y1, cv=TRUE)
plot(y1,y2-y1, xlab = 'replicate A',ylab = 'replicate B - replicate A', main = 'Before removing bias')
lines(ss, col=2)
abline(h=0, col=4,lwd=2)
### calculating the bias
bias <- c()
for(i in 1:length(y1)){
  bias[i] <- as.numeric(predict(ss, y1[i])$y)
}
y2.unbiased <- y2-bias  
### replace the foldchange of original data set with unbiased values
B$FoldChange <- y2.unbiased
A$FoldChange <- y1
### combine two data sets
gene.processed <- rbind(A,B) # NOTE: the foldchange of gene.processed
                             # is log(fold change)


# EDA 2 part
## distributions of gene expression 
### join tables together
gene.join <- merge(gene.processed, tag, by = 'Code')
gene.join.A <- merge(A, tag, by = 'Code') %>%
  select(Gene,Medium,Time,Temp,Rep, FoldChange)
gene.join.B <- merge(B, tag, by = 'Code') %>%
  select(Gene,Medium,Time,Temp,Rep, FoldChange)
### visualization of distributions
plot.A <- gene.join.A %>%
  ggplot(aes(FoldChange, colour = Medium)) + 
  geom_freqpoly(binwidth = 1) + 
  facet_grid(Temp ~ Time)
print(plot.A+labs(title= "Replicate: A"))

plot.B <- gene.join.B %>%
  ggplot(aes(FoldChange, colour = Medium)) + 
  geom_freqpoly(binwidth = 1) + 
  facet_grid(Temp ~ Time)
  
print(plot.B+labs(title= "Replicate: B"))


## Testing for Normality
hist(gene.processed$FoldChange,probability=T, main="Histogram of log(fold change)",xlab="log(fold change)", ylim = c(0,0.7), xlim = c(-6,6))
lines(density(gene.processed$FoldChange),col=2)

## Correlation across samples
### reshape data to wide format
gene.wide <- pivot_wider(gene.processed,id_cols = Gene, names_from = Code, values_from = FoldChange)
### calculate correlations
samples.corr <- gene.wide %>% 
  select(-Gene) %>% 
  cor(method = "pearson")
corrplot(samples.corr, 
         method = "circle", is.corr = FALSE) 


# Q1 
set.seed(20)
## K-means
### setup
gene.kmeans <- as.matrix(gene.wide[,-1]) %>%
  t() %>%
  scale() %>%
  t()
kmeans.combined <- matrix(nrow = nrow(gene.kmeans),ncol = (ncol(gene.kmeans)/2))
for(i in 1:18){
  kmeans.combined[,i] <- (gene.kmeans[,i]+gene.kmeans[,18+i])/2
}
### clustering
cluster.gene <- kmeans(kmeans.combined, 4, nstart = 20)
gene.kmeans <- data.table(kmeans.combined)
gene.kmeans$Gene <- gene.wide$Gene
gene.kmeans$cluster <- as.factor(cluster.gene$cluster)
names(gene.kmeans) <- c(unique(gene$OtherCon),'Gene','cluster')
### summary
kmeans.summary <- gene.kmeans %>%
  group_by(cluster) %>%
  summarise(mean.RC1 = mean(RC1), mean.RC2 = mean(RC2), mean.RC3 = mean(RC3), mean.RC4 = mean(RC4),
            mean.RH1 = mean(RH1),mean.RH2 = mean(RH2),mean.RH3 = mean(RH3),mean.RH4 = mean(RH4),
            mean.YC0 = mean(YC0),mean.YC1 = mean(YC1),mean.YC2 = mean(YC2),mean.YC3 = mean(YC3),
            mean.YC4 = mean(YC4),mean.YH1 = mean(YH1),mean.YH2 = mean(YH2),mean.YH3 = mean(YH3),
            mean.YH4 = mean(YH4),mean.YM9 = mean(YM9))

###Number of members in each cluster
table(gene.kmeans$cluster)

### plot
fviz_cluster(cluster.gene, data = kmeans.combined)

### choose the best number of clusters: The elbow method
fviz_nbclust(kmeans.combined, kmeans, method = "wss", k.max = 20) + theme_minimal() + ggtitle("the Elbow Method")

## Hierarchical Agglomerative
### creating distance matrix
d <- dist(kmeans.combined, method = "euclidean") 
### clustering
fit <- hclust(d, method="ward.D")

cut.tree <- cut(as.dendrogram(fit), h=40)
plot(cut.tree$upper, main="Upper tree of cut at h=40")

#cut tree and vasualization
abline(h = 300, col = "RED", lwd = 2)
sub.group <- cutree(fit, k=4) # cut tree into 4 clusters
rect.hclust(fit, k=4, border = 2:5)
#Number of members in each cluster
table(sub.group)
#the optimal number of clusters
fviz_nbclust(kmeans.combined, FUN = hcut, method = "wss")


# Q2
pca.combined <- gene.kmeans[,-c('cluster','Gene')]
## perform pca
pca <- prcomp(pca.combined)
## variance explained by PCs
pca.tidy <- tidy(pca, matrix = "pcs")
## Scree Plot
tidy(pca, matrix = "pcs") %>% 
  ggplot(aes(x = factor(PC))) +
  geom_col(aes(y = percent)) +
  geom_line(aes(y = cumulative, group = 1)) + 
  geom_point(aes(y = cumulative)) +
  labs(x = "Principal component", y = "Fraction variance explained")
## the influence of conditions on PCs
(top_conditions <- pca$rotation %>% 
  as_tibble(rownames = "conditions") %>%
  pivot_longer(cols = -conditions, names_to = 'PC', values_to = 'value') %>%
  filter(PC %in% c('PC1', 'PC2')) %>%
  group_by(PC) %>%
  arrange(desc(abs(value))) %>%
  slice(1:5) %>%
  pull(conditions) %>%
  unique())

con_loadings <- pca$rotation %>% 
  as_tibble(rownames = "conditions") %>% 
  filter(conditions %in% top_conditions)

(loadings_plot <- ggplot(con_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(data = con_loadings, 
            aes(x = PC1, y = PC2, label = conditions),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)))


# Q3
## Combine equivalent motifs
#setup 
pos <- c('A','T','C','G')
neg <- c('T','A','G','C')
kmers.name = colnames(fivemers.filtered)
fivemers.combined <- matrix(nrow = nrow(fivemers.filtered),ncol = 513)

kmers.name.combined <- vector(length = 512)
#combine these columns(kmers) matched
strreverse <- function(x){ 
  strsplit(x, NULL) %>% 
    lapply(rev) %>% 
    sapply(paste, collapse="") 
}
  

for (i in 2:513) {
  code = unlist(strsplit(kmers.name[i],split = ""))
  temp = character(5)
  for (j in 1:5) {
    ind <- which(code[j]==pos)
    temp[j] <- neg[ind]
  }
  temp = strreverse(paste(temp,collapse = '')) 
  match.ind = which(kmers.name == temp)
  fivemers.combined[,i] = (fivemers.filtered[[i]]+fivemers.filtered[[match.ind]])/2
  kmers.name.combined[i-1] = paste(kmers.name[i],temp,collapse = ',')
}
fivemers.combined <- data.table(fivemers.combined)
fivemers.combined[,1] <- fivemers.filtered[[1]]

names(fivemers.combined) <- c('Gene',kmers.name.combined)


## Data table
gene.names <- gene.processed$Gene
fivemers.filtered <- fivemers[fivemers$Gene %in% gene.names,] 
fivemers.grouped <- merge(fivemers.combined,gene.kmeans[,c(19,20)], by = 'Gene')

##employ multinomial logistic lasso regression
yy <- fivemers.grouped$cluster
lasso.multi <- cv.glmnet(x,yy, family='multinomial')
fit_lasso_multi <- glmnet(x,yy, lambda = lasso.multi$lambda.min,family = 'multinomial')

#plot
plot(lasso.multi)
par(mfrow = c(2,2))
plot(lasso.multi$glmnet.fit)
par(mfrow = c(1,1))

#over-or-under-represented
non1 <- colnames(x)[which(fit_lasso_multi$beta[[1]]!=0)]
non2 <- colnames(x)[which(fit_lasso_multi$beta[[2]]!=0)]
non3 <- colnames(x)[which(fit_lasso_multi$beta[[3]]!=0)]
non4 <- colnames(x)[which(fit_lasso_multi$beta[[4]]!=0)]

non.df <- data.frame(beta1 = vector(length = ncol(x)), beta2 = vector(length = ncol(x)),
                   beta3 = vector(length = ncol(x)), beta4 = vector(length = ncol(x)))

for(i in 1:4){
  non.df[,i] <- fit_lasso_multi$beta[[i]][,1]
}
#return the number of motifs whose coefficients are non-zero for different number of responses.
sum(rowSums(non.df != 0) == 1)
#for each cluster, calculate the total number of under-or-over-represented motifs
under.over <- function(i){
  ind.temp <- non.df[,i] !=0 
  ind.temp.2 <- rowSums(non.df[,-i]) == 0 
  ind.temp.under <- non.df[,i] > 0
  ind.temp.over <- non.df[,i] <0
  return(c(sum(ind.temp+ind.temp.2 == 2),sum(ind.temp.under+ind.temp.2 == 2),
           sum(ind.temp.over+ind.temp.2 == 2)))
}
for(i in 1:4){
  print(paste('cluster ',i,' :', under.over(i)[1],'(total), ',under.over(i)[2],
              '(under), ',under.over(i)[3],'(over)',sep = ''))
}

#subset of predictors
length(unique(c(non1,non2,non3,non4)))

# Q4

## produce 6-mer count tables from fasta-format file

### r kmer_functions
make_kmer_list <- function(k=2L, alphabet = c("A","C","G","T")) {
  
  list(alphabet) %>%
    
    base::rep(times = k) %>%
    
    purrr::set_names(LETTERS[1:k]) %>%
    
    base::expand.grid() %>%
    
    tidyr::unite("kmer",LETTERS[1:k],sep="") %>%
    
    dplyr::pull(kmer)
  
}

count_kmer_one <- function(string,kmer_list,k=2L) {
  
  stopifnot( nchar(kmer_list[1]) == k )
  
  
  
  kCount <- rep(0L,length(kmer_list)) %>% purrr::set_names(kmer_list)
  
  kminus1 <- k - 1L
  
  for(i in 1:( nchar(string) - kminus1 ) ) {
    
    kmer <- stringr::str_sub(string, start=i, end = i + kminus1)
    
    if(kmer %in% kmer_list) {
      
      kCount[kmer] <- kCount[kmer] + 1L
      
    }
    
  }
  
  return(kCount)
  
}


count_kmer_tibble <- function(stringset,kmer_list,k=2L) {
  
  lapply(stringset,count_kmer_one,kmer_list=kmer_list,k=k) %>%
    
    lapply(as.list) %>%
    
    dplyr::bind_rows()
  
}

### Count 6-mers for Cryptococcus promoters
promoterseqs <- Biostrings::readDNAStringSet('H99_allorfs_p500.fasta')



promoterseqids <- names(promoterseqs) %>%
  
  stringr::str_extract(pattern="\\w+")



all_6mers <- make_kmer_list(k=6L)



promoter_6mer_countsonly <- promoterseqs %>%
  
  as.character() %>%
  
  count_kmer_tibble(kmer_list = all_6mers, k=6L)



promoter_6mer_counts <- bind_cols(tibble(Gene=promoterseqids),
                                  
                                  promoter_6mer_countsonly)


## combine equivalent motif
#setup 
sixmers.filtered <- promoter_6mer_counts[promoter_6mer_counts$Gene %in% gene.names,] 
sixmers.name = colnames(sixmers.filtered)

sixmers.combined <- matrix(nrow = nrow(sixmers.filtered),ncol = 2049)

sixmers.name.combined <- vector(length = 2048)
#combine these columns(6mers) matched
for (i in 2:2049) {
  code = unlist(strsplit(sixmers.name[i],split = ""))
  temp = character(6)
  for (j in 1:6) {
    ind <- which(code[j]==pos)
    temp[j] <- neg[ind]
  }
  temp = strreverse(paste(temp,collapse = '')) 
  match.ind = which(sixmers.name == temp)
  sixmers.combined[,i] = (sixmers.filtered[[i]]+sixmers.filtered[[match.ind]])/2
  sixmers.name.combined[i-1] = paste(sixmers.name[i],temp,collapse = ',')
}
sixmers.combined <- data.table(sixmers.combined)
sixmers.combined[,1] <- sixmers.filtered[[1]]

names(sixmers.combined) <- c('Gene',sixmers.name.combined)
sixmers.grouped <- merge(sixmers.combined,gene.kmeans[,c(19,20)], by = 'Gene')

## Lasso regression on 6mers data
yyy <- sixmers.grouped$cluster
xx <- as.matrix(sixmers.combined[,-1])
lambda.seq <- seq(exp(-6),exp(-2),by=0.001)
lasso.multi.4 <- cv.glmnet(xx,yyy, lambda = lambda.seq,family='multinomial')
fit_lasso_multi.4 <- glmnet(xx,yyy, lambda = lasso.multi.4$lambda.min,family = 'multinomial')

#plot
plot(lasso.multi.4)
par(mfrow = c(2,2))
plot(lasso.multi$glmnet.fit)
par(mfrow=c(1,1))
non1.4 <- colnames(x)[which(fit_lasso_multi.4$beta[[1]]!=0)]
non2.4 <- colnames(x)[which(fit_lasso_multi.4$beta[[2]]!=0)]
non3.4 <- colnames(x)[which(fit_lasso_multi.4$beta[[3]]!=0)]
non4.4 <- colnames(x)[which(fit_lasso_multi.4$beta[[4]]!=0)]

non.df.4 <- data.frame(beta1 = vector(length = ncol(xx)), beta2 = vector(length = ncol(xx)),
                   beta3 = vector(length = ncol(xx)), beta4 = vector(length = ncol(xx)))

for(i in 1:4){
  non.df.4[,i] <- fit_lasso_multi.4$beta[[i]][,1]
}

sum(rowSums(non.df.4 != 0) == 1)