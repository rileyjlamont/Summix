library(dplyr)
mergeframe <- fread(paste0(my_directory, 'SNP_dataframe_f87.text.gz'), header = TRUE)
mergeframe=data.frame(mergeframe)

#create AFR network 
afr_sum <- summix(mergeframe, reference = c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri"), observed = "AFafr")
results_afr <- t(afr_sum[5:54])

ancs <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")
conts <- c("AFR", "AFR", "EUR", "SAS", "EAS", "EAS", "EUR", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "EUR", "EUR", "EUR", "SAS", "AFR", "EAS", "EAS", "EUR", "EUR", "SAS", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "AFR", "IAM", "EAS", "EAS", "AFR", "EAS", "EUR", "EAS", "IAM", "SAS", "EUR", "EAS", "SAS", "EUR", "EAS", "EAS", "EUR", "EAS", "EAS", "AFR", "AFR")
sizes_afr = as.vector(as.numeric(results_afr))
dat_afr <- cbind.data.frame("props" =sizes_afr, "ancs" = rownames(results_afr), "conts" = conts)
rownames(dat_afr) <- ancs

#nonzero
non_zero_sizes_afr = dat_afr[dat_afr$props>=.0000001,]
sumafr_AFR <- round(100*sum(non_zero_sizes_afr[non_zero_sizes_afr$conts == 'AFR',]$props), 2)
sumafr_EUR <- round(100*sum(non_zero_sizes_afr[non_zero_sizes_afr$conts == 'EUR',]$props), 2)
sumafr_EAS <- round(100*sum(non_zero_sizes_afr[non_zero_sizes_afr$conts == 'EAS',]$props), 2)
sumafr_SAS <- round(100*sum(non_zero_sizes_afr[non_zero_sizes_afr$conts == 'SAS',]$props), 2)
sumafr_IAM <- round(100*sum(non_zero_sizes_afr[non_zero_sizes_afr$conts == 'IAM',]$props), 2)
contafr_sum = cbind.data.frame('prop'= rbind(sumafr_AFR, sumafr_EAS, sumafr_EUR, sumafr_IAM, sumafr_SAS), 'cont'= rbind('AFR', 'EAS', 'EUR', 'IAM', 'SAS'))

#filter for fine_scale ancestries with at least 1% for network plot                               
dat_afr_2 = dat_afr[dat_afr$props>=.01,]

afr_2 <- cut_fst %>%
  dplyr::filter(row.names(cut_fst) %in% rownames(dat_afr_2))
afr_2 = data.frame(afr_2)
fst_filt_afr <- afr_2 %>% dplyr::select(rownames(dat_afr_2))                 
fst_filt_afr <- as.matrix(fst_filt_afr)

write.csv(dat_afr_2, file=paste("~/Downloads/", "afr_fs_dat", ".csv", sep=""))   #export as csv


#create EAS network 
eas_sum <- summix(mergeframe, reference = c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri"), observed = "AFeas")
results_eas <- t(eas_sum[5:54])

ancs <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")
conts <- c("AFR", "AFR", "EUR", "SAS", "EAS", "EAS", "EUR", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "EUR", "EUR", "EUR", "SAS", "AFR", "EAS", "EAS", "EUR", "EUR", "SAS", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "AFR", "IAM", "EAS", "EAS", "AFR", "EAS", "EUR", "EAS", "IAM", "SAS", "EUR", "EAS", "SAS", "EUR", "EAS", "EAS", "EUR", "EAS", "EAS", "AFR", "AFR")
sizes_eas = as.vector(as.numeric(results_eas))
dat_eas <- cbind.data.frame("props" =sizes_eas, "ancs" = rownames(results_eas), "conts" = conts)
rownames(dat_eas) <- ancs

non_zero_sizes_eas = dat_eas[dat_eas$props>=.0000001,]
sumeas_AFR <- round(100*sum(non_zero_sizes_eas[non_zero_sizes_eas$conts == 'AFR',]$props), 2)
sumeas_EUR <- round(100*sum(non_zero_sizes_eas[non_zero_sizes_eas$conts == 'EUR',]$props), 2)
sumeas_EAS <- round(100*sum(non_zero_sizes_eas[non_zero_sizes_eas$conts == 'EAS',]$props), 2)
sumeas_SAS <- round(100*sum(non_zero_sizes_eas[non_zero_sizes_eas$conts == 'SAS',]$props), 2)
sumeas_IAM <- round(100*sum(non_zero_sizes_eas[non_zero_sizes_eas$conts == 'IAM',]$props), 2)
conteas_sum = cbind.data.frame('prop'= rbind(sumeas_AFR, sumeas_EAS, sumeas_EUR, sumeas_IAM, sumeas_SAS), 'cont'= rbind('AFR', 'EAS', 'EUR', 'IAM', 'SAS'))
#filter for fine_scale ancestries with at least 1% for network plot                               
dat_eas_2 = dat_eas[dat_eas$props>=.01,]

eas_2 <- cut_fst %>%
  dplyr::filter(row.names(cut_fst) %in% rownames(dat_eas_2))
eas_2 = data.frame(eas_2)
fst_filt_eas <- eas_2 %>% dplyr::select(rownames(dat_eas_2))                 
fst_filt_eas <- as.matrix(fst_filt_eas)

write.csv(dat_eas_2, file=paste("~/Downloads/", "eas_fs_dat", ".csv", sep=""))   #export as csv


#create NFE
eur_sum <- summix(mergeframe, reference = c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri"), observed = "AFnfe")
results_eur <- t(eur_sum[5:54])

ancs <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")
conts <- c("AFR", "AFR", "EUR", "SAS", "EAS", "EAS", "EUR", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "EUR", "EUR", "EUR", "SAS", "AFR", "EAS", "EAS", "EUR", "EUR", "SAS", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "AFR", "IAM", "EAS", "EAS", "AFR", "EAS", "EUR", "EAS", "IAM", "SAS", "EUR", "EAS", "SAS", "EUR", "EAS", "EAS", "EUR", "EAS", "EAS", "AFR", "AFR")
sizes_eur = as.vector(as.numeric(results_eur))
dat_eur <- cbind.data.frame("props" =sizes_eur, "ancs" = rownames(results_eur), "conts" = conts)
rownames(dat_eur) <- ancs

non_zero_sizes_eur = dat_eur[dat_eur$props>=.0000001,]
sumeur_AFR <- round(100*sum(non_zero_sizes_eur[non_zero_sizes_eur$conts == 'AFR',]$props), 2)
sumeur_EUR <- round(100*sum(non_zero_sizes_eur[non_zero_sizes_eur$conts == 'EUR',]$props), 2)
sumeur_EAS <- round(100*sum(non_zero_sizes_eur[non_zero_sizes_eur$conts == 'EAS',]$props), 2)
sumeur_SAS <- round(100*sum(non_zero_sizes_eur[non_zero_sizes_eur$conts == 'SAS',]$props), 2)
sumeur_IAM <- round(100*sum(non_zero_sizes_eur[non_zero_sizes_eur$conts == 'IAM',]$props), 2)
conteur_sum = cbind.data.frame('prop'= rbind(sumeur_AFR, sumeur_EAS, sumeur_EUR, sumeur_IAM, sumeur_SAS), 'cont'= rbind('AFR', 'EAS', 'EUR', 'IAM', 'SAS'))

dat_eur_2 = dat_eur[dat_eur$props>=.01,]

eur_2 <- cut_fst %>%
  dplyr::filter(row.names(cut_fst) %in% rownames(dat_eur_2))
eur_2 = data.frame(eur_2)
fst_filt_eur <- eur_2 %>% dplyr::select(rownames(dat_eur_2))                 
fst_filt_eur <- as.matrix(fst_filt_eur)

write.csv(dat_eur_2, file=paste("~/Downloads/", "eur_fs_dat", ".csv", sep=""))   #export as csv


#create AMR
amr_sum <- summix(mergeframe, reference = c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri"), observed = "AFamr")
results_amr <- t(amr_sum[5:54])

ancs <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")
conts <- c("AFR", "AFR", "EUR", "SAS", "EAS", "EAS", "EUR", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "EUR", "EUR", "EUR", "SAS", "AFR", "EAS", "EAS", "EUR", "EUR", "SAS", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "AFR", "IAM", "EAS", "EAS", "AFR", "EAS", "EUR", "EAS", "IAM", "SAS", "EUR", "EAS", "SAS", "EUR", "EAS", "EAS", "EUR", "EAS", "EAS", "AFR", "AFR")
sizes_amr = as.vector(as.numeric(results_amr))
dat_amr <- cbind.data.frame("props" =sizes_amr, "ancs" = rownames(results_amr), "conts" = conts)
rownames(dat_amr) <- ancs

non_zero_sizes_amr = dat_amr[dat_amr$props>=.0000001,]
sumamr_AFR <- round(100*sum(non_zero_sizes_amr[non_zero_sizes_amr$conts == 'AFR',]$props), 2)
sumamr_EUR <- round(100*sum(non_zero_sizes_amr[non_zero_sizes_amr$conts == 'EUR',]$props), 2)
sumamr_EAS <- round(100*sum(non_zero_sizes_amr[non_zero_sizes_amr$conts == 'EAS',]$props), 2)
sumamr_SAS <- round(100*sum(non_zero_sizes_amr[non_zero_sizes_amr$conts == 'SAS',]$props), 2)
sumamr_AMR <- round(100*sum(non_zero_sizes_amr[non_zero_sizes_amr$conts == 'IAM',]$props), 2)
contamr_sum = cbind.data.frame('prop'= rbind(sumamr_AFR, sumamr_EAS, sumamr_EUR, sumamr_AMR, sumamr_SAS), 'cont'= rbind('AFR', 'EAS', 'EUR', 'IAM', 'SAS'))


dat_amr_2 = dat_amr[dat_amr$props>=.01,]

amr_2 <- cut_fst %>%
  dplyr::filter(row.names(cut_fst) %in% rownames(dat_amr_2))
amr_2 = data.frame(amr_2)
fst_filt_amr <- amr_2 %>% dplyr::select(rownames(dat_amr_2))                 
fst_filt_amr <- as.matrix(fst_filt_amr)

write.csv(dat_iam_2, file=paste("~/Downloads/", "iam_fs_dat", ".csv", sep=""))   #export as csv


#create sas
sas_sum <- summix(mergeframe, reference = c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri"), observed = "AFsas")
results_sas <- t(sas_sum[5:54])

ancs <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")
conts <- c("AFR", "AFR", "EUR", "SAS", "EAS", "EAS", "EUR", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "EUR", "EUR", "EUR", "SAS", "AFR", "EAS", "EAS", "EUR", "EUR", "SAS", "EAS", "EAS", "IAM", "EAS", "EAS", "AFR", "AFR", "IAM", "EAS", "EAS", "AFR", "EAS", "EUR", "EAS", "IAM", "SAS", "EUR", "EAS", "SAS", "EUR", "EAS", "EAS", "EUR", "EAS", "EAS", "AFR", "AFR")
sizes_sas = as.vector(as.numeric(results_sas))
dat_sas <- cbind.data.frame("props" =sizes_sas, "ancs" = rownames(results_sas), "conts" = conts)
rownames(dat_sas) <- ancs

non_zero_sizes_sas = dat_sas[dat_sas$props>=.0000001,]
sumsas_AFR <- round(100*sum(non_zero_sizes_sas[non_zero_sizes_sas$conts == 'AFR',]$props), 2)
sumsas_EUR <- round(100*sum(non_zero_sizes_sas[non_zero_sizes_sas$conts == 'EUR',]$props), 2)
sumsas_EAS <- round(100*sum(non_zero_sizes_sas[non_zero_sizes_sas$conts == 'EAS',]$props), 2)
sumsas_SAS <- round(100*sum(non_zero_sizes_sas[non_zero_sizes_sas$conts == 'SAS',]$props), 2)
sumsas_IAM <- round(100*sum(non_zero_sizes_sas[non_zero_sizes_sas$conts == 'IAM',]$props), 2)
contsas_sum = cbind.data.frame('prop'= rbind(sumsas_AFR, sumsas_EAS, sumsas_EUR, sumsas_IAM, sumsas_SAS), 'cont'= rbind('AFR', 'EAS', 'EUR', 'IAM', 'SAS'))


dat_sas_2 = dat_sas[dat_sas$props>=.01,]

sas_2 <- cut_fst %>%
  dplyr::filter(row.names(cut_fst) %in% rownames(dat_sas_2))
sas_2 = data.frame(sas_2)
fst_filt_sas <- sas_2 %>% dplyr::select(rownames(dat_sas_2))                 
fst_filt_sas <- as.matrix(fst_filt_sas)

write.csv(dat_sas_2, file=paste("~/Downloads/", "sas_fs_dat", ".csv", sep=""))   #export as csv


##SET network plot threshold for fst between fine-scale ancestries 
fst_threshold = .01

#filter fst dataframe for appropriate fst threshold 
cut_fst <- as.data.frame(apply(fst_dat, 2, function(x) ifelse(x > fst_threshold, 0, x)))


##### STATIC GRAPHS#install igraph
install.packages("igraph")
library(igraph)
library(RColorBrewer)
library(ggplot2)

############### GRAPH 1
my_directory = '~/Downloads/'

#read in data
cont_nodes <- read.csv(paste0(my_directory,"cont_nodes.csv"), header=T, as.is=T)
cont_links <- read.csv(paste0(my_directory,"cont_edges.csv"), header=T, as.is=T)
cont_links <- cont_links[order(cont_links$from, cont_links$to),]
rownames(cont_links) <- NULL

cont_net <- graph_from_data_frame(d=cont_links, vertices=cont_nodes, directed=T) 

#set attributes
V(cont_net)$size <- V(cont_net)$anc.size*30
V(cont_net)$label.color <- "black"
E(cont_net)$width <- E(cont_net)$weight/1.2
pal <- c(AFR = "#66C2A5", EAS = "#FC8D62", EUR ="#8DA0CB", IAM = "#E78AC3", SAS = "#A6D854")

#pal <- brewer.pal(length(unique(V(cont_net)$anc)), "Set2")

#plot
plot(cont_net, edge.arrow.mode=0, vertex.label = cont_nodes$anc, vertex.color = pal[as.numeric(as.factor(vertex_attr(cont_net, "anc")))])


################ GRAPH 2

barchart_dat <- read.csv(paste0(my_directory, "data_cont_barchart.csv"))
ggplot(barchart_dat, aes(fill=reference, y=value, x=ancestry)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette="Set2") + xlab("Continental Ancestry") +
  ylab("Percentage") 

################ GRAPH 3

#code for AFR, EAS, NFE, AMR, SAS networks


#AFR data
afr_node <- read.csv(paste0(my_directory,"afr_fs_dat.csv"), header=T, as.is=T)
afr_link <- read.csv(paste0(my_directory, "afr_link_dat.csv"), header=T, as.is=T)
afr_link <- afr_link[order(afr_link$from, afr_link$to),]
rownames(afr_link) <- NULL

afr_net <- graph_from_data_frame(d=afr_link, vertices=afr_node, directed=T) 

V(afr_net)$color = pal[V(afr_net)$conts]
V(afr_net)$size <- V(afr_net)$props*100
V(afr_net)$label.color <- "black"
E(afr_net)$width <- 10- 1000*E(afr_net)$weight

#plot
plot(afr_net, edge.arrow.mode=0, vertex.label = afr_node$X,vertex.label.dist=2.6)
ig2ggplot(afr_net)
#EAS data
eas_node <- read.csv(paste0(my_directory,"eas_fs_dat.csv"), header=T, as.is=T)
eas_link <- read.csv(paste0(my_directory, "eas_link_dat.csv"), header=T, as.is=T)
eas_link <- eas_link[order(eas_link$from, eas_link$to),]
rownames(eas_link) <- NULL

eas_net <- graph_from_data_frame(d=eas_link, vertices=eas_node, directed=T) 

V(eas_net)$color = pal[V(eas_net)$conts]
V(eas_net)$size <- V(eas_net)$props*100
V(eas_net)$label.color <- "black"
E(eas_net)$width <- 10- 1000*E(eas_net)$weight

plot(eas_net, edge.arrow.mode=0, vertex.label = eas_node$X,vertex.label.dist=2.6)

#EUR data
eur_node <- read.csv(paste0(my_directory,"eur_fs_dat.csv"), header=T, as.is=T)
eur_link <- read.csv(paste0(my_directory, "eur_link_dat.csv"), header=T, as.is=T)
eur_link <- eur_link[order(eur_link$from, eur_link$to),]
rownames(eur_link) <- NULL

eur_net <- graph_from_data_frame(d=eur_link, vertices=eur_node, directed=T) 

V(eur_net)$color = pal[V(eur_net)$conts]
V(eur_net)$size <- V(eur_net)$props*100
V(eur_net)$label.color <- "black"
E(eur_net)$width <- 10- 1000*E(eur_net)$weight

plot(eur_net, edge.arrow.mode=0, vertex.label = eur_node$X,vertex.label.dist=2.6)

#IAM data  
iam_node <- read.csv(paste0(my_directory,"iam_fs_dat.csv"), header=T, as.is=T)
iam_link <- read.csv(paste0(my_directory, "iam_link_dat.csv"), header=T, as.is=T)
iam_link <- iam_link[order(iam_link$from, iam_link$to),]
rownames(iam_link) <- NULL

iam_net <- graph_from_data_frame(d=iam_link, vertices=iam_node, directed=T) 

V(iam_net)$color = pal[V(iam_net)$conts]
V(iam_net)$size <- V(iam_net)$props*100
V(iam_net)$label.color <- "black"
E(iam_net)$width <- 10- 1000*E(iam_net)$weight

plot(iam_net, edge.arrow.mode=0, vertex.label = iam_node$X,vertex.label.dist=2.6)

#SAS data
sas_node <- read.csv(paste0(my_directory,"sas_fs_dat.csv"), header=T, as.is=T)
sas_link <- read.csv(paste0(my_directory, "sas_link_dat.csv"), header=T, as.is=T)
sas_link <- sas_link[order(sas_link$from, sas_link$to),]
rownames(sas_link) <- NULL

sas_net <- graph_from_data_frame(d=sas_link, vertices=sas_node, directed=T) 

#test before adding anc_slider into app code
sas_link_trim <- sas_link %>%
  dplyr::filter(weight <= .007)

sas_node_trim <- sas_node %>%
  dplyr::filter(props >= .1)

res <- subset(sas_link_trim, !(to %in% sas_node_trim$X))
df = sas_link_trim

for (i in 1:nrow(res)){
  name_from = res[i,1]
  name_to = res[i, 2]
  
  #remove name_from from both columns
  df0 <- df[!grepl(name_from, df$from),]
  df1 <- df0[!grepl(name_from, df0$to),]
  #remove name_to from both columns
  df2 <- df1[!grepl(name_to, df1$from),]
  df3 <- df2[!grepl(name_to, df2$to),]
  
  df = df3
}

#combine

#remove duplicates
sas_link_trim_2 <- sas_link_trim_2[!duplicated(sas_link_trim_2),]

# 
# datalist = vector("list")
# 
# for (i in 1:nrow(sas_node_trim)){
#   sas_link_trim_2_from = subset(sas_link_trim, from == sas_node_trim[i, 1], 
#                               select = c(from, to, weight))
#   sas_link_trim_2_to = subset(sas_link_trim, to == sas_node_trim[i, 1],
#                               select = c(from, to, weight))
#   datalist[[i]] <- rbind(sas_link_trim_2_from, sas_link_trim_2_to)
# }
# 
sas_link_trim_2 = do.call(rbind, datalist)

#remove duplicates
sas_link_trim_2 <- sas_link_trim_2[!duplicated(sas_link_trim_2),]

sas_net <- graph_from_data_frame(d=sas_link_trim_2, vertices=sas_node_trim, directed=T) 


V(sas_net)$color = pal[V(sas_net)$conts]
V(sas_net)$size <- V(sas_net)$props*100
V(sas_net)$label.color <- "black"
E(sas_net)$width <- 10- 1000*E(sas_net)$weight

plot(sas_net, edge.arrow.mode=0, vertex.label = sas_node$X,vertex.label.dist=2.6)







#generic to use in interactive net
#select input
anc_node <- input_node
anc_link <- input_link
anc_link <- anc_link[order(anc_link$from, anc_link$to),]
rownames(anc_link) <- NULL

anc_net <- graph_from_data_frame(d=anc_link, vertices=anc_node, directed=T) 

V(anc_net)$color = pal[V(anc_net)$conts]
V(anc_net)$size <- V(anc_net)$props*100
V(anc_net)$label.color <- "black"
E(anc_net)$width <- 10- 1000*E(anc_net)$weight

#plot
plot(anc_net, edge.arrow.mode=0, vertex.label = anc_node$X,vertex.label.dist=2.6)
























