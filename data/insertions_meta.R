
library(tidyverse)



plasmid_sites <- readRDS('data/plasmid_data/out1/sites.rds')
dog_sites <- readRDS('data/Aradhana_data/output_10Q/sites.rds')

all_samples <-c("GTSP4175","GTSP4176","GTSP4181","GTSP4182","UninfectedControl-210422-MN","NoTemplateControl-210422-MN")
sample_type <-c("plasmid","plasmid","dog","dog","control","control")
count_ins_dog <- sapply(all_samples,function(x){sum(dog_sites$sample==x)})
count_ins_plasmid <- sapply(all_samples,function(x){sum(plasmid_sites$sample==x)})
sample_desc <- as.data.frame(t(rbind(sample_type,count_ins_dog,count_ins_plasmid)))



all_sites <- rbind(plasmid_sites,dog_sites)

kk <- all_sites %>% mutate(insetion_type=as.factor(ifelse(startsWith(as.character(seqnames), 'chr'),'dog','plasmid')),
                           readID.list=NULL,
                           ltrRepSeq=NULL,
                           ltrRepSeq2=NULL)
all_sites_clean <- group_by(kk,sample,replicate) %>% mutate(relAbund=estAbund/sum(estAbund)) %>% ungroup() %>% arrange(sample,replicate)


all_sites_summ <- all_sites_clean %>% group_by(sample,replicate) %>% 
  summarise(plasmid_relAbund=sum(relAbund*ifelse(insetion_type=="plasmid",1,0)),dog_relAbund=sum(relAbund*ifelse(insetion_type=="dog",1,0))) %>%
  mutate(dog_sd=sd(dog_relAbund),plasmid_sd=sd(plasmid_relAbund)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0))
  ungroup()


png(height = 4, width = 6,units = 'in', res=300, file = 'data/plasmidplot.png')
ggplot(all_sites_summ, aes(x=sample, y=plasmid_relAbund)) + 
  geom_boxplot()
dev.off()

long_all_sites_summ <- all_sites_summ %>% gather(insertion_type,relAbund, plasmid_relAbund,dog_relAbund) %>%
  mutate(sd=ifelse(insertion_type=='plasmid_relAbund',plasmid_sd,dog_sd),plasmid_sd=NULL,dog_sd=NULL)

png(height = 4, width = 6,units = 'in', res=300, file = 'data/plasmid_bar.png')
ggplot(long_all_sites_summ, aes(fill=insertion_type, y=relAbund, x=sample)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar( aes(x=sample,ymin=relAbund-sd, ymax=relAbund+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)
dev.off()

png(height = 4, width = 6,units = 'in', res=300, file = 'data/plasmidplot_2.png')
ggplot(long_all_sites_summ, aes(x=sample, y=relAbund,fill=insertion_type)) + 
  geom_boxplot()
dev.off()

names(sample_type) <- all_samples
plasmid_sites %>% transmute(sample=sample,sample_type=sample_type[sample],replicate=replicate,position=paste0(position,strand),estAbund=estAbund) %>% arrange(position)

