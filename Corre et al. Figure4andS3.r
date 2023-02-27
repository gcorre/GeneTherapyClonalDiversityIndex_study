## Corre Guillaume
## GENETHON
## Sept 2022

## Figure 4

# Load libraries --------------------------------------------------------------------



library(tidyverse); 
library(ggpubr);
library(sonicLength)




# Load dataset of published papers --------------------------------------------------


files = list.files(path = "manuscript.V2/data/figure4_VISAresults_corre.et.al.2022/", pattern = "ISCollapsed2frag.bed", recursive = T,full.names = T)

my <- list()


for(file in files){
  x = read.delim(file,header = T,sep = "\t")
  my[[str_match(basename(file),pattern = "VISAs([[:alnum:]]+)")[,2]]] <- x
}

my.df <- bind_rows(my,.id = "sample")

rm(x);rm(files);rm(file)



my.df <- my.df %>% 
  mutate(id = paste(sample,X.chr,start,strand,sep = "_"))


# SonicLength estimation --------------------------------------------------


fit <- estAbund(locations = my.df$id,my.df$length)

# add sonic length quantification to data.frame
my.df.sum <- my.df %>% 
  group_by(id) %>% 
  summarise(nFrags = n(),
            nReads = sum(readCount)) %>%
  left_join(data.frame(sonic = (fit$theta)) %>%
              rownames_to_column("id"))

rm(fit)





# Attribute IS to known clones (Corre et al. 2022)


clones = data.frame(clone = c("KS10","KS39","KS39","KS40","KS40","KS40"),
                    IS = c("chr17_1719890_+",
                           "chr22_41230375_-",
                           "chr4_188092607_-",
                           "chr8_144812321_+" ,
                           "chr5_179151158_+",
                           "chr2_74442451_+"))


# Annotate Is to clones and calculate relative abundance of IS per sample

my.df.sum2 <- my.df.sum  %>%
  separate(id,into = c("sample","X.chr","start","strand"),sep = "_",remove = T) %>%
  mutate(IS = paste(X.chr,start,strand,sep="_")) %>% 
  left_join(clones, by = "IS") %>%
  mutate(cloneID = ifelse(is.na(clone), IS, paste(clone))) %>%
  group_by(sample) %>%
  mutate(prop = nFrags / sum(nFrags),
         propSonic = sonic / sum(sonic),
         mix = recode(sample, E = "F", F = "E")) %>% 
  ungroup %>%
  select(-sample)






# Group is with low abundance together for ploting

dat_plot = my.df.sum2 %>% 
  group_by(mix) %>%
  mutate(IS = case_when(propSonic>0.01~ IS,
                        TRUE ~ "LowAbundance"))%>%
  group_by(IS,mix) %>%
  summarise(propSonic = sum(propSonic)) %>%
  mutate(IS = factor(IS),
         IS = relevel(IS, "LowAbundance"))











# calculate diversity indices -------------------------------------------------------------


div = my.df.sum2 %>% 
  group_by(mix) %>% 
  summarise(nIS = n_distinct(IS),
            cells = sum(sonic),
            ShannonSonic = -sum((sonic/sum(sonic)) * log(sonic/sum(sonic))),
            pielou = ShannonSonic / log(nIS),
            simpson = vegan::diversity(sonic,index = "simpson"),
            gini = ineq::Gini(sonic)) %>%
  filter(mix != "B1") %>% 
  mutate(mix = recode(mix, B2 = "B"))



# calculate chao1

chao_data = my.df.sum2 %>% pivot_wider(id_cols = mix,names_from = IS,values_from = nFrags,values_fill = 0) %>% column_to_rownames("mix")
vegan::estimateR(x = chao_data)


theme_set(theme_bw() + theme(text = element_text(color = "black",size = 11), 
                             axis.text = element_text(color = "black",size=11),
                             axis.title = element_text(color = "black",size=11)))

x11(width = 3.4, height = 3);

ggplot(dat_plot %>% 
                             filter(mix != "B1") %>%
                             mutate(mix = recode(mix, B2 = "B")), 
                           aes(mix, y= propSonic, fill = IS)) + 
  geom_col(show.legend = F, col = "black") +
  scale_fill_manual(values = c("grey", colorRamps::matlab.like2(12)))+
  xlab(NULL) +
  ylab("Relative abundance") +
  theme(legend.position = "bottom") +
  
  #add pielou index
  geom_point(inherit.aes = F, data = div ,
             aes(mix,y=pielou, shape = "Pielou")) +
  geom_path(inherit.aes = F, data = div ,
            aes(mix,y=pielou,group = 1),lwd=.4) + 
  
  # add simpson index
  geom_point(inherit.aes = F, data = div,
             aes(mix,y=simpson, shape = "Simpson"))+ 
  geom_path(inherit.aes = F, data = div ,
            aes(mix,y=simpson, group = 1),lwd=.4) +
  
  
  scale_shape_manual(name="Diversity index",values = c(19,17)) +
  scale_y_continuous(labels = scales::percent_format(),sec.axis = sec_axis(trans = ~.,name = "Diversity")) + 
  guides(fill = "none", color = "none")


# make plots

ggsave(filename = "manuscript.V2/graphs/Figure_4.png", device = "png", width = 3.4, height = 3,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_4.svg", device = "svg",  width = 3.4, height = 3,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_4.ps", device = "ps",  width = 3.4, height = 3,units = "in",dpi = 72)








# Figure S3 -------------------------------------------------------------------------

# Make a simulation from real polyclonal cells abundance ------------------

# remove KS clones IS

my_polyclonal <- my.df %>% 
  select(id,length,readCount) %>% 
  separate(id, into = c("mix","IS"),sep = "_", extra = "merge") %>% 
  filter(mix == "A") %>% ## select only sample A
  select(-mix) %>%
  group_by(IS) %>%
  summarise(nFrags = n())

poly_count = sum(my_polyclonal$nFrags)



for(ndom in c(1,2,3,4,5)){

mix_fake = list()
for(p in c(seq(0.10,0.90,0.10))){
  
  dom = (p * poly_count) / (1-p)

  fake = data.frame(IS = paste("spike",1:ndom,sep = "_"),
                    clone = paste("spike",1:ndom,sep = "_"),
                    nFrags = rep(dom/ndom,times = ndom))
  
  mix_fake[[paste(p)]] <- bind_rows(my_polyclonal,fake)
  
}

mix_fake = bind_rows(mix_fake, .id = "dominance")%>%
  mutate(dominance = as.numeric(as.character(dominance))*100)


mix_fake_summary = mix_fake %>%
  group_by(dominance) %>%
  summarise(IS = n_distinct(IS),
            shannon = vegan::diversity(nFrags,index = "shannon"),
            pielou = shannon / log(IS),
            simpson = vegan::diversity(nFrags, index = "simpson"))




theme_set(theme_bw() + theme(text = element_text(color = "black",size = 11), 
                             axis.text = element_text(color = "black",size=11),
                             axis.title = element_text(color = "black",size=11)))



p1 = ggplot(mix_fake_summary, aes(factor(dominance))) + 
  geom_path(aes(y = simpson, group = 1,linetype = "simpson"))+
  geom_point(aes(y = simpson, group = 1))+
  geom_path(aes(y = pielou, group = 1, linetype = "pielou")) +
  geom_point(aes(y = pielou, group = 1))+
  scale_linetype_discrete(name = "Diversity index") +
  ylim(c(0,1)) +
  scale_x_discrete(breaks=seq(0,100,10)) +
  labs(y = "Diversity", x = "Dominance (%)")





p2 = ggplot(mix_fake %>%   
              group_by(clone, dominance) %>% 
              summarise(count = sum(nFrags)), 
            aes(factor(dominance), count, fill = clone)) +
  geom_col(show.legend = F, position = 'fill', col = "black")+
  scale_x_discrete(breaks = seq(0,100,10)) +
  scale_y_continuous() +
  labs(x = "Dominance (%)",y =" Abundance",title = paste(ndom,"dominant clones",sep=" ")) + 
  theme(text = element_text(color = "black",size = 11), 
        axis.text = element_text(color = "black",size=11),
        axis.title = element_text(color = "black",size=11))

#ploti <- p1 + patchwork::inset_element(p2,left = 0.01,bottom = 0.01, right = .50,top = .5, clip = T,on_top = T)

ploti <-  p2| p1

assign(paste("plot_",ndom,"clones",sep=""), ploti)
}


x11(width = 12,height = 14)

(plot_1clones / plot_2clones / plot_3clones / plot_4clones / plot_5clones)  +  patchwork::plot_layout(guides = "collect") +plot_annotation(tag_levels = "A")



ggsave(filename = paste("manuscript.V2/graphs/Figure_S3",".png",sep = ""), device = "png", width = 8, height = 10,units = "in",dpi = 72)
ggsave(filename = paste("manuscript.V2/graphs/Figure_S3",".svg",sep = ""), device = "svg", width = 8, height = 10,units = "in",dpi = 72)
ggsave(filename = paste("manuscript.V2/graphs/Figure_S3",".ps",sep = ""), device = "ps", width = 8, height = 10,units = "in",dpi = 72)


