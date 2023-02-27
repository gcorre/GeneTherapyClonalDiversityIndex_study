## Corre Guillaume
## GENETHON
## Fev.2023

## Script for the simulation of polyclonal population of cells with increasing dominance over time
## Steps : 
##  - set the starting relative abundance of clones @ t0 using a known distribution (uniform, lognormal, gaussian, ...) 
##  - attribute a proliferation rate to each clone. Give a single clone a higher proliferation rate to mimick advantage.
##  - Simulate clonal evolution over n pseudo-timepoints. The proliferation rate is modified at each time point with a small gaussian noise
##  - At each step, determine the new abundance of each clone and store it in a matrix
##  - Simulate different sampling using the relative abundance at each time point as sampling probabilities
##  - estimate diversity values from each samples and make some cool plots


## THIS SCRIPT IS USED TO PRODUCE FIGURES 2,S1 AND S2


# Libraries ---------------------------------------------------------------

  library(vegan);
  library(ineq);
  library(tidyverse)

  theme_set(ggpubr::theme_pubr(base_size = 10)+theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")))
# generate a polyclonale population ---------------------------------------


  Nclones <- 1000
  
  nameClones <- replicate(Nclones, paste(sample(LETTERS,size = 6,replace = T),collapse = ""))
  
  any(duplicated(nameClones)) # check that all clones have a different name



# Set the relative abundance of each clone at T0

  #prop <- rpois(n = Nclones,lambda = .5) + 1
  #prop <- prop / sum(prop)
  
  prop <- 1/Nclones


  
  # set the initial proliferation rate
  prol_rate <- rnorm(Nclones,mean = 1,sd = 0.1)#0.1
  
  #prol_rate <- rlnorm(Nclones,meanlog = 0,sdlog = .3)
  
  # set a higher value for the last clone to mimick advantage (choose the last clone for reproducibility and easier data managment for ploting)
  prop_dominant <- case_when(Nclones >=10000 ~ 4, between(Nclones,1000,9999)~ 3.5, between(Nclones,1,999)~ 2.5)
  prol_rate[Nclones] <- mean(prol_rate)*prop_dominant
  



# Simulation --------------------------------------------------------------



  


## make a matrix with relative proportion of clones (rows) at each time point (columns)
  time = 20
  mat <- matrix(nrow = Nclones,ncol = time)
  
  mat[,1] <- prop # proportion at T0
  
  
  for(i in 2:(time)){
    prop <- prop + (prop * prol_rate)
    prop <- prop / sum(prop)
    
    prol_rate <- prol_rate + rnorm(Nclones,mean = 0,sd = 0.1) # add a small noise for next day
    prol_rate[which(prol_rate<0)] <- mean(prol_rate) # set negative proliferation rate to 0
    
    mat[,i] <- prop 
  }
  
  dimnames(mat)<- list(nameClones,1:time) 
  


  rm(list = c("prop","prol_rate","i"))
  

# Population analysis -----------------------------------------------------

  
  apply(mat,MARGIN = 2, sum) # check that total abundance is 1
  

  
  # compute population "true" diversity  
  
  pop_shannon = apply(mat,MARGIN = 2, function(x){
    -sum(x*log(x))
  }
  )
  
  
  
  pop_pielou = pop_shannon / log(Nclones)
  
  pop_simpson = apply(mat,MARGIN = 2, function(x){
    sum(x*x)
  }
  )
  
  pop_gini =  apply(mat, MARGIN = 2, function(x) {
    Gini(x)
  })
  
  
  pop_uc50 <- apply(mat, MARGIN = 2, function(x){
    y = cumsum(sort(x))
    length(which(y>0.5))
    
  })
  
  
  pop_diversity <- data.frame(time = 1:time, 
                              shannon = pop_shannon,
                              pielou = pop_pielou,
                              gini = pop_gini,
                              simpson = pop_simpson,
                              UC50 = pop_uc50)
  
  rm(pop_shannon); rm(pop_pielou);rm(pop_gini);rm(pop_simpson)
  





# sample different number of cells ----------------------------------------
  
  
  
s = Nclones * c(0.01,0.1,0.5,1,2,5,10)
rep = 20 # number of replicate sampling



res_array <- list("clones" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "shannon" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "shannonbc" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "pielou" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "uc50" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "gini" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "simpson" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "shannon_HT" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "S.chao" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)) ,
                  "S.ACE" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)),
                  "dominant_prop" = array(dim = c(length(s),(time),rep),dimnames = list(s, 1:time, 1:rep)))


clones_names_sampled = list()

# make the sampling and compute diversity values or intermediates

for(i in 1:rep){  # for each replicate
  cat(i,"\n")
  for(x in 1:length(s)){ # for each number of cells to sample
    for(y in 1:time){ # for each time point
      smp <- sample(nameClones, size = s[x], replace = T, prob = mat[,y])
      tb = table(smp)
      
      if(y%in%c(1,5,10,15,20)){
        clones_names_sampled[[paste(y,"time")]][[paste(s[x],"cells")]][[paste("rep",i)]] <- tb
      }
      
      uc50 = length(which(cumsum(sort(tb))>(s[x]/2)))
      gini = ineq(tb,type = "Gini")
      f1 = length(which(tb == 1))
      Cov = 1 - (f1/s[x])
      S.chao1 = estimateR(tb)["S.chao1"]
      S.ACE = estimateR(tb)["S.ACE"]
      pi <- tb/s[x]
      prop_dominant <- pi[nameClones[Nclones]]
      if(is.na(prop_dominant)) prop_dominant=0
      simpson = sum(pi*pi)
      shannon <-  -sum((pi * log(pi)))#/ log(length(pi))
      H.MillerMadow = shannon + ((length(pi)-1)/(2*s[x]))
      pielou = shannon / log(length(pi))
      shannon_HT = -sum((pi * log(pi))/ (1-(1-pi)^s[x]))
      bcshannon = -sum((Cov*pi * log(Cov*pi))/(1-(1-Cov*pi)^s[x])) # chao & shen 2003 # coverage and Horvitz-Thompson correction

      res_array$clones[x,y,i] <- length(tb)
      res_array$dominant_prop[x,y,i] <- prop_dominant
      res_array$shannon[x,y,i] <- shannon
      res_array$shannonbc[x,y,i] <- bcshannon 
      res_array$pielou[x,y,i] <- pielou
      res_array$uc50[x,y,i] <- uc50
      res_array$gini[x,y,i] <- gini
      res_array$simpson[x,y,i] <- simpson
      res_array$shannon_HT[x,y,i] <- shannon_HT
      res_array$S.chao[x,y,i] <- S.chao1
      res_array$S.ACE[x,y,i] <- S.ACE
      res_array$simpson[x,y,i] <- simpson
    }
  }
}


# summarize the results array over replicates for each parameter (time and sample size)
 
  avg_res = lapply(res_array, function(y) {
    apply(y,MARGIN = c(1,2), function(x) {mean(x,na.rm = T)})})
  
  sd_res = lapply(res_array, function(y) {
    apply(y,MARGIN = c(1,2), function(x) {sd(x,na.rm = T)})})
  


# convert to a dataframe
  z = lapply(avg_res, function(x) {x %>% 
      data.frame %>%
      rownames_to_column('Cells') %>% 
      gather(key = "pseudotime",value = "mean", -Cells)}) %>%
    bind_rows(.id = "variable") %>%
    left_join(lapply(sd_res, function(x) {x %>% 
        data.frame %>%
        rownames_to_column('Cells') %>% 
        gather(key = "pseudotime",value = "sd", -Cells)}) %>%
          bind_rows(.id = "variable")
    )
  
  z <- z %>% mutate(pseudotime = str_remove(pseudotime,"X"))
  
  

  
  
  
  
  
  
  
  
  # Graphics ----------------------------------------------------------------

q <- ggplot(mat[Nclones,] %>% data.frame(time = 1:time, abundance = .), 
            aes(time, abundance)) +
  geom_col(col = "white",fill = "black") + 
  ggpubr::theme_pubr(base_size = 12)+ 
  theme(axis.ticks = element_line(colour =  "black")) +
  labs(y = "Rel. abundance", x = "pseudotime") + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(legend.position = "none") 
  
  
  # top10 <- mat %>% 
  #   data.frame %>%
  #   rownames_to_column("clone") %>%
  #   pivot_longer(cols = starts_with("X"),names_to = "pseudotime", values_to = "abundance") %>%
  #   group_by(pseudotime) %>%
  #   slice_max(order_by = abundance,n = 10,with_ties = F) %>% arrange(desc(abundance),.by_group = T) %>%
  #   mutate(pseudotime = factor(pseudotime, levels = paste("X",1:20,sep="")), rank = row_number()) %>% arrange(rank)
  # 
  # top10 %>% pivot_wider(names_from = "rank", values_from = "clone",-abundance) %>% arrange(pseudotime)
  # 
  # ggplot(top10, aes(pseudotime, y=abundance, fill = factor(rank,levels = 1:10))) + geom_col(col = "black")
  # 

  
 p0 <-  ggplot(data = z %>% filter(variable=="dominant_prop"),
         aes(factor(pseudotime, levels = paste(1:time)),
             mean,
             col=factor(Cells,levels = s), 
             group = factor(Cells))) + 
    geom_path() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
    scale_x_discrete(breaks = c(1,5,10,15,20))+
    theme(axis.ticks = element_line(colour =  "black"))+
    labs(x= "pseudotime", col = "sample size", y = "Relative abundance") +
    scale_color_viridis_d() +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_path(data = mat[Nclones,] %>% data.frame(time = 1:time, abundance = .),aes(time, abundance), inherit.aes = F, lty = 2, lwd = 1, col ="black")

  
  

idx = "clones"
p1 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                             mean,
                                             col = factor(Cells,levels = s), 
                                             group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "Richness") +
  scale_color_viridis_d() +
  geom_hline(yintercept = Nclones, lty = 2, lwd = 1)



idx = "shannon"
p2 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "Shannon index (H)") +
  scale_color_viridis_d()  +
  geom_path(data=pop_diversity,aes(time,shannon),inherit.aes = F, lwd = 1, lty = 2)



idx = "pielou"
p3 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "Pielou index (J)") +
  scale_color_viridis_d() +
  geom_path(data=pop_diversity,aes(time,pielou),inherit.aes = F, lwd = 1, lty = 2)


idx = "uc50"
p4 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "UC50") +
  scale_color_viridis_d()  +
  geom_path(data=pop_diversity,aes(time,UC50),inherit.aes = F, lwd = 1, lty = 2)


idx = "gini"
p5 <-ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                             mean,
                                             col = factor(Cells,levels = s), 
                                             group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "Gini index") +
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  scale_color_viridis_d() +
  geom_path(data=pop_diversity,aes(time,gini),inherit.aes = F, lwd = 1, lty = 2)




idx = "shannonbc"
p6 <-ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                             mean,
                                             col = factor(Cells,levels = s), 
                                             group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "shannon corrected", y = "Shannon corrected index") +
  scale_color_viridis_d() +
  geom_path(data=pop_diversity,aes(time,shannon),inherit.aes = F, lwd = 1, lty = 2)






idx = "S.chao"
p7 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Chao1", y = "Richness (Chao1)") +
  scale_color_viridis_d()  +
  geom_hline(yintercept = Nclones, lty = 2, lwd = 1)






idx = "S.ACE"
p8 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "ACE", y = "Richness (ACE)") +
  scale_color_viridis_d()  +
  geom_hline(yintercept = Nclones, lty = 2, lwd = 1)



idx = "simpson"
p9 <- ggplot(z %>% filter(variable==idx), aes(factor(pseudotime, levels = paste(1:time)), 
                                              1-mean,
                                              col = factor(Cells,levels = s), 
                                              group = factor(Cells))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = 1-(mean - sd), ymax = 1-(mean + sd)), width = .4) +
  #scale_x_log10() + 
  scale_x_discrete(breaks = c(1,5,10,15,20))+
  theme(axis.ticks = element_line(colour =  "black"))+
  labs(x= "pseudotime", col = "Sample size", y = "Simpson index (D)") +
  scale_color_viridis_d()  +
  geom_path(data=pop_diversity,aes(time,1-simpson),inherit.aes = F, lwd = 1, lty = 2)


legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + 
    guides(color = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom")
)


#x11()
cowplot::plot_grid(p0+ theme(legend.position="none"),
                   p1+ theme(legend.position="none"),
                   p7+ theme(legend.position="none"),
                   #p8+ theme(legend.position="none"),
                   p2+ theme(legend.position="none"),
                   p5+ theme(legend.position="none"),
                   p4+ theme(legend.position="none"),
                   p6+ theme(legend.position="none"),
                   p3+ theme(legend.position="none"),
                   p9+ theme(legend.position="none"),
                   NULL,
                   legend,
                   NULL,
                   ncol = 3, axis = "lr", align = "hv",labels = LETTERS[1:9])



ggsave(filename = "manuscript.V2/graphs/Figure_2.png", device = "png", width = 16, height = 16,units = "cm",dpi = 300)
ggsave(filename = "manuscript.V2/graphs/Figure_2.svg", device = "svg", width = 16, height = 16,units = "cm",dpi = 300)
ggsave(filename = "manuscript.V2/graphs/Figure_2.ps", device = "ps", width = 16, height = 16,units = "cm",dpi = 300)







## 3D visualization : 
mat = z %>% filter(variable =="gini") %>%
  pivot_wider(names_from = "pseudotime", values_from = "mean",-c(variable,sd)) %>%
  column_to_rownames("Cells")%>% as.matrix


zlim <- range(mat)
zlen <- zlim[2] - zlim[1] + 1
zlen=30
colorlut <- hcl.colors(zlen)

col <- colorlut[ scales::rescale(mat,to = c(1,30)) ]




rgl::open3d()
rgl::persp3d(x=(as.numeric(rownames(mat))),
             y=as.numeric(colnames(mat)),
             z=-mat, 
             color = col,xlab = "", ylab = " ", zlab = "", box= F, axes=T)
rgl::surface3d(x=as.numeric(rownames(mat)),
               y=as.numeric(colnames(mat)),
               z=-mat, back = "lines")




## coupon collector problem

# cc <- z %>% filter(variable =="clones", pseudotime %in% c(1,5,10,15,20))
# 
# ggplot(cc, aes(as.numeric(Cells),mean, col = pseudotime)) + geom_path() + geom_point() + labs(x = "Cells", y = "Clones") + geom_vline(xintercept = Nclones * log(Nclones) + 0.58*Nclones) + geom_errorbar(aes(ymin = mean-2*sd, ymax = mean +2*sd))
# 
# 

