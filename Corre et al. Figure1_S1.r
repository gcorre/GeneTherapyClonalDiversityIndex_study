## Corre Guillaume
## GENETHON
## Feb 2023

## Figure 1

# Load libraries --------------------------------------------------------------------



library(tidyverse); 
library(ggpubr);




# Load dataset of published papers --------------------------------------------------


dt <- read.delim("manuscript.V2/data/corre et al._publishedData_figure1.csv", sep=";")


# convert time to numeric, change log base if necessary
dt <- dt %>%
  mutate(TIME_PAPER = as.numeric(str_remove_all(TIME_PAPER,pattern = "[A-Za-z]+")),
         PATIENT = as.numeric(str_extract(PATIENT,pattern = "[0-9]+")),
         TIME = case_when(UNIT == "days" ~ TIME_PAPER / 30.5,
                          UNIT == "year" ~ TIME_PAPER * 12,
                          TRUE ~ TIME_PAPER),
         SHANNON = ifelse(BASE_LOG == "exp(1)",yes =  SHANNON_DIGITALIZED, no = SHANNON_DIGITALIZED / log2(exp(1))))



# select PBMC and myeloid cell types | filter trials of interest

dt_sub = dt %>% 
  filter(CELL_POP%in% c("PBMC","","MIXED","PB MYELOID","MYELOID"), 
         TRIAL != "MORRIS et AL. 2017", !str_detect(TRIAL, "ACEIN-BEY"))


# Make plots  plot -----------------------------------------------------------
theme_set(theme_bw() + theme(text = element_text(color = "black",size = 11), 
                             axis.text = element_text(color = "black",size=11),
                             axis.title = element_text(color = "black",size=11)))

# plot of pielou index 

a <- ggplot(dt_sub,
            aes(TIME, SHANNON / log(TOTAL_IS), 
                col = factor(PATIENT), group = paste(PATIENT, CELL_POP))) +
  facet_wrap(~TRIAL, scales = "free_x",ncol = 1) + 
  geom_point(size=1) +
  geom_line(lwd = 0.6) + 
  labs(y = "Pielou index",
       x = NULL,
       col = "Patient") +
  lims(y = c(0,1)) +
  #geom_text_repel(aes(label = CLONE),show.legend = F,segment.colour = "black",force = 100) +
  geom_point(inherit.aes = F,data = dt %>% filter(LEUKEMIA!=""), aes(TIME, SHANNON/log(TOTAL_IS)),  cex = 5, pch = 1, col = "red")+
  scale_color_viridis_d()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0.5, lty=2)




# plot of Shannon index 

b <- ggplot(dt_sub ,
            aes(TIME, SHANNON , 
                col = factor(PATIENT), group = paste(PATIENT, CELL_POP))) +
  facet_wrap(~TRIAL,ncol = 1, scale = "free_x") + 
  geom_point(size=1) +
  geom_line(lwd = 0.6) + 
  labs(y = "Shannon index",
       x = "Time (Months)",
       col = "Patient") +
  #geom_text_repel(aes(label = CLONE),show.legend = F,segment.colour = "black",force = 100) +
  geom_point(inherit.aes = F,data = dt %>% filter(LEUKEMIA!=""), aes(TIME, SHANNON), cex = 5, pch = 1, col = "red")+
  scale_colour_viridis_d()+
  theme(legend.position = "none")



# plot of IS number 

c <- ggplot(dt_sub,
            aes(TIME, (TOTAL_IS+1) , 
                col = factor(PATIENT), group = paste(PATIENT, CELL_POP))) +
  facet_wrap(~TRIAL,ncol = 1, scale = "free_x") + 
  geom_point(size = 1) +
  geom_line(lwd=0.6) + 
  labs(y = "Unique IS\n",
       x = NULL,
       col = "Patient") +
  #geom_text_repel(aes(label = CLONE),show.legend = F,segment.colour = "black",force = 100) +
  geom_point(inherit.aes = F,data = dt %>% filter(LEUKEMIA!=""), aes(TIME, (TOTAL_IS)), cex = 5, pch = 1, col="red") +
  scale_colour_viridis_d() +
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "none")


x11(width = 7, height = 7);
ggarrange(c, b, a, ncol = 3, common.legend = T, labels = "AUTO",legend = "bottom",align = 'hv' )


ggsave(filename = "manuscript.V2/graphs/Figure_1.png", device = "png", width = 7, height = 7,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_1.svg", device = "svg",  width = 7, height = 7,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_1.ps", device = "ps",  width = 7, height = 7,units = "in",dpi = 72)





##### END of script #####


# add logistical regression

dt_sub_logistic_keep <- dt_sub %>% filter(LEUKEMIA != "") %>% distinct(TIME,TRIAL) %>% mutate(TIME_LL = TIME-2,
                                                                                              TIME_UL = TIME+2) %>% select(-TIME) %>% distinct()


dt_sub_logistic <- dt_sub %>% 
  #left_join(dt_sub_logistic_keep) %>% 
  filter(TIME>1) %>%
  mutate(pielou = SHANNON / log(TOTAL_IS),
         clonal_dominance = case_when(LEUKEMIA !="" ~ 1,
                                      (TRIAL == "BRAUN ET AL.2014" & PATIENT==9 & between(TIME,30,40)) ~ 1,
                                        TRUE ~ 0)) %>%
  distinct(TRIAL, PATIENT, TIME,pielou, clonal_dominance)



ggplot(dt_sub_logistic, aes(factor(clonal_dominance), pielou)) + geom_point(position = position_jitter(width = .1))


model <- glm( clonal_dominance ~ pielou, data = dt_sub_logistic, family = binomial())
summary(model)

pred = data.frame(x=seq(0,1,0.01),
                  pred=predict.glm(model,newdata = data.frame(pielou=seq(0,1,.01)),type = "response",se.fit = T))


UC50 = -coefficients(model)[1]/coefficients(model)[2]


ggplot(dt_sub_logistic, aes(pielou,clonal_dominance)) +
  geom_point()+
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T) +
  geom_vline(xintercept = UC50, lty = 2)

