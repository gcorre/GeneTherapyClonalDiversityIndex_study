## Corre Guillaume
## GENETHON
## Sept 2022

## Figure 1

# Load libraries --------------------------------------------------------------------



library(tidyverse); 
library(data.table)
library(ggpubr);



theme_set(theme_bw() + theme(text = element_text(color = "black",size = 11), 
                             axis.text = element_text(color = "black",size=11),
                             axis.title = element_text(color = "black",size=11)))



# Load dataset of published papers --------------------------------------------------

files = list.files(path = "manuscript.V2/data/figure5_VISA_bushman/",pattern = "RData_summary", full.names = T)

sum = lapply(files, fread)

sum = bind_rows(sum)



fig5_data = sum  %>%
  mutate(patient = recode(patient, 
                          pFR01 = 2, 
                          pFR03 = 4, 
                          pFR04 = 5, 
                          pFR05 = 7,
                          pWAS_UK02 = 1,
                          pWAS_UK04 = 6,
                          pWAS_UK05 = 9, 
                          pWAS_UK07 = 11,
                          pWAS_UK11 = 8)) %>% 
  filter(cellType=="PBMC") %>% 
  mutate(time = str_extract(timePoint,pattern = "[0-9]+"), 
         time_num = ifelse(str_starts(timePoint,pattern = "Y"),yes = as.numeric(time)*12,no = as.numeric(time))) %>% 
  select(time_num, patient,Pielou=pielou2, Simpson=simpson) %>%
  pivot_longer(cols = Pielou:Simpson)





indices_plot <- ggplot(fig5_data,
                       aes(group = patient,time_num, value, col = factor(patient))) +
  geom_point()  +
  geom_line() +
  scale_y_continuous(limits = c(0,1))+
  facet_grid(~name) +
  #theme_bw() + 
  scale_color_viridis_d()+
  labs(x = "Time post gene therapy (months)", y = "Index value", col = "patient") +
  theme(legend.position = "none")


corr_plot <- ggplot(fig5_data %>% pivot_wider(id_cols = c(time_num, patient), names_from = "name", values_from = "value", values_fn = mean),
                    aes(Pielou, Simpson, col = factor(patient))) +
  geom_point() +
  scale_x_continuous(limits = c(0.5,1))+
  scale_y_continuous(limits = c(0.5,1))+
  scale_color_viridis_d() + 
  theme_bw() +
  labs(col = "patient")+
  theme(legend.text = element_text(size=10),legend.key.width= unit(5,units="points")) +
  guides(color = guide_legend(ncol = 2)) + coord_fixed()

indices_plot

ggsave(filename = "manuscript.V2/graphs/Figure_5.png", device = "png", width = 3.4, height = 2,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_5.svg", device = "svg", width = 3.4, height = 2,units = "in",dpi = 72)
ggsave(filename = "manuscript.V2/graphs/Figure_5.ps", device = "ps",  width = 3.4, height = 2,units = "in",dpi = 72)


