## Corre Guillaume
## GENETHON
## Sept 2022

## Figure 3

# Load libraries --------------------------------------------------------------------



library(vegan); 
library(ineq);
library(rgl)





clones = seq(10,1000,20)
dominance = seq(0,1,.02)


res = list()

for(c in clones){
  for(d in dominance){
    
    if(d == 0){
      contrib = c(rep((1-d)/(c),times = c))
    } else {
      contrib = c(rep((1-d)/(c-1),times = c-1),d)
    }
    
    sh = vegan::diversity(contrib,index = "shannon",base = exp(1))
    pielou = sh / log(c)
    simpson = diversity(contrib, index = "simpson")
    gini = ineq(contrib,type = "Gini")
    
    res[[paste(c,d)]] <- data.frame(clones=c, dominance=d, Shannon = sh, Pielou = pielou, Simpson = simpson, Gini = gini)
  }
}
res = bind_rows(res)




## 3D visualization : 


mat = res %>% select(clones,dominance, Simpson) %>%
  pivot_wider(names_from = "dominance", values_from = "Simpson") %>%
  column_to_rownames("clones")%>% as.matrix


zlim <- range(mat)
zlen <- zlim[2] - zlim[1] + 1
zlen=30
colorlut <- hcl.colors(zlen)

col <- colorlut[ scales::rescale(mat,to = c(1,30)) ]




rgl::open3d()
rgl::clear3d()
rgl::persp3d(x=(as.numeric(rownames(mat))),
             y=as.numeric(colnames(mat)),
             z=(-mat), 
             #color = col,
             color = "grey",
             xlab = "Clones",
             lit = F,
             ylab = "Dominance", 
             zlab = "Simpson", xlim = c(0,1000),
             box= F, axes=T)
rgl::surface3d(x=as.numeric(rownames(mat)),
               y=as.numeric(colnames(mat)),
               z=(-mat), back = "lines")




rgl::rgl.snapshot(filename = "manuscript.V2/graphs/Figure_3C.png")
