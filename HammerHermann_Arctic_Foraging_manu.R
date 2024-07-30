





# Top ---------------------------------------------------
set.seed(42)

setwd("") #identify location for data files

library(tidyverse)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(grid)
library(RColorBrewer)
library(readr)
library(gridExtra)
library(lubridate)
library(reshape2)
library(MuMIn)
library(vegan)
library(BiodiversityR)
library(rphylopic)
library(png)
library(RVAideMemoire)
library(usedist)
library(indicspecies)
library(ggpubr)
library(rgl)
library(raster)


theme_set(theme_bw(base_size=45))

'%notin%'<-Negate('%in%')

details<-read_csv("all_dietcontents_final_v2.csv")
summaries<-read_csv("all_dietsummaries_final.csv")
dietMat<-read_csv("all_diet_PAdietMat_final_v2.csv")




# Data Summaries ----------------------------------------------------------

dataSummaries<-summaries%>%
  mutate(family=ifelse(Species=="Char","Arctic char","Sculpin"))%>%
  group_by(family)%>%
  summarise(N=n(),
            Mass=mean(Mass_g,na.rm=T),
            Mass_SD=sd(Mass_g,na.rm=T),
            TL=mean(TL_mm,na.rm=T),
            TL_SD=sd(TL_mm,na.rm=T),
            HSI=mean(100*(Liver_Mass_g/(Mass_g-Liver_Mass_g)),na.rm=T),
            HSI_SD=sd(100*(Liver_Mass_g/(Mass_g-Liver_Mass_g)),na.rm=T),
            GSI=mean(100*(Gonad_Mass_g/(Mass_g-Gonad_Mass_g)),na.rm=T),
            GSI_SD=sd(100*(Gonad_Mass_g/(Mass_g-Gonad_Mass_g)),na.rm=T),
            fultonK=mean(10^5*Mass_g/TL_mm^3,na.rm=T),
            fultonK_SD=sd(10^5*Mass_g/TL_mm^3,na.rm=T))

yearSummaries<-summaries%>%
  mutate(family=ifelse(Species=="Char","Arctic char","Sculpin"))%>%
  group_by(family,Year)%>%
  summarise(N=n(),
            Mass=mean(Mass_g,na.rm=T),
            Mass_SD=sd(Mass_g,na.rm=T),
            TL=mean(TL_mm,na.rm=T),
            TL_SD=sd(TL_mm,na.rm=T),
            HSI=mean(100*(Liver_Mass_g/(Mass_g-Liver_Mass_g)),na.rm=T),
            HSI_SD=sd(100*(Liver_Mass_g/(Mass_g-Liver_Mass_g)),na.rm=T),
            GSI=mean(100*(Gonad_Mass_g/(Mass_g-Gonad_Mass_g)),na.rm=T),
            GSI_SD=sd(100*(Gonad_Mass_g/(Mass_g-Gonad_Mass_g)),na.rm=T),
            fultonK=mean(10^5*Mass_g/TL_mm^3,na.rm=T),
            fultonK_SD=sd(10^5*Mass_g/TL_mm^3,na.rm=T))

empty<-left_join(details,summaries)%>%
  mutate(species2=ifelse(Species=="Arctic char","Arctic char","Sculpin"))%>%
  group_by(species2)%>%
  mutate(N=n_distinct(Fish_ID))%>%
  group_by(species2,Year)%>%
  mutate(N_year=n_distinct(Fish_ID))%>%
  filter(Category=="Empty" | Diet_Mass_g ==0)%>%
  group_by(species2)%>%
  mutate(p=n_distinct(Fish_ID)/N*100)%>%
  group_by(species2,Year)%>%
  mutate(p_year=n_distinct(Fish_ID)/N_year*100)%>%
  dplyr::select(species2,Year,N,N_year,p,p_year)%>%distinct()

# Analyses--frequency ----------------------------------------------------------

#### Frequency of Occurrence

sumFOO<-filter(dietMat)%>%
  dplyr::select(matches("^prey"))%>%
  mutate(prey_totalFish=ifelse(prey_Fish==1 | prey_Sand.lance==1 | prey_Arctic.cod==1,1,0),
         prey_totalAmph=ifelse(prey_Amphipod==1 | prey_Gammaracanthus.sp.==1 | prey_Gammarus.sp.==1 | prey_Onisimus.sp.==1 | prey_Themisto.sp.==1,1,0))%>%
  t()%>%as.data.frame()%>%
  reframe(foo=rowSums(.)/nrow(dietMat),
            N=rowSums(.),
            prey=substr(row.names(.),6,100))


#Cumulative across years
SculpinFOO<-filter(dietMat,Species %in% c("Slimy","Fourhorn"))%>%
  dplyr::select(matches("^prey"))%>%
  mutate(prey_totalFish=ifelse(prey_Fish==1 | prey_Sand.lance==1 | prey_Arctic.cod==1,1,0),
         prey_totalAmph=ifelse(prey_Amphipod==1 | prey_Gammaracanthus.sp.==1 | prey_Gammarus.sp.==1 | prey_Onisimus.sp.==1 | prey_Themisto.sp.==1,1,0))%>%
  colSums()/nrow(dietMat%>%filter(Species %in% c("Slimy","Fourhorn")))

CharFOO<-filter(dietMat,Species=="Arctic char")%>%
  dplyr::select(matches("^prey"))%>%
  mutate(prey_totalFish=ifelse(prey_Fish==1 | prey_Sand.lance==1 | prey_Arctic.cod==1,1,0),
         prey_totalAmph=ifelse(prey_Amphipod==1 | prey_Gammaracanthus.sp.==1 | prey_Gammarus.sp.==1 | prey_Onisimus.sp.==1 | prey_Themisto.sp.==1,1,0))%>%
  colSums()/nrow(dietMat%>%filter(Species=="Arctic char"))

foo<-bind_rows(data.frame(foo=CharFOO,species="Arctic char",
                          prey=c(gsub("prey_","",colnames(dplyr::select(dietMat,matches("^prey")))),"total.Fish","total.Amphipod")),
               data.frame(foo=SculpinFOO,species="Sculpin",
                          prey=c(gsub("prey_","",colnames(dplyr::select(dietMat,matches("^prey")))),"total.Fish","total.Amphipod")))%>%
  mutate(prey=gsub("\\.","\n",prey))%>%
  group_by(prey)%>%
  mutate(diff=max(foo)-min(foo))

ggplot(foo,aes(prey,foo,fill=species))+
  geom_col(position="dodge")

#Separated by years
presenceYears<-dietMat%>%
  mutate(species=ifelse(Species=="Arctic char","Arctic char","Sculpin"))%>%
  dplyr::select(species,Year,matches("^prey"))%>%
  mutate(prey_totalFish=ifelse(prey_Fish==1 | prey_Sand.lance==1 | prey_Arctic.cod==1,1,0),
         prey_totalAmphipod=ifelse(prey_Amphipod==1 | prey_Gammaracanthus.sp.==1 | prey_Gammarus.sp.==1 | prey_Onisimus.sp.==1 | prey_Themisto.sp.==1,1,0))%>%
  group_by(species,Year)%>%
  summarise_each(funs=sum)
totalYears<-dietMat%>%
  mutate(species=ifelse(Species=="Arctic char","Arctic char","Sculpin"))%>%
  group_by(species,Year)%>%
  summarise(N=n())
fooYears<-left_join(presenceYears,totalYears)%>%
  pivot_longer(cols=matches("^prey"),names_to="prey",values_to="noo")%>%
  ungroup()%>%
  mutate(prey=gsub("^prey_","",prey),
         prey=gsub("\\.([A-z])"," \\1",prey),
         prey=case_when(prey=="Fish"~"Unid. Fish",
                        prey=="Amphipod"~"Unid. Amphipod",
                        grepl("total",prey)~gsub("total","Total ",prey),
                        TRUE~prey),
         foo=noo/N,
         Year=as.character(Year))
#Ordering the prey items so like are by like
fooYears$preyF<-factor(gsub("sp\\.","spp\\.",fooYears$prey),
                       levels=c("Total Fish","Arctic cod","Sand lance","Unid. Fish",
                                              "Total Amphipod","Gammarus spp.","Gammaracanthus spp.","Onisimus spp.","Themisto spp.","Unid. Amphipod",
                                              "Krill","Mysid","Copepod","Chironomid","Jellyfish","Sea Angel","Miscellaneous Invert",
                                              "Miscellaneous","Algae","Undigestible","Digested"))


#Emphasizing species differences
ggplot(fooYears,aes(preyF,foo,fill=Year))+
  geom_col(position="dodge",color="black")+
  scale_fill_viridis_d(name="Species")+
  scale_y_continuous(name="Frequency of Occurrence",limits=c(0,1),expand=expansion(add=0))+
  scale_x_discrete(name="Diet Item")+
  theme(axis.text.x=element_text(angle=20,hjust=1,vjust=1.12),
        axis.title.x=element_text(vjust=5))+
  facet_wrap(~species,nrow=2)


#Pulling from RPhylopic but they have been lost, need to update
amphipod<-rasterGrob(readPNG("silhouettes/amphipod.png",native=T))
chironomid<-rasterGrob(readPNG("silhouettes/chironomid.png",native=T))
cod<-rasterGrob(readPNG("silhouettes/cod.png",native=T))
copepod<-rasterGrob(readPNG("silhouettes/copepod.png",native=T))
jellyfish<-rasterGrob(readPNG("silhouettes/jellyfish.png",native=T))
krill<-rasterGrob(readPNG("silhouettes/krill.png",native=T))
mysid<-rasterGrob(readPNG("silhouettes/mysid.png",native=T))
sandlance<-rasterGrob(readPNG("silhouettes/sand_lance.png",native=T))
seaAngel<-rasterGrob(readPNG("silhouettes/sea_angel.png",native=T))

#Emphasizing inter-annual variation
ggplot(data=fooYears)+
  geom_col(aes(preyF,foo,fill=species),
           position="dodge",color="black",width=0.5)+
  annotation_custom(cod,xmin=1.5,xmax=2.5,ymin=0.8,ymax=1)+
  annotation_custom(sandlance,xmin=2.5,xmax=3.5,ymin=0.8,ymax=1)+
  annotation_custom(amphipod,xmin=7,xmax=8,ymin=0.8,ymax=1)+
  annotation_custom(krill,xmin=10.5,xmax=11.5,ymin=0.8,ymax=1)+
  annotation_custom(mysid,xmin=11.5,xmax=12.5,ymin=0.8,ymax=1)+
  annotation_custom(copepod,xmin=12.5,xmax=13.5,ymin=0.8,ymax=1)+
  annotation_custom(chironomid,xmin=13.5,xmax=14.5,ymin=0.8,ymax=1)+
  annotation_custom(jellyfish,xmin=14.6,xmax=15.6,ymin=0.8,ymax=1)+
  annotation_custom(seaAngel,xmin=15.5,xmax=16.5,ymin=0.8,ymax=1)+
  geom_vline(aes(xintercept=4.5),lty=3,color="black",size=1.5)+
  geom_vline(aes(xintercept=10.5),lty=3,color="black",size=1.5)+
  geom_vline(aes(xintercept=17.5),lty=3,color="black",size=1.5)+
  geom_segment(aes(x=4.62,xend=4.62,y=0.92,yend=0.98),size=1.5)+
  geom_segment(aes(x=10.38,xend=10.38,y=0.92,yend=0.98),size=1.5)+
  geom_segment(aes(x=4.6,xend=6.9,y=0.98,yend=0.98),size=1.5)+
  geom_segment(aes(x=8.2,xend=10.4,y=0.98,yend=0.98),size=1.5)+
  scale_fill_viridis_d(name="Consumer")+
  scale_y_continuous(name="Frequency of Occurrence",
                     limits=c(0,1),expand=expansion(add=0),
                     labels=scales::percent_format())+
  scale_x_discrete(name="Diet Item")+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.12),axis.text.y=element_text(size=30),
        axis.title.x=element_text(vjust=5),strip.text=element_text(size=30),
        panel.grid.major.x=element_line(size=20),legend.margin=margin(0,0,0,-33))+
  facet_grid(Year~.)



# NMDS -------------------------------------------

#Duplicating with ALL the counts
allCount.Mat<-as.matrix(dplyr::select(dietMat,matches("^prey")))
allEnv.Mat<-(dplyr::select(dietMat,-matches("^prey")))

#Percent zeroes--number of zeroes times dim of matrix
(sum(allCount.Mat==0)/(nrow(allCount.Mat)*ncol(allCount.Mat)))*100
#Definitely a sparse dataset then, typical of community data, need to try with an NMDS and not a PCA


#Using the whole dataset--but can use the one without minor species (Jellyfish, Sea Angel)
set.seed(42)

#For the all diets NMDS
original.dist<-vegdist(allCount.Mat)
stress_values<-numeric(6)
r2<-numeric(6)

for (n in 1:6) {
  nmds.resu <- metaMDS(allCount.Mat, k=n, distance = "bray", try=250, autotransform=F)
  stress_values[n]<-nmds.resu$stress
  nmds.scores<-vegan::scores(nmds.resu)$sites
  nmds.dist<-dist(nmds.scores)
  r2[n]<-summary(lm(original.dist~nmds.dist))[[8]]
}
plot(stress_values, xlab="Number of axes", ylab="Stress",type="b")
abline(h=0.2,col="red")

View(stress_values) 

#Go back and create the output for the 3 dimensions NMDS
count_NMDS<-metaMDS(allCount.Mat, distance = "bray", k = 3, try=250, autotransform=F)
r2<-summary(lm(original.dist~dist(vegan::scores(count_NMDS)$sites)))[[8]]
actualStress<-count_NMDS$stress
stressplot(count_NMDS) #This is the visual of stress, the divergence of observed and ordinated distance. It's random, that's good

#Writing it out to use in PC-ORD for comparison and for the axes R2 values
#temp<-cbind(allEnv.Mat[,1],allCount.Mat)
#write.csv(temp,"/Users/nh1087/Documents/NMDS_Matrix.csv",row.names = F)




# Printing NMDS -----------------------------------------------------------

#Print the species scores and sample scores
NMDS_species<-as.data.frame(count_NMDS$species)
NMDS_scores<-as.data.frame(count_NMDS$points)

speciesPal<-scales::viridis_pal(option="G",begin=0.05,end=0.95)(2)
yearPal<-scales::viridis_pal(option="H",begin=0.05,end=0.95)(3)

allEnv.Mat<-allEnv.Mat%>%
  mutate(family=ifelse(Species=="Arctic char","Arctic char","Sculpin"),
         famCol=ifelse(family=="Arctic char",speciesPal[1],speciesPal[2]),
         yearCol=case_when(Year==2017~yearPal[1],
                           Year==2018~yearPal[2],
                           Year==2019~yearPal[3]))


dgg_ordiplot <- function(ord, groups, scaling = 1, choices = c(1,2), kind = c("sd", "se", "ehull"), conf=NULL, show.groups="all", ellipse = TRUE, label = FALSE, hull = FALSE, spiders = FALSE, pt.size = 3, plot=TRUE) {
  groups <- as.factor(groups)
  if (show.groups[1]=="all") {
    show.groups <- as.vector(levels(groups))
  }
  
  # Get site coordinates to plot.
  df_ord <- vegan::scores(ord, display = "sites", scaling=scaling, choices=choices)
  axis.labels <- ord_labels(ord)[choices]
  df_ord <- data.frame(x=df_ord[ , 1], y=df_ord[ , 2], Group=groups)
  
  # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
  
  # Get parameters from the ordiellipse function.
  if (is.null(conf)) {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", label = label)
  } else {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", conf = conf, label = label)
  }
  
  # Get points to plot for the ellipses.
  df_ellipse <- data.frame()
  for(g in show.groups) {
    df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                             vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
  }
  colnames(df_ellipse) <- c("x", "y", "Group")
  df_ellipse <- df_ellipse[ , c(3,1,2)]
  
  # Make data frame for hulls.
  rslt.hull <- vegan::ordihull(ord, groups = groups, scaling = scaling, choices = choices, show.groups = show.groups, draw = "none")
  df_hull <- data.frame()
  df_temp <- data.frame()
  for (g in show.groups) {
    x <- rslt.hull[[g]][ , 1]
    y <- rslt.hull[[g]][ , 2]
    Group <- rep(g, length(x))
    df_temp <- data.frame(Group = Group, x=x, y=y)
    df_hull <- rbind(df_hull, df_temp)
  }
  
  # Make a data frame for the spiders.
  df_spiders <- df_ord
  df_spiders$cntr.x <- NA
  df_spiders$cntr.y <- NA
  for (g in show.groups) {
    df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
  }
  df_spiders <- df_spiders[ , c(3,4,5,1,2)]
  df_spiders <- df_spiders[order(df_spiders$Group), ]
  df_spiders <- df_spiders[df_spiders$Group %in% show.groups, ]
  
  # Make basic ggplot with ellipses.
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  plt <- ggplot2::ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = pt.size) +
    xlab(xlab) + ylab(ylab)
  
  # Add ellipses.
  if (ellipse == TRUE) {
    plt <- plt + geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add labels.
  if (label == TRUE) {
    plt <- plt + geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE)
  }
  
  # Add hulls.
  if (hull == TRUE) {
    plt <- plt + geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add spiders.
  if (spiders == TRUE) {
    plt <- plt + geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE)
  }
  
  plt <- plt + coord_fixed(ratio=1)
  
  # Plot?
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, df_mean.ord=df_mean.ord, df_ellipse=df_ellipse, df_hull=df_hull, df_spiders=df_spiders, plot=plt))
}


speciesCentroids <- aggregate(NMDS_scores, by=list(allEnv.Mat$family),mean)
yearCentoids <- aggregate(NMDS_scores, by=list(allEnv.Mat$Year),mean)
spyrCentoids <- aggregate(NMDS_scores, by=list(allEnv.Mat$family,allEnv.Mat$Year),mean)%>%
  mutate(lab=paste0(ifelse(Group.1=="Sculpin","S","AC"),substr(Group.2,3,4)))


axes12<-ggplot()+
  geom_polygon(data=NMDS_scores%>%group_by(family=allEnv.Mat$family,Year=as.character(allEnv.Mat$Year))%>%slice(chull(MDS1,MDS2)),
               aes(x=MDS1,y=MDS2,group=paste(family,Year),color=family,lty=Year),alpha=0.1,lwd=1.5)+
  geom_point(data=NMDS_scores,aes(MDS1,MDS2,fill=allEnv.Mat$family,shape=as.character(allEnv.Mat$Year)),color="black",size=6)+
  geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=gsub(" sp"," spp",gsub("\\.([A-z])"," \\1",gsub("prey_","",row.names(NMDS_species))))),size=9)+
  scale_fill_viridis_d(option="E",name="Consumer")+
  scale_color_viridis_d(option="E",name="Consumer",guide="none")+
  scale_shape_manual(values=c(21,22,23),name="Year")+
  scale_linetype_manual(values=c(1,5,3))+
  guides(fill=guide_legend(override.aes = list(shape=21,size=20)),
         linetype=guide_legend(override.aes=list(color="black",size=20,fill="transparent")))+
  theme(legend.spacing = unit(0,"cm"))+
  geom_point(data=spyrCentoids,aes(MDS1,MDS2,fill=Group.1,shape=as.character(Group.2)),
             color="red",size=16,stroke=2,alpha=0.7,show.legend=F)+
  geom_label(data=spyrCentoids,aes(MDS1,MDS2,label=lab),
             color="red",size=8,alpha=0.7,show.legend=F)+
  xlab("Axis 1")+ylab("Axis 2")


axes13<-ggplot()+
  geom_polygon(data=NMDS_scores%>%group_by(family=allEnv.Mat$family,Year=as.character(allEnv.Mat$Year))%>%slice(chull(MDS1,MDS3)),
               aes(x=MDS1,y=MDS3,group=paste(family,Year),color=family,lty=Year),alpha=0.1,lwd=1.5)+
  geom_point(data=NMDS_scores,aes(MDS1,MDS3,fill=allEnv.Mat$family,shape=as.character(allEnv.Mat$Year)),color="black",size=6)+
  geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS3))+
  geom_text(data=NMDS_species,aes(x=MDS1,y=MDS3,label=gsub(" sp"," spp",gsub("\\.([A-z])"," \\1",gsub("prey_","",row.names(NMDS_species))))),size=9)+
  scale_fill_viridis_d(option="E",name="Consumer")+
  scale_color_viridis_d(option="E",name="Consumer",guide="none")+
  scale_shape_manual(values=c(21,22,23),name="Year")+
  scale_linetype_manual(values=c(1,5,3))+
  guides(fill=guide_legend(override.aes = list(shape=21,size=20)),
         linetype=guide_legend(override.aes=list(color="black",size=20,fill="transparent")))+
  #geom_text(aes(x=-1.6,y=-1.6,label=paste0("Stress: ",round(actualStress*100,2))),size=6,hjust=0)+
  geom_point(data=spyrCentoids,aes(MDS1,MDS3,fill=Group.1,shape=as.character(Group.2)),
             color="red",size=16,stroke=2,alpha=0.7,show.legend=F)+
  geom_label(data=spyrCentoids,aes(MDS1,MDS3,label=lab),
             color="red",size=8,alpha=0.7,show.legend=F)+
  xlab("Axis 1")+ylab("Axis 3")

axes23<-ggplot()+
  geom_polygon(data=NMDS_scores%>%group_by(family=allEnv.Mat$family,Year=as.character(allEnv.Mat$Year))%>%slice(chull(MDS2,MDS3)),
               aes(x=MDS2,y=MDS3,group=paste(family,Year),color=family,lty=Year),alpha=0.1,lwd=1.5)+
  geom_point(data=NMDS_scores,aes(MDS2,MDS3,fill=allEnv.Mat$family,shape=as.character(allEnv.Mat$Year)),color="black",size=6)+
  geom_segment(data=NMDS_species,aes(x=0,xend=MDS2,y=0,yend=MDS3))+
  geom_text(data=NMDS_species,aes(x=MDS2,y=MDS3,label=gsub(" sp"," spp",gsub("\\.([A-z])"," \\1",gsub("prey_","",row.names(NMDS_species))))),size=9)+
  scale_fill_viridis_d(option="E",name="Consumer")+
  scale_color_viridis_d(option="E",name="Consumer",guide="none")+
  scale_shape_manual(values=c(21,22,23),name="Year")+
  scale_linetype_manual(values=c(1,5,3))+
  scale_y_continuous(breaks=c(-1,0,1))+
  guides(fill=guide_legend(override.aes = list(shape=21,size=20)),
         linetype=guide_legend(override.aes=list(color="black",size=20,fill="transparent")))+
  geom_text(aes(x=-2.55,y=-1.2,label=paste0("Stress: ",round(actualStress*100,2))),size=15,hjust=0)+
  geom_point(data=spyrCentoids,aes(MDS2,MDS3,fill=Group.1,shape=as.character(Group.2)),
             color="red",size=16,stroke=2,alpha=0.7,show.legend=F)+
  geom_label(data=spyrCentoids,aes(MDS2,MDS3,label=lab),
             color="red",size=8,alpha=0.7,show.legend=F)+
  xlab("Axis 2")+ylab("Axis 3")

ggarrange(axes12,axes13,axes23,common.legend = T,legend="top",ncol=1)


interleave <- function(v1, v2) as.vector(rbind(v1,v2))

plot3d(NMDS_scores$MDS1,NMDS_scores$MDS2,NMDS_scores$MDS3,
       size=8,col=allEnv.Mat$famCol)
segments3d(interleave(0, NMDS_species$MDS1),
           interleave(0, NMDS_species$MDS2),
           interleave(0, NMDS_species$MDS3),
           alpha = 0.4, col = "blue")
text3d(NMDS_species$MDS1,
       NMDS_species$MDS2,
       NMDS_species$MDS3,
       texts=gsub("prey_","",row.names(NMDS_species)),add=T)
play3d(spin3d())





# Indicator Species Analysis ----------------------------------------------




#ISA to see what diet items might be associated with the different years that make them different from each other
ISA_year = multipatt(as.data.frame(allCount.Mat), allEnv.Mat$Year,
                     func = "IndVal.g", duleg=TRUE, control = how(nperm=4999))

#What species?
summary(ISA_year, indvalcomp = T) 
ISA_year$str

#Extract them, with the year they're significant for
ISA_year_df<-dplyr::select(subset(ISA_year$sign, p.value<0.05),index,stat,p.value)
ISA_year_df$Species<-rownames(ISA_year_df)
ISA_year_df$group<-colnames(ISA_year$B)[ISA_year_df$index]
ISA_year_df

fullISA_year<-as.data.frame(ISA_year$str)%>%
  mutate(Species=rownames(ISA_year$str))%>%
  pivot_longer(cols=colnames(ISA_year$str),names_to="group",values_to="stat")


#ISA to see what diet items might be associated with the different time periods
ISA_species = multipatt(as.data.frame(allCount.Mat), allEnv.Mat$family,
                        func = "IndVal.g", duleg=TRUE, control = how(nperm=4999))

#What species?
summary(ISA_species,indvalcomp = T) 
ISA_species$str

#Extract them, with the year they're significant for
ISA_species_df<-dplyr::select(subset(ISA_species$sign, p.value<0.05),index,stat,p.value)
ISA_species_df$Species<-rownames(ISA_species_df)
ISA_species_df$group<-colnames(ISA_species$B)[ISA_species_df$index]
ISA_species_df

fullISA_species<-as.data.frame(ISA_species$str)%>%
  mutate(Species=rownames(ISA_species$str))%>%
  pivot_longer(cols=colnames(ISA_species$str),names_to="group",values_to="stat")

## Skipping ##
#ISA to see what diet items might be associated with the different time periods
ISA_spyr = multipatt(as.data.frame(allCount.Mat), paste(allEnv.Mat$family,allEnv.Mat$Year),
                     func = "IndVal.g", duleg=TRUE, control = how(nperm=4999))

#What species?
summary(ISA_spyr,indvalcomp = T) 
ISA_spyr$str

#Extract them, with the year they're significant for
ISA_spyr_df<-dplyr::select(subset(ISA_spyr$sign, p.value<0.05),index,stat,p.value)
ISA_spyr_df$Species<-rownames(ISA_spyr_df)
ISA_spyr_df$group<-colnames(ISA_spyr$B)[ISA_spyr_df$index]
ISA_spyr_df

fullISA_spyr<-as.data.frame(ISA_spyr$str)%>%
  mutate(Species=rownames(ISA_spyr$str))%>%
  pivot_longer(cols=colnames(ISA_spyr$str),names_to="group",values_to="stat")


fullISA<-bind_rows(fullISA_species,fullISA_year)%>% #ADD IN spyr IF NEEDED
  mutate(Species=gsub("prey_","",Species),
         Species=gsub("\\.([A-z])"," \\1",Species),
         Species=factor(Species,levels=c("Arctic cod","Sand lance","Fish",
                                         "Gammarus sp.","Gammaracanthus sp.","Onisimus sp.","Themisto sp.","Amphipod",
                                         "Krill","Mysid","Copepod","Chironomid","Jellyfish","Sea Angel","Miscellaneous Invert",
                                         "Miscellaneous","Algae","Undigestible","Digested")),
         group=factor(group,levels=c("Arctic char","Sculpin","2017","2018","2019")))

significantISA<-bind_rows(ISA_species_df,ISA_year_df)%>%
  mutate(Species=gsub("prey_","",Species),
         Species=gsub("\\.([A-z])"," \\1",Species),
         Species=factor(Species,levels=c("Arctic cod","Sand lance","Fish",
                                         "Gammarus sp.","Gammaracanthus sp.","Onisimus sp.","Themisto sp.","Amphipod",
                                         "Krill","Mysid","Copepod","Chironomid","Jellyfish","Sea Angel","Miscellaneous Invert",
                                         "Miscellaneous","Algae","Undigestible","Digested")),
         group=factor(group,levels=c("Arctic char","Sculpin","2017","2018","2019")))
fullISA<-left_join(fullISA,significantISA)%>%
  mutate(p.value.sig=case_when(p.value<0.001~"***",
                               p.value<0.01 & p.value>=0.001~"**",
                               p.value<=0.05 & p.value>=0.01~"*"))
fullISA%>%
  group_by(group)%>%
  summarise(meanISA=mean(stat),
            medianISA=median(stat),
            seISA=sd(stat)/n_distinct(Species),
            minISA=min(stat),max=max(stat))

legendMap<-data.frame(group=unique(fullISA$group),
                      xpos=c(1,2,0,0,0,1,1,1,2,2,2),
                      ypos=c(0,0,1,2,3,1,2,3,1,2,3))
legendGrob<-ggplotGrob(ggplot(legendMap)+
  geom_point(aes(xpos,ypos,fill=group),shape=22,size=15,show.legend=F)+
  geom_text(aes(0.5,-1.6,label="Arctic char",angle=45),size=7)+
  geom_text(aes(1.75,-1.25,label="Sculpin",angle=45),size=7)+
  geom_text(aes(-1,1.1,label="2017"),size=7)+
  geom_text(aes(-1,2.1,label="2018"),size=7)+
  geom_text(aes(-1,3.1,label="2019"),size=7)+
  scale_fill_manual(values=c(hcl(h=0,c=200,l=33),hcl(h=0,c=100,l=10),hcl(h=0,c=100,l=50),hcl(h=0,c=100,l=100),
                             hcl(h=256,c=200,l=33),hcl(h=256,c=100,l=10),hcl(h=256,c=100,l=50),hcl(h=256,c=100,l=100),
                             hcl(h=0,c=0,l=10),hcl(h=0,c=0,l=50),hcl(h=0,c=0,l=100)))+
  ylim(-5,5)+xlim(-5,5)+
  theme_void())
  

ggplot(fullISA)+
  geom_col(aes(x=Species,y=stat,fill=group),size=0.05,color="black",position=position_dodge(width=1),show.legend=F)+
  geom_text(aes(x=Species,y=stat+0.02,label=p.value,group=group),color="black",position=position_dodge(width=1),vjust=0.6)+
  #scale_fill_brewer(palette="Paired",name="Association")+
  scale_fill_manual(values=c(hcl(h=0,c=200,l=33),hcl(h=0,c=100,l=10),hcl(h=0,c=100,l=50),hcl(h=0,c=100,l=100),
                             hcl(h=256,c=200,l=33),hcl(h=256,c=100,l=10),hcl(h=256,c=100,l=50),hcl(h=256,c=100,l=100),
                             hcl(h=0,c=0,l=10),hcl(h=0,c=0,l=50),hcl(h=0,c=0,l=100)))+
  scale_x_discrete(name="Diet Item")+
  scale_y_continuous(name="ISA Stat (A*B)",limits=c(0,1),expand=expansion(0))+
  coord_flip()+
  theme(legend.position=c(0.8,0.75),legend.background=element_rect(color="black"),
        plot.margin=margin(15,20,15,15,unit="pt"),panel.grid.major.y=element_blank())+
  annotation_custom(legendGrob,xmin=7.5,xmax=20,ymin=-0.15,ymax=1.4)


#Making it as a pretty table
tableISA<-fullISA%>%
  mutate(stat=ifelse(stat==0,"",round(stat,4)),
         stat.sig=gsub("NA","",paste0(stat,p.value.sig)))%>%
  pivot_wider(id_cols="Species",names_from="group",values_from="stat.sig")%>%
  arrange(Species)

library(reactable)
reactable(tableISA,
          defaultColDef = colDef(
            style = function(value) {
             if (grepl("\\*",value)) {
                 color <- "black"
                 weight <- 600
             } else {
                 color <- "grey50"
                 weight <- 300
             }
               list(color = color, fontWeight = weight)
             },
            minWidth = 90,align="center"),
          columns = list(
            Species = colDef(minWidth=155,align="left")),
          
          defaultPageSize = 20,highlight = T,striped=T)



# PerMANOVA and post-hocs ----------------------------------------------------------

#PERMANOVA for the interaction of Year and fortnight 
set.seed(42)
adonis(original.dist~as.character(Year)*family,
       data=allEnv.Mat,permutations=1000,method="bray")

#significant pairwise test of the time period
pairwise.perm.manova(original.dist,paste(allEnv.Mat$Year,allEnv.Mat$family),nperm=1000)

fdisp<-betadisper(original.dist,allEnv.Mat$Year)
fdisp
permutest(fdisp)

library(usedist)
dist_multi_centroids(original.dist,allEnv.Mat$Year)


ydisp<-betadisper(original.dist,allEnv.Mat$family)
ydisp
permutest(ydisp)
dist_multi_centroids(original.dist,allEnv.Mat$family)



sydisp<-betadisper(original.dist,paste(allEnv.Mat$family,allEnv.Mat$Year,sep="."))
sydisp
permutest(sydisp)
centroidDists<-dist_multi_centroids(original.dist,paste(allEnv.Mat$family,allEnv.Mat$Year))
centroidDists
permutest(centroidDists)

charCentroids<-c(0.5168236,0.6386732,0.4974933)
sculpCentroids<-c(0.4926337,0.2285743,0.3625625)
sameYearCentroids<-c(0.2207148,0.3005605,0.5710684)
interSpeciesCentroids<-c(sameYearCentroids,0.4323695,0.6286381,0.5417022,0.6648193,0.3000332,0.3775107)
mean(charCentroids)
mean(sculpCentroids)
mean(c(charCentroids,sculpCentroids))
mean(sameYearCentroids)
mean(interSpeciesCentroids)

centroidDists%>%
  as.matrix()%>%as.data.frame()%>%
  mutate(species1=labels(centroidDists))%>%
  pivot_longer(cols=labels(centroidDists),names_to = "species2", values_to = "dist")%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=dist))+
  geom_text(aes(species1,species2,label=round(dist,digits=2)),size=8)+
  geom_rect(aes(xmin=0.5,ymin=0.5,xmax=3.5,ymax=3.5),fill="transparent",color="black",size=2)+
  geom_rect(aes(xmin=3.5,ymin=0.5,xmax=6.5,ymax=3.5),fill="transparent",color="black",size=2)+
  geom_rect(aes(xmin=3.5,ymin=3.5,xmax=6.5,ymax=6.5),fill="transparent",color="black",size=2)+
  geom_rect(aes(xmin=3.5,ymin=0.5,xmax=4.5,ymax=1.5),fill="transparent",color="red",size=2)+
  geom_rect(aes(xmin=4.5,ymin=1.5,xmax=5.5,ymax=2.5),fill="transparent",color="red",size=2)+
  geom_rect(aes(xmin=5.5,ymin=2.5,xmax=6.5,ymax=3.5),fill="transparent",color="red",size=2)+
  geom_rect(aes(xmin=0.5,ymin=3.5,xmax=1.5,ymax=4.5),fill="transparent",color="red",size=2)+
  geom_rect(aes(xmin=1.5,ymin=4.5,xmax=2.5,ymax=5.5),fill="transparent",color="red",size=2)+
  geom_rect(aes(xmin=2.5,ymin=5.5,xmax=3.5,ymax=6.5),fill="transparent",color="red",size=2)+
  scale_fill_viridis_c(option="B",name="Centroid\ndistance",na.value="white")+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
#Note that the centroids between char are further apart on average (0.55) than the sculpin are (0.36) 
  #or the two species are from each other in the same year (0.36) or even generally (0.44)

bind_cols(group=sydisp[["group"]],distance=sydisp[["distances"]])%>%
  separate(group,into=c("Species","Year"),sep="\\.")%>%
  ggplot()+
  geom_boxplot(aes(x=Year,y=distance,fill=Species),outlier.shape=NA,size=1.5,show.legend = F)+
  geom_point(aes(x=Year,y=distance,fill=Species,color=Species),
             shape=21,size=5,position=position_jitterdodge(jitter.width=0.1,seed=42))+
  scale_fill_viridis_d(option="E",name="Consumer")+
  scale_color_manual(values=c("white","black"),guide="none")+
  scale_y_continuous(name="Distance to Centroid")
  


# Analyses accumulation -----------------------------------------------------------

##prey accumulation curves using dietmat dataset

sculppaa <- subset(dietMat, Species == 'Fourhorn' | Species == 'Slimy')
charpaa <- subset(dietMat, Species == 'Arctic char')

##needs to be dataframe for accumcomp to work
sculpenv <- data.frame(sculppaa[,1:3])
charenv <- data.frame(charpaa[,1:3])

sculpcom <- data.frame(sculppaa[,13:31])
charcom <- data.frame(charpaa[,13:31])

##creation of data must be done separately as far as I know
##sculpin
Accum.sculp <- accumcomp(sculpcom, y=sculpenv, factor = 'Year',
                     method ='exact', conditioned = F, plotit = F)
Accum.sculp
accum.longsculp <- accumcomp.long(Accum.sculp, ci = NA, label.freq = 5)
head(accum.longsculp)

##test for sculpin
sculpenv2 <- sculpenv
sculpenv2$test <- 1
Accum.sculp2 <- accumcomp(sculpcom, y = sculpenv2, factor = 'test', method = 'exact', conditioned = F, plotit = F)

accum.longsculp2 <- accumcomp.long(Accum.sculp2, ci = NA, label.freq = 5)
accum.longsculp2$Species <- 'Sculpin'

##char
Accum.char <- accumcomp(charcom, y=charenv, factor = 'Year',
                         method ='exact', conditioned = F, plotit = F)
Accum.char
accum.longchar <- accumcomp.long(Accum.char, ci = NA, label.freq = 5)
head(accum.longchar)

##test for char
charenv2 <- charenv
charenv2$test <- 1
Accum.char2 <- accumcomp(charcom, y = charenv2, factor = 'test', method = 'exact', conditioned = F, plotit = F)

accum.longchar2 <- accumcomp.long(Accum.char2, ci = NA, label.freq = 5)
accum.longchar2$Species <- 'Arctic char'


accum.longchar$Species <- 'Arctic char' ##creating a column for species so we can bind together
accum.longsculp$Species <- 'Sculpin'


accum.long <- rbind(accum.longchar, accum.longsculp)
accum.long2 <- rbind(accum.longchar2, accum.longsculp2)

##plot
bind_rows(accum.long,accum.long2%>%mutate(Grouping="Cumulative"))
PAAplot <- ggplot(data = accum.long, aes(x = Sites, y = Richness,ymax = UPR, ymin = LWR))+
  facet_wrap(~Species, ncol = 1)+
  geom_line(aes(color = Grouping), size = 2,show.legend = F)+
  geom_point(data=subset(accum.long, labelit==TRUE),
             aes(colour=Grouping),size = 5)+
  scale_color_manual(values = c('black', 'grey50', 'grey80','black'))+
  geom_line(data = accum.long2, aes(x = Sites, y = Richness,color="Cumulative"), lty = 2, size = 2)+
  labs(x = 'Number of Samples', y = 'Species Richness', color = 'Year')+
  guides(color=guide_legend(override.aes=list(linetype=c(1,1,1,2))))



bind_rows(accum.long%>%mutate(level="a"),accum.long2%>%mutate(level="Cumulative"))%>%
  ggplot(aes(x = Sites, y = Richness,ymax = UPR, ymin = LWR))+
  facet_wrap(~Species, ncol = 1)+
  geom_line(aes(color = Grouping,linetype=level), size = 2)+
  geom_point(data=subset(accum.long, labelit==TRUE),
             aes(colour=Grouping),size = 5)+
  scale_color_manual(values = c('black', 'grey50', 'grey80','black'),labels=c("Cumulative","2017","2018","2019"))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,19,19,19),linetype=c(3,1,1,1))),
         linetype="none")+
  labs(x = 'Number of Samples', y = 'Species Richness', color = 'Year')+
  theme(legend.key.width = unit(1.7,"cm"))

##only do this once
# tiff('/Users/larshammer/Hammer_Hermann/Figures/PreyAccumulationCurve.tiff', width = 6.5, height = 5.5, res = 300, units = 'in')

PAAplot

dev.off()

# relative consumption GLM work using summaries dataset ---------------------------------

library(car)
library(MuMIn)
library(lme4)
library(lmerTest)
##removing unneccesary columns and unusable data
relcon <- summaries[,1:9] ## unneccessary columns
relcon <- subset(relcon, !is.na(Diet_Mass_g)) ## removing fish with unknown diet masses
relcon <- subset(relcon, !is.na(Mass_g)) ## removing fish with unknown masses
n_distinct(relcon$Fish_ID) ##127 useable fish for relcon analysis

##relative consumption column
relcon$relative_consumption <- (relcon$Diet_Mass_g/relcon$Mass_g)*100

relcon$Date_pos <- as.POSIXct(relcon$Date, format = '%m/%d/%y')
relcon$month <- factor(month(relcon$Date_pos)) ##make sure month is a factor
relcon$Year <- factor(relcon$Year) ## make sure year is a factor
relcon$species2 <- factor(ifelse(relcon$Species == 'Slimy' | relcon$Species == 'Fourhorn', 'Sculpin', 'Arctic char'))


relcon %>% 
  group_by(species2, Year) %>% 
  summarize(mean = mean(relative_consumption),
            sd = sd(relative_consumption))

relcon %>% group_by(species2) %>% 
  summarize(min = min(relative_consumption), max = max(relative_consumption),
            mean = mean(relative_consumption), se = sd(relative_consumption)/sqrt(n()),
            sd = sd(relative_consumption))

##running a hurdle model

##first we need to determine factors affecting whether we feed or not
relcon$feed <- ifelse(relcon$relative_consumption > 0, 1, 0)

table(relcon$feed, relcon$species2) ##looks right at least for char
relcon%>%ungroup()%>%group_by(species2)%>%mutate(N=n())%>%group_by(species2,feed)%>%reframe(p=n()/N*100)%>%distinct()

##don't bother with mass because species are so different in weight already
bingelm1 <- glm(feed~species2 + Year + month, data = relcon, na.action = 'na.fail', family = 'binomial')
vif(bingelm1)
dredge(bingelm1) ##species is the only thing that comes out might be best as a table

##now removing those who didn't feed
relcon2 <- subset(relcon, feed == 1)

bingelm2 <- glm(relative_consumption ~ species2 + Year + month, data = relcon2, 
                na.action = 'na.fail', family = 'Gamma')
vif(bingelm2)
dredge(bingelm2)


relconplot <- ggplot(relcon2)+
  geom_boxplot(aes(x = Year, y = relative_consumption), outlier.shape = NA)+
  geom_point(aes(x = Year, y = relative_consumption), position = position_jitter(width = .2))+
  facet_grid(~species2)+
  labs(y = 'Relative Consumption (%)')

##only do this once
# tiff('/users/larshammer/Hammer_Hermann/Figures/RelativeConsumption_year.tiff', width = 6.5, height = 5.5, res = 300, units = 'in')

relconplot

dev.off()

##looking at mass
ggplot(relcon)+
  geom_point(aes(x = Mass_g, y = relative_consumption), position = position_jitter(width = .2))+
  facet_grid(~species2, scales = "free_x")+
  labs(y = 'Relative Consumption (%)')

##looking at month
ggplot(relcon)+
  geom_boxplot(aes(x = month, y = relative_consumption), outlier.shape = NA)+
  geom_point(aes(x = month, y = relative_consumption), position = position_jitter(width = .2))+
  facet_grid(~species2)+
  labs(y = 'Relative Consumption (%)')




# Binge-feeding, just a simple assessment from other species -----------------------------------------------------------

#Model formulas for calculating mass and temperature dependent Cmax
Cmax<-function(df,W) {
  Cmax=df$CA*W^df$CB
}
Ft<-function(df,T) {
  G1=(1/(df$CTO-df$CQ))*log((0.98*(1-df$CK1))/(df$CK1*0.02))
  G2=(1/(df$CTL-df$CTM))*log((0.98*(1-df$CK4))/(df$CK4*0.02))
  L1=exp(G1*(T-df$CQ))
  L2=exp(G2*(df$CTL-T))
  Ka=(df$CK1*L1)/(1+df$CK1*(L1-1))
  Kb=(df$CK4*L2)/(1+df$CK4*(L2-1))
  F=Ka*Kb
}
#Data (Deslauriers et al. 2017 provide the dataframe but the values from Mesa and Moss)
bioenerg<-read_csv("Parameters_official.csv")


#Size range for Arctic char and Sculpins
summaries%>%
  mutate(family=ifelse(Species=="Char","Arctic char","Sculpin"))%>%
  group_by(family)%>%
  summarise(Mass=range(Mass_g,na.rm=T),
            mean=mean(Mass_g,na.rm=T))

#Estimated Cmax for those max and min sizes
print(100*Ft(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),c(-2.5,-2.5,11,11))*Cmax(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),c(8.2,6600)))

print(100*Ft(filter(bioenerg,Sci_Name=="Cottus asper"),c(-2.5,-2.5,11,11))*Cmax(filter(bioenerg,Sci_Name=="Cottus asper"),c(0.8,400)))


#For the mean size at a maximum temperature
print(100*Ft(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),11)*Cmax(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),2609))

print(100*Ft(filter(bioenerg,Sci_Name=="Cottus asper"),11)*Cmax(filter(bioenerg,Sci_Name=="Cottus asper"),78.9))

#All individual values
print(100*Ft(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),11)*Cmax(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),
                                                                          filter(summaries,Species=="Char")$Mass_g))
print(100*Ft(filter(bioenerg,Sci_Name=="Cottus asper"),11)*Cmax(filter(bioenerg,Sci_Name=="Cottus asper"),
                                                                filter(summaries,Species!="Char")$Mass_g))


#Comparing with the actual relcon for each individual
relcon%>%
  mutate(Cmax=ifelse(species2=="Sculpin",
                     100*Ft(filter(bioenerg,Sci_Name=="Cottus asper"),11)*Cmax(filter(bioenerg,Sci_Name=="Cottus asper"),Mass_g),
                     100*Ft(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),11)*Cmax(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),Mass_g)))%>%
  ggplot()+
  geom_point(aes(x=relative_consumption,y=Cmax,fill=species2),shape=21,size=5)+
  geom_abline(aes(slope=1,intercept=0))+
  scale_x_continuous(limits=c(0,11))+
  scale_y_continuous(limits=c(0,11))+
  scale_fill_viridis_d(option="E")+
  theme(legend.position = c(0.5,0.8))
  
relcon%>%
  mutate(Cmax=ifelse(species2=="Sculpin",
                     100*Ft(filter(bioenerg,Sci_Name=="Cottus asper"),11)*Cmax(filter(bioenerg,Sci_Name=="Cottus asper"),Mass_g),
                     100*Ft(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),11)*Cmax(filter(bioenerg,Sci_Name=="Salvelinus confluentus"),Mass_g)),
         exceed=ifelse(Cmax<=relative_consumption,"y","n"))%>%
  group_by(species2)%>%
  summarise(N=n(),
            N_exceed=sum(exceed=="y"),
            p_exceed=N_exceed/N)
(8+39)/(52+75)
