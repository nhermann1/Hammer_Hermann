
#Most of this is just data prep and won't be included in actual analyses

rm(list=ls())

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/UNH Research/Data/CSV_Files/")
set.seed(42)

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
library(vegan)
library(MuMIn)



#Colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#To use, add command for scale_fill_manual(values=cbPalette)

theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')


#Functions to use later
cat4<-function(a,b) { #a is the category 3 with 16 levels and b is the category 4 with 8 (each a column)
  if (a=="Gammarus sp."|a=="Onisimus sp."|a=="Gammaracanthus sp."|a=="Amphipod") {
    b="Amphipod"
  } else if (a=="Krill"|a=="Mysid") {
    b="Krill/Mysid"
  } else if (a=="Jellyfish"|a=="Sea Angel") {
    b="Jellyfish/Sea Angel"
  } else if (a=="Fish"|a=="Fish Bits") {
    b="Fish"
  } else if (a=="Copepod") {
    b="Copepod"
  } else if (a=="Miscellaneous Invert"|a=="Chironomid") {
    b="Miscellaneous Invert"
  } else if (a=="Algae"|a=="Undigestible") {
    b="Undigestible"
  } else {
    b="Digested"
  }
}
cat5<-function (a,b) { #a is category 4 with 8 and b is category 5 with 6
  if (a=="Amphipod") {
    b="Amphipod"
  } else if (a=="Krill/Mysid") {
    b="Krill/Mysid"
  } else if (a=="Fish") {
    b="Fish"
  } else if (a=="Copepod"|a=="Miscellaneous Invert"|a=="Jellyfish/Sea Angel") {
    b="Miscellaneous Invert"
  } else  if (a=="Undigestible") {
    b="Undigestible"
  } else {
    b="Digested"
  }
}



#Niche breadth
B<-function(p) {
  a<-apply(p,2,sum)/sum(p)
  b<-(1/(sum(a^2))-1)/(length(a)-1)
  return(b)
}

myboot <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  list(std.dev=std.dev, levins=mean(levin))   
}
#Usage: myboot(df,999)
myboot_reps <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  return(levin)   
}




# Data Input --------------------------------------------------------------

#Looking at the Diet stomachs
#2019
stomachs2019<-read.csv("2019_Sculpin_Dissections.csv",stringsAsFactors = F)

sculpStomachs2019<-stomachs2019[1:17,]%>%
  mutate(Diet_Mass=Stomach_Mass_Full-Stomach_Mass_Empty,
         Consump=Diet_Mass/Mass*100)
sculpStomachs2019$Year<-"2019"




#2018
dissections2018<-subset(read.csv("2018_Fish_Data.csv",stringsAsFactors = F),grepl("SSM[0-9]+",Fish.ID))
colnames(dissections2018)[colnames(dissections2018)=="Fish.ID"]<-"FishID"
stomachs2018<-subset(read_csv("2018_StomachMass.csv", col_types = cols(StomachContentMass = col_number(),
                                                                       StomachMass = col_number())),
                     grepl("SSM[0-9]+",FishID))
stomachs2018<-merge(dissections2018,stomachs2018,by.x="FishID",by.y="FishID")
stomachs2018<-dplyr::select(stomachs2018,FishID,Date,Species,TL,FL,Mass,Sex.,Stomach.Weight,StomachContentMass,StomachMass,GutFullness,Comments.x,Comments.y,Liver.Weight,Heart.Weight,Gonad.Weight)
stomachs2018<-stomachs2018[!is.na(stomachs2018$StomachMass),]
stomachs2018[is.na(stomachs2018)]<-0
stomachs2018$Stomach.Weight<-as.numeric(stomachs2018$Stomach.Weight)
stomachs2018$Year<-"2018"

#Conserved mass because the diets were weighed and picked through in different places/times
stomachs2018<-stomachs2018%>%
  mutate(Conserved_Mass=(StomachContentMass+StomachMass)/Stomach.Weight*100, #How much water was lost
         diffDietMass=Stomach.Weight-StomachMass) #what is the diet based on the full stomach-empty stomach

#Make diet equal to the average of the measured mass (SCM) and the mass from the difference of the empty stomach (dDM)
stomachs2018$Diet_Mass<-(stomachs2018$StomachContentMass+stomachs2018$diffDietMass)/2
#If it was empty, make it zero
stomachs2018$Diet_Mass<-ifelse(grepl("^empty",stomachs2018$Comments.y),0,stomachs2018$Diet_Mass)

#Find the relative consumption for each individual (1 didn't get fish weighed, so missing consumption)
stomachs2018<-stomachs2018%>%
  mutate(diffConsump=diffDietMass/as.numeric(Mass)*100,
         Consump=Diet_Mass/as.numeric(Mass)*100)

stomachs2018$TL<-as.numeric(stomachs2018$TL)
stomachs2018$Mass<-as.numeric(stomachs2018$Mass)
stomachs2018$Sex.<-gsub("\\?","",stomachs2018$Sex.)
stomachs2018$Sex.<-gsub("unkown","Unk",stomachs2018$Sex.) #Remove from the data if it's unknown

sculpStomachs2018<-dplyr::select(stomachs2018,
                                 Fish_ID=FishID,Date,Species,TL_mm=TL,Mass_g=Mass,Sex=Sex.,
                                 Diet_Mass_g=StomachContentMass,GutFullness,Year,
                                 Liver_Mass_g=Liver.Weight,Gonad_Mass_g=Gonad.Weight,Heart_Mass_g=Heart.Weight)%>%
  mutate(Species=ifelse(Species=="slimy","Slimy","Fourhorn"),
         Liver_Mass_g=as.numeric(Liver_Mass_g),
         Gonad_Mass_g=as.numeric(Gonad_Mass_g),
         Heart_Mass_g=as.numeric(Heart_Mass_g),
         Date=dmy(Date))


#2017, but with rows for each fish
sculpStomachs2017<-read.csv("2017_Just_Sculpin.csv",stringsAsFactors = F)%>%
  filter(Gut.fullness..1.5.>0)%>%
  mutate(Year="2017",
         Diet_Mass_g=NA,
         Heart_Mass_g=NA)

sculpStomachs2017<-dplyr::select(sculpStomachs2017,
                                   Fish_ID=Sample.ID,Date,Species=Common,TL_mm=Total.L..cm.,Mass_g=Wet.mass..g.,
                                   Sex,Diet_Mass_g,GutFullness=Gut.fullness..1.5.,Year,
                                   Liver_Mass_g=Liver.mass..g.,Gonad_Mass_g=Gonad.mass..g.,Heart_Mass_g)%>%
  mutate(Species=ifelse(Species=="Slimy sculpin","Slimy","Fourhorn"),
         Liver_Mass_g=as.numeric(Liver_Mass_g),
         Gonad_Mass_g=as.numeric(Gonad_Mass_g),
         Date=dmy(Date))


sculpStomachs2019<-sculpStomachs2019%>%
  mutate(GutFullness=NA)%>%
  dplyr::select(Fish_ID,Date,Species,TL_mm=TL,Mass_g=Mass,
                Sex,Diet_Mass_g=Diet_Mass,GutFullness,Year,Liver_Mass_g=Liver_Mass,Gonad_Mass_g=Gonad_Mass,Heart_Mass_g=Heart_Mass)%>%
  mutate(Species=ifelse(Species=="Slimy","Slimy","Fourhorn"),
         Date=dmy(Date))

sculpStomachs<-bind_rows(sculpStomachs2017,sculpStomachs2018,sculpStomachs2019)%>%
  mutate(GutFullness=as.numeric(substr(GutFullness,1,1)),
         Sex=ifelse(Sex=="Unk","Unknown",Sex),
         Sex=factor(Sex,levels=c("F","M","I","Unknown")))

#write.csv(sculpStomachs,"allSculpin_dietSummaries_final.csv",row.names = F)















# Read-in Data ------------------------------------------------------------

sculpStomachs<-read_csv("Hammer_Hermann/allSculpin_dietSummaries_final.csv")



detail2017stomachs <- read_csv("R CSVs/2017_fullStomachData.csv", 
                               col_types = cols(Count = col_number(), 
                                                Date = col_date(format = "%d/%m/%Y"), 
                                                Fyke_Num = col_number(), 
                                                GutFullness = col_character(), 
                                                Mass = col_number(), 
                                                Mass_Fish = col_number(), 
                                                Mass_Gonad = col_number(), 
                                                Mass_Liver = col_number(), 
                                                SL_Fish = col_number(), 
                                                TL_Fish = col_number(), 
                                                Year = col_character()))


detail2018stomachs <- read_csv("2018_StomachContents.csv")%>%
  mutate(Fish_ID=gsub("A","a",FishID),
         Fish_ID=gsub("B","b",Fish_ID),
         Count=as.numeric(EstimatedNumber),
         Category=case_when(grepl("gammarus",Taxa,ignore.case=T)~"Gammarus sp.",
                            grepl("onisimus",Taxa,ignore.case=T)~"Onisimus sp.",
                            grepl("mysid",Taxa,ignore.case=T)~"Mysid",
                            grepl("krill",Taxa,ignore.case=T)~"Krill",
                            grepl("amphipod",Taxa,ignore.case=T)~"Amphipod",
                            Taxa=="Unidentified"~"Digested",
                            Taxa=="Copepod"~"Miscellaneous",
                            TRUE~Taxa),
         Year="2018")%>%
  filter(grepl("SSM",FishID))%>%
  dplyr::select(Fish_ID,Year,Content=Taxa,Category,Count,Mass_g=TotalMass,Comments)%>%
  full_join(sculpStomachs[which(sculpStomachs$Year==2018),],by=c("Fish_ID","Year"))


detail2019stomachs <- read_csv("R CSVs/2019_fullStomachData.csv",
                               col_types = cols(Date = col_date(format = "%d/%m/%Y"),
                                                GutFullness = col_character(), 
                                                Mass_Fish = col_number(), 
                                                Mass_Gonad = col_number(), 
                                                Mass_Liver = col_number(), 
                                                SL_Fish = col_number(), 
                                                TL_Fish = col_number(), 
                                                Year = col_character()))


#Making the different categories for all the diets
detail2017stomachs$Category4<-sapply(detail2017stomachs$Category,cat4)
detail2017stomachs$Category5<-sapply(detail2017stomachs$Category4,cat5)
#And 2018
detail2018stomachs$Category4<-sapply(detail2018stomachs$Category,cat4)
detail2018stomachs$Category5<-sapply(detail2018stomachs$Category4,cat5)
#And 2019
detail2019stomachs$Category4<-sapply(detail2019stomachs$Category,cat4)
detail2019stomachs$Category5<-sapply(detail2019stomachs$Category4,cat5)






# Species Accumulation ----------------------------------------------------
detailAll<-bind_rows(dplyr::select(detail2017stomachs,Fish_ID,Year,Date,Content,Category,Count,Mass),
                     dplyr::select(detail2018stomachs,Fish_ID,Year,Date,Content,Category,Count,Mass=Mass_g.x),
                     dplyr::select(detail2019stomachs,Fish_ID,Year,Date,Content,Category,Count,Mass))%>%
  mutate(Content=tolower(Content),
         Category=ifelse(grepl("Fish",Category),"Fish",Category))%>%
  #filter(Category!="Empty")%>%
  full_join(sculpStomachs[,which(colnames(sculpStomachs)!="Date")])%>%
  mutate(Content=case_when(is.na(Content) & Diet_Mass_g>0~"Unidentified",
                           is.na(Content) & Diet_Mass_g==0~"Empty",
                           TRUE~Content),
         Category=case_when(is.na(Category) & Diet_Mass_g>0~"Unidentified",
                            is.na(Category) & Diet_Mass_g==0~"Empty",
                            TRUE~Category))

#write.csv(detailAll,"Hammer_Hermann/allSculpin_dietDetails_Contents_final.csv",row.names = F)

toPA<-function(x) {
  x<-as.numeric(gsub("[A-z \\.]+",1,x))
  x<-ifelse(is.na(x),0,1)
}

dietsWide<-pivot_wider(detailAll,id_cols=c("Fish_ID","Year","Date","Species","TL_mm","Mass_g","Diet_Mass_g",
                                           "GutFullness","Sex","Liver_Mass_g","Gonad_Mass_g","Heart_Mass_g"),
                       names_from="Category",names_prefix="prey_",
                       values_from = "Category",values_fn=min)%>%
  mutate_at(vars(matches("prey")),toPA)

#write.csv(dietsWide,"Hammer_Hermann/allSculpin_PAdietMat.csv",row.names=F)


#Remove empty and unidentified (that means they weren't analyzed, not the same as digested)
dietsWide_nEmpty<-filter(dietsWide,prey_Empty==0)%>%dplyr::select(-prey_Empty)
dietsWide_final<-filter(dietsWide_nEmpty,prey_Unidentified==0)%>%dplyr::select(-prey_Unidentified)

#write.csv(dietsWide_final,"Hammer_Hermann/dietSculpin_PAdietMat_final.csv",row.names=F)






# Editing the full datasets (with Char) -----------------------------------

rm(list=ls())
setwd("../../Hammer_Hermann/Data/") #Nate's file structure
setwd("") #Lars' file structure


details<-read_csv("all_dietcontents_final.csv")%>%
  mutate(Date=mdy(Date),
         Year=as.character(Year),
         GutFullness=as.factor(GutFullness),
         Category=case_when(grepl("jelly",Content)~"Jellyfish",
                            Content=="Sea Angel"~"Sea Angel",
                            grepl("copepod",Content)|grepl("lice",Content)~"Copepod",
                            grepl("flea",Content)~"Amphipod",
                            Category=="Unidentified" & Species=="Arctic char"~"Digested",
                            TRUE~Category))
#write.csv(details,"all_dietcontents_final_v2.csv",row.names=F)
summaries<-read_csv("all_dietsummaries_final.csv")%>%
  mutate(Date=mdy(Date),
         Year=as.character(Year),
         GutFullness=as.factor(GutFullness))

dietMat<-read_csv("all_diet_PAdietMat_final.csv")%>%
  mutate(Date=mdy(Date),
         Year=as.character(Year),
         GutFullness=as.factor(GutFullness))
dietMat2<-pivot_wider(details,id_cols=c("Fish_ID","Year","Date","Species","TL_mm","Mass_g","Diet_Mass_g",
                                            "GutFullness","Sex","Liver_Mass_g","Gonad_Mass_g","Heart_Mass_g"),
                        names_from="Category",names_prefix="prey_",
                        values_from = "Category",values_fn=min)%>%
  mutate_at(vars(matches("prey")),toPA)%>%
  filter(prey_Empty==0)%>%dplyr::select(-prey_Empty)%>%
  filter(prey_Unidentified==0)%>%dplyr::select(-prey_Unidentified)
colnames(dietMat2)<-gsub(" ",".",colnames(dietMat2))
#write.csv(dietMat2,"all_diet_PAdietMat_final_v2.csv",row.names=F)


#How many do we have for the major groups
table(summaries$Species,summaries$Year)
#How many do we have relative consumption for
table(filter(summaries,!is.na(Diet_Mass_g)&!is.na(Mass_g))$Species,
      filter(summaries,!is.na(Diet_Mass_g)&!is.na(Mass_g))$Year)







# Actual data ---------------------------------------------------
set.seed(42)


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
library(vegan)
library(MuMIn)



#Colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#To use, add command for scale_fill_manual(values=cbPalette)

theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

details<-read_csv("Data/all_dietcontents_final_v2.csv")
summaries<-read_csv("Data/all_dietsummaries_final.csv")
dietMat<-read_csv("Data/all_diet_PAdietMat_final_v2.csv")




# Analyses--Nate ----------------------------------------------------------

#### Frequency of Occurrence

#Cumulative across years
SculpinFOO<-filter(dietMat,Species %in% c("Slimy","Fourhorn"))%>%
  dplyr::select(matches("^prey"))%>%
  colSums()/nrow(dietMat%>%filter(Species %in% c("Slimy","Fourhorn")))

CharFOO<-filter(dietMat,Species=="Arctic char")%>%
  dplyr::select(matches("^prey"))%>%
  colSums()/nrow(dietMat%>%filter(Species=="Arctic char"))

foo<-bind_rows(data.frame(foo=CharFOO,species="Arctic char",prey=gsub("prey_","",colnames(dplyr::select(dietMat,matches("^prey"))))),
               data.frame(foo=SculpinFOO,species="Sculpin",prey=gsub("prey_","",colnames(dplyr::select(dietMat,matches("^prey"))))))%>%
  mutate(prey=gsub("\\.","\n",prey))

ggplot(foo,aes(prey,foo,fill=species))+
  geom_col(position="dodge")

#Separated by years
presenceYears<-dietMat%>%
  mutate(species=ifelse(Species=="Arctic char","Arctic char","Sculpin"))%>%
  dplyr::select(species,Year,matches("^prey"))%>%
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
         foo=noo/N,
         Year=as.character(Year))

#Emphasizing species differences
ggplot(fooYears,aes(prey,foo,fill=Year))+
  geom_col(position="dodge",color="black")+
  scale_fill_viridis_d(name="Species")+
  scale_y_continuous(name="Frequency of Occurrence",limits=c(0,1),expand=expansion(add=0))+
  scale_x_discrete(name="Diet Item")+
  theme(axis.text.x=element_text(angle=20,hjust=1,vjust=1.12),
        axis.title.x=element_text(vjust=5))+
  facet_wrap(~species,nrow=2)

#Emphasizing inter-annual variation
ggplot(fooYears,aes(prey,foo,fill=species))+
  geom_col(position="dodge",color="black")+
  scale_fill_viridis_d(name="Species")+
  scale_y_continuous(name="Frequency of Occurrence",limits=c(0,1),expand=expansion(add=0))+
  scale_x_discrete(name="Diet Item")+
  theme(axis.text.x=element_text(angle=20,hjust=1,vjust=1.12),
        axis.title.x=element_text(vjust=5))+
  facet_wrap(~Year,nrow=3)
