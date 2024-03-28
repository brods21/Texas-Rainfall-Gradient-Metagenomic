
# this is the data setup from the eggnog_KEGG_microtraitKOs_relabun_withcovariates_72623.r

# 1-19-24
rm(list=ls())


library(mvabund)
library(tidyverse)
library(lme4)
library(emmeans)
library(vegan)
library(ggrepel)
library(gridExtra)
library (edgeR)
library(lterpalettefinder)
library(visreg)
library(MuMIn)
library(car)
library(dplyr)
library(ggthemes)
library(RRPP)


setwd("~/Documents/Texas gradient")

data <-read.csv( "proteins_eggnog_nobbmap.csv") 
#sitedata <-read.csv( "Texas site characteristics.csv")
trt <- read.csv("Texas_Metadata.csv")
timedata <- read.csv("TXGrad_RAPID_sitedata.csv")

enzymes <- read.csv("TXRAPID_Enz_forCaitlin.csv")
co2 <- read.csv("TXGrad_AvgHighCO2_2024_01_09.csv")

microtrait <- read.csv("~/Documents/Texas gradient/rules_microtrait_2023-05-16.csv")
#sporegene <- read.csv("Lietal2022_MolEcol_S2_SporulationGenes.csv")


bothmeta <- merge (trt, timedata, by = c ("Site", "SamplingPt"))

ggplot(bothmeta, aes ( x = MAP, y = PercC, group = Site)) + 
  geom_point() + 
  theme_light() + 
  facet_grid(~SamplingPt) # percent C across the gradient

bothmeta2 <- bothmeta # this will just be for getting the standardized MAP for initial
# models for soil properties,
bothmeta2$MAP.std <- ( bothmeta2$MAP - mean(bothmeta2$MAP) ) / (2 * sd ( bothmeta2$MAP))

ggplot(bothmeta, aes ( x = MAP, y = PercC )) +  # combining the two wettest sites 
  stat_summary(aes (y = PercC) , geom = "line", fun.y = "mean", col = "red") + 
  stat_summary(aes (y = PercN * 10) , geom = "line", fun.y = "mean", col = "green") + 
  theme_light() 

ggplot(bothmeta, aes ( x = MAP, y = pH, group = Site, color = SamplingPt)) +  # pH
  geom_point() + 
  #geom_smooth(method = "lm") +
  theme_light() 


hist(bothmeta$PercMoist)
ggplot(data = bothmeta, aes (x = PercMoist, fill = SamplingPt)) + 
  theme_light() +
  geom_histogram( alpha = 0.5) #+ 
#facet_grid(~SamplingPt)


names(bothmeta)
dim(bothmeta)

# how does soil moisture vary with MAP at each sampling point?
ggplot(bothmeta, aes ( x = MAP, y = PercMoist , group = Site)) +  # combining the two wettest sites 
  geom_point() + 
  theme_light()  +
  facet_grid(~SamplingPt)

summary(lm(PercMoist~MAP.std*SamplingPt, bothmeta2))
# This shows that percMoist varies with SamplingPoint 
summary(lm(PercMoist~MAP, bothmeta)) # marg sig
summary(aov(PercMoist~SamplingPt, bothmeta)) # this is highly significant
TukeyHSD(aov(PercMoist~SamplingPt, bothmeta))
bothmeta %>% group_by(SamplingPt) %>% summarise(Moist = mean(PercMoist))
#SamplingPt  Moist
#<chr>       <dbl>
#  1 Fall 2015    23.6
#2 Spring 2016  14.7
#3 Summ 2016    12.5

(14.7+12.5) / 2  # 13.6 is average for later two months.
(23.6- 14.7) / 23.6 * 100


summary(lm(PercMoist~MAP.std, bothmeta2) )  # marginally significant. 

summary(lm(PercC~MAP.std*SamplingPt, bothmeta2)) # nothing.
summary(lm(PercN~MAP.std*SamplingPt, bothmeta2)) # nothing
summary(lm(PercSand~MAP.std*SamplingPt, bothmeta2)) # nothing
summary(lm(MBC_mg_g~MAP.std*SamplingPt, bothmeta2)) # SamplingPoint
summary(lm(DOC_mg_g~MAP.std*SamplingPt, bothmeta2)) # SamplingPoint
summary(lm(pH~MAP.std*SamplingPt, bothmeta2)) # MAP significant !! .
range(bothmeta$pH)

plot(bothmeta2$DOC_mg_g, bothmeta2$MBC_mg_g) # not really correlated. ? 

summary(lm(PercC~MAP, bothmeta))
summary(lm(PercN~MAP, bothmeta))
summary(lm(PercSand~MAP, bothmeta))
summary(lm(MBC_mg_g~MAP, bothmeta))
summary(lm(pH~MAP, bothmeta))

summary(aov(PercC~SamplingPt, bothmeta))
summary(aov(PercN~SamplingPt, bothmeta))
summary(aov(PercSand~SamplingPt, bothmeta))
summary(aov(MBC_mg_g~SamplingPt, bothmeta))
TukeyHSD(aov(MBC_mg_g~SamplingPt, bothmeta)) # fall diff than sspring 
bothmeta %>% group_by(SamplingPt) %>% summarise(MBC = mean(MBC_mg_g))
# very high in spring for some reason....

# FOR PAPER DID SEPARATE ANOVA AND LM.




# get rid of extra useless output rows from eggnog # !grepl("RTB",TrackingPixel)

proteins.init <- data %>% 
  dplyr::filter(str_detect(query, "TX_SGMG")) %>%  # got almost 50% more annotations 
  dplyr::filter(str_detect(query, "mnt", negate = TRUE ))


# get proper soil sample name by grepping from query name
proteins.init$sample <- gsub(pattern = "_S[0-9].*", replacement = "", proteins.init$query)

dim(proteins.init %>% filter(str_detect(eggNOG_OGs, "Bacteria"))) # 284409 
dim(proteins.init %>% filter(str_detect(eggNOG_OGs, "Archaea"))) #61142
dim(proteins.init %>% filter(str_detect(eggNOG_OGs, "Fungi"))) # 3533
284409 + 61142 + 3533 # so some still unaccounted for....


# get rid of fungi
proteins.init.nofun <- proteins.init %>% 
  filter(str_detect(eggNOG_OGs, "Bacteria") | str_detect(eggNOG_OGs, "Archaea"))
dim(proteins.init.nofun) # around 4000 had fungi in it....

dim(proteins.init %>% filter(str_detect(eggNOG_OGs, "Bacteria")))


# merge with metadata
proteins.merged.init <- merge( proteins.init.nofun, trt , by.x = "sample", by.y = "Sample.Name")
proteins.merged <- merge ( proteins.merged.init, timedata, by =c("Site", "SamplingPt"))
length(unique(proteins.merged$PercMoist))


# standardize predictor variables ....
names (proteins.merged)
proteins.merged$MAP.std <- ( proteins.merged$MAP - mean(proteins.merged$MAP) ) / (2 * sd ( proteins.merged$MAP))
proteins.merged$PercMoist.std <- ( proteins.merged$PercMoist - mean(proteins.merged$PercMoist) ) / (2 * sd ( proteins.merged$PercMoist))
proteins.merged$PercClay.std <- ( proteins.merged$PercClay - mean(proteins.merged$PercClay) ) / (2 * sd ( proteins.merged$PercClay))
proteins.merged$PercC.std <- ( proteins.merged$PercC - mean(proteins.merged$PercC) ) / (2 * sd ( proteins.merged$PercC))
proteins.merged$PercN.std <- ( proteins.merged$PercN - mean(proteins.merged$PercN) ) / (2 * sd ( proteins.merged$PercN))
proteins.merged$pH.std <- ( proteins.merged$pH - mean(proteins.merged$pH) ) / (2 * sd ( proteins.merged$pH))
proteins.merged$MBC.std <- ( proteins.merged$MBC_mg_g - mean(proteins.merged$MBC_mg_g) ) / (2 * sd ( proteins.merged$MBC_mg_g))
proteins.merged$DOC.std <- ( proteins.merged$DOC_mg_g - mean(proteins.merged$DOC_mg_g) ) / (2 * sd ( proteins.merged$DOC_mg_g))

# REMOVE 26!!!
proteins <- proteins.merged %>% 
  filter(sample != "TX_SGMG_26" )
# make KO column that gets rid of the ko: lead-in, makes it possible to match w other stuff
proteins$KO <-  gsub(pattern = "ko:", replacement = "", proteins$KEGG_ko)
proteins$KO

#make every KO its own row !!!!
proteins.uniqueKO <- separate_rows(proteins, KO)
dim (proteins.uniqueKO)

unique(proteins.uniqueKO$KO)
table(nchar(proteins.uniqueKO$KO))
table(nchar(proteins$KO))


names(proteins)


unique(proteins$Preferred_name)

length(unique(proteins$KEGG_Module)) #close to 1100 unique .... 900 with BBmap
length(unique(proteins$KEGG_Pathway)) # >2000
# KEGG pathway N metabolism = ko00910 


unique(proteins$KEGG_Module)

proteins %>% 
  group_by(sample, MAP, SamplingPt, Site, Year) %>% 
  dplyr::count() %>% 
  print(n=57) # saample 42 - wettest  sum 16, still an outlier here

length(grep("Cell wall", proteins$Description))
length(grep("cell wall", proteins$Description))
length(grep("Peptidoglycan", proteins$Description)) #177
length(grep("peptidoglycan", proteins$Description)) # 451
length(grep("ko00550", proteins$KEGG_Pathway)) 
length(grep("03708", proteins$KO))  #transcriptional regulator of stress and heat shock response..
# guess not there?
length(grep("03708", proteins$eggNOG_OGs)) # not here either...

length(grep("01846", proteins$KO))  # sigma factors...
# also I searched for sigma factors in KEGG
# and they're under biofilm production so we already got it in stresss....

length(grep("03110", proteins$BRITE)) # 2810 in BRITE - chaperones.


nrow(proteins[proteins$CAZy!="-",]) 
# 2951 rows have a CAZy designation, pretty few but we want these
nrow(proteins[proteins$CAZy!="-" & proteins$KEGG_Module != "-",])
# 458 have both a cazy and a kegg module designation
nrow(proteins[proteins$CAZy!="-" & proteins$KEGG_Pathway != "-",])
# 2283 have both a cazy and a kegg module designation
nrow(proteins[proteins$CAZy!="-" & proteins$KEGG_ko != "-",])
# ALL of them have a KEGG_ko!!!! so if we select rows with a KO we'd also get these cazy rows
nrow(proteins[proteins$CAZy!="-" & proteins$KO != "-",])
# ALL of them have a KO!!!!!!


nrow(proteins[proteins$PFAMs!="-",]) #300106
nrow(proteins[proteins$KO!="-",])  # 162388 rows have a KO
nrow(proteins[proteins$KEGG_Pathway!="-",])  # only 95982 have a pathway
nrow(proteins[proteins$KEGG_Pathway!="-" & proteins$KO != "-",])
# this means that: everything that has a KEGG pathway, has a kegg KO!!!!!
# so it might be OKAY to use pathways to give designations to KOs,
# but USE KOs as my site x species matrix !!!!!!!!!!!!!!!
# BASED ON THIS : I AM THINKING.......
# I  will do this
# "Pathway" will be used to assign functions to rows
# thaat will be grouped into "species" by KO

# DID NOT END UP USING THESE



names(proteins)
names (microtrait)
microtrait$microtrait_hmm.dbxref_kegg
proteins$KO
dim(proteins) # 1] 345236     49
# Do we want only the KOs that were in the microtrait databses? I do not think so.
#proteins_trait <-  fuzzy_left_join (proteins, microtrait, by = c("KO" = "microtrait_hmm.dbxref_kegg") ,  match_fun = stringr::str_detect)
proteins_trait <-  merge (proteins.uniqueKO, microtrait, by.x = "KO", by.y = "microtrait_hmm.dbxref_kegg" , all.x  =TRUE)
# is part of the problem that different microtrait rules match different parts of "compound KOs"?
# that might make multiple rows for the same KO
# YES I THINK THIS ISS THE PROBLEM

349750  - 345236 # 4000 have multiple matches. what to do ?

dim(proteins_trait) # [1] 349750     67 

names(proteins_trait)
#View(proteins_trait)
table (proteins_trait$KO)
table (proteins_trait$YAS ) # very few for resource use.
# when original join used: 
#Resource Acquisition         Resource Use     Stress Tolerance 
#14882                  150                 2836 
# when fuzzy join used: 
#Resource Acquisition         Resource Use     Stress Tolerance 
#19535                  212                 3532 



proteins.cazy <- proteins_trait %>% 
  filter (CAZy != "-")
nrow(proteins.cazy)
table (proteins.cazy$microtrait_trait.name1)
proteins.cazy %>% 
  filter(!is.na(microtrait_trait.name1)) %>% 
  nrow()





table (proteins_trait$microtrait_trait.name1)

# make stress tolerance, and stress tolerance general, the same category
proteins_trait$microtrait_trait.name1[proteins_trait$microtrait_trait.name1 == "Stress Tolerance:General"] <- "Stress Tolerance"
table (proteins_trait$microtrait_trait.name1)

table (proteins_trait$microtrait_trait.name3)

proteins_trait_DessicStress <- proteins_trait %>% 
  filter( microtrait_trait.name1 == "Stress Tolerance:Specific:desiccation/osmotic/salt stress")  

table(proteins_trait_DessicStress$microtrait_trait.name2) # 122 observations: accumulation of compat ssolutes

table(proteins_trait_DessicStress$microtrait_trait.name3)

proteins.ko_init <- proteins_trait %>% 
  filter(KO != "-" & KO != "")
dim (proteins.ko_init)
# filter - must have at least 10 0ccurrences
proteins.ko <- proteins.ko_init %>% 
  group_by(KO) %>% 
  filter(n() >= 10)
dim ( proteins.ko)


names(proteins.ko)
dim(proteins.ko %>% filter(str_detect(eggNOG_OGs, "Bacteria"))) #
dim(proteins.ko %>% filter(str_detect(eggNOG_OGs, "Archaea"))) 
dim(proteins.ko %>% filter(str_detect(eggNOG_OGs, "Fungi"))) # 

unique(proteins.ko$microtrait_trait.name1)

proteins.ko %>% 
  filter (microtrait_trait.name1 == "Stress Tolerance") %>% 
  distinct(Description) %>% 
  print(n=200)

export_ko_description <- proteins.ko %>% 
  drop_na(microtrait_trait.name1) %>% 
  group_by (microtrait_trait.name1, KO) %>% 
  distinct(Description) %>% 
  print(n=17000)

#write.csv( export_ko_description, "Exported KOs.csv")

length(unique(proteins.ko$KO)) # 3028 KOs
proteins.ko.yas <- proteins.ko %>% filter(!is.na(YAS))
length(unique(proteins.ko.yas$KO)) # only  581 unique proteins in this dataset that are YAS
#444 if filter out rares

# so I think applying function based on a higher KEGG designation will not be a problem here
# since no "overlapping"

#write.csv( proteins.ko, "Protein annotations with microtrait separaterows 6-30-23.csv")

#proteins.ko_onlyKOandfun <- proteins.ko[,c("KO", "fun")]
#listKOs.init <- as.data.frame(unique(proteins.ko$KO))
#names(listKOs.init) <- "KO"
#listKOs <- inner_join(listKOs.init, proteins.ko_onlyKOandfun, by = "KO")

listKOs <- proteins.ko %>% 
  select(c(KO, YAS , microtrait_trait.name1 )) %>% 
  group_by(KO) %>% filter(row_number(microtrait_trait.name1) == 1)
listKOs #

dim(listKOs )

#write.csv( listKOs, "List of YAS KOs.csv")
# this does not include any KOs that do not have a function attached to them@
# but I think this should be okay.
# THIS IS HOW I WILL RETURN FUNCTIONAL CATEGORIES TO THE "OTU" TABLE OR WHATEVER.
# will probably make a dataframe from the OTUs
# then join the fun designation based on the KO.





# so if we want to do NORMALIZATION: 

# # get counts for eack KO  microtrait_trait.name1
data_counts <- proteins.ko %>% 
  group_by(sample, MAP, PercMoist, SamplingPt, Site, Year,   KO , YAS, microtrait_trait.name1,
           MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std, MBC.std, DOC.std,
           PercC, PercN, pH, MBC_mg_g, DOC_mg_g) %>% 
  dplyr::count()
data_counts <- data_counts %>% 
  dplyr::rename("counts" = "n")

names (data_counts)

total_counts <- proteins.ko %>% 
  group_by(sample, MAP, SamplingPt, PercMoist, Site, Year) %>% 
  dplyr::count()

total_counts <- total_counts %>% 
  dplyr::rename("totalsample" = "n")
total_counts # 42 has 10,000 - this is the wettest site in sum 2016


data_counts_withtotal <- merge(data_counts, total_counts, by = c("sample", "MAP", "PercMoist", "Site", "SamplingPt", "Year") )
data_counts_withtotal$relabun<- data_counts_withtotal$counts /  data_counts_withtotal$totalsample * 100
data_counts_withtotal$relabun
# this will be used to compare with TMM  below
str(data_counts_withtotal)
# 69142 when not doing microtrait. 
# this means something is awry, I should not really have more here!

data_counts_withtotal %>% 
  group_by(sample,KO) %>% filter(n() > 1) %>%
  ungroup 

#write.csv(data_counts_withtotal, file = "Texas gene abundances.csv")

# SAME THING BUT WITH CAZY
cazy_counts <- proteins.ko %>% 
  group_by(sample, MAP, PercMoist, SamplingPt, Site, Year,   CAZy ,
           MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std) %>% 
  dplyr::count()
cazy_counts <- cazy_counts %>% 
  dplyr::rename("counts" = "n")




# what is the distribution of total genes in each sample pre norm.

ggplot(data_counts, aes(x = MAP, y= counts, group = Site)) +
  stat_summary(fun.y = "sum", geom = "point") +
  ylab( "Total Raw Counts") + 
  facet_grid(~SamplingPt) # outlier plot 60 

ggplot(data_counts, aes(x = MAP, y= counts, group = Site)) +
  stat_summary(fun.y = "sum", geom = "point")

data_counts_withtotal[!is.na(data_counts_withtotal$YAS),]

data_counts_withtotal_sel <-  data_counts_withtotal %>% 
  select(c(sample, MAP, PercMoist, MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std, Site, SamplingPt, KO, relabun ))
str(data_counts_withtotal_sel$KO)
str(data_counts_withtotal_sel$relabun)
sitexsp_wmeta <- data_counts_withtotal_sel  %>%  pivot_wider(names_from = KO, values_from = relabun ,values_fill = 0)

head(sitexsp_wmeta) 6
dim(sitexsp_wmeta) 

data_counts_withtotal_sel %>%
  dplyr::group_by(sample, MAP, PercMoist, MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std,
                  Site, SamplingPt, KO) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)  

#sitexsp <- sitexsp_wmeta[,6:9092]
sitexsp <- sitexsp_wmeta[,12:3039]  #2997 if not unique KOs separaate rows
coldata <- sitexsp_wmeta[,1:11]# this is basically metadata
rownames(coldata) <- sitexsp_wmeta$sample
dim(coldata)
names(sitexsp)
names(coldata)

