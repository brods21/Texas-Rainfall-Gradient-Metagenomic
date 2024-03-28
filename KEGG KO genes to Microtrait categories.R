# This code uses the supplemental tables from the MicroTrait tool in Karaoz et al. 2022
    # to match up KEGG KO assignments to functional categories. 
# this is modeled after code from Jennifer Jones. 

rm(list=ls())
library(tidyverse)
setwd("~/Downloads")

#Combining microtrait rules to traits databases ----
#Supplementary table 5 
#This file links the microtrait rules and the traits. The rules have multiple gene names (hmm name). 
#This file also links some microtrait rules to substrates. The link between substrates and traits are in the next file. 

# ST5 matches the Microtrait rules (if this gene) to trait assignment
ST5.microtrait<-read.csv("microtrait_st5_rulestotraits.csv", na.strings=c("","NA"))
names(ST5.microtrait) #microtrait_rule.substrate
sum(is.na(ST5.microtrait$microtrait_trait.name1))

#Supplementary Table 6
#ST6 has links substrates to traits. 
#Note: there are multiple substrates that correspond to rules and traits in ST5, so the number of rows in ST6 don't line up with the number of rows with missing data in ST5. 
ST6.microtrait<-read.csv("microtrait_st6_substratetotrait.csv")
names(ST6.microtrait) #microtrait_substrate.name 
ST6.microtrait <- ST6.microtrait %>%
  select(-microtrait_substrate.subclass1,-microtrait_substrate.subclass2,-microtrait_substrate.subclass3)

# Merging ST5 and ST6. 
  # this links up microtrait riules, microtrait substrate, microtrait trait
ST5.microtrait.merged<-merge(ST5.microtrait, ST6.microtrait, by.x= "microtrait_rule.substrate", by.y="microtrait_substrate.name", all.x=T)
names(ST5.microtrait.merged)

#Cleaning up the data frame
ST5.microtrait.merged<-ST5.microtrait.merged %>%
  mutate(microtrait_trait.name1=ifelse(is.na(microtrait_trait.name1.x),microtrait_trait.name1.y,microtrait_trait.name1.x),
         microtrait_trait.name2=ifelse(is.na(microtrait_trait.name2.x),microtrait_trait.name2.y,microtrait_trait.name2.x),
         microtrait_trait.name3=ifelse(is.na(microtrait_trait.name3.x),microtrait_trait.name3.y,microtrait_trait.name3.x)) %>%
  select(-microtrait_trait.name1.x,-microtrait_trait.name2.x,
         -microtrait_trait.name3.x,-microtrait_trait.name1.y,-microtrait_trait.name2.y,-microtrait_trait.name3.y)
names(ST5.microtrait.merged)
sum(is.na(ST5.microtrait.merged$microtrait_trait.name1))
#Note: some of the substrates in ST5 have more than one substrate listed, and these don't merge. 
  #There are also some substrates that aren't listed in the substrate table (ST6), so we have to go through the table manually. We lose 215 rules because of it. 

#write.csv(ST5.microtrait.merged,file="microtrait_rules-to-traits_full.csv")


# Here, stopped and made changes in Excel !!!!!
#Changes in excel: 
#For all the substrates that included a substrate that was used for a resource acquisition trait, I added the Resource Acquisition trait.
#For osmolytes, I added Stress Tolerance trait. 
#most of the unclassified substrates that are left are elements (that aren't obvious nutrients like N, P, and K)
# did not do any outside research on these substrates so far


#Combining microtrait gene names/ko # and traits databases ----
rules.microtrait<-read.csv("microtrait_rules-to-traits_full 2023-05-16.csv", na.strings=c("","NA"))
names(rules.microtrait)

#Supplementary Table 2
#I edited the microtrait supplementary table 2, by adding ko numbers that were associated with GH groups (CAZy database) to the table. 
#I only used ko numbers that were associated with individual GH groups. 
#If ko numbers were associated with more than 1 CAZy group, then I did not include them. 
#I also checked the ko numbers with the existing ko numbers in supp. table 2. 
  #If there were any overlapping ko numbers, I left the original ko number in table 2 and didn't add to it. 
# this should be a list of ko numbers and their associated gene name or "hmm name". 
genes.microtrait<-read.csv("ST2.microtrait_hmms wGH 2023-05-16.csv", na.strings=c("","NA"))
names(genes.microtrait)
genes.microtrait$microtrait_hmm.name<-paste("'",genes.microtrait$microtrait_hmm.name,"'",sep="") # You need to add the quotes to the beginning and end of the hmm name because some of the hmm names are parts of other hmm names. Adding the quotes prevents multiple joining when joining with the microtrait rules database.  

# Merging the rules with the gene names 
#names(rules.microtrait)
rules.microtrait <- rules.microtrait %>%
  select(-microtrait_trait.version) %>%
  mutate(microtrait_hmm.name=microtraitrule.rule) 
# The column that you join with needs to be the same name. So I renamed microtraitrule.rule column
library(fuzzyjoin)

#When joining the microtrait database with the ko number, there end up being duplicate ko numbers
rules.microtrait1 <- rules.microtrait %>%
  fuzzy_join(genes.microtrait, match_fun = str_detect, by = "microtrait_hmm.name", mode = "full") 

rules.microtrait1 <- rules.microtrait1 %>%
  distinct(microtrait_hmm.dbxref_kegg, .keep_all = TRUE) %>%
  separate(microtrait_trait.name1,c("YAS","YAS.2","YAS.3"),sep=":", remove=FALSE)
names(rules.microtrait1)#2280
# There are quite a few hmms names that don't correspond to a microtrait trait
sum(is.na(rules.microtrait1$microtrait_hmm.name.x)) #1048 out of 2280 total hmm names


#write.csv(rules.microtrait1,"~/Documents/Texas gradient/rules_microtrait_2023-05-16.csv")
