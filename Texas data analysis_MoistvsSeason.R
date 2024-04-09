
# this is the analysis from the eggnog_KEGG_microtraitKOs_relabun_withcovariates_72623.r

rm(list=ls())

library(mvabund)
library(tidyverse)
library(lme4)
library(emmeans)
library(vegan)
library(ggrepel)
library(gridExtra)
library(lterpalettefinder)
library(visreg)
library(car)
library(dplyr)
library(ggthemes)
library(RRPP) # just updated !!!



setwd("~/Documents/Texas gradient")
data_counts_withtotal <- read.csv( "Texas gene abundances.csv")
#data <-read.csv( "proteins_eggnog_nobbmap.csv")  #this was what i waas originally using 
#sitedata <-read.csv( "Texas site characteristics.csv")
trt <- read.csv("Texas_Metadata.csv")
timedata <- read.csv("TXGrad_RAPID_sitedata.csv")
#microtrait <- read.csv("~/Documents/Texas gradient/rules_microtrait_2023-05-16.csv")
#sporegene <- read.csv("Lietal2022_MolEcol_S2_SporulationGenes.csv")



enzymes <- read.csv("TXRAPID_Enz_forCaitlin.csv")
co2 <- read.csv("TXGrad_AvgHighCO2_2024_01_09.csv")


MAPcolors <- c("#3B99B1FF",  "#35A4ABFF", "#3BA8A7FF", "#48AEA2FF", "#81BB95FF",
               "#A2C194FF","#B3C58FFF", "#CBC988FF",  "#EAC728FF", "#E9B41FFF" ,
               "#E9AB1CFF", "#E79C15FF",  "#E78D08FF" , "#E78200FF", "#E87600FF" , 
               "#EA6800FF" , "#EC5C00FF" , "#F04105FF", "#F5191CFF")



# subest with only cols im interested in for sitexsp matrix
data_counts_withtotal_sel <-  data_counts_withtotal %>% 
  select(c(sample, MAP, PercMoist, MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std, Site, SamplingPt, KO, relabun ))
str(data_counts_withtotal_sel$KO)
str(data_counts_withtotal_sel$relabun)
sitexsp_wmeta <- data_counts_withtotal_sel  %>%  pivot_wider(names_from = KO, values_from = relabun ,values_fill = 0)

head(sitexsp_wmeta) # genes start at col 12
dim(sitexsp_wmeta)  # 3039 total



sitexsp <- sitexsp_wmeta[,12:3039]  #2997 if not unique KOs separaate rows
coldata <- sitexsp_wmeta[,1:11]# this is basically metadata
rownames(coldata) <- sitexsp_wmeta$sample
dim(coldata)
names(sitexsp)
names(coldata)








# I removed the Differential gene abundance code 1-19-24, still in the
    # eggnog_KEGG_microtraitKOs_relabun_withcovariates_72623.r file 




########################################
# environmental models
#######################################

bothmeta <- merge (trt, timedata, by = c ("Site", "SamplingPt"))
bothmeta2 <- bothmeta
  # all 60, but remember that there are two sites that have the same MAP !!!
bothmeta2$MAP.std <- ( bothmeta2$MAP - mean(bothmeta2$MAP) ) / (2 * sd ( bothmeta2$MAP))
levels(bothmeta$SamplingPt)
str(bothmeta$SamplingPt)
# six comparisons.
percC.mod <- lm(PercC~MAP.std*SamplingPt, bothmeta2) # N/S
percN.mod <- lm(PercN~MAP.std*SamplingPt, bothmeta2) # N/S
sand.mod <- lm(PercSand~MAP.std*SamplingPt, bothmeta2) # N/S
mbc.mod <- lm(MBC_mg_g~MAP.std*SamplingPt, bothmeta2) # SamplingPoint
doc.mod <- lm(DOC_mg_g~MAP.std*SamplingPt, bothmeta2) # SamplingPoint
ph.mod <- lm(pH~MAP.std*SamplingPt, bothmeta2) # MAP significant !! .
moist.mod <- lm(PercMoist~MAP.std*SamplingPt, bothmeta2) # SamplingPoint


# this is fine order of terms does not matter
  # use this do not use anova
  # do not change stats tables just adjust alpha to 0.008
  # these are T and P values for table s2
summary(moist.mod) 
summary(percC.mod)
summary(percN.mod)
summary(sand.mod)
summary(mbc.mod)
summary(doc.mod)
summary(ph.mod)

bothmeta2 %>% 
  group_by( MAP, Site) %>% 
  summarize(pH = mean (pH),
            C = mean (PercC),
            N = mean (PercN), 
            Sand = mean(PercSand),
            Clay = mean(PercClay))



ggplot(data_counts_withtotal, aes( x = MAP, y= PercMoist, group = Site)) + 
  geom_point() + 
  facet_grid(~SamplingPt) # relationship between MAP and percmoist

ggplot(data_counts_withtotal, aes( x = SamplingPt, y= MBC_mg_g)) + 
  geom_boxplot() # sampling point and MBC - high in SPRING

ggplot(data_counts_withtotal, aes( x = SamplingPt, y= DOC_mg_g)) + 
  geom_boxplot() # sampling point and DOC - high in SUMMER


#######################################
# Make data frames for stress and resource and both
  # both overall and fun
#######################################

# stress: 
stress <- data_counts_withtotal  %>% filter(YAS == "Stress Tolerance")

stress$label <- str_to_title (sub('.*Specific:', '', stress$microtrait_trait.name1) ) 
table (stress$label )
sum(stress$counts) 
stress_sum <- stress %>%  # sum over all subcaategories 
  group_by(sample, MAP, Site, SamplingPt, PercMoist,
           MAP.std, PercMoist.std, MBC.std, DOC.std, PercN.std, pH.std) %>% 
  summarise(relabun = sum(relabun))

stress_sum_fun <- stress %>%  # sum for each sub-category
  #filter(microtrait_trait.name1 != "Stress Tolerance:Specific:oxygen limitation") %>%  # IDK why I did this !!!
  group_by(sample, MAP, MAP.std, PercMoist, PercMoist.std, MBC.std, DOC.std, SamplingPt, microtrait_trait.name1,label, Site) %>% 
  summarise(relabun = sum(relabun)) # this low cat of ox lim will be removed anyway,,,

# resource 
resource <- data_counts_withtotal  %>% filter(YAS == "Resource Acquisition")
resource$label <- as.character (str_to_title (sub('.*:.*:', '', resource$microtrait_trait.name1) ) )

resource_sum <- resource %>%  # sum over all subcategories 
  group_by(sample, MAP, SamplingPt, Site, PercMoist, 
           MAP.std, PercMoist.std, MBC.std, DOC.std, PercN.std, pH.std) %>% 
  summarise(relabun = sum(relabun))

resource_sum_fun <- resource %>%  # sum for each sub-category
  group_by(sample, MAP,  MAP.std, PercMoist, PercMoist.std, MBC.std, DOC.std, SamplingPt, microtrait_trait.name1,label, Site) %>% 
  summarise(relabun = sum(relabun))


# both together
names(stress_sum_fun)
names(resource_sum_fun)
both_sum_fun <- rbind (stress_sum_fun, resource_sum_fun )
#MOVED TO bELOW to include the zeroes counts !!!!!!!!!!!!!!!!!!
# do the same for supplemental  stress counts ..... 

# make this wide
table(both_sum_fun$label)
57*.9 # in at least 90% of samples 
tt <- table(both_sum_fun$label)
both_sum_fun_onlyabundantlabel <- subset(both_sum_fun, label %in% names(tt[tt > 50]))

length(unique(both_sum_fun_onlyabundantlabel$label) ) # 18 total categories
table(both_sum_fun_onlyabundantlabel$label)

both_sum_fun_wide <- both_sum_fun_onlyabundantlabel %>%  # make wide 
  ungroup %>% 
  select(-c(microtrait_trait.name1)) %>%
  pivot_wider( names_from = label, values_from = relabun , values_fill = 0)
names(both_sum_fun_wide)

categories_matrix <- both_sum_fun_wide %>%  select(`Stress Tolerance`:`Vitamin Transport`) %>% as.matrix()



#######################################
# Stress line graphs 
#######################################

head ( data_counts_withtotal )


summary(lm(relabun ~ pH.std*SamplingPt, stress_sum)) # not sig wih pH
summary(lm(relabun ~ pH.std, stress_sum)) # not sig with pH.

# MAP fig
stress_plot <- ggplot(stress_sum, aes(x=MAP, y=relabun)) + # higher with higher MAP early on!
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_point() + 
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) + 
  geom_smooth(method = "lm") +
  theme_classic() + 
  ylab("Percent of total genes") + 
  facet_grid(~SamplingPt , scales="free")

stress_plot

# soil moisture fig
stress_plot_moist <- ggplot(stress_sum, aes(x=PercMoist, y=relabun)) + # higher with higher MAP early on!
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_point() + 
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) + 
  geom_smooth(data = subset (stress_sum, SamplingPt=="Fall 2015"), method = "lm", color = "grey20") +
  #ggthemes::theme_clean() + 
  theme_classic(base_size = 13) +  # changed from theme_clean on 1-10-24
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent stress tolerance genes") + 
  xlab("Soil moisture (%)") + 
  facet_grid(~SamplingPt )

stress_plot_moist # 

#ggsave("Stress_microtrait_moist_1102024.png", plot = stress_plot_moist, device = "png",
#   width = 8, height = 3.2, dpi = 300)


hist(stress_sum$relabun) # fine to treat as normal methinks

stress_sum %>% filter (relabun > 3) # BREECO MAP 850 might be driving this. 



stress.model.nofun.season <- lm(relabun~MAP.std* SamplingPt, stress_sum)
summary(stress.model.nofun.season) # season matters, overall R2 0.1327
emmeans(stress.model.nofun.season, pairwise ~SamplingPt)
# fall marg sig from spring, sig from summer. 
ggplot(stress.model.nofun.season, aes(x = SamplingPt, y = relabun)) + 
  geom_boxplot() # BREECO MAP 850 might be driving this. in fall


stresss_boxplot_season <- ggplot(stress.model.nofun.season, aes(x = SamplingPt, y = relabun)) + 
  geom_boxplot()  + 
  ylab("Percent stress tolerance genes") +  
  theme_classic()

stresss_boxplot_season 

#ggsave("Stress_boxplot season.png", plot = stresss_boxplot_season, device = "png",
#   width = 3.5, height = 3.5, dpi = 300)


# follow up if sampling point sig: 
stress.model.nofun.moist <- lm(relabun~ MAP.std*PercMoist.std, stress_sum)
summary(stress.model.nofun.moist) # 
# marg sig interaction between MAP*moist.
# adjusted R2 = 0.05585


# driven by only one site ????

stress_sum_noout <- stress_sum %>%  filter (relabun < 3)
stress.model.nofun.season.noout <- lm(relabun~MAP.std* SamplingPt, stress_sum_noout)
summary(stress.model.nofun.season.noout) # marg
emmeans(stress.model.nofun.season.noout, pairwise ~SamplingPt)
# hmmm but nothing pairwisse .... 
ggplot(stress.model.nofun.season.noout, aes(x = SamplingPt, y = relabun)) + 
  geom_boxplot() # this looks like nothing . 

stress.model.nofun.moist.noout <- lm(relabun~MAP.std* PercMoist.std , stress_sum_noout)
summary(stress.model.nofun.moist.noout) # # nothing honey
  # that means marg interaction driven by that one weird point. 


interactions::interact_plot(model = stress.model.nofun.moist, 
                            modx = MAP.std, pred = PercMoist.std,
                            vary.lty = FALSE ,
                            colors = c("lightblue", "blue", "darkblue"))
# suggests positive relationship @ high MAP, neg at low MAP?

# this is for percmoist * map marginal interaction 
ggplot(stress_sum, aes (x = PercMoist , y= relabun, color = as.character(MAP)))  + 
  geom_line() +
  geom_point() +
  scale_color_manual(values = rev(MAPcolors)) + 
  ylab("Rel. abun of resource genes") + 
  xlab("Soil Moisture (%)") + 
  theme_classic(base_size = 20) 



#######################################
# Stress line graphs by function 
#######################################

#### MANYLM in mvabund #############

table(stress_sum_fun$label)
57*.9
tt <- table(stress_sum_fun$label)
stress_sum_fun_onlyabundantlabel <- subset(stress_sum_fun, label %in% names(tt[tt > 50]))
  # subset so haave to be in 90% of samples!!!!!!!! 
table(stress_sum_fun_onlyabundantlabel$label)

stress_sum_fun_wide <- stress_sum_fun_onlyabundantlabel %>% 
  ungroup %>% 
  select(-c(microtrait_trait.name1)) %>%
  pivot_wider( names_from = label, values_from = relabun , values_fill = 0)
names(stress_sum_fun_wide)
head(stress_sum_fun_wide)



manylm_stress_season <- manylm(as.matrix(stress_sum_fun_wide[,10:14]) ~ 
                          MAP.std*SamplingPt, data = stress_sum_fun_wide)
manylm_stress_season
#plot(manylm_stress)
summary.manylm(manylm_stress_season, nBoot = 9999 , rep.seed = T) # this says nothing matters

summary.manylm(manylm_stress_season, p.uni="adjusted", nBoot = 9999, rep.seed = T)
# Hooper's R-squared: 0.05793
# nothing matters, not even close.



manylm_stress_moist <- manylm(as.matrix(stress_sum_fun_wide[,10:14]) ~ 
                                 MAP.std*PercMoist.std, data = stress_sum_fun_wide)
manylm_stress_moist
#plot(manylm_stress)
summary.manylm(manylm_stress_moist, nBoot = 9999 , rep.seed = T) # nothing 
summary.manylm(manylm_stress_moist, p.uni="adjusted", nBoot = 9999, rep.seed = T)
# Hooper's R-squared: 0.06589 









# make sure the line graphs have zeroes!!!

stress_sum_fun_backtolong <- stress_sum_fun_wide %>% 
  pivot_longer( cols = "Stress Tolerance":"Temperature:low" ,
                names_to = "label",
                values_to = "relabun")


# Thiss iss where color = moisture?
stress_fun_plot_moistcol <- ggplot( data= stress_sum_fun_backtolong #, subset(stress_sum_fun,( label =="Ph Stress" ) )
                                    , aes(x=MAP, y=relabun, 
                                          color = PercMoist)) + # higher with higher MAP early on!
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_smooth(method = "lm", color = "grey", linetype = "dashed") +
  geom_point( size = 2) + 
  xlab("MAP") + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  ylab("Relative Abundance ") + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  #theme(legend.position = "none") + 
  #facet_grid(rows =vars(SamplingPt), cols =vars(label ))
  facet_wrap(~label , nrow = 2, scales= "fixed")

stress_fun_plot_moistcol

ggsave("StressFunctions_KEGGmicrotrait_moistcol_3-24.png", plot = stress_fun_plot_moistcol, device = "png",
      width = 8, height = 5, dpi = 300)





#######################################
# Resource line graphs 
#######################################


# do pH as a driver. include in supplemental
summary(lm (relabun~pH.std*SamplingPt, resource_sum)) # nothing honey
summary(lm (relabun~pH.std, resource_sum)) # 0.0591, R2 = 0.046




resource_plot <- ggplot(resource_sum, aes(x=MAP, y=relabun)) + # higher with higher MAP early on!
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_point() + 
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) + 
  geom_smooth(data = subset (resource_sum, SamplingPt=="Spring 2016"), method = "lm", color = "grey20") +
  theme_classic(base_size = 13) +  # changed from theme_clean on 1-10-24
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent resource acquisition genes") + 
  xlab("MAP") + 
  facet_grid(~SamplingPt )

resource_plot # 

#ggsave("Resource_microtrait_MAP_1102024.png", plot = resource_plot , device = "png",
#      width = 8, height = 3.2, dpi = 300)


resource_plot_mapxmoist_a <- ggplot(resource_sum, aes(x=MAP, y=relabun, col = PercMoist)) + 
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_point() + 
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) + 
  #geom_smooth(data = subset (resource_sum, SamplingPt=="Spring 2016"), method = "lm", color = "grey20") +
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  theme_classic(base_size = 13) +  # changed from theme_clean on 1-10-24
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent resource acquisition genes") + 
  xlab("MAP")

#ggsave("Resource_microtrait_MAPxMoist_A.png", plot = resource_plot_mapxmoist_a , device = "png",
#      width = 4.5, height = 4, dpi = 300)


hist(resource_sum$relabun) #this is good

resource.model.nofun.season <- lm(relabun~MAP.std*SamplingPt, resource_sum)
summary(resource.model.nofun.season) # nothinig.
# R2 0.08793

resource.model.nofun.moist <- lm(relabun~MAP.std*PercMoist.std, resource_sum)
summary(resource.model.nofun.moist) # MAP*moist interaction!
# R2 0.2137

interactions::interact_plot(model = resource.model.nofun.moist, 
                            modx = MAP.std, pred = PercMoist.std,
                            vary.lty = FALSE ,
                            colors = c("lightblue", "blue", "darkblue"))
  # suggests relabun high at high moist low MAP, and vice versa

# this is for percmoist * map interaction 
ggplot(resource_sum, aes (x = PercMoist , y= relabun, color = as.character(MAP)))  + 
  geom_line() +
  geom_point() +
  scale_color_manual(values = rev(MAPcolors)) + 
  ylab("Rel. abun of resource genes") + 
  xlab("Soil Moisture (%)") + 
  theme_classic(base_size = 20) 

ggplot(resource_sum, aes (x = PercMoist , y= relabun, color = as.character(MAP)))  + 
  geom_smooth(method = lm, se = F) + # not sure it's appropriate to draw line through 
  geom_point() +
  scale_color_manual(values = rev(MAPcolors)) + 
  ylab("Rel. abun of resource genes") + 
  xlab("Soil Moisture (%)") + 
  theme_classic(base_size = 20) 



# STATS and FIGURE for how MOIST relationship changes with MAP category?


# make df
resource_sum_wetdry <-resource_sum  %>% ungroup(.) %>% mutate(wetdry_MAP = case_when(
  resource_sum$MAP <= 600 ~ "Dry", 
  resource_sum$MAP > 600 & resource_sum$MAP < 800 ~ "Mid", 
  resource_sum$MAP >= 800 ~ "Wet"
))

# stats

resource.model.mapcat <- lm(relabun~wetdry_MAP*PercMoist.std, resource_sum_wetdry)
summary(resource.model.mapcat) # yes this is sig!
joint_tests(resource.model.mapcat, by= "wetdry_MAP") # only sig for dry??? 

# figure - lines for each MAP category
resource_plot_mapxmoist_b<- ggplot(resource_sum_wetdry, aes (x = PercMoist , y= relabun, color = wetdry_MAP, linetype = wetdry_MAP))  + 
  geom_smooth(data = subset (resource_sum_wetdry, wetdry_MAP == "Dry"), 
              method = lm, linetype = "solid") + # line and SE cuz sig
  geom_smooth(data = subset (resource_sum_wetdry, wetdry_MAP != "Dry"), 
              method = lm, se = F ,linetype = "dashed") + # no line and no se,
  geom_point() +
  scale_color_manual(values = c("#F04105FF", "#E9B41FFF",  "#3B99B1FF")) + 
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) + 
  theme_classic(base_size = 12) +
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent resource acquisition genes") + 
  xlab("Soil Moisture (%)") 
  
 

#ggsave("Resource_microtrait_MAPxMoist_C.png", plot = resource_plot_mapxmoist_b , device = "png",
#      width = 5, height = 4, dpi = 300)



# make df

hist(resource_sum$PercMoist)
resource_sum_wetdrycurrent <-resource_sum  %>% ungroup(.) %>% mutate(wetdry_moist = case_when(
  resource_sum$PercMoist <= 10 ~ "Dry", 
  resource_sum$PercMoist > 10 & resource_sum$PercMoist < 20 ~ "Mid", 
  resource_sum$PercMoist >= 20 ~ "Wet"
))


resource.model.moistcat <- lm(relabun~MAP.std*wetdry_moist, resource_sum_wetdrycurrent)
summary(resource.model.moistcat) # MAP, interaction sig.
joint_tests(resource.model.moistcat, by= "wetdry_moist") # only sig for dry??? 
# MAP only matters when it is currently dry, is the other interpretation??????????


# maybe resource figure could be:
# 1: MAP-abundance, colored by percmoist
  # then MAP categories : ie most sensitive at dry site
  # then moist categories: ie legacy matters when dry. these both work!!



# figure - lines for each MAP category
resource_plot_mapxmoist_c <- ggplot(resource_sum_wetdrycurrent, aes (x = MAP , y= relabun, color = wetdry_moist, linetype = wetdry_moist))  + 
  geom_point() +
  geom_smooth(data = subset (resource_sum_wetdrycurrent, wetdry_moist == "Dry"), 
              method = lm, linetype = "solid") + # line and SE cuz sig
  geom_smooth(data = subset (resource_sum_wetdrycurrent, wetdry_moist != "Dry"), 
              method = lm, se = F ,linetype = "dashed") + # no line and no se,
  scale_color_manual(values = c("#F04105FF", "#E9B41FFF",  "#3B99B1FF")) + 
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) + 
  theme_classic(base_size = 12) +
  theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent resource acquisition genes") + 
  xlab("MAP")
 # theme_classic(base_size = 20) 


 # ggsave("Resource_microtrait_MAPxMoist_D.png", plot = resource_plot_mapxmoist_c , device = "png",
 #        width = 5, height = 4, dpi = 300)




#######################################
# Resource line graphs by function 
#######################################


#MOVED TO bELOW to include the zeroes counts !!!!!!!!!!!!!!!!!!
# do the same for supplemental  stress counts ..... 


#### MANYLM in mvabund #############

length(unique(resource_sum_fun$label) )
table(resource_sum_fun$label)
57*.9
tt <- table(resource_sum_fun$label)
resource_sum_fun_onlyabundantlabel <- subset(resource_sum_fun, label %in% names(tt[tt > 50]))

length(unique(resource_sum_fun_onlyabundantlabel$label) )
table(resource_sum_fun_onlyabundantlabel$label)

resource_sum_fun_wide <- resource_sum_fun_onlyabundantlabel %>% 
  ungroup %>% 
  select(-c(microtrait_trait.name1)) %>%
  pivot_wider( names_from = label, values_from = relabun , values_fill = 0)
names(resource_sum_fun_wide)
names(resource_sum_fun_wide[,10:22])

# maybe do not do it this way, since the overall model nothing was sig here. 
manylm_resource_season <- manylm(as.matrix(resource_sum_fun_wide[,10:22]) ~ 
                            MAP.std*SamplingPt, data = resource_sum_fun_wide)
summary(manylm_resource_season, nBoot = 9999 , rep.seed = T)
summary.manylm(manylm_resource_season, p.uni="adjusted", nBoot = 9999 , rep.seed = T) 
# Vitamin transport, Complex Carbohydrate Depolymerizatio - nope ......
# Marg Carboxylate Transport MAP*season



# these are the results I will show. 
manylm_resource_moist <- manylm(as.matrix(resource_sum_fun_wide[,10:22]) ~ 
                                   MAP.std*PercMoist, data = resource_sum_fun_wide)
summary(manylm_resource_moist, nBoot = 9999 , rep.seed = T)
summary.manylm(manylm_resource_moist, p.uni="adjusted", nBoot = 9999 , rep.seed = T) 
# MAP: Free Amino Acids Transport  , Carbohydrate Transport marg
# same for MAP*moist.

# free amino acid transport - MAPxMoist (no samplingPt) 
# carbohydrate transport - MAPxMoist (no samplingPt) 
# Complex Carbohydrate Depolymerization - MAP

AAmod <- lm(`Free Amino Acids Transport`~ MAP.std*PercMoist.std, resource_sum_fun_wide)
summary(AAmod)
library(interactions)
interactions::interact_plot(model = AAmod,  # main way
                            pred = MAP.std, modx = PercMoist.std) 
interactions::interact_plot(model = AAmod,  #other way
                            modx = MAP.std, pred = PercMoist.std,
                            vary.lty = FALSE ,
                            colors = c("lightblue", "blue", "darkblue"))
visreg2d(AAmod, "MAP.std", "PercMoist.std")

# third way - Christine suggestion in 


AA_mapxmoist_plot <- ggplot(resource_sum_fun_wide, aes (x = PercMoist , y= `Free Amino Acids Transport`, color = as.character(MAP)))  + 
  geom_line() +
  scale_color_manual(values = rev(MAPcolors)) + 
  xlab("Soil Moisture (%)") + 
  theme_classic(base_size = 20)
AA_mapxmoist_plot

#ggsave("AATransport_interaction_1-9-2024.png", plot = AA_mapxmoist_plot , device = "png",
#    width = 7, height = 5, dpi = 300)


Ctransmod <- lm(`Carbohydrate Transport`~ MAP.std*PercMoist.std, resource_sum_fun_wide)
summary(Ctransmod)
library(interactions)
interactions::interact_plot(model = Ctransmod, 
                            pred = MAP.std, modx = PercMoist.std)
interactions::interact_plot(model = Ctransmod, 
                            modx = MAP.std, pred = PercMoist.std,
                            vary.lty = FALSE ,
                            colors = c("lightblue", "blue", "darkblue"))
visreg2d(Ctransmod, "MAP.std", "PercMoist.std")

CarbTrans_mapxmoist_plot <- ggplot(resource_sum_fun_wide, aes (x = PercMoist , y= `Carbohydrate Transport`, color = as.character(MAP)))  + 
  geom_line() +
  scale_color_manual(values = rev(MAPcolors)) + 
  xlab("Soil Moisture (%)") + 
  theme_classic(base_size = 20)
CarbTrans_mapxmoist_plot 

#ggsave("CarbTransport_interaction_1-9-2024.png", plot = CarbTrans_mapxmoist_plot , device = "png",
#       width = 7, height = 5, dpi = 300)




names ( resource_sum_fun_wide)
resource_sum_fun_backtolong <- resource_sum_fun_wide %>% 
  pivot_longer( cols = "Resource Acquisition":"Vitamin Transport" ,
                names_to = "label",
                values_to = "relabun")


resource_fun_plot_moistcol_haszeroes_sel <- ggplot( data= subset(resource_sum_fun_backtolong, 
                                                             
                                                               ( label =="Carbohydrate Transport" ) |
                                                               
                                                               ( label =="Free Amino Acids Transport" ) )
                                                , aes(x=MAP, y=relabun, 
                                                      colour = PercMoist)) + # higher with higher MAP early on!
  #stat_summary(fun.y = "mean", geom = "point") +
  geom_point( size = 2) + 
  xlab("MAP") + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  ylab("Relative Abundance ") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  facet_wrap(~factor(label))
  #theme(legend.position = "none") + 
  #facet_wrap( ~factor(label,
   #                   levels = c(
   #                              "Carbohydrate Transport",
                                 
   #                              "Free Amino Acids Transport")), scales = "free_y" ) # not - nodif between sampling points!

resource_fun_plot_moistcol_haszeroes_se;

#ggsave("ResourceFunctions_KEGGmicrotrait_moistcol_11024_all.png", plot = resource_fun_plot_moistcol_haszeroes_sel, device = "png",
#      width = 6, height = 4.8, dpi = 300)


resource_fun_plot_moistcol_haszeroes <- ggplot( data= resource_sum_fun_backtolong,
                                                aes(x=MAP, y=relabun )) + # higher with higher MAP early on!
  geom_point( size = 2, aes ( colour = PercMoist)) + 
  geom_smooth(data = subset (resource_sum_fun_backtolong, label =="Free Amino Acids Transport"), 
              method = lm, linetype = "solid", color = "black") + # line and SE cuz sig
  geom_smooth(data = subset (resource_sum_fun_backtolong, label =="Carbohydrate Transport"), 
              method = lm, linetype = "solid", color = "black") + # line and SE cuz sig
  geom_smooth(data = subset (resource_sum_fun_backtolong, label !="Carbohydrate Transport" & 
                               label  != "Free Amino Acids Transport"), 
              method = lm, linetype = "dashed", color = "grey") + # line and SE cuz sig
  xlab("MAP") + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  ylab("Relative Abundance ") + 
  theme_classic(base_size = 9) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  facet_wrap(~factor(label), scales = "free_y")

resource_fun_plot_moistcol_haszeroes

ggsave("ResourceFunctions_KEGGmicrotrait_moistcol_3_24_all.png", plot = resource_fun_plot_moistcol_haszeroes, device = "png",
      width = 8, height = 7, dpi = 300)




#######################################
# look at individual KOs in free AA transport and carbohydrate transport ... ?
#######################################

names(data_counts_withtotal)
names(resource)

AA_carb <- resource %>% 
  filter (label == "Free Amino Acids Transport" |
            label == "Carbohydrate Transport")
table (AA_carb$label)
table (AA_carb$KO)
names (AA_carb)

AA_carb_wide <- AA_carb %>% 
  select(-c(X,label, microtrait_trait.name1, YAS,counts, totalsample)) %>% 
  pivot_wider (names_from = KO, values_from = relabun, values_fill = 0)

manylm.aa.carb <- manylm (as.matrix(AA_carb_wide[,20:89]) ~MAP.std*PercMoist.std, data = AA_carb_wide )
summary( manylm.aa.carb) # map matters. interaction matters 
summary.manylm( manylm.aa.carb, p.uni="adjusted", nBoot = 9999 ) 
  # K16210, this is not interesting 

##############################################
# Resource AND STRESS line graphs by function 
# if we were to put them in the same model. 
############################################

# both_sum_fun


#### MANYLM in mvabund #############

# what do I need:
# both_sum_fun_wide has categories and relevant metadata
# categories_matrix is subset matrix of just the categories from this df
  # (eg site x species matrix but gene cats.)

manylm_both_season <- manylm(categories_matrix ~ 
                        MAP.std*SamplingPt, data = both_sum_fun_wide)

summary.manylm(manylm_both_season)  # Hooper R2 = 0.1001 , marg sig MAP
summary.manylm(manylm_both_season, p.uni = "adjusted", nBoot = 9999, rep.seed = T) # nothing is sig this way,

manylm_both_moist <- manylm(categories_matrix ~ 
                        MAP.std*PercMoist.std, data = both_sum_fun_wide)

summary.manylm(manylm_both_moist)  # Hooper R2 = 0.08898 , MAP sig
summary.manylm(manylm_both_moist, p.uni = "adjusted", nBoot = 9999, rep.seed = T) 
  # MAP*Moist : Free AA trans, carb marg, same as when separate models



########################################################################
# stress and resource categories as predictors of CO2
########################################################################
names(both_sum_fun_wide)
dim(both_sum_fun_wide)
both_funs_co2  <- merge (both_sum_fun_wide, co2, by = c( "Site", "SamplingPt"))
dim(both_funs_co2)
names(both_funs_co2)
dim(categories_matrix)


# first off: does it vary by gradient
co2_gradient_season <- lm(CO2_AVG~MAP.std  *SamplingPt, data = both_funs_co2)
summary(co2_gradient_season)  # yes varies by gradient and sampling point, no interaction. 

co2_gradient_season <- lm(CO2_AVG~MAP.std  *SamplingPt, data = both_funs_co2)
summary(co2_gradient_season)  # yes varies by gradient and sampling point, no interaction. 

co2_gradient_moist <- lm(CO2_AVG~MAP.std  *PercMoist.std, data = both_funs_co2)
summary(co2_gradient_moist)  # yes varies by gradient but not sampling point, no interaction. 

# now linking genes and function!

# make gene cats matrix - same as before but making sure in same order ....
genecats_co2 <- both_funs_co2  %>% select (`Stress Tolerance`:`Vitamin Transport`) %>% as.matrix()
dim (genecats_co2) # 55


co2predict_MAP <- lm.rrpp(CO2_AVG~MAP.std +SamplingPt + PercMoist.std,  # R2  = 0.16
                          data = both_funs_co2, print.progress = FALSE, SS.type = "III",
                          turbo = FALSE, verbose = TRUE, iter = 9999)
co2predict_genecats <- lm.rrpp(CO2_AVG~genecats_co2   , 
                               data = both_funs_co2, print.progress = FALSE, SS.type = "III",
                               turbo = FALSE, verbose = TRUE, iter = 9999) # R2 = 0.47

co2predict_genecats_MAP <- lm.rrpp(CO2_AVG~genecats_co2 +MAP.std + SamplingPt+ PercMoist.std , 
                                   data = both_funs_co2, print.progress = FALSE, SS.type = "III",
                                   turbo = FALSE, verbose = TRUE, iter = 9999) #R2 = 0.57


summary(co2predict_MAP) # p<0.001, R2 = 0.359 ( 40 w moist)
summary(co2predict_genecats) # p=0.042, R2 = 0.525 . p value is significant.. NOT when no ssamplingpt
summary(co2predict_genecats_MAP) # p=0.005, R2 = 0.626 (6335 with moist)


anova(co2predict_MAP) #MAP and samplingpt sig (percmoist marg)
anova(co2predict_genecats) # this IS NOT SIGNIFICANT when using type III :(
anova(co2predict_genecats_MAP) # matrix not sig, but explains lots of variation... sad

model.comparison(co2predict_MAP, co2predict_genecats, co2predict_genecats_MAP, type ="logLik")
# aic suggests, here, that the MAP only model is best. 

model.comparison(co2predict_MAP,  co2predict_genecats_MAP, type ="logLik")
# aic suggests, here, that the MAP only model is best. 


anova(co2predict_MAP,  co2predict_genecats_MAP)
# this suggests not valuable to include gene cats, possibly??
coef(co2predict_genecats, test = TRUE) #

coef(co2predict_genecats_MAP, test = TRUE) #
# stress tolerance, carboxylate transport, 
# marg - vitamin, organophos trans, ion trans,


ggplot (both_funs_co2, aes (x = `Stress Tolerance`, y = CO2_AVG)) + 
  geom_smooth(method = lm) + 
  geom_point() 

ggplot (both_funs_co2, aes (x = `Carboxylate Transport`, y = CO2_AVG)) + 
  geom_smooth(method = lm) +
  geom_point() 


########################################################################
# stress and resource categories as predictors of ENZYMES!!!!!!
########################################################################
names(both_sum_fun_wide)
dim(both_sum_fun_wide)
both_funs_enz  <- merge (both_sum_fun_wide, enzymes, by = c( "Site", "SamplingPt"))
dim(both_funs_enz)
names(both_funs_enz)
 # has one fewer than CO2 matrix becayse 


# remove the row that had outliers for several enzymes
both_funs_enz_noout <- both_funs_enz %>%  #enzymes_meta %>% 
  filter(CBH < 2000)

# make gene cats matrix 
genecats_enz <- both_funs_enz_noout %>% select (`Stress Tolerance`:`Vitamin Transport`) %>% as.matrix()
dim (genecats_enz) # 55


hist(both_funs_enz_noout$CBH)  # transform?
hist(both_funs_enz_noout$AG) # transform?
hist(both_funs_enz_noout$BG) # transform?
hist(both_funs_enz_noout$BX) # transform?

hist(log(both_funs_enz_noout$CBH))  # better
hist(log(both_funs_enz_noout$AG) )# better
hist(log(both_funs_enz_noout$BG) )# better
hist(log(both_funs_enz_noout$BX)) # better

names(both_funs_enz_noout)

both_funs_enz_noout$CBH_log <- log(both_funs_enz_noout$CBH) # always log enzyme assays
both_funs_enz_noout$BG_log <- log(both_funs_enz_noout$BG)
both_funs_enz_noout$AG_log <- log(both_funs_enz_noout$AG)
both_funs_enz_noout$NAG_log <- log(both_funs_enz_noout$NAG)
both_funs_enz_noout$BX_log <- log(both_funs_enz_noout$BX)
both_funs_enz_noout$AP_log <- log(both_funs_enz_noout$AP)

both_funs_enz_noout$total = log(
  both_funs_enz_noout$CBH +
    both_funs_enz_noout$BG+
    both_funs_enz_noout$AG +
    both_funs_enz_noout$NAG +
    both_funs_enz_noout$BX  + # calculate total enzyme abundance - 
    both_funs_enz_noout$AP  # measure of microbial investment in resource acquisiton
)




names(both_funs_enz_noout)


# first off: does it vary across the gradient.
enz_gradient_season <- lm(total~MAP.std +SamplingPt, both_funs_enz_noout)
summary(enz_gradient_season) # marg sig effect of MAP, highly sig if take out interaction. 

enz_gradient_moist <- lm(total~MAP.std *PercMoist.std, both_funs_enz_noout)
summary(enz_gradient_moist) # MAP and Moist both sig, no interaction. 


# with total enzyme abundance as response variable. 

enz_predict_MAP <- lm.rrpp(total~MAP.std +  SamplingPt + PercMoist.std ,  # R2  = 0.23
                           data = both_funs_enz_noout, print.progress = FALSE, SS.type = "III", # just MAP
                           turbo = FALSE, verbose = TRUE, iter = 9999)
enz_predict_genecats <- lm.rrpp(total~genecats_enz ,  # enzymes
                                data = both_funs_enz_noout, print.progress = FALSE, SS.type = "III",
                                turbo = FALSE, verbose = TRUE, iter = 9999) # R2 = 0.60
enz_predict_genecats_MAP <- lm.rrpp(total~genecats_enz + MAP.std + SamplingPt + PercMoist.std, 
                                    data = both_funs_enz_noout, print.progress = FALSE, SS.type = "III", # enz + MAP
                                    turbo = FALSE, verbose = TRUE, iter = 9999) #R2 = 0.57
# what does it all mean

summary(enz_predict_MAP) # p=0.004, R2 = 0.235
summary(enz_predict_genecats) # p=0.006, R2 = 0.604 #still sig (p=0.009) if remove SamplingPt
summary(enz_predict_genecats_MAP) # p=0.004, R2 = 0.638
# this suggests including genes good

anova(enz_predict_MAP, type = "III") # map significant
anova(enz_predict_genecats, type = "III") # genes significant , even without 
anova(enz_predict_genecats_MAP, type = "III") # matrix significant , moreso than MAP

model.comparison( enz_predict_MAP,enz_predict_genecats, enz_predict_genecats_MAP, type = c( "logLik"))
# what tf does this mean.
# aic - enzymes + MAP best!!!!!!
  # if remove sampling point this is no longer true....
model.comparison(enz_predict_MAP,  enz_predict_genecats_MAP, type = "logLik")
  # when comparing only these two models enzymes help!
anova(enz_predict_MAP,  enz_predict_genecats_MAP) # p = 0.036 compaared to MAP only model.
coef(enz_predict_genecats_MAP, test = TRUE) #
# hm this iss confusing cuz none appear to be sig......
  #. but very large effect sizes!!!!!

names(both_funs_enz_noout)
ggplot (both_funs_enz_noout, aes (x = `Stress Tolerance`, y = total)) +  # marg sig, large effect size
  geom_point() 

ggplot (both_funs_enz_noout, aes (x = `Desiccation/Osmotic/Salt Stress`, y = total)) +  # huge effect size
  geom_point() # but this iss nothinig lmao

ggplot (both_funs_enz_noout, aes (x = `N Compound Transport`, y = total)) +  # huge effect size
  geom_point() 

ggplot (both_funs_enz_noout, aes (x = `Carboxylate Transport`, y = total)) +  # huge effect size
  geom_point() 

ggplot (both_funs_enz_noout, aes (x = MAP, y = total)) + 
  geom_point() + 
  facet_grid (~SamplingPt)


# with individual enzymes abundance as response variable. 
enz_all_predict_genecats <- lm.rrpp(as.matrix(both_funs_enz_noout[,34:39])
                                    ~as.matrix(both_funs_enz_noout[,8:25])+SamplingPt , 
                                    data = both_funs_enz_noout, print.progress = FALSE,
                                    turbo = FALSE, verbose = TRUE) # R2 = 0.60

enz_all_predict_MAP <- lm.rrpp(as.matrix(both_funs_enz_noout[,34:39])
                               ~MAP.std +SamplingPt,  # R2  = 0.16
                               data = both_funs_enz_noout, print.progress = FALSE,
                               turbo = FALSE, verbose = TRUE)
enz_all_predict_genecats_MAP <- lm.rrpp(as.matrix(both_funs_enz_noout[,34:39])
                                        ~as.matrix(both_funs_enz_noout[,8:25])+MAP.std +SamplingPt, 
                                        data = both_funs_enz_noout, print.progress = FALSE,
                                        turbo = FALSE, verbose = TRUE) #R2 = 0.57
# what does it all mean
summary(enz_all_predict_MAP) # p=0.001, R2 = 0.23
summary(enz_all_predict_genecats) # p=0.002, R2 = 0.58
summary(enz_all_predict_genecats_MAP) # p=0.001, R2 = 0.62
# this suggests including genes good


anova(enz_all_predict_genecats) # 0.005 when SamplingPt included. how do we get which ones matter?

anova(enz_all_predict_MAP,  enz_all_predict_genecats_MAP) # p = 0.031 compaared to MAP only model.

model.comparison(enz_all_predict_MAP, enz_all_predict_genecats, enz_all_predict_genecats_MAP, type = "logLik")

model.comparison(enz_all_predict_MAP, enz_all_predict_genecats_MAP, type = "logLik")
coef(enz_all_predict_genecats_MAP, test = TRUE) #
# hm this iss confusing cuz none appear to be sig......







# what figure do I want to show with this???????



########################################################################
# stress and resource categories as predictors of MASS_SPECIFIC ENZYMES!!!!!!
########################################################################


names(timedata)

both_funs_enz_noout_MBC <- merge (both_funs_enz_noout, timedata, by = c( "Site", "SamplingPt"))

head(both_funs_enz_noout_MBC)


both_funs_enz_noout_MBC$CBH_mbc <- (both_funs_enz_noout_MBC$CBH) / both_funs_enz_noout_MBC$MBC_mg_g
both_funs_enz_noout_MBC$BG_mbc <- (both_funs_enz_noout_MBC$BG)/ both_funs_enz_noout_MBC$MBC_mg_g
both_funs_enz_noout_MBC$AG_mbc <- (both_funs_enz_noout_MBC$AG)/ both_funs_enz_noout_MBC$MBC_mg_g
both_funs_enz_noout_MBC$NAG_mbc <- (both_funs_enz_noout_MBC$NAG)/ both_funs_enz_noout_MBC$MBC_mg_g
both_funs_enz_noout_MBC$BX_mbc <- (both_funs_enz_noout_MBC$BX)/ both_funs_enz_noout_MBC$MBC_mg_g
both_funs_enz_noout_MBC$AP_mbc <- (both_funs_enz_noout_MBC$AP)/ both_funs_enz_noout_MBC$MBC_mg_g


both_funs_enz_noout_MBC$CBH_log <- log(both_funs_enz_noout_MBC$CBH_mbc) # always log enzyme assays
both_funs_enz_noout_MBC$BG_log <- log(both_funs_enz_noout_MBC$BG_mbc)
both_funs_enz_noout_MBC$AG_log <- log(both_funs_enz_noout_MBC$AG_mbc)
both_funs_enz_noout_MBC$NAG_log <- log(both_funs_enz_noout_MBC$NAG_mbc)
both_funs_enz_noout_MBC$BX_log <- log(both_funs_enz_noout_MBC$BX_mbc)
both_funs_enz_noout_MBC$AP_log <- log(both_funs_enz_noout_MBC$AP_mbc)

both_funs_enz_noout_MBC$total_mbc = log(
  both_funs_enz_noout_MBC$CBH_mbc +
    both_funs_enz_noout_MBC$BG_mbc+
    both_funs_enz_noout_MBC$AG_mbc +
    both_funs_enz_noout_MBC$NAG_mbc +
    both_funs_enz_noout_MBC$BX_mbc  + # calculate total enzyme abundance - 
    both_funs_enz_noout_MBC$AP_mbc  # measure of microbial investment in resource acquisiton
)




names(both_funs_enz_noout_MBC)


summary(lm(total_mbc~MAP.std*PercMoist.std*SamplingPt, both_funs_enz_noout_MBC))
# interactions....


# with total enzyme abundance as response variable. 
enz_predict_genecats_mbc <- lm.rrpp(total_mbc~as.matrix(both_funs_enz_noout_MBC[,8:25]) +SamplingPt,  # enzymes
                                data = both_funs_enz_noout_MBC, print.progress = FALSE,
                                turbo = FALSE, verbose = TRUE, iter = 9999) # R2 = 0.60

enz_predict_MAP_mbc <- lm.rrpp(total_mbc~MAP.std + SamplingPt,  # R2  = 0.16
                           data = both_funs_enz_noout_MBC, print.progress = FALSE, # just MAP
                           turbo = FALSE, verbose = TRUE, iter = 9999)
enz_predict_genecats_MAP_mbc <- lm.rrpp(total_mbc~as.matrix(both_funs_enz_noout_MBC[,8:25]) + MAP.std + SamplingPt, 
                                    data = both_funs_enz_noout_MBC, print.progress = FALSE, # enz + MAP
                                    turbo = FALSE, verbose = TRUE, iter = 9999) #R2 = 0.57
# what does it all mean
summary(enz_predict_MAP_mbc) # p=0.002, R2 = 0.24
summary(enz_predict_genecats_mbc) # p=0.02, R2 = 0.57
summary(enz_predict_genecats_MAP_mbc) # p=0.03, R2 = 0.57
# this suggests including genes good


anova(enz_predict_genecats_MAP_mbc) # only sampling time matterss, which makes ssenses
  # MBC varied with strongly with sampling time so the denominator in our calcs is very different for season

model.comparison(enz_predict_MAP_mbc, enz_predict_genecats_mbc, enz_predict_genecats_MAP_mbc, type = c( "logLik"))
# MAP model besst

anova(enz_predict_MAP_mbc,  enz_predict_genecats_MAP_mbc) # p = 0.2 compaared to MAP only model.
coef(enz_predict_genecats_MAP_mbc, test = TRUE) #
# hm this iss confusing cuz none appear to be sig......

names(both_funs_enz_noout_MBC)
ggplot (both_funs_enz_noout_MBC, aes (x = `Ph Stress`, y = total_mbc)) + 
  geom_point() + 
  facet_grid (~SamplingPt) # nothing honey

ggplot (both_funs_enz_noout_MBC, aes (x = `Carbohydrate Transport`, y = total_mbc)) + 
  geom_point() + 
  facet_grid (~SamplingPt) # nothing honey

ggplot (both_funs_enz_noout_MBC, aes (x = MAP, y = total)) + 
  geom_point() + 
  facet_grid (~SamplingPt)

########################################################################
# correlation stress and resource
########################################################################

dim(stress_sum)
dim(resource_sum)

names(stress_sum)
names(resource_sum)

stressandres <- merge(resource_sum, stress_sum, by = c("sample","MAP", "SamplingPt","Site","PercMoist","MAP.std",      
                                                       "PercMoist.std"  , "PercN.std", "pH.std" , "MBC.std", "DOC.std"))
names(stressandres)
stressandres <- stressandres  %>% 
  dplyr::rename(   "resource_relabun" = "relabun.x",    "stress_relabun" = "relabun.y")
names(stressandres)

ggplot(stressandres, aes(x = resource_relabun, y = stress_relabun, color = MAP)) +
  geom_point() + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + # this is where we took other colors from I think
  theme_classic() 

stress_vs_resource <- ggplot(stressandres, aes(x = resource_relabun, y = stress_relabun, color = MAP)) +
  geom_point() + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + # this is where we took other colors from I think
  theme_classic(base_size = 14) + 
  facet_wrap(~SamplingPt)

#ggsave("stressvsresource_11024.png", plot = stress_vs_resource, device = "png",
#    width = 8, height = 3, dpi = 300)

stress_vs_resource 

summary(lm(resource_relabun~stress_relabun, data = stressandres)) # nope
summary(lm(resource_relabun~stress_relabun*SamplingPt, data = stressandres)) # nay lady.

cor.test(stressandres$stress_relabun, stressandres$resource_relabun) # nope

# correlation by group
library(rstatix) 
stressandres %>% 
  group_by ( SamplingPt) %>% 
  cor_test(stress_relabun, resource_relabun, method = "pearson") # nay lady
citation("rstatix") 



########################################################################
# NMDSs and permanovas
########################################################################

# first - total.
dim(sitexsp_wmeta) #3039
head(sitexsp_wmeta)


#permanova_all <- adonis2(sitexsp~ sitexsp_wmeta$MAP.std* sitexsp_wmeta$PercMoist.std * sitexsp_wmeta$SamplingPt )
#permanova_all
set.seed(42)
permanova_all_season <- adonis2(sitexsp~ sitexsp_wmeta$MAP.std*  sitexsp_wmeta$SamplingPt , iter = 9999)
permanova_all_season # map matters thats it :)
#sitexsp_wmeta$MAP.std                           1   0.3461 0.03679 2.1236  0.002 **
#sitexsp_wmeta$SamplingPt                        2   0.3594 0.03821 1.1028  0.181   
#sitexsp_wmeta$MAP.std:sitexsp_wmeta$SamplingPt  2   0.3911 0.04157 1.1998  0.096 . 

set.seed(42)
permanova_all_moist <- adonis2(sitexsp~ sitexsp_wmeta$MAP.std* sitexsp_wmeta$PercMoist.std , iter = 9999 )
permanova_all_moist # MAP matters thats it :)
#sitexsp_wmeta$MAP.std                              1   0.3461 0.03679 2.0996  0.002 **
#sitexsp_wmeta$PercMoist.std                        1   0.1426 0.01516 0.8654  0.730   
#sitexsp_wmeta$MAP.std:sitexsp_wmeta$PercMoist.std  1   0.1833 0.01948 1.1118  0.226 
# these are both done and in results table 3-25



nmds_all <- metaMDS(sitexsp, distance="bray",  trymax= 500 )

ggplot_nmds <- as.data.frame(nmds_all[["points"]])
ggplot_nmds$MAP = sitexsp_wmeta$MAP
ggplot_nmds$SamplingPt = sitexsp_wmeta$SamplingPt


# MAP matters but not naything else!!!!!!!!!!!!!!! so just have one NMDS


mycol.rev <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
mycol = rev(mycol.rev)

nmds <- ggplot(ggplot_nmds, aes(x = MDS1, y = MDS2)) +  # do not need treatament in here :()) +
  geom_point(aes (color =  MAP, shape = SamplingPt), size = 3, alpha= 0.7) + 
  #scale_color_gradientn(colours = mycol) + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + # this is where we took other colors from I think
  #scale_color_viridis (direction = -1) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  theme_classic(base_size = 16)     
nmds

#ggsave("overallNMDS_microtrait_11024.png", plot = nmds, device = "png",
#     width = 5.5, height = 4, dpi = 300)





# categories_matrix 


# NEXT  - only on 18 categories coonsidered in 
dim(both_funs_co2) #57x27
head(both_funs_co2)
names(both_funs_co2)

set.seed(42)
permanova_funcats_season<- adonis2(categories_matrix~MAP.std*SamplingPt,  
                            data = both_sum_fun_wide, iter = 9999) # yes MAP matters. 
permanova_funcats_season # MAP

set.seed(42)
permanova_funcats_moist<- adonis2(categories_matrix~MAP.std*PercMoist.std,  
                                   data = both_sum_fun_wide, iter = 9999) # yes MAP matters. 
permanova_funcats_moist # MAP*Moissst..... 



set.seed(42069)
nmds_funcats <- metaMDS(as.matrix(both_funs_co2[,8:25]), distance="bray",  trymax= 500 )
nmds_funcats
ggplot_funcats_nmds <- as.data.frame(nmds_funcats[["points"]])
ggplot_funcats_nmds$MAP = both_funs_co2$MAP
ggplot_funcats_nmds$SamplingPt = both_funs_co2$SamplingPt

categories.nmds.scores <- as.data.frame(scores(nmds_funcats, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
categories.nmds.scores$species <- rownames(categories.nmds.scores)  # create a column of species, from the rownames of species.scores
categories.nmds.scores$Category <- gsub("Transport" ,"trans.", categories.nmds.scores$species)
categories.nmds.scores$Category <- gsub("Carbohydrate" ,"carb.", categories.nmds.scores$Category)
categories.nmds.scores$Category <- gsub("Temperature" ,"Temp", categories.nmds.scores$Category)
categories.nmds.scores$Category <- gsub("Depolymerization" ,"depoly.", categories.nmds.scores$Category)
categories.nmds.scores$Category <- gsub("Free Amino Acids" ,"A.A.", categories.nmds.scores$Category)
categories.nmds.scores$Category <- gsub("Ph Stress" ,"pH stress", categories.nmds.scores$Category)
categories.nmds.scores$Category <- gsub("Desiccation/Osmotic/Salt Stress" ,"Desiccation stress", categories.nmds.scores$Category)
categories.nmds.scores$Category 

mycol.rev <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
mycol = rev(mycol.rev)

nmds_funcats <- ggplot(ggplot_funcats_nmds, aes(x = MDS1, y = MDS2)) +  # do not need treatament in here :()) +
  geom_point(aes (color =  MAP, shape = SamplingPt), size = 3, alpha= 0.7) + 
  #scale_color_gradientn(colours = mycol) + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + # this is where we took other colors from I think
  #scale_color_viridis (direction = -1) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  geom_text(data=categories.nmds.scores,aes(x=NMDS1,y=NMDS2,label=Category),alpha=0.3, size =3.5) +  # add the species labels
  theme_classic(base_size = 16)     
nmds_funcats

#ggsave("funcatsNMDS_microtrait_11024.png", plot = nmds_funcats, device = "png",
#    width = 5.5, height = 4, dpi = 300)



##### STRESS
head(stress)
length(unique(stress$KO))

stress_sel <- stress %>% select(c(KO, sample, counts, MAP.std, MAP, PercMoist.std, PercMoist, Site, Year, SamplingPt))
stress_sitexsp_wmeta <- stress_sel  %>%  pivot_wider(names_from = KO, values_from = counts , values_fill = 0)
head(stress_sitexsp_wmeta)
dim(stress_sitexsp_wmeta)
dim (stress_sitexsp_wmeta[,9:71])
stress_sitexsp_wmeta[,71]
stress_sitexsp_wmeta[,1:10] # so wan

set.seed(42)
permanova.stress.season <- adonis2(stress_sitexsp_wmeta[,9:71]~ stress_sitexsp_wmeta$MAP.std *stress_sitexsp_wmeta$SamplingPt,
                                   iter = 9999)
permanova.stress.season # marg MAP

set.seed(42)
permanova.stress.moist <- adonis2(stress_sitexsp_wmeta[,9:71]~ stress_sitexsp_wmeta$MAP.std *stress_sitexsp_wmeta$PercMoist.std,
                                  iter = 9999)
permanova.stress.moist # also marg MAP

nmds_stress <- metaMDS(stress_sitexsp_wmeta[,9:71], distance="bray" , trymax = 1000)

ggplot_nmds_stress <- as.data.frame(nmds_stress[["points"]])
ggplot_nmds_stress$MAP = stress_sitexsp_wmeta$MAP
ggplot_nmds_stress$SamplingPt = stress_sitexsp_wmeta$SamplingPt

nmds_stress <- ggplot(ggplot_nmds_stress, aes(x = MDS1, y = MDS2)) +  # do not need treatament in here :()) +
  geom_point(aes (color =  MAP, shape = SamplingPt), size = 3, alpha = 0.7) + 
  #scale_color_viridis (direction = -1) + 
  #scale_color_gradientn(colours = mycol) + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  theme_classic(base_size = 16)     
nmds_stress

#ggsave("stressNMDS_microtrait_11024.png", plot = nmds_stress, device = "png",
#      width = 5.5, height = 4, dpi = 300)


nmds_stress_axis2 <- ggplot(ggplot_nmds_stress, aes(x = MAP, y = MDS2)) +  # do not need treatament in here :()) +
  geom_point(aes ( shape = SamplingPt), size = 3) + 
  scale_color_viridis (direction = -1) + 
  geom_smooth(method = lm) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  theme_classic()

nmds_stress_axis2

summary(lm (MDS2~MAP, data = ggplot_nmds_stress))
# very strong relationship with NMDS axis 2!

#ggsave("stressNMDS_axis2_ 5-5-23.png", plot = nmds_stress_axis2, device = "png",
#      width = 4, height = 3, dpi = 300)




##### RESOURCE
head(resource)

resource_sel <- resource %>% select(c(KO, sample, counts, MAP,MAP.std, PercMoist, PercMoist.std, Site, Year, SamplingPt))
resource_sitexsp_wmeta <- resource_sel  %>%  pivot_wider(names_from = KO, values_from = counts , values_fill = 0)
head(resource_sitexsp_wmeta)
dim(resource_sitexsp_wmeta)
resource_sitexsp_wmeta[,324]
resource_sitexsp_wmeta[,1:10] # so wan
dim (resource_sitexsp_wmeta[,9:324])

set.seed(42)
permanova.resource.season <- adonis2(resource_sitexsp_wmeta[,9:324]~ resource_sitexsp_wmeta$MAP.std*resource_sitexsp_wmeta$SamplingPt,
                                     iter = 9999)
permanova.resource.season # MAP

set.seed(42)
permanova.resource.moist <- adonis2(resource_sitexsp_wmeta[,9:324]~ resource_sitexsp_wmeta$MAP.std*resource_sitexsp_wmeta$PercMoist.std, 
                                    iter = 9999)
permanova.resource.moist


set.seed(42069)
nmds_resource <- metaMDS(resource_sitexsp_wmeta[,9:324], distance="bray" , trymax = 500)

ggplot_nmds_resource <- as.data.frame(nmds_resource[["points"]])
ggplot_nmds_resource$MAP = resource_sitexsp_wmeta$MAP
ggplot_nmds_resource$SamplingPt = resource_sitexsp_wmeta$SamplingPt

nmds_resource <- ggplot(ggplot_nmds_resource, aes(x = MDS1, y = MDS2)) +  # do not need treatament in here :()) +
  geom_point(aes (color =  MAP, shape = SamplingPt), size = 3) + 
  #scale_color_viridis (direction = -1) + 
  #scale_color_gradientn(colours = mycol) + 
  paletteer::scale_colour_paletteer_c("grDevices::Zissou 1", direction = -1) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  theme_classic(base_size = 16)     
nmds_resource

#ggsave("resourceNMDS_microtrait_11024.png", plot = nmds_resource, device = "png",
#      width = 5.5, height = 4, dpi = 300)
#


nmds_resource_axis1 <- ggplot(ggplot_nmds_resource, aes(x = MAP, y = MDS1)) +  # do not need treatament in here :()) +
  geom_point(aes ( shape = SamplingPt), size = 3) + 
  scale_color_viridis (direction = -1) + 
  #geom_smooth(method = lm) + 
  #stat_ellipse() + 
  #scale_color_manual(values = c("#A12919","#CE5732", "#979461" , "#9bafba", "#687F96") )+
  #scale_shape_manual(values=c(15:19))+
  theme_classic()     

#ggsave("resourceNMDS_axis1_ 5-5-23.png", plot = nmds_resource_axis1, device = "png",
#       width = 4, height = 3, dpi = 300)

summary(lm (MDS1 ~MAP, data = ggplot_nmds_resource) )
# not really.




########################################################
# specific genes in categories we care about
########################################################

# drought stress genes 
unique(stress$microtrait_trait.name1)
droughtgenes <- stress %>% 
  filter(microtrait_trait.name1 == "Stress Tolerance:Specific:desiccation/osmotic/salt stress" )

droughtgenes$KO
table(droughtgenes$KO)


# c degradation genes 
unique(resource$microtrait_trait.name1)
cdeggenes <- resource %>% 
  filter(microtrait_trait.name1 == "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization" )
table(cdeggenes$KO)

# transport genes 
unique(resource$microtrait_trait.name1)
transgenes <- resource %>% 
  filter( str_detect (microtrait_trait.name1, "ransport" ) )
table(transgenes$KO)

dim(droughtgenes)
dim(cdeggenes)
dim(transgenes)
names(transgenes)
goi <- rbind (droughtgenes, cdeggenes, transgenes)

tt <- table(goi$label)
goi_abun <- subset(goi, label %in% names(tt[tt > 29]))

names(goi_abun)
unique(goi_abun$KO)

# need to make into wide format then run manylm
goi_abun_wide <- pivot_wider(goi_abun, names_from = KO, values_from = relabun ,values_fill = 0)
names(goi_abun_wide)
manyglm_goi <- manylm(as.matrix(goi_abun_wide[,18:220])~MAP.std*SamplingPt,  data = goi_abun_wide)
summary(manyglm_goi)
anova.manylm(manyglm_plhyum_stress, p.uni="adjusted")










#




#




#

########################################################
# let's try taxonomy again
########################################################

head (proteins.ko$eggNOG_OGs)
hist (lengths(strsplit(proteins.ko$eggNOG_OGs, split = ",")  ) )

# what if the firsts step was to take out everything before first "root
length(unique(proteins.ko$KO))
length(unique(resource$KO))
proteins.ko$eggnog_OG_noroot <- sapply(strsplit(proteins.ko$eggNOG_OGs, "root,"), tail, 1)
#proteins.ko$eggnog_OG_noroot <- sapply(strsplit(proteins.ko$eggNOG_OGs, "*.Archaea,"), tail, 1)
proteins.ko$eggnog_OG_noroot
proteins.ko.tax.init <- proteins.ko %>%    separate(eggnog_OG_noroot, c("Domain", "Phylum", "Class", "Order", "Delete"), ",")
proteins.ko.tax.init$eggNOG_OGs <-  proteins.ko$eggNOG_OGs
proteins.ko.tax.init$Domain
#proteins.ko.tax.init %>%  mutate (across (Domain:Order,str_replace_all ( regex("\\W+"), " ")) )  
# this is notworking but I do not remember what it did
proteins.ko.tax.init$Domain
proteins.ko.tax.init$Domain <- gsub(".*\\|","",proteins.ko.tax.init$Domain) # HELL YES !!!!!!!!
table ( proteins.ko.tax.init$Domain) # hm what's going on with Root
# i think | was a special character so needed the \\ to show it was actually a thing to match !!!
proteins.ko.tax.init$Phylum <- gsub(".*\\|","",proteins.ko.tax.init$Phylum)
table (proteins.ko.tax.init$Phylum) # Still somehow some root, some Bacteria, some Archaea
# somehow sub im the value from the class column
table(is.na (proteins.ko.tax.init$Phylum) )
length(proteins.ko.tax.init$Phylum) # 158467
table(is.na(proteins.ko.tax.init$Phylum) ) # 15839 are NA
# i believe it was about 700 that got taken out for only having bacteria or archaea below
(15839 + 690 + 60) /  158467 * 100 # 10% got taken out
proteins.ko.tax.init$Class <- gsub(".*\\|","",proteins.ko.tax.init$Class)
proteins.ko.tax.init$Order <- gsub(".*\\|","",proteins.ko.tax.init$Order)

# replace Phylum designation with class if Bacteria or Archaea !!!!!


proteins.ko.tax.init2 <-proteins.ko.tax.init %>% select (-c (Delete))
proteins.ko.tax.init3 <- proteins.ko.tax.init2[!is.na(proteins.ko.tax.init2$Phylum),]
proteins.ko.tax <- proteins.ko.tax.init3[proteins.ko.tax.init3$Phylum!="root",]
dim (proteins.ko.tax.init2)
dim (proteins.ko.tax)


table(is.na(proteins.ko.tax$Domain)) # all have domain
table(is.na(proteins.ko.tax$Phylum)) # all now have Phylum
table(is.na(proteins.ko.tax$Class))  # more than half haave classs
table(is.na(proteins.ko.tax$Order)) # most no order



# COME BACK HERE TO ALTER CATEGORIES
# DO I WANT THIS TO BE TWO CATEGORIES OR THREE CATEGORIES???
proteins.ko.tax <-proteins.ko.tax  %>% ungroup(.) %>% mutate(wetdry_MAP = case_when(
  proteins.ko.tax$MAP <= 600 ~ "Dry", 
  proteins.ko.tax$MAP > 600 & proteins.ko.tax$MAP < 800 ~ "Mid", 
  proteins.ko.tax$MAP >= 800 ~ "Wet"
))


hist(bothmeta$PercMoist, breaks = 30)
proteins.ko.tax <-proteins.ko.tax  %>% ungroup(.) %>% mutate(wetdry_moist = case_when(
  proteins.ko.tax$PercMoist <= 14 ~ "Dry", 
  proteins.ko.tax$PercMoist > 14  & proteins.ko.tax$PercMoist < 23  ~ "Mid", 
  proteins.ko.tax$PercMoist >= 23 ~ "Wet"
))

# SAME FOR SOIL MOISTURE?

table ( proteins.ko.tax$ wetdry_MAP) # hm, more wet than dry.
table ( proteins.ko.tax$ wetdry_moist)

table ( proteins.ko.tax$sample)
table ( proteins.ko.tax$Phylum)
# Before: 7194 Bacteria, 1626 Archaea
# check out what's going on here
#proteins.ko.tax <- proteins.ko.tax %>% 
#  mutate(Phylum = ifelse(Phylum == 'Bacteria' | Phylum == 'Archaea', Class, Phylum))
proteins.ko.tax <- proteins.ko.tax %>% 
  mutate(Phylum = ifelse(Phylum == 'Bacteria', Class, Phylum))
proteins.ko.tax <- proteins.ko.tax %>% 
  mutate(Phylum = ifelse(Phylum == 'Archaea', Class, Phylum))
table ( proteins.ko.tax$Phylum)
# # Before: 7194 Bacteria, 1626 Archaea
# NOW: 692 Bacteria, 60 Archaea



##########################################
# stress!
##########################################
names(proteins.ko.tax)
table(proteins.ko.tax$YAS)
proteins.ko.tax.stress <- proteins.ko.tax %>%  filter (YAS == "Stress Tolerance" & 
                                                         Phylum != "Bacteria" & 
                                                         Phylum != "Archaea" & 
                                                         Phylum != "unclassified Bacteria")
table(is.na(proteins.ko.tax.stress$Phylum)) # almost 90% got phylum
table(is.na(proteins.ko.tax.stress$Order)) # most NO order
table(proteins.ko.tax.stress$Phylum)
# what to take out? NA, Archaea, Bacteria? 


# WHAT DO WE WANT OUR DATASET TO LOOK LIKE?
# proportions for each sample, or like overall in each cat?
# started out the second way but think first way is necessary for stats.
# so for bubble plot, should I take the average of these averages ?


proteins.phylum.stress.totalcounts <- proteins.ko.tax.stress %>% 
  group_by(sample, MAP, MAP.std, PercMoist, PercMoist.std, wetdry_moist, SamplingPt) %>% 
  filter(!is.na(Phylum)) %>% 
  dplyr::count() %>% 
  dplyr::rename ("totalphycounts" = "n") # this is total counts of reads with a phylum in stress, 
# per sample
head(proteins.phylum.stress.totalcounts)

proteins.phylum.stress.counts <- proteins.ko.tax.stress %>% 
  group_by(sample, MAP, MAP.std, PercMoist, PercMoist.std, wetdry_moist, SamplingPt, Phylum) %>% 
  filter(!is.na(Phylum)) %>% # this is counts of each phylum in eaach sample for stress.
  dplyr::count() %>% 
  dplyr::rename ("counts" = "n")

proteins.phylum.stress <- merge (proteins.phylum.stress.counts, proteins.phylum.stress.totalcounts,
                                 by = c ("sample", "wetdry_moist", "MAP", "MAP.std", "PercMoist", "PercMoist.std","SamplingPt"))
proteins.phylum.stress$prop <- proteins.phylum.stress$counts / proteins.phylum.stress$totalphycounts *100
# THIS ISS PERCENT OF TOTAL KOS THAT HAVE PHYLUM DESIGNATION IN EACH PHYLUM, FOR EACH SAMPLE
head(proteins.phylum.stress)
#View(proteins.phylum.stress)

#CHANGE THIS BACK AND FORTH IF DOING WETDRY_MOIST or PERCMOIST CONTINUOUS !!!!!!!!!!


proteins.phylum.stress %>% 
  group_by(Phylum, wetdry_moist, SamplingPt) %>% 
  summarize(mean = mean (prop)) %>% 
  View()

proteins.phylum.stress %>% 
  group_by(sample) %>% summarize(check = sum (prop)) # all sum to 100 for every sample !

range ( proteins.phylum.stress$prop)









# manyglm

# FIRST get rid of  counts column thats fucking us up, need to get one row for each sample
proteins.phylum.stress.sel <- proteins.phylum.stress %>% 
  select(-c (counts))
names(proteins.phylum.stress.sel)
# then, go to wide pivot_wider(names_from = KO, values_from = relabun ,values_fill = 0)
table(proteins.phylum.stress.sel$Phylum)
# subset only abundant phyla 

tt = table(proteins.phylum.stress.sel$Phylum)
tt
proteins.phylum.stress_onlyabun <- subset(proteins.phylum.stress.sel, Phylum %in% names(tt[tt > 28]))


proteins.phylum.stress_wide <- pivot_wider(proteins.phylum.stress_onlyabun, names_from = Phylum, values_from = prop ,values_fill = 0)
names(proteins.phylum.stress_wide)
dim(proteins.phylum.stress_wide) # this is not averaged over anything !!!! 
# so we should be good to toggle back and forth between 
manyglm_phylum_stress <- manylm(as.matrix(proteins.phylum.stress_wide[,9:16])~PercMoist.std*SamplingPt,  data = proteins.phylum.stress_wide)
summary(manyglm_phylum_stress) # R2 only 0.02
anova.manylm(manyglm_phylum_stress, p.uni="adjusted")
# categories: acidobacteria, cyanobacteria
# NOTHING if in linear model with PercMoist.std


phylum_acido_moist <- lm (Acidobacteria ~PercMoist.std*SamplingPt, data = proteins.phylum.stress_wide)
summary(phylum_acido_moist )
Anova(phylum_acido_moist, type = "III")  # nothing honey.....



# THESE GRAPHS SHOULD INCLUDE ZEROES!!!! which it does. good.
ggplot(data = proteins.phylum.stress_wide,
       aes(x = MAP, y = Acidobacteria)  ) + 
  stat_summary(fun = "mean", geom = "point", alpha = 0.8) + 
  geom_smooth(method = "lm") + 
  theme_light() + 
  facet_grid(cols = vars(SamplingPt)) # this includes zeroess, better. 
# thiss shows just now many zeroes we are dealing with here.
# seems like the pattern is real??




ggplot(data = proteins.phylum.stress_wide,
       aes(x = MAP, y = Cyanobacteria)  ) + 
  stat_summary(fun = "mean", geom = "point", alpha = 0.8) + 
  geom_smooth(method = "lm") + 
  theme_light() + 
  facet_grid(cols = vars(SamplingPt)) # this includes zeroes.

#this feels like it's largely driven by one site at one sampling point :(


ggplot(data = proteins.phylum.stress_wide,
       aes(x = wetdry_moist, y = Acidobacteria)  ) + 
  geom_boxplot() + 
  geom_smooth(method = "lm") + 
  theme_light() + 
  facet_grid(cols = vars(SamplingPt))
ggplot(data = proteins.phylum.stress_wide,
       aes(x = wetdry_moist, y = Cyanobacteria)  ) + 
  geom_boxplot() + 
  geom_smooth(method = "lm") + 
  theme_light() + 
  facet_grid(cols = vars(SamplingPt))




# bubble graphs that INCLUDE ZEROES
# so make 

names(proteins.phylum.stress_wide)

proteins.phylum.stress_backtolong <- proteins.phylum.stress_wide %>% 
  pivot_longer( cols = "Actinobacteria":"Cyanobacteria" ,
                names_to = "Phylum",
                values_to = "prop")

proteins.phylum.stress.avg <- proteins.phylum.stress_backtolong %>% 
  group_by( wetdry_moist, SamplingPt, Phylum) %>% 
  summarize(prop = mean (prop))

proteins.phylum.stress.avg.persample <- proteins.phylum.stress_backtolong %>% 
  group_by( PercMoist.std, sample, SamplingPt, Phylum) %>% 
  summarize(prop = mean (prop))

hist(proteins.phylum.stress.avg$prop)

# bubble plot -- this might be able to work with the full MAP dataset.... but how to handle multiple site issue?
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")

bubble_stress= ggplot(data=proteins.phylum.stress.avg, 
                      aes(x = wetdry_moist, y = forcats::fct_rev(Phylum))) + 
  geom_point(aes(size = prop, fill = Phylum), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(0.25,17), breaks = c(1,5,10,25,50, 80)) + 
  labs( x= "", y = "", size = "Proportion", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = NULL) +  
  #scale_fill_manual(values = colours, guide = FALSE) + 
  facet_grid (~SamplingPt)

bubble_stress

##ggsave("Bubble plot stress microtrait 6-14-23.png", plot = bubble_stress, device = "png",
#      width = 9.5, height = 8, dpi = 300)




##########################################
# resource!
##########################################

table (proteins.ko.tax$YAS)
proteins.ko.tax.resource <- proteins.ko.tax %>%  filter (YAS == "Resource Acquisition" & 
                                                           Phylum != "Bacteria" & 
                                                           Phylum != "Archaea" & 
                                                           Phylum != "unclassified Bacteria")
table(is.na(proteins.ko.tax.resource$Phylum)) # almost 90% got phylum
table(is.na(proteins.ko.tax.resource$Order)) # most NO order
table(proteins.ko.tax.resource$Phylum)
# what to take out? NA, Archaea, Bacteria? 


# WHAT DO WE WANT OUR DATASET TO LOOK LIKE?
# proportions for each sample, or like overall in each cat?
# started out the second way but think first way is necessary for stats.
# so for bubble plot, should I take the average of these averages ?


proteins.phylum.resource.totalcounts <- proteins.ko.tax.resource %>% 
  group_by(sample, wetdry_MAP, SamplingPt, MAP, MAP.std, PercMoist, PercMoist.std) %>% 
  filter(!is.na(Phylum)) %>% 
  dplyr::count() %>% 
  dplyr::rename ("totalphycounts" = "n") # this is total counts of reads with a phylum in resource, 
# per sample

proteins.phylum.resource.counts <- proteins.ko.tax.resource %>% 
  group_by(sample, wetdry_MAP, SamplingPt, Phylum, MAP, MAP.std, PercMoist, PercMoist.std)  %>% 
  filter(!is.na(Phylum)) %>% 
  dplyr::count() %>% 
  dplyr::rename ("counts" = "n")

proteins.phylum.resource <- merge (proteins.phylum.resource.counts, proteins.phylum.resource.totalcounts,
                                   by = c ("sample", "wetdry_MAP", "SamplingPt", "MAP", "MAP.std", "PercMoist", "PercMoist.std"))
proteins.phylum.resource$prop <- proteins.phylum.resource$counts / proteins.phylum.resource$totalphycounts * 100

range ( proteins.phylum.resource$prop)






# FIRST get rid of  counts column thats fucking us up, need to get one row for each sample
proteins.phylum.resource.sel <- proteins.phylum.resource %>% 
  select(-c (counts))
names(proteins.phylum.resource.sel)



table(proteins.phylum.resource.sel$Phylum)
# subset only abundant phyla 
57*.75

tt = table(proteins.phylum.resource.sel$Phylum)
proteins.phylum.resource_onlyabun <- subset(proteins.phylum.resource.sel, Phylum %in% names(tt[tt > 28]))

proteins.phylum.resource_wide <- pivot_wider(proteins.phylum.resource_onlyabun, names_from = Phylum, values_from = prop ,values_fill = 0)
names(proteins.phylum.resource_wide)
dim(proteins.phylum.resource_wide)
#nothing was sig 
manyglm_phylum_resource <- manylm(as.matrix(proteins.phylum.resource_wide[,9:20])~MAP.std*SamplingPt,  data = proteins.phylum.resource_wide)
summary(manyglm_phylum_resource) # R2 better here.
anova.manylm(manyglm_phylum_resource, p.uni="adjusted")
# Thaumarchaeota, cyanobacteria, is the only one significant. when map is linear
# when categories, cyanobacteria only.


phylum_resource_cyano <- lm (Cyanobacteria~MAP.std*SamplingPt, proteins.phylum.resource_wide)
phylum_resource_cyano.slopes <- emtrends(phylum_resource_cyano,~ SamplingPt, var="MAP.std")
summary(phylum_resource_cyano.slopes , infer=c(TRUE,TRUE)) 
# sooooo in fall 2015 negative trend - relative abundnace of 



ggplot(data = proteins.phylum.resource_wide,
       aes(x = MAP, y = Cyanobacteria)  ) + 
  stat_summary(fun = "mean", geom = "point", alpha = 0.8) + 
  geom_smooth(method = "lm") + 
  theme_light() + 
  facet_grid(cols = vars(SamplingPt)) # this includes zeroes.



# bubble graphs that INCLUDE ZEROES
# so make 

names(proteins.phylum.resource_wide)

proteins.phylum.resource_backtolong <- proteins.phylum.resource_wide %>% 
  pivot_longer( cols = "Acidobacteria":"Thaumarchaeota" ,
                names_to = "Phylum",
                values_to = "prop")

proteins.phylum.resource.avg <- proteins.phylum.resource_backtolong %>% 
  group_by( wetdry_MAP, SamplingPt, Phylum) %>% 
  summarize(prop = mean (prop))

proteins.phylum.resource.avg.persample <- proteins.phylum.resource_backtolong %>% 
  group_by( MAP.std, sample, SamplingPt, Phylum) %>% 
  summarize(prop = mean (prop))

hist(proteins.phylum.resource.avg$prop)

# bubble plot -- this might be able to work with the full MAP dataset.... but how to handle multiple site issue?
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")

bubble_resource= ggplot(data=proteins.phylum.resource.avg, 
                        aes(x = wetdry_MAP, y = forcats::fct_rev(Phylum))) + 
  geom_point(aes(size = prop, fill = Phylum), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(0.25,17), breaks = c(1,5,10,25,50, 80)) + 
  labs( x= "", y = "", size = "Proportion", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = NULL) +  
  #scale_fill_manual(values = colours, guide = FALSE) + 
  facet_grid (~SamplingPt)

bubble_resource

##ggsave("Bubble plot resource microtrait 6-14-23.png", plot = bubble_resource, device = "png",
#      width = 9.5, height = 8, dpi = 300)












###############################################################
# Cazy??
###############################################################
proteins.ko.cazy <- proteins.ko %>% 
  filter(CAZy != "-")

dim(proteins.ko.cazy)
head(proteins.ko.cazy)
unique(proteins.ko.cazy$CAZy)


proteins.ko.GH <- proteins.ko %>% 
  filter(str_detect(CAZy, "GH"))

dim(proteins.ko.GH)
head(proteins.ko.GH)
unique(proteins.ko.GH$CAZy)



# get into long format
data_counts_GH <- proteins.ko.GH %>% 
  group_by(sample, MAP, PercMoist, SamplingPt, Site, Year,   KO , YAS, microtrait_trait.name1,
           MAP.std, PercMoist.std, PercClay.std, PercC.std, PercN.std, pH.std) %>% 
  dplyr::count()
data_counts_GH <- data_counts_GH %>% 
  dplyr::rename("counts" = "n")

names (data_counts_GH)

total_counts_GH <- proteins.ko.GH %>% 
  group_by(sample, MAP, SamplingPt, PercMoist, Site, Year) %>% 
  dplyr::count()

total_counts_GH <- total_counts_GH %>% 
  dplyr::rename("totalsample" = "n")
total_counts # 42 has 10,000 - this is the wettest site in sum 2016


data_counts_withtotal_GH <- merge(data_counts_GH, total_counts_GH, by = c("sample", "MAP", "PercMoist", "Site", "SamplingPt", "Year") )
data_counts_withtotal_GH$relabun<- data_counts_withtotal_GH$counts /  data_counts_withtotal_GH$totalsample * 100
data_counts_withtotal_GH$relabun
# this will be used to compare with TMM  below
str(data_counts_withtotal_GH)
dim(data_counts_withtotal_GH)



# enzyme data
enzymes <- read.csv("TXRAPID_Enz_forCaitlin.csv")
head(enzymes)
#enzymes$SamplingPt  <- paste("Season", "Year", sep = " ")
enzymes_meta <- merge ( enzymes, trt, by = c ("Site", "SamplingPt") )


enzymes_meta[enzymes_meta$CBH > 2000,] # LBJWFC 2016  is a universal stinker?
head(enzymes_meta)


enzymes_meta_noout <- enzymes_meta %>% 
  filter(CBH < 2000)

ggplot (data = enzymes_meta_noout, aes(x = MAP , y = BG)) + 
  geom_point() + 
  theme_light() + 
  facet_wrap(~SamplingPt)

summary(lm(BG~MAP, enzymes_meta_noout)) # marg

ggplot (data = enzymes_meta_noout, aes(x = MAP , y = CBH)) + 
  geom_point() + 
  theme_light() + 
  facet_wrap(~SamplingPt)

summary(lm(CBH~MAP, enzymes_meta_noout))  # nope

#α-glucosidase 
ggplot (data = enzymes_meta_noout, aes(x = MAP , y = AG)) + 
  geom_point() + 
  theme_light() + 
  facet_wrap(~SamplingPt)

summary(lm(AG~MAP, enzymes_meta_noout))  # nope

# both glucosidases
enzymes_meta_noout$gluc <- enzymes_meta_noout$BG + enzymes_meta_noout$AG

ggplot (data = enzymes_meta_noout, aes(x = MAP , y = gluc)) + 
  geom_point() + 
  theme_light() + 
  facet_wrap(~SamplingPt)

summary(lm(gluc~MAP, enzymes_meta_noout))  # nope



