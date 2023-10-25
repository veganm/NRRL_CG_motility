# Script for processing data from NRRL synthetic community evolution run 2
pacman::p_load(tidyverse, cowplot, vegan, FactoMineR, readxl, ggpubr, MASS, RColorBrewer,
               readxl, growthrates, growthcurver)

load(".Rdata")

save.image(".Rdata")

# Global variables
xTextSize<-14

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               FIGURE 1
# Read in multispecies data from batch worm digests

NRRLMulti2<-read_excel("ExperimentalEvolutionPlateData.xlsx", sheet="MultispeciesTimeSeries")
names(NRRLMulti2)

#grab the count data for each of the 7 bacteria
MultiLog<-decostand(NRRLMulti2[,5:11], method="log") 
MStandLog<-decostand(MultiLog, method="standardize") 
MStandLog$Pass<-NRRLMulti2$Pass 
pcaM<-PCA(MStandLog, quali.sup=8) 


# Extract the first two components
pcaM.top2<-data.frame(PC1=pcaM$ind$coord[,1],
                      PC2=pcaM$ind$coord[,2],
                      Pass=as.factor(NRRLMulti2$Pass))
pcaM$eig

#colvec1<-c("black", "orangered", "orange1", "gold", "yellowgreen", "green3", "blue", "lilac", "darkorchid4") 
colvec1<-c("#000000", "#FF4500", "#FFA500", "#FFD700", "#9ACD32", "#00CD00", "#0000FF", "#EE82EE", "#7D26CD") 
pPCAM<-pcaM.top2 %>%
  ggplot(aes(x=PC1,y=PC2,color=Pass)) +
  geom_point(size=2) +
  #theme_classic() +
  theme(axis.text=element_text(size=xTextSize), 
        axis.title=element_text(size=xTextSize),
        legend.text=element_text(size=xTextSize),
        plot.title=element_text(hjust=0.5, size=xTextSize+2),
        #legend.title = element_text(size=myTextSize),
        legend.position = "none") + 
  labs(title="Worm communities", y="PC2 (18.6%)", x="PC1 (31.4%)")
pPCAM + scale_color_manual(values=colvec1)
#pPCAM + scale_color_viridis_d(option="inferno", begin=0.9, end=0)
ggsave("pNRRL2_PCAMultiWorms2.svg", width=5, height=4, units="in", dpi=400)
ggsave("pNRRL2_PCAMultiWorms2.png", width=5, height=4, units="in", dpi=400)

# make stacked bar plots
NRRLMulti2$Total<-rowSums(NRRLMulti2[,5:11])
fAA<-NRRLMulti2$AA/NRRLMulti2$Total 
fMO<-NRRLMulti2$MO/NRRLMulti2$Total 
fRE<-NRRLMulti2$RE/NRRLMulti2$Total 
fBS<-NRRLMulti2$BS/NRRLMulti2$Total 
fOA<-NRRLMulti2$OA/NRRLMulti2$Total 
fCX<-NRRLMulti2$CX/NRRLMulti2$Total 
fSX<-NRRLMulti2$SX/NRRLMulti2$Total 

values<-c(fAA,fMO,fRE,fBS,fOA,fCX,fSX) 
mylen<-length(fAA)
bact<-c(rep("AA",mylen), rep("MO",mylen), rep("RE",mylen), rep("BS",mylen), rep("OA",mylen), rep("CX",mylen), rep("SX",mylen)) 
MultiF<-data.frame(Pass=NRRLMulti2$Pass, 
                   ID=NRRLMulti2$ID, 
                   Set=NRRLMulti2$Set, 
                   Rep=NRRLMulti2$Rep, 
                   bact, values)
MultiF<-as_tibble(MultiF)
glimpse(MultiF)

#plot out stacked bars form worm-associated batch digest communities
#pMultiF_allpasses<-MultiF %>%
pMultiF<-MultiF %>%
  filter(Pass=="1" | Pass=="5" | Pass=="10") %>% # we don't need to see ALL of these surely
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple')) + 
  theme(
    axis.text=element_text(size=xTextSize+2), 
    axis.title=element_text(size=xTextSize+2),
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_blank(),
    strip.text=element_text(size=xTextSize+4),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize+2),
    #legend.position = "none")  
    )+
  labs(y="Relative Abundance", x="", title="Worm+")+
  facet_wrap(~Pass)
#pMultiF_allpasses
pMultiF
#ggsave("pNRRL2_StackBarWorms2.png", width=4, height=8, units="in", dpi=300)
#ggsave("pNRRL2_StackBarWorms2_wide.png", width=6, height=4, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's do the same for the no-worm data 
MultiNoWorm<-read.table("NRRLNoWormMulti.txt", header=TRUE) 
names(MultiNoWorm)
MultiNoWorm<-as_tibble(MultiNoWorm)
MultiNoWorm$Set<-substr(MultiNoWorm$Rep,1,1)
MultiNoWorm$ID<-MultiNoWorm$Rep
MultiNoWorm$Rep<-as.numeric(substr(MultiNoWorm$ID,2,2))
MultiNoWorm<-MultiNoWorm %>% relocate(ID, Set)
MultiNoWorm

# PCA
MultiNWLog<-decostand(MultiNoWorm[,5:11], method="log") 
MNWStandLog<-decostand(MultiNWLog, method="standardize") 
MNWStandLog$Gen<-MultiNoWorm$Pass 
pcaMNW<-PCA(MNWStandLog, quali.sup=8) 

pcaMNW.top2<-data.frame(PC1=pcaMNW$ind$coord[,1],
                      PC2=pcaMNW$ind$coord[,2],
                      Pass=as.factor(MultiNoWorm$Pass))
pcaMNW$eig
#dim1: 36.46%; dim2=16.09% 

# Plot out the first two PCs
#colvec2<-c("black", "red", "orangered", "orange1", "gold", "yellowgreen", "darkgreen", "blue", "purple", "darkorchid2") 
colvec2<-c("#000000", "#FF0000", "#FF4500", "#FFA500", "#FFD700", "#9ACD32", "#00CD00", "#0000FF", "#EE82EE", "#7D26CD") 

pPCAMNW<-pcaMNW.top2 %>%
  ggplot(aes(x=PC1,y=PC2,color=Pass)) +
  geom_point(size=2) +
  #theme_classic() +
  theme(axis.text=element_text(size=xTextSize), 
        axis.title=element_text(size=xTextSize),
        legend.text=element_text(size=xTextSize),
        plot.title=element_text(hjust=0.5, size=xTextSize+2),
        legend.title = element_text(size=xTextSize)) + 
  labs(title="No worms", y="PC2 (16.1%)", x="PC1 (36.5%)")
pPCAMNW + scale_color_manual(values=colvec2)
ggsave("pPCAMultiNoWorms.svg", width=5, height=4, units="in", dpi=300)
ggsave("pPCAMultiNoWorms.png", width=5, height=4, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~
# Ordinate both data sets together?

#grab the count data for each of the 7 bacteria
#MultiLog<-decostand(NRRLMulti2[,5:11], method="log") 
#MStandLog<-decostand(MultiLog, method="standardize") 
#MStandLog$Pass<-NRRLMulti2$Pass 

NRRLMulti2$Worms<-"Worm+"
MultiNoWorm$Worms<-"Worm-"
MultiLogMerge<-rbind(MultiLog, MultiNWLog) # combine
MultiLogMergeStand<-decostand(MultiLogMerge, method="standardize")
MultiLogMergeStand$Pass<-c(NRRLMulti2$Pass, MultiNoWorm$Pass)
MultiLogMergeStand$Worms<-c(NRRLMulti2$Worms, MultiNoWorm$Worms)

# sanity check structure
str(MultiLogMergeStand)
View(MultiLogMergeStand)

# joint PCA
pcaMNW_Merge<-PCA(MultiLogMergeStand, quali.sup=c(8,9)) 

pcaMNW_Merge$eig
#dim1: 48.87%; dim2=21.42% 

pcaMNW_Merge.top2<-data.frame(PC1=pcaMNW_Merge$ind$coord[,1],
                        PC2=pcaMNW_Merge$ind$coord[,2],
                        Pass=as.factor(MultiLogMergeStand$Pass),
                        Worms=as.factor(MultiLogMergeStand$Worms))


# Plot out the first two PCs
#colvec2<-c("black", "red", "orangered", "orange1", "gold", "yellowgreen", "darkgreen", "blue", "purple", "darkorchid2") 
#colvec2<-c("#000000", "#FF0000", "#FF4500", "#FFA500", "#FFD700", "#9ACD32", "#00CD00", "#0000FF", "#EE82EE", "#7D26CD") 

pPCAMNW_Merge<-pcaMNW_Merge.top2 %>%
  ggplot(aes(x=PC1,y=PC2,color=Pass)) +
  geom_point(size=2) +
  #theme_classic() +
  scale_color_viridis_d()+
  theme(axis.text=element_text(size=xTextSize+2), 
        axis.title=element_text(size=xTextSize+2),
        legend.text=element_text(size=xTextSize+2),
        strip.text=element_text(size=xTextSize+4),
        #plot.title=element_text(hjust=0.5, size=myTextSize+2),
        legend.title = element_text(size=xTextSize+2)) + 
  labs(y="PC2 (21.4%)", x="PC1 (48.9%)")+
  facet_wrap(~Worms)
pPCAMNW_Merge

# Back to just no-worm data
# make the stacked barplots
fNAA<-MultiNoWorm$AA/MultiNoWorm$Total 
fNMO<-MultiNoWorm$MO/MultiNoWorm$Total 
fNRE<-MultiNoWorm$RE/MultiNoWorm$Total 
fNBS<-MultiNoWorm$BS/MultiNoWorm$Total 
fNOA<-MultiNoWorm$OA/MultiNoWorm$Total 
fNCX<-MultiNoWorm$CX/MultiNoWorm$Total 
fNSX<-MultiNoWorm$SX/MultiNoWorm$Total 

values<-c(fNAA, fNMO, fNRE, fNBS, fNOA, fNCX, fNSX)
mylen2<- length(fNAA)
bact<-c(rep("AA",mylen2), rep("MO",mylen2), rep("RE",mylen2), rep("BS",mylen2), rep("OA",mylen2), rep("CX",mylen2), rep("SX",mylen2))  
MultiF_NoWorm<-data.frame(Pass=MultiNoWorm$Pass, 
                   ID=MultiNoWorm$ID, 
                   Set=MultiNoWorm$Set, 
                   Rep=MultiNoWorm$Rep, 
                   bact, values)
MultiF_NoWorm<-as_tibble(MultiF_NoWorm)
glimpse(MultiF_NoWorm)

# plot out no-worm communities barplots
pMultiF_NW<-MultiF_NoWorm %>%
  filter(Pass=="1" | Pass=="5" | Pass=="10") %>% # we don't need to see ALL of these surely
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple')) + 
  theme(
    axis.text=element_text(size=xTextSize+2), 
    axis.title=element_text(size=xTextSize+2),
    axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    strip.text=element_text(size=xTextSize+4),
    #legend.title=element_blank(),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    #legend.text=element_text(size=xTextSize),
    legend.position = "none")  +
  labs(y="Relative Abundance", x="", title="Worm-")+
  facet_wrap(~Pass)
pMultiF_NW
ggsave("pNRRL2_StackBarNoWorms.svg", width=4, height=8, units="in", dpi=300)
ggsave("pNRRL2_StackBarNoWorms_wide.png", width=6, height=4, units="in", dpi=300)

#Plot out together
#pMultiPlotsTop<-plot_grid(pPCAMNM  + scale_color_manual(values=colvec1), 
#                          pPCAM  + scale_color_manual(values=colvec2),
#          pMultiF_NW, pMultiF,
#          rel_heights=c(0.9,1), rel_widths=c(1, 1.2), labels="AUTO")
#pMultiPlotsTop
#ggsave("pNRRL2_BatchDigests_PCA_Stackbar.png", width=10, height=8, units="in", dpi=400)

# And with the collective ordination
pMultiPlotsTop_Merge_BC<-plot_grid(pMultiF_NW, pMultiF,
                          rel_widths=c(1, 1.2), labels=c("B", "C"))
pMultiPlotsTop_Merge<-plot_grid(pPCAMNW_Merge, pMultiPlotsTop_Merge_BC,
                                ncol=1, rel_heights=c(0.9,1), labels=c("A", ""))
pMultiPlotsTop_Merge


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Multispecies worm digests: When did alt morphs appear in these data?
pVariantEmergenceGrid<-NRRLMulti2 %>%
  rename("Chryseobacterium" = "varCX",
         "Microbacterium"="varMO") %>%
  pivot_longer(cols=Chryseobacterium:Microbacterium, names_to="Variant", values_to="Presence") %>%
  ggplot(aes(x=factor(Pass), y=ID, fill=Presence))+
  geom_tile(color="white", lwd=0.5)+
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=xTextSize-5), 
        axis.text.x=element_text(size=xTextSize), 
        axis.title=element_text(size=xTextSize),
        plot.title=element_text(hjust=0.5, size=xTextSize) 
  )+
  labs(y="Community", x="Pass", title="Emergence of Morphological Variants")+
  facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceGrid
#ggsave("pVariantEmergenceGrid.png", width=10, height=4, units="in", dpi=400)
  
# Plot entire Figure 1
# with the joint ordination
#plot_grid(pMultiPlotsTop_Merge, pVariantEmergenceGrid, ncol=1, rel_heights = c(1.8,1), labels=c("", "D"))
#ggsave("pNRRL2_Multispecies_AllPassSummary_Merge.png", width=10, height=12, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~
# let's heatmap the variant emergence instead
# Note that values of -1 indicate cases where variants were not counted
# but were observed on plates (passes 6-7) 
# and/or retrieved from glycerol stocks (8+)

NRRLMulti2_Variants<-read_excel("122121EvolutionNRRL2.xlsx", sheet="MultispeciesVariants")
pVariantEmergenceHeat<-NRRLMulti2_Variants %>%
  ggplot(aes(x=factor(Pass), y=ID, fill=fAlt))+
  geom_tile(lwd=0.5)+
  #scale_fill_viridis_c() +
  #scale_fill_gradient2(low="blue", mid="grey95", high="red3")+
  scale_fill_gradient2(low="blue", mid="ivory1", high="red3")+
  theme_classic()+
  theme(axis.text.y=element_text(size=xTextSize-5), 
        axis.text.x=element_text(size=xTextSize), 
        axis.title=element_text(size=xTextSize),
        plot.title=element_text(hjust=0.5, size=xTextSize),
        legend.text = element_text(size=xTextSize),
        legend.title = element_text(size=xTextSize),
        panel.background = element_rect(fill = 'grey80', color = 'grey80')
  )+
  labs(y="Community", x="Pass", title="Emergence of Morphological Variants", fill="fAlt")+
  facet_wrap(~Species, ncol=3)+
  theme(strip.text=element_text(size=xTextSize))
  #theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceHeat

# Plot Fig 1 with the joint ordination
plot_grid(pMultiPlotsTop_Merge, pVariantEmergenceHeat, ncol=1, rel_heights = c(2,1), labels=c("", "D"))
ggsave("pFig1_Multispecies_AllPassSummary_Merge_Heat_ivory1_grey80.png", width=10, height=12, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~
# Any cases where we see both MO and CG alternate morphs in the same community?
glimpse(NRRLMulti2)
NRRLMulti2<-NRRLMulti2 %>%
  mutate(varTotal = varMO + varCX)
  
max(NRRLMulti2$varTotal) #nope
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's generate diversity estimates
H<-NRRLMulti2 %>%
  dplyr::select(AA:SX) %>%
  diversity()
NRRLMulti2$H<-H
glimpse(NRRLMulti2)

D_inv<-NRRLMulti2 %>%
  dplyr::select(AA:SX) %>%
  diversity(index="invsimpson")
NRRLMulti2$D_inv<-D_inv
glimpse(NRRLMulti2)

S<-NRRLMulti2 %>%
  dplyr::select(AA:SX) %>%
  specnumber()
NRRLMulti2$S<-S
glimpse(NRRLMulti2)

pVariantEmergenceRichness_CX<-NRRLMulti2 %>%
  ggplot(aes(x=factor(Pass), y=S, color=factor(varCX)))+
  geom_jitter(width=0.15)+
  theme_classic()+
#  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Richness", x="Pass", title="Emergence of Morphological Variants: Chryseobacterium")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceRichness_CX

# Filter down to CG only
pVariantEmergenceRichness_CG<-NRRLMulti2 %>%
  filter(Set=="A" | Set=="B" | Set=="E" | Set=="F"| 
           Set=="I" | Set=="J") %>%
  ggplot(aes(x=factor(Pass), y=S, color=factor(varCX)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Richness", x="Pass", title="Emergence of Morphological Variants: C. gleum")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceRichness_CG

# Filter down to CG only and plot Shannon
pVariantEmergenceH_CG<-NRRLMulti2 %>%
  filter(Set=="A" | Set=="B" | Set=="E" | Set=="F"| 
           Set=="I" | Set=="J") %>%
  ggplot(aes(x=factor(Pass), y=H, color=factor(varCX)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Shannon H'", x="Pass", title="C. gleum")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceH_CG

# Filter down to CG only and plot Simpson
pVariantEmergenceIS_CG<-NRRLMulti2 %>%
  filter(Set=="A" | Set=="B" | Set=="E" | Set=="F"| 
           Set=="I" | Set=="J") %>%
  ggplot(aes(x=factor(Pass), y=D_inv, color=factor(varCX)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Inverse Simpson", x="Pass", title="C. gleum")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceIS_CG

# how about MO?
pVariantEmergenceRichness_MO<-NRRLMulti2 %>%
  ggplot(aes(x=factor(Pass), y=S, color=factor(varMO)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Richness", x="Pass", title="Emergence of Morphological Variants: M. oxydans")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceRichness_MO

pVariantEmergenceH_MO<-NRRLMulti2 %>%
  ggplot(aes(x=factor(Pass), y=H, color=factor(varMO)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Shannon H'", x="Pass", title="M. oxydans")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceH_MO

pVariantEmergenceIS_MO<-NRRLMulti2 %>%
  ggplot(aes(x=factor(Pass), y=D_inv, color=factor(varMO)))+
  geom_jitter(width=0.15)+
  theme_classic()+
  #  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize), 
        axis.text.x=element_text(size=myTextSize), 
        axis.title=element_text(size=myTextSize),
        plot.title=element_text(hjust=0.5, size=myTextSize) 
  )+
  labs(y="Inverse Shannon", x="Pass", title="M. oxydans")+
  #facet_wrap(~Variant, ncol=1)+
  theme(strip.text=element_text(size=myTextSize, face="italic"))
pVariantEmergenceIS_MO

plot_grid(pVariantEmergenceRichness_CG, pVariantEmergenceRichness_MO,
          pVariantEmergenceH_CG, pVariantEmergenceH_MO, 
          pVariantEmergenceIS_CG, pVariantEmergenceIS_MO,
          ncol=2, align="h")
ggsave("pvariantEmergence_CG_MO_AllIndices.png", height=6, width=8, units="in", dpi=300)

# calculations - median diversity
summarizeVariantEmergenceD<-NRRLMulti2 %>%
  group_by(Pass) %>%
  summarize(median_H=median(H),
            median_D=median(D_inv),
            median_S=median(S))
summarizeVariantEmergenceD

# plot out
summarizeVariantEmergenceD %>%
  pivot_longer(!Pass, names_to="Measure", values_to="median_diversity") %>%
  ggplot(aes(x=factor(Pass), y=median_diversity, color=Measure)) +
  geom_point()+
  theme(legend.position = "none",
        axis.text.y=element_text(size=myTextSize+2), 
        axis.text.x=element_text(size=myTextSize+2), 
        axis.title=element_text(size=myTextSize+2),
        #plot.title=element_text(hjust=0.5, size=myTextSize)
        )+
  labs(y="Median Diversity", x="Pass", title="")
  





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         Pairwise in vitro Chryseobacterium 
# The first run (20220504) has day 21; this was dropped in later runs

NRRL2_CG_Morph_Compete<-read_excel("NRRL2CXMorphologyCompetition.xlsx", sheet="summary")
names(NRRL2_CG_Morph_Compete)
NRRL2_CG_Morph_Compete<-as_tibble(NRRL2_CG_Morph_Compete, stringsAsFactors="TRUE")

glimpse(NRRL2_CG_Morph_Compete)
View(NRRL2_CG_Morph_Compete)

# max proportion alternate
# when total counts sufficient
max(NRRL2_CG_Morph_Compete$fAlt[NRRL2_CG_Morph_Compete$cAlt>5]) #99.45%


# rename ancestor as ANC
NRRL2_CG_Morph_Compete$strain2<-replace(NRRL2_CG_Morph_Compete$strain2,
                                        NRRL2_CG_Morph_Compete$strain2=="CG0",
                                        "ANC") 

# plot out raw proportion fAlt data
# first reorder levels
unique(NRRL2_CG_Morph_Compete$strain1)
NRRL2_CG_Morph_Compete$strain1 <-factor(NRRL2_CG_Morph_Compete$strain1,
                                        levels=c("A2a6", "A2a10", "F2a6", "F2a10"))
unique(NRRL2_CG_Morph_Compete$strain2)
NRRL2_CG_Morph_Compete$strain2 <-factor(NRRL2_CG_Morph_Compete$strain2,
                                        levels=c("ANC", "A1o6", "A1o10", "A2o6", "A2o10",
                                                 "F1o6", "F1o10", "F2o6", "F2o10"))

# Take out the day that appears in only one set
NRRL2_CG_Morph_Compete<-NRRL2_CG_Morph_Compete  %>%
  filter(day!=21)

pNRRL2_MorphCompeteInVitro_fAlt<-NRRL2_CG_Morph_Compete %>%
  ggplot(aes(x=day, y=fAlt, color=strain1)) +
  geom_line(aes(linetype=factor(rep)))+
  theme_bw()+
  scale_colour_discrete()+
  theme(
    axis.text=element_text(size=xTextSize), 
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize),
    legend.title =element_text(size=xTextSize),
    plot.margin = margin(r=40)
    )  +
  labs(y="Fraction Alternate Morph", x="Day", title="In Vitro Competition",
       linetype="Isolate", color="Alt Morph")+
  facet_wrap(~strain2)
pNRRL2_MorphCompeteInVitro_fAlt


# Condensed plot of competition data
pNRRL2_CG_Morph_Compete_d14_Same_small<-NRRL2_CG_Morph_Compete_d14 %>%
  ggplot(aes(x=strain2, y=fAlt, color=factor(SamePair), pch=strain1)) +
  geom_jitter(height=0, width=0.1, size=2)+
  theme_bw()+
  #scale_color_viridis_d(option="inferno", begin=0.4, end=0.9, direction=-1)+
  scale_colour_discrete()+
  theme(
    axis.text.x = element_text(size=xTextSize, angle = 45, vjust = 1, hjust=1), 
    axis.text.y = element_text(size=xTextSize), 
    axis.title.y=element_text(size=xTextSize+2),
    #axis.text.x=element_blank(),
    #legend.title=element_blank(),
    #legend.position = "none",
    legend.title=element_text(size=xTextSize),
    strip.text = element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize),
    plot.margin = margin(r=40)
    )  +
  labs(y="fAlt", x="", 
       color="Matched", pch="Alt Morph",
       title="In Vitro Competition, Day 14")
  #facet_wrap(~strain2)
pNRRL2_CG_Morph_Compete_d14_Same_small









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Compare with growth parameter data. First we have to re-name "CG0" to "WT":
idx<-which(NRRL2_CG_Morph_Compete$strain2=="CG0")
NRRL2_CG_Morph_Compete$strain2[idx]<-"WT"
glimpse(NRRL2_CG_Morph_Compete)

#Now we can take the various deltas. Let's set aside the day 14 data from the pairwise comparisons. We'll add some columns to hold the deltas, then fill them using the strain1 and strain2 labels to find the right data in the other data frames.
NRRL2_CG_Morph_Compete_14<-NRRL2_CG_Morph_Compete %>%
  filter(day==14) %>%
  add_column(delta_r=NA, delta_CFU=NA)
glimpse(NRRL2_CG_Morph_Compete_14)

#Now to grab the deltas. First check out the objects - CFU
str(CG_CFU)

#and growth rates
str(allplates_manyspline_coef)

#Make the deltas (we'll do Alt-Ori).
my_dim<-dim(NRRL2_CG_Morph_Compete_14)
for (i in 1:my_dim[1]){
  strain_a<-NRRL2_CG_Morph_Compete_14$strain1[i]
  strain_o<-NRRL2_CG_Morph_Compete_14$strain2[i]
  NRRL2_CG_Morph_Compete_14$delta_r[i]<-median(allplates_manyspline_coef$mumax[allplates_manyspline_coef$Strain==strain_a &
                                                                                 allplates_manyspline_coef$Rep==NRRL2_CG_Morph_Compete_14$rep[i]]) - 
                                           median(allplates_manyspline_coef$mumax[allplates_manyspline_coef$Strain==strain_o &
                                                                                 allplates_manyspline_coef$Rep==NRRL2_CG_Morph_Compete_14$rep[i]])
  NRRL2_CG_Morph_Compete_14$delta_CFU<-median(CG_CFU$logCFU[CG_CFU$ID==strain_a &
                                              CG_CFU$Replicate==NRRL2_CG_Morph_Compete_14$rep[i]]) - 
                                           median(CG_CFU$logCFU[CG_CFU$ID==strain_o &
                                                  CG_CFU$Replicate==NRRL2_CG_Morph_Compete_14$rep[i]])
}

glimpse(NRRL2_CG_Morph_Compete_14)


#Linear fits. Growth rate is not significant
lm_deltagrowth<-lm(delta_r~fAlt, data=NRRL2_CG_Morph_Compete_14)
summary(lm_deltagrowth)

#Nor is delta log(CFU)
lm_deltaCFU<-lm(delta_CFU~fAlt, data=NRRL2_CG_Morph_Compete_14)
summary(lm_deltaCFU)


#Plot out:
pCompeteDeltaGrowth_byORI<-NRRL2_CG_Morph_Compete_14 %>%
  ggplot(aes(x=fAlt, y=delta_r, color=strain2))+
  geom_jitter(width=0.1)+
  theme_light()+
    theme(
    axis.text=element_text(size=xTextSize), 
    axis.title=element_text(size=xTextSize+2),
    #axis.text.x=element_blank(),
    legend.title=element_blank(),
    #legend.position = "none",
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="Change in Growth Rate", x="Fraction Alt Morph", title="")
pCompeteDeltaGrowth_byORI
```
 Color by alt morph?
```{r}
pCompeteDeltaGrowth_byALT<-NRRL2_CG_Morph_Compete_14 %>%
  ggplot(aes(x=fAlt, y=delta_r, color=strain1))+
  geom_jitter(width=0.1)+
  theme_light()+
    theme(
    axis.text=element_text(size=xTextSize), 
    axis.title=element_text(size=xTextSize+2),
    #axis.text.x=element_blank(),
    legend.title=element_blank(),
    #legend.position = "none",
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="Change in Growth Rate", x="Fraction Alt Morph", title="")
pCompeteDeltaGrowth_byALT
```
 

And plot vs CFU?
```{r}
pCompeteDeltaCFU_byORI<-NRRL2_CG_Morph_Compete_14 %>%
  ggplot(aes(x=fAlt, y=delta_CFU, color=strain2))+
  geom_jitter(width=0.1)+
  theme_light()+
    theme(
    axis.text=element_text(size=xTextSize), 
    axis.title=element_text(size=xTextSize+2),
    #axis.text.x=element_blank(),
    legend.title=element_blank(),
    #legend.position = "none",
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="Change in log(Max CFU)", x="Fraction Alt Morph", title="")
pCompeteDeltaCFU_byORI
```
again, color by alt morph?
```{r}
pCompeteDeltaCFU_byALT<-NRRL2_CG_Morph_Compete_14 %>%
  ggplot(aes(x=fAlt, y=delta_CFU, color=strain1))+
  geom_jitter(width=0.1)+
  theme_light()+
    theme(
    axis.text=element_text(size=xTextSize), 
    axis.title=element_text(size=xTextSize+2),
    #axis.text.x=element_blank(),
    legend.title=element_blank(),
    #legend.position = "none",
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="Change in log(Max CFU)", x="Fraction Alt Morph", title="")
pCompeteDeltaCFU_byALT
```
And together:
```{r}
plot_grid(pCompeteDeltaGrowth_byORI, pCompeteDeltaCFU_byORI,
          pCompeteDeltaGrowth_byALT, pCompeteDeltaCFU_byALT,
          ncol=2, labels="AUTO")
ggsave("pCompeteDeltas.png", width=10, height=8, units="in", dpi=300)
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Pairwise competition with worms on plates
CG_Compete_PlusWorm<-read_excel("NRRL2CXMorphologyCompetition.xlsx", sheet="20230731PlusWorm")
names(CG_Compete_PlusWorm)
CG_Compete_PlusWorm<-as_tibble(CG_Compete_PlusWorm, stringsAsFactors="TRUE")

glimpse(CG_Compete_PlusWorm)
View(CG_Compete_PlusWorm)

# max proportion alternate
max(CG_Compete_PlusWorm$fAlt, na.rm=TRUE) #82.6%


# plot out raw proportion fAlt data
# first reorder levels
unique(CG_Compete_PlusWorm$alt)
CG_Compete_PlusWorm$alt <-factor(CG_Compete_PlusWorm$alt,
                                        levels=c("A2a10", "F2a10"))
unique(CG_Compete_PlusWorm$ori)
CG_Compete_PlusWorm$ori <-factor(CG_Compete_PlusWorm$ori,
                                        levels=c("ANC", "A1o10", "A2o10",
                                                 "F1o10", "F2o10"))
# plot
CG_Compete_PlusWorm %>%
  ggplot(aes(x=day, y=fAlt, color=alt)) +
  geom_line(aes(linetype=factor(isolate)))+
  theme_bw()+
  scale_colour_discrete()+
  theme(
    axis.text=element_text(size=myTextSize), 
    axis.title=element_text(size=myTextSize+2),
    strip.text=element_text(size=myTextSize),
    plot.title=element_text(hjust=0.5, size=myTextSize+2, face="bold"), 
    legend.text=element_text(size=myTextSize))  +
  labs(y="Fraction Alternate Morph", x="Day", title="+Worms",
       linetype="Isolate", color="Alt Morph")+
  facet_grid(vars(ori), vars(dAlt)) # by original morph and starting alt dilution (1:10 or 1:1)

# maybe a difference map?
temp1<-CG_Compete_PlusWorm %>%
  dplyr::select(day, isolate, ori, alt, condition, dAlt, fAlt) %>%
  pivot_wider(names_from = day, values_from=fAlt, 
              names_prefix="day") %>%
  #mutate(delta7=day7-day0, delta14=day14-day7) %>% # get the deltas between time points
  #dplyr::select(isolate, condition, dAlt, day0, day14, delta7, delta14) %>%
  dplyr::select(isolate, ori, alt, condition, dAlt, day0, day7) %>%
  rename(f1=day0, f2=day7)
temp2<-CG_Compete_PlusWorm %>%
  dplyr::select(day, isolate, ori, alt, condition, dAlt, fAlt) %>%
  pivot_wider(names_from = day, values_from=fAlt, 
              names_prefix="day") %>%
  #mutate(delta7=day7-day0, delta14=day14-day7) %>% # get the deltas between time points
  #dplyr::select(isolate, condition, dAlt, day0, day14, delta7, delta14) %>%
  dplyr::select(isolate, ori, alt, condition, dAlt, day7, day14) %>%
  rename(f1=day7, f2=day14)

CG_Compete_PlusWorm_Diff<-rbind(temp1, temp2)
rm(temp1, temp2)

CG_Compete_PlusWorm_Diff$Community<-substr(CG_Compete_PlusWorm_Diff$alt, 1,1)

pCG_Compete_PlusWorm_Diff<-CG_Compete_PlusWorm_Diff %>%
  ggplot(aes(x=f1, y=f2, color=ori))+
  geom_point(size=2)+
  theme_bw()+
  geom_abline(slope=1, intercept=0)+
  #scale_colour_viridis_d(begin=0.92, end=0)+
  scale_color_brewer(palette="Dark2")+
  theme(
    axis.text.x=element_text(size=xTextSize, angle = 45, vjust = 1, hjust=1),
    axis.text.y=element_text(size=xTextSize),
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="Final fAlt", x="Initial fAlt",
       title="+Worms", color="")+
  facet_wrap(~alt)
pCG_Compete_PlusWorm_Diff

#~~~~~~~~~~~~~~
# FAlt data from chunks cut out of plates passaged in quadruplicate
# from day 14 of the experiment above,
# Diameters taken at 3 days,
# CFUs taken after 7 further days

CG_Compete_WormDiameter<-read_excel("NRRL2CXMorphologyCompetition.xlsx", sheet="20230817WormDiameter")
names(CG_Compete_WormDiameter)
CG_Compete_WormDiameter<-as_tibble(CG_Compete_WormDiameter, stringsAsFactors="TRUE")
glimpse(CG_Compete_WormDiameter)

# reorder
unique(CG_Compete_WormDiameter$Ori)
CG_Compete_WormDiameter$Ori <-factor(CG_Compete_WormDiameter$Ori,
                                  levels=c("ANC", "A1o10", "A2o10",
                                           "F1o10", "F2o10"))


pCG_Compete_WormDiameter<-CG_Compete_WormDiameter %>%
  ggplot(aes(x=Ori, y=Diameter, color=Alt))+
  geom_jitter(width=0.2, size=2)+
  theme_bw()+
  #scale_colour_viridis_d(begin=0.92, end=0)+
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("chartreuse3", "darkorchid2")) +
  theme(
    axis.text.x=element_text(size=xTextSize-2, angle = 45, vjust = 1, hjust=1),
    axis.text.y=element_text(size=xTextSize),
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize),
    #legend.position = "top"
    #legend.position=c(0.8, 0.8)
    )  +
  labs(title="+Worms", y="Diameter (cm)", x="", color="")
pCG_Compete_WormDiameter

# and the CFU data
CG_Compete_WormChunk<-read_excel("NRRL2CXMorphologyCompetition.xlsx", sheet="20230814WormChunk")
names(CG_Compete_WormChunk)
CG_Compete_WormChunk<-as_tibble(CG_Compete_WormChunk, stringsAsFactors="TRUE")
glimpse(CG_Compete_WormChunk)

# first reorder levels
unique(CG_Compete_WormChunk$Ori)
CG_Compete_WormChunk$Ori <-factor(CG_Compete_WormChunk$Ori,
                                        levels=c("ANC", "A1o10", "A2o10",
                                                 "F1o10", "F2o10"))

pCG_Compete_WormChunk<-CG_Compete_WormChunk %>%
  ggplot(aes(x=Chunk, y=fAlt, pch=Chunk)) +
  geom_jitter(width=0.1, size=2, aes(color=Alt))+
  theme_bw()+
  ylim(-0.1, 1.1)+
  #scale_colour_discrete()+
  scale_color_manual(values=c("chartreuse3", "darkorchid2")) +
  theme(
    axis.text.x=element_text(size=xTextSize, angle = 45, vjust = 1, hjust=1),
    #axis.text.x=element_blank(),
    axis.text.y=element_text(size=xTextSize),
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y="fAlt", x="", title="+Worms",
       color="", pch="")+
  facet_wrap(~Ori, nrow=1)+
  stat_compare_means(method="wilcox.test", 
                     label.y=1.05,
                     aes(label = paste0("p =", ..p.format..)))
pCG_Compete_WormChunk


#~~~~~~~~~~~~~~~~~~~~~~~
# Are day-3 expansion diameter and over-representation of Alt at plate edge correlated?
# First, let's get over-representation quantified

CG_Compete_WormChunk_OR<-CG_Compete_WormChunk %>%  # pivot out so fAlt at center and edge are accessible
  dplyr::select(Ori, Alt, Well, Isolate, techrep, Chunk, fAlt) %>%
  pivot_wider(names_from = Chunk, values_from = fAlt)

CG_Compete_WormChunk_OR<-CG_Compete_WormChunk_OR %>% # get the difference (OR)
  mutate(delta_fAlt=Edge-Center, # as absolute
         delta_fAlt_normalized=(Edge-Center)/Edge) #and relative

# quick check
CG_Compete_WormChunk_OR %>%
  ggplot(aes(x=Ori, y=delta_fAlt, color=Alt))+
  geom_jitter(width=0.2)
CG_Compete_WormChunk_OR %>% # yikes
  ggplot(aes(x=Ori, y=delta_fAlt_normalized, color=Alt))+
  geom_jitter(width=0.2)

# now the expansion diameters at day 3
# The data frames are already in the same order, so this is easy
CG_Compete_WormChunk_OR$Diameter<-CG_Compete_WormDiameter$Diameter

# quick check
CG_Compete_WormChunk_OR %>%
  ggplot(aes(x=Diameter, y=delta_fAlt, color=Ori))+
  geom_jitter(width=0.2)

# clearly this is complex. Separate by pair type
CG_Compete_WormChunk_OR$PairType<-"Matched"
idx<-which(CG_Compete_WormChunk_OR$Ori=="A1o10" |
             CG_Compete_WormChunk_OR$Ori=="F1o10")
CG_Compete_WormChunk_OR$PairType[idx]<-"Non-Matched"
idx<-which(CG_Compete_WormChunk_OR$Ori=="ANC")
CG_Compete_WormChunk_OR$PairType[idx]<-"Ancestor"

pCG_Compete_WormChunk_OR_slopes<-CG_Compete_WormChunk_OR %>%
  ggplot(aes(x=Diameter, y=delta_fAlt, color=Alt))+
  geom_jitter(width=0.2)+
  ylim(-0.5, 1)+
  theme_bw()+
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("chartreuse3", "darkorchid2")) +
  geom_smooth(method="lm", formula = y ~ x, fill="NA")+
  theme(
    axis.text.x=element_text(size=xTextSize, vjust = 1, hjust=1),
    #axis.text.x=element_blank(),
    axis.text.y=element_text(size=xTextSize),
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y=expression(fAlt[Edge]-fAlt[Center]), x="Diameter (cm)", 
       title="",
       color="")+
  facet_wrap(~PairType)
pCG_Compete_WormChunk_OR_slopes
#ggsave("pCG_Compete_WormChunk_OR_slopes.png", width=8, height=3, units="in", dpi=400)

pCG_Compete_WormChunk_fEdge_slopes<-CG_Compete_WormChunk_OR %>%
  ggplot(aes(x=Diameter, y=Edge, color=Alt))+
  geom_jitter(width=0.2)+
  ylim(-0.1, 1.1)+
  theme_bw()+
  scale_color_manual(values=c("chartreuse3", "darkorchid2")) +
  #scale_color_brewer(palette="Dark2")+
  geom_smooth(method="lm", formula = y ~ x, fill="NA")+
  theme(
    axis.text.x=element_text(size=xTextSize, vjust = 1, hjust=1),
    #axis.text.x=element_blank(),
    axis.text.y=element_text(size=xTextSize),
    axis.title=element_text(size=xTextSize+2),
    strip.text=element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5, size=xTextSize+2, face="bold"), 
    legend.text=element_text(size=xTextSize))  +
  labs(y=expression(fAlt[Edge]), x="Diameter (cm)", 
       title="",
       color="")+
  facet_wrap(~PairType)
pCG_Compete_WormChunk_fEdge_slopes
#ggsave("pCG_Compete_WormChunk_fEdge_slopes.png", width=8, height=3, units="in", dpi=400)

# together
plot_grid(pCG_Compete_WormChunk_fEdge_slopes,
          pCG_Compete_WormChunk_OR_slopes + ylim(-0.2, 1.1),
          ncol=1, labels="AUTO")
ggsave("pCG_Compete_WormChunk_slopes.png", width=8, height=6, units="in", dpi=400)


#### FIGURE 2

# version with the large grid of the 14-day in vitro data
#pCG_Compete_AB<-plot_grid(pNRRL2_MorphCompeteInVitro_fAlt,
#          pNRRL2_CG_Morph_Compete_d14_Same,
#          nrow=1,
#          #rel_widths = c(1.2, 1), 
#          labels="AUTO")
#pCG_Compete_CD<-plot_grid(pCG_Compete_PlusWorm_Diff, 
#                          pCG_Compete_WormDiameter,
#                          pCG_Compete_WormChunk,
#                          nrow=1,
#                          rel_widths=c(0.6, 0.4, 1), 
#                          labels=c("C", "D", "E"))
#pCG_Compete<-plot_grid(pCG_Compete_AB,
#                       pCG_Compete_CD,
#                       ncol=1,
#                       rel_heights = c(1.5, 1),
#                       labels="none")

# Version with the condensed in vitro day 14 data

pCG_Compete_AB<-plot_grid(pNRRL2_MorphCompeteInVitro_fAlt,
                          pNRRL2_CG_Morph_Compete_d14_Same_small,
                          ncol=1,
                          rel_heights = c(1.8, 1),
                          labels="AUTO")
pCG_Compete_CD<-plot_grid(pCG_Compete_PlusWorm_Diff, 
                          pCG_Compete_WormDiameter,
                          nrow=1,
                          rel_widths=c(1, 0.7), 
                          labels=c("C", "D"))
pCG_Compete_CDE<-plot_grid(pCG_Compete_CD,
                           pCG_Compete_WormChunk,
                           ncol=1,
                           labels=c("", "E"))

pCG_Compete<-plot_grid(pCG_Compete_AB,
                       pCG_Compete_CDE,
                       ncol=2,
                       #rel_heights = c(1.8, 0.8, 0.9, 1),
                       labels=c("", ""))

pCG_Compete
ggsave("Fig2_CG_Compete.png", width=14, height=7, units="in", dpi=400)

