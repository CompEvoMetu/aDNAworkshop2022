if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
​
if(!require(tidyverse)){
    install.packages("tidyverse")
    library(tidyverse)
}
​
#if(!require(ggpubr)){
#  install.packages("ggpubr")
#  library(ggpubr)
#}
​
# Loading READ results ----------------------------------------------------
p0_vals <- read.table("./meansP0_AncientDNA_normalized", header = T, sep = " ")
snpCount <- read.table("./SNPCount.READ_output_ordered", header = T, sep = "\t")
​
p0_vals <- merge(p0_vals, snpCount, by = "PairIndividuals", all.x = T)
​
str(p0_vals)
​
# NORMALIZATION: MEDIAN ---------------------------------------------------
p0_vals <- p0_vals %>%
  filter(SNPcount >= 5000) %>%
  mutate(median_p0 = median(NonNormalizedP0)) %>%
  separate(col = PairIndividuals, into = c("Pair1", "Pair2"), sep = "\\|", remove = F)
​
# Constants for kinship relatedness
SecondUpCons <- 0.90625
FirstUpCons <- 0.8125
IdentUpCons <- 0.625
​
#Multiplying medians with kinship constants for estimation of kinship status
p0_vals$SecondUp <- p0_vals$median_p0 * SecondUpCons
p0_vals$FirstUp <- p0_vals$median_p0 * FirstUpCons
p0_vals$IdentUp <- p0_vals$median_p0 * IdentUpCons
​
#Relationship treshold determination
p0_vals <- p0_vals %>%
  mutate(sd = sd(NonNormalizedP0)) %>%
  mutate(n = n()) %>%
  mutate(standart_error = (qnorm(0.975)*sd/sqrt(n))) %>%
  mutate(SecondUp_max = (median_p0*SecondUpCons) + standart_error,
         SecondUp_min = (median_p0*SecondUpCons) - standart_error,
         FirstUp_max = (median_p0*FirstUpCons) + standart_error,
         FirstUp_min = (median_p0*FirstUpCons) - standart_error,
         IdentUp_max = (median_p0*IdentUpCons) + standart_error,
         IdentUp_min = (median_p0*IdentUpCons) - standart_error) %>%
  mutate(degree =
          ifelse(NonNormalizedP0 > SecondUp_max, "Unrelated",
           ifelse(NonNormalizedP0 < SecondUp_max & NonNormalizedP0 > SecondUp_min, "NA",
            ifelse(NonNormalizedP0 < SecondUp_min & NonNormalizedP0 > FirstUp_max, "SecondDegree",
             ifelse(NonNormalizedP0 < FirstUp_max & NonNormalizedP0 > FirstUp_min, "FD_or_SD",
              ifelse(NonNormalizedP0 < FirstUp_min & NonNormalizedP0 > IdentUp_max, "FirstDegree",
               ifelse(NonNormalizedP0 < IdentUp_max & NonNormalizedP0 > IdentUp_min, "FD_or_IT",
                ifelse(NonNormalizedP0 < IdentUp_min, "IdenticalTwin", "error")))))))) %>%
  filter(SecondUp_max > SecondUp_min & SecondUp_min > FirstUp_max & FirstUp_max > FirstUp_min &
           FirstUp_min > IdentUp_max & IdentUp_max > IdentUp_min) %>%
  arrange(NonNormalizedP0)
​
p0_vals_brief <- p0_vals %>%
  select(PairIndividuals, degree, SNPcount, NonNormalizedP0, median_p0)
​
#Plotting
p_median <- p0_vals %>%
  ggplot(aes(reorder(PairIndividuals,NonNormalizedP0),NonNormalizedP0))+
  annotate('ribbon', x = c(-Inf, Inf), 
           ymin = unique(p0_vals$IdentUp_min), ymax = unique(p0_vals$IdentUp_max), 
           alpha = 0.4, fill = 'grey')+
  annotate('ribbon', x = c(-Inf, Inf), 
           ymin = unique(p0_vals$FirstUp_min), ymax = unique(p0_vals$FirstUp_max), 
           alpha = 0.4, fill = 'grey')+
  annotate('ribbon', x = c(-Inf, Inf), 
           ymin = unique(p0_vals$SecondUp_min), ymax = unique(p0_vals$SecondUp_max), 
           alpha = 0.4, fill = 'grey')+
  geom_point() +
  geom_errorbar(aes(ymin = as.numeric(as.character(NonNormalizedP0)) - 
                      as.numeric(as.character(NonNormalizedStandardError)),
                    ymax = as.numeric(as.character(NonNormalizedP0)) + 
                      as.numeric(as.character(NonNormalizedStandardError)), width = 0.01))+
  geom_hline(yintercept = p0_vals$median_p0*SecondUpCons, color = "black") +
  geom_hline(yintercept = p0_vals$median_p0*FirstUpCons, linetype = "dashed", color = "black") +
  geom_hline(yintercept = p0_vals$median_p0*IdentUpCons, linetype = "dotted", color = "black") +
  
  #Identical Twin Area
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1, y = 0, yend = unique(IdentUp_min) - 0.001),
               color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = unique(IdentUp_min) - 0.001, yend = unique(IdentUp_min)- 0.001), color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = 0, yend = 0), color = "grey70")+
  annotate("text", x = nrow(p0_vals)-5 , y = unique(p0_vals$IdentUp_min)/2, 
           label = "Identical Twin\nArea", angle = 0, fontface = 3, size = 5 , alpha = .5)+
  #First Degree Area
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1, 
                   y = unique(IdentUp_max)+ 0.001, yend = unique(FirstUp_min) - 0.001),color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = unique(IdentUp_max)+ 0.001, yend = unique(IdentUp_max)+ 0.001), color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = unique(FirstUp_min) - 0.001, yend = unique(FirstUp_min) - 0.001), color = "grey70")+
  annotate("text", x = nrow(p0_vals)-5 , y = (unique(p0_vals$FirstUp_min) + unique(p0_vals$IdentUp_max))/2, 
           label = "First Degree Area", angle = 0, fontface = 3, size = 5 , alpha = .5)+
  #Second Degree Area
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1, 
                   y = unique(FirstUp_max)+ 0.001, yend = unique(SecondUp_min) - 0.001),color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = unique(FirstUp_max)+ 0.001, yend = unique(FirstUp_max)+ 0.001), color = "grey70")+
  geom_segment(aes(x = nrow(p0_vals)-1, xend = nrow(p0_vals)-1.5, 
                   y = unique(SecondUp_min) - 0.001, yend = unique(SecondUp_min) - 0.001), color = "grey70")+
  annotate("text", x = nrow(p0_vals)-5 , y = (unique(p0_vals$FirstUp_max) + unique(p0_vals$SecondUp_min))/2, 
           label = "Second Degree Area", angle = 0, fontface = 3, size = 5 , alpha = .5)+
  ylab("Pairwise mismatch value (P0)")+
  ylim(c(0, .25))+
  labs(title = "READ Results - Normalization method: median")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12));p_median
​
ggsave(p_median, filename = "./READ_plot_snpfilter_median.pdf", device = "pdf", 
       height = 7, width = 15, units = "in")
write.table(p0_vals, "./READ_results_snpfilter_median", row.names = F, sep = "\t", quote = F)
write.table(p0_vals_brief, "./READ_results_snpfilter_median_brief", row.names = F, sep = "\t", quote = F)
​
# NORMALIZATION: MEAN -----------------------------------------------------
​
# Loading READ results ----------------------------------------------------
p0_vals <- read.table("./meansP0_AncientDNA_normalized"...
