library(tidyverse)
library(readxl)
library(ggExtra)
library(scales)
library(viridis)
library(RColorBrewer)


data <- read_excel("Supplementary_Table.xlsx")



# PANEL A ######################################################################
data %>% 
  ggplot(aes(x=log10(Cont_N50), 
             y=log10(Scaff_N50),
             size=Asm_Size,
             shape = LongReads_tec,
             color=Status,
             fill = Status,
             label=ToLID)) + 
  geom_point(aes(fill = Status), alpha=0.8, stroke = 1.1) +
  scale_shape_manual(values = c(1, 21)) +
  scale_size_continuous(range = c(4, 16),breaks = c(1e+08, 5e+08, 1e+09, 2e+09, 3e+09), name = "Asm Size (Gbp)") +
  scale_color_manual(values = c("Curated"= "#e01b24cc", "Pre-curation" = "#1c71d8cc", "Non-final" = "#2ec27ecc"),
                     breaks = c("Curated", "Pre-curation", "Non-final")) +
  scale_fill_manual(values = c("Curated"= "#e01b24cc", "Pre-curation" = "#1c71d8cc", "Non-final" = "#2ec27ecc"),
                    breaks = c("Curated", "Pre-curation", "Non-final")) +
  ylab(expression(log[10]*"(Scaffold N50)")) +
  xlab(expression(log[10]*"(Contig N50)")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) +
  geom_text(data = subset(data, log10(Scaff_N50) < 7 | log10(Cont_N50) < 6), 
            aes(x = log10(Cont_N50), y = log10(Scaff_N50), label = ToLID),
            size = 3.75, show.legend = FALSE, nudge_y = 0.02, vjust = -2.5, color="black") +
  geom_hline(yintercept = 7, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 6, linetype = "dashed", color = "black")



  
# PANEL B ######################################################################
data2 <- data %>% filter(!grepl("Non-final", Status))

#ONT
data2 %>% 
  ggplot(aes(x=log10(Gaps_perGbp+1), 
             y=QV,
             size=Asm_Size,
             color= Completeness,
             label=ToLID)) + 
  geom_point(alpha=0.8) +
  scale_size_continuous(range = c(4, 16),breaks = c(1e+08, 5e+08, 1e+09, 2e+09, 3e+09), name = "Asm Size (Gbp)") +
  scale_color_gradient(low = "#9B111E", high = "#0F52BA")+
  ylab("QV") +
  xlab("log10(Gaps/Gbp)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black") +
  geom_text(data = subset(data2, QV < 40 | log10(Gaps_perGbp+1) > 3 | Completeness < 90), 
            aes(x = log10(Gaps_perGbp+1), y = QV, label = ToLID),
            size = 3.75, show.legend = FALSE, nudge_y = 0.02, vjust = -2.5, color="black")

#HiFi: (low = "#F4A460", high = "#228B22")




# PANEL C ######################################################################

data_lolli <- data %>% filter(!grepl("Non-final", Status))
data_lolli <- data_lolli %>% arrange(Taxa_1)

data_lolli <- data_lolli %>% mutate(BUSCOs_ancient=BUSCOs_ancient*100)
data_lolli <- data_lolli %>% mutate(BUSCOs_euka=BUSCOs_euka*100)
data_lolli <- data_lolli %>% mutate(BUSCOd_ancient=BUSCOd_ancient*100)
data_lolli <- data_lolli %>% mutate(BUSCOd_euka=BUSCOd_euka*100)
data_lolli <- data_lolli %>% mutate(BUSCOc_ancient=BUSCOs_ancient+BUSCOd_ancient)
data_lolli <- data_lolli %>% mutate(BUSCOc_euka=BUSCOs_euka+BUSCOd_euka)

new_data <- data_lolli %>% select(ToLID, Taxa_1, BUSCOc_ancient, BUSCOc_euka)

new_data <- new_data %>%
  rowwise() %>% 
  mutate( mymean = mean(c(BUSCOc_ancient,BUSCOc_euka) )) %>% 
  mutate(ToLID=factor(ToLID, ToLID))



ggplot(new_data) +
  geom_segment(aes(x = ToLID, xend = ToLID, y = BUSCOc_ancient, yend = BUSCOc_euka), color = "grey") +
  geom_point(aes(x = ToLID, y = BUSCOc_euka), color = "#2e3436", size = 3) +
  geom_point(aes(x = ToLID, y = BUSCOc_ancient, color = Taxa_1), size = 3) +
  coord_flip() +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.grid.major.x = element_line(color = "grey", linewidth = 0.1, linetype = "dashed")) +
  theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.2, linetype = "dotted")) +
  scale_color_manual(values = c("Viridiplantae" = "#4e9a06", "Vertebrata" = "#3465a4", "Phaeophyceae" = "#f57900", "Nematomorpha" = "#c17d11", "Porifera" = "#c17d11", "Insecta" = "#cc0000", "Fungi" = "#edd400")) +
  labs(x = "Assembly", y = "BUSCO complete") +
  theme(axis.text.x = element_text(size=14), axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size=12), axis.title = element_text(size = 15)) +
  theme(legend.position = "none")

  


# SUPPL ########################################################################

data_bar <- data2[-6,]

data_bar <- data_bar %>% select(ToLID,Taxa_1,HiC_tec,PCR_Dup,NoDup) %>%
  filter(!grepl("None", HiC_tec)) %>%
  mutate(PCR_Dup=as.numeric(PCR_Dup)*100) %>%
  mutate(NoDup=as.numeric(NoDup)*100)


# Order data by BUSCOdb_ancient
ordered_data <- data_bar %>%
  arrange(Taxa_1) %>%
  mutate(ToLID = factor(ToLID, levels = ToLID))

# Reshape the data to a long format
long_data <- tidyr::gather(ordered_data, key = "variable", value = "value", PCR_Dup, NoDup)


ggplot(long_data, aes(x = ToLID, y = value, fill = interaction(variable, HiC_tec))) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("ToLID ordered by Taxa") +
  ylab("% mapped reads") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "grey", linewidth = 0.1, linetype = "dashed")
  ) +
  scale_fill_manual(values = c("PCR_Dup.OmniC" = "#ce5c00", "NoDup.OmniC" = "#fcaf3e",
                               "PCR_Dup.Arima" = "#a40000", "NoDup.Arima"= "#ef2929"))

