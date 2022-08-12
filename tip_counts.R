library(tidyverse)
setwd("~/Desktop/caitlin")
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

ogs = c("OG5_142622", "OG5_131553", "OG5_126569", "OG5_127274", "OG5_135474", "OG5_126573", "OG5_128667", "OG5_128568")
tips <- read_csv("tip_counts.csv") %>% 
  mutate(taxa = X1) %>%
  select(-X1) %>%
  arrange(taxa)

tips <- tips %>%
  pivot_longer(cols=colnames(tips)[startsWith(colnames(tips), "OG5")], names_to="OG5", values_to="count") %>%
  select(-colnames(tips)[startsWith(colnames(tips), "X")]) %>%
  replace_na(replace=list(count = 0)) %>%
  separate(OG5, into=c("OG5", "dup"), sep=10, fill="left") %>%
  select(-dup) %>%
  group_by(taxa, OG5) %>%
  summarize(count = mean(count)) %>%
  separate(taxa, into=c('major_clade', 'minor_clade', 'species'), sep="_")



counts_by_major = tips %>%
  group_by(OG5, major_clade) %>%
  summarize(total_taxa = n_distinct(species),
            n_taxa = n_distinct(species[count>0]),
            n_tips = sum(count)) %>%
  mutate(major_clade = factor(major_clade, levels=c("Ba", "Za", "Op",
                                                    "Am", "Ex", "EE",
                                                    "Pl", "Sr"))) 

taxon_counts = counts_by_major %>%
  group_by(major_clade) %>%
  summarize(total_taxa=mean(total_taxa)) %>%
  .$total_taxa

gene_list = read_csv("ogs_genes.csv")[,c(1,14,19)]
colnames(gene_list) <- c("gene_family", "OG5", "source")
gene_list <- gene_list %>%
  drop_na() %>%
  separate(OG5, into=as.character(rep(c(1:11), each=1)), sep=", ") %>%
  filter(2 != "NONE") %>%
  pivot_longer(cols=c(2:12), names_to="drop", values_to="OG5") %>%
  drop_na() 

gene_list <- gene_list %>%
  mutate(OG5 = gsub(",","",gene_list$OG5)) %>%
  filter(OG5 != "NONE") %>%
  select(-drop)

counts_by_major <- counts_by_major %>%
  left_join(gene_list, by = "OG5") %>%
  group_by(OG5, major_clade) %>%
  summarize(
            total_taxa = mean(total_taxa),
            n_taxa = mean(n_taxa),
            n_tips = mean(n_tips),
            gene_family = paste(Modes(gene_family), collapse=", "),
            source = paste(Modes(source), collapse=", "),
            paralogness = if_else(n_taxa > 0, round(n_tips/n_taxa, 1), 0))
counts_by_major$gene_family[counts_by_major$OG5=="OG5_126569"] <- "H3"
counts_by_major$gene_family[counts_by_major$OG5=="OG5_126573"] <- "H4"
counts_by_major$source[counts_by_major$OG5=="OG5_126569"] <- "+"
counts_by_major$source[counts_by_major$OG5=="OG5_126573"] <- "+"

og_list = counts_by_major %>%
  group_by(OG5) %>%
  summarize(n()) %>%
  .$OG5

gene_family_list = counts_by_major %>%
  arrange(OG5) %>%
  group_by(OG5) %>%
  summarize(gene_family = paste(Modes(gene_family), collapse=", ")) %>%
  .$gene_family

source_list = counts_by_major %>%
  arrange(OG5) %>%
  group_by(OG5) %>%
  summarize(source = paste(Modes(source), collapse=", ")) %>%
  .$source

colors = c("#FFFFFF", RColorBrewer::brewer.pal(9, "YlOrRd"))
#values = c(0,5,10,20,35,50,75,100,200,300)
values = c(0,1,2,3,4,5,6,7,8,9)
ggplot(data = counts_by_major, aes(x = OG5, y = major_clade, fill=paralogness)) +
  geom_tile() +
  geom_text(aes(label=round(paralogness,0)), size=2.5) +
  theme_bw() +
  scale_fill_gradientn(colors=colors, values=values/9) +
  scale_x_discrete(labels=paste0(og_list, ": ", gene_family_list, " (", source_list, ")")) +
  scale_y_discrete(labels=paste0(levels(counts_by_major$major_clade), " (N=",taxon_counts, ")")) +
  labs(x = "Gene Family", y = "Major Clade Code", fill = "Avg. Paralog-ness") +
  theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1, size=9, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        legend.text = element_text(size=10, color="black"),
        legend.title= element_text(size=14, color="black"),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("tips_heatmap_paralog.png", plot=last_plot(), device="png", dpi=600,
       width=55, height=30, units='cm')


colors = c("#FFFFFF", RColorBrewer::brewer.pal(9, "GnBu"))
values = c(0,5,10,15,20,25,30,35,40,45,50)
ggplot(data = counts_by_major, aes(x = OG5, y = major_clade, fill=(n_taxa/total_taxa)*100)) +
  geom_tile() +
  theme_bw() +
  geom_text(aes(label=n_taxa), size=2.5) +
  theme_bw() +
  scale_fill_gradientn(colors=colors, values=values/50) +
  scale_x_discrete(labels=paste0(og_list, ": ", gene_family_list, " (", source_list, ")")) +
  scale_y_discrete(labels=paste0(levels(counts_by_major$major_clade), " (N=",taxon_counts, ")")) +
  labs(x = "Gene Family", y = "Major Clade Code", fill = "Proportion of Taxa (%)") +
  theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1, size=9, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        legend.text = element_text(size=10, color="black"),
        legend.title= element_text(size=14, color="black"),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("taxa_heatmap_proportional.png", plot=last_plot(), device="png", dpi=600,
       width=55, height=30, units='cm')





# preguidance tips
all_preguidance <- data.frame()
for (i in c(1:length(og_list))) {
  if (file.exists(paste0("/Users/caitlintimmons/Desktop/Katzlab/Meiosis_Genes/preguidance_tips/", og_list[i], "_preguidancetips.csv"))) {
    tips <- read.delim(paste0("/Users/caitlintimmons/Desktop/Katzlab/Meiosis_Genes/preguidance_tips/", og_list[i], "_preguidancetips.csv"), sep=",", header=TRUE) %>%
      select(2)
    colnames(tips) <- "taxon"
    tips <- tips %>%
      separate(taxon, into=c("drop", "minor_clade", "taxon"), sep="_") %>%
      select(-drop) %>%
      right_join(all_taxa, by = "taxon") %>%
      mutate("OG5" = og_list[i])
    
    all_preguidance <- rbind(all_preguidance, tips)
  }
}

all_taxa <- read.delim("/Users/caitlintimmons/Desktop/Katzlab/Meiosis_Genes/outgroup_taxa_Meiosis", header=FALSE, col.names="code") %>%
  separate(code, into=c("major_clade", "minor_clade", "taxon")) %>%
  select(-minor_clade)

all_preguidance <- all_preguidance %>%
  mutate(count = if_else(is.na(minor_clade), 0, 1)) %>%
  group_by(OG5, major_clade, minor_clade, taxon) %>%
  summarize(n_tips = sum(count))

counts_by_major = all_preguidance %>%
  group_by(OG5, major_clade) %>%
  summarize(total_taxa = n_distinct(taxon),
            n_taxa = n_distinct(taxon[n_tips>0]),
            n_tips = sum(n_tips)) %>%
  mutate(major_clade = factor(major_clade, levels=c("Ba", "Za", "Op",
                                                    "Am", "Ex", "EE",
                                                    "Pl", "Sr"))) 

taxon_counts = counts_by_major %>%
  group_by(major_clade) %>%
  summarize(total_taxa=mean(total_taxa)) %>%
  .$total_taxa

counts_by_major <- counts_by_major %>%
  left_join(gene_list, by = "OG5") %>%
  group_by(OG5, major_clade) %>%
  summarize(
    total_taxa = mean(total_taxa),
    n_taxa = mean(n_taxa),
    n_tips = mean(n_tips),
    gene_family = paste(Modes(gene_family), collapse=", "),
    source = paste(Modes(source), collapse=", "))
counts_by_major$gene_family[counts_by_major$OG5=="OG5_126569"] <- "H3"
counts_by_major$gene_family[counts_by_major$OG5=="OG5_126573"] <- "H4"
counts_by_major$source[counts_by_major$OG5=="OG5_126569"] <- "+"
counts_by_major$source[counts_by_major$OG5=="OG5_126573"] <- "+"

og_list = counts_by_major %>%
  group_by(OG5) %>%
  summarize(n()) %>%
  .$OG5

gene_family_list = counts_by_major %>%
  arrange(OG5) %>%
  group_by(OG5) %>%
  summarize(gene_family = paste(Modes(gene_family), collapse=", ")) %>%
  .$gene_family

source_list = counts_by_major %>%
  arrange(OG5) %>%
  group_by(OG5) %>%
  summarize(source = paste(Modes(source), collapse=", ")) %>%
  .$source

colors = c("#FFFFFF", RColorBrewer::brewer.pal(9, "YlOrRd"))
values = c(0,5,10,20,35,50,75,100,200,300)
ggplot(data = counts_by_major, aes(x = OG5, y = major_clade, fill=n_tips)) +
  geom_tile() +
  geom_text(aes(label=n_tips), size=3) +
  theme_bw() +
  scale_fill_gradientn(colors=colors, values=values/300) +
  scale_x_discrete(labels=paste0(og_list, ": \n", gene_family_list, " (", source_list, ")")) +
  scale_y_discrete(labels=paste0(levels(counts_by_major$major_clade), " (N=",taxon_counts, ")")) +
  labs(x = "Gene Family", y = "Major Clade Code", fill = "Number of Tips",
       title="Pre-Guidance") +
  theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1, size=10, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        legend.text = element_text(size=10, color="black"),
        legend.title= element_text(size=14, color="black"),
        legend.position="top")

ggsave("~/Desktop/Katzlab/Meiosis_Genes/tips_heatmap_preguidance.png", plot=last_plot(), device="png", dpi=600,
       width=20, height=25, units='cm')


colors = c("#FFFFFF", RColorBrewer::brewer.pal(9, "GnBu"))
values = c(0,5,10,15,20,25,30,40,60,80)
ggplot(data = counts_by_major, aes(x = OG5, y = major_clade, fill=n_taxa)) +
  geom_tile() +
  geom_text(aes(label=n_taxa)) +
  theme_bw() +
  geom_text(aes(label=n_taxa), size=3) +
  theme_bw() +
  scale_fill_gradientn(colors=colors, values=values/80) +
  scale_x_discrete(labels=paste0(og_list, ": \n", gene_family_list, " (", source_list, ")")) +
  scale_y_discrete(labels=paste0(levels(counts_by_major$major_clade), " (N=",taxon_counts, ")")) +
  labs(x = "Gene Family", y = "Major Clade Code", fill = "Number of Taxa",
       title="Pre-Guidance") +
  theme(axis.text.x = element_text(angle=90, vjust=0, hjust=1, size=10, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        legend.text = element_text(size=10, color="black"),
        legend.title= element_text(size=14, color="black"),
        legend.position="top")

ggsave("~/Desktop/Katzlab/Meiosis_Genes/taxa_heatmap_preguidance.png", plot=last_plot(), device="png", dpi=600,
       width=20, height=25, units='cm')
