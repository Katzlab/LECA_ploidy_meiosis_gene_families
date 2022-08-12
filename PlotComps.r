library(tidyverse)
library(FactoMineR)
library(factoextra)
library(data.table)
library(ggrepel)
library(grid)
library(gridExtra)
library(gt)
library(dplyr)

# change to your working directory
setwd('~/Desktop/caitlin/jf')

gc3 <- data.frame(read.csv('gc3_master.csv')) %>%
  filter(Type != "Type") %>% # remove excess header rows
  mutate(taxon = paste(substr(Seq, 1, 4), substr(Seq,6,10), sep = ''))
  
#gc3 <- gc3[sample(10000),]

taxa = unique(gc3$taxon)

gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)

diffs <- gc3 %>%
  group_by(taxon) %>%
  summarise(diff = max(GC3.Degen) - min(GC3.Degen))

gc3 <- gc3 %>%
  merge(diffs, by = 'taxon')

gc3 <- gc3[order(gc3$diff, decreasing = TRUE),]

gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

enc_null <- data.frame(read_tsv('ENc.Null.tsv'))
 
gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)
gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

gc3_am = gc3 %>%
  filter(startsWith(taxon, "Am"))
gc3_ba = gc3 %>%
  filter(startsWith(taxon, "Ba"))
gc3_ee = gc3 %>%
  filter(startsWith(taxon, "EE"))
gc3_ex = gc3 %>%
  filter(startsWith(taxon, "Ex"))
gc3_pl = gc3 %>%
  filter(startsWith(taxon, "Pl"))
gc3_op = gc3 %>%
  filter(startsWith(taxon, "Op"))
gc3_sr = gc3 %>%
  filter(startsWith(taxon, "Sr"))
gc3_za = gc3 %>%
  filter(startsWith(taxon, "Za"))

# colors for major clades (same as in phylotol)
ba = '#000000'
za = '#808080'
sr = '#7b2516'
op = '#12aaff' 
pl = '#006300' 
ex = '#ffa100'
ee = '#ff6288'
am = '#aa00ff'

# change data input to plotting code for each clade, and change point color to corresponding clade code. 
# example below is for EE clade
gc3_plot <- ggplot(gc3_ee, aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_6Fold))) + # change data here
  geom_point(size = .4, color=ee) + # change color here
  scale_color_manual(values = c('black', 'red')) +
  geom_line(data = enc_null, aes(GC3, ENc)) +
  theme_classic() +
  labs(x = 'GC3 Degen', y = 'ObsWrightENc_6Fold') +
  facet_wrap(~factor(taxon, levels = unique(taxon)))
        
gc3_plot

# saving plot output; change filename and path as necessary
ggsave(filename="../EE_jf.png", plot=last_plot(), device = "png", width=20, height=20, units='cm', dpi=600)
