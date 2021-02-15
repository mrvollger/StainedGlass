#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ncolors  = 11 # max value is 11

read_bedpe <- function(f){
  dfa = fread(f, col.names = c("q", "q_st","q_en","r","r_st","r_en", "perID_by_events") )
  #df = dfa[q == "chr8" & r == "chr8"]
  # color by quartiles of the data
  bot = floor(min(df$perID_by_events)); top = 100
  breaks = unique( c(quantile(df$perID_by_events, probs = seq(0, 1, by = 1/ncolors))) )
  labels = seq(length(breaks)-1)
  #print(length(labels))
  #print(length(breaks))
  df$discrete <- cut(df$perID_by_events, 
                     breaks=breaks,
                     labels = labels,
                     include.lowest=TRUE)
  window=max(df$q_en-df$q_st)
  df$first_pos = df$q_st/window
  df$second_pos = df$r_st/window
  return(df)
}
f="results/chm13.draft.v1.0.aln.bed"
f="results/chr1cen.10000.aln.bed"
df = read_bedpe(f); df

# get the lowest 0.1% of the data so we can not plot it
bot = quantile(df$perID_by_events, probs=0.001)[[1]]
p1 = ggplot(df, aes(df$perID_by_events, fill = discrete)) + geom_histogram( bins = 300) +
  theme_cowplot() + 
  scale_fill_brewer(palette = "Spectral", direction = -1) + theme(legend.position = "none") + 
  coord_cartesian(xlim = c(bot, 100))
p2 = ggplot(df) +
  geom_raster(aes(first_pos,second_pos, fill = discrete)) + 
  theme_cowplot() + scale_fill_brewer(palette = "Spectral", direction = -1) + theme(legend.position = "none", strip.text.x = element_blank()) + 
  coord_fixed(ratio = 1)
p3 = plot_grid(p1, p2, ncol=1,rel_heights = c(1,4));p3

ggsave(paste0(f, ".pdf"), plot=p3, height = 9, width = 9)


display.brewer.pal(10,"Spectral")
brewer.pal(10,"Spectral")
