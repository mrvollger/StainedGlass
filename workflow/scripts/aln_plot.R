#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(glue)
require("argparse")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ncolors  = 11 # max value is 11
########################################################################################################
make_scale = function(vals){
  comma(vals/1e6)
}

make_k = function(vals){
  comma(vals/1e3)
}

get_colors = function(sdf){
  bot = floor(min(sdf$perID_by_events)); top = 100
  breaks = unique( c(quantile(sdf$perID_by_events, probs = seq(0, 1, by = 1/ncolors))) )
  labels = seq(length(breaks)-1)
  return( cut(sdf$perID_by_events, breaks=breaks, labels = labels, include.lowest=TRUE)  )
}

read_bedpe <- function(all.files){
  l <- lapply(all.files, fread, sep="\t")
  df <- rbindlist( l )
  
  df$discrete = get_colors(df)
  if("#query_name" %in% colnames(df)){
    df$q = df$`#query_name`
    df$q_st = df$query_start
    df$q_en = df$query_end
    df$r = df$reference_name
    df$r_st = df$reference_start
    df$r_en = df$reference_end
  } 
  window=max(df$q_en-df$q_st)
  df$first_pos = df$q_st/window
  df$second_pos = df$r_st/window
  return(df)
}

# get the lowest 0.1% of the data so we can not plot it
make_hist = function(sdf){
  bot = quantile(sdf$perID_by_events, probs=0.001)[[1]]
  p = ggplot(data=sdf, aes(perID_by_events, fill = discrete)) + geom_histogram( bins = 300) +
    theme_cowplot() + 
    scale_fill_brewer(palette = "Spectral", direction = -1) + theme(legend.position = "none") + 
    scale_y_continuous(labels=make_k) +
    coord_cartesian(xlim = c(bot, 100))+
    xlab("% identity")+ylab("# of alignments (thousands)")
  p
}

diamond <- function( row ){
  #side_length, x, y, color) {
  #print(row)
  side_length = as.numeric(row["window"])
  x=as.numeric(row["w"])
  y=as.numeric(row["z"])
  
  base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) * sqrt(2) / 2
  trans <- (base * side_length) + c(x,y)
  df = as.data.frame(t(trans))
  colnames(df) = c("w","z")
  df$discrete = as.numeric(row["discrete"])
  df$group = as.numeric(row["group"])
  df
}

make_tri = function(sdf, rname=""){
  #sdf = read_bedpe("~/Desktop/repos/StainedGlass/results/chr8.2000.10000.bed")[`#query_name` == "chr8" & reference_name == "chr8"]
  sdf$w = (sdf$first_pos  + sdf$second_pos) 
  sdf$z = -sdf$first_pos + sdf$second_pos
  window = max(sdf$q_en - sdf$q_st)
  scale =  max(sdf$q_st)/max(sdf$w) 
  sdf$window = max(sdf$q_en - sdf$q_st)/scale
  sdf$group = seq(nrow(sdf))
  df.d = rbindlist(apply(sdf, 1, diamond ))
  #df.d %>% summarise(max(w-z))
  #ggplot(sdf) +
    #geom_tile(aes(x = w*scale, y = z*window , fill = discrete)) + 
  ggplot(df.d)+
    geom_polygon(aes(x = w*scale, y = z*window , group=group, fill = factor(discrete))) + 
    theme_cowplot() + 
    scale_fill_brewer(palette = "Spectral", direction = -1) + 
    scale_x_continuous(labels=make_scale) +
    scale_y_continuous(labels=make_scale, limits = c(0,NA)) +
    #coord_fixed(ratio = 1) +
    xlab("Genomic position (Mbp)") + ylab("") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) + 
    ggtitle(rname)
}
#sdf
#rbindlist(apply(head(sdf), 1, diamond ))
#t = df[r==df$r[1] & q_en < 1e6 & r_en < 1e6]
#t = df[r == df$r[1]]
#make_tri(t, rname="a")
#make_dot(t, rname="a")

make_dot = function(sdf, rname=""){
  max = max(sdf$q_en, sdf$r_en)
  window = max(sdf$query_end - sdf$query_start)
  ggplot(sdf) +
    geom_tile(aes(x = q_st, y = r_st, fill = discrete, height=window, width=window)) + 
    theme_cowplot() + 
    scale_fill_brewer(palette = "Spectral", direction = -1) + 
    theme(legend.position = "none") + 
    scale_x_continuous(labels=make_scale, limits = c(0,max)) +
    scale_y_continuous(labels=make_scale, limits = c(0,max)) +
    coord_fixed(ratio = 1) +
    facet_grid( r ~ q )+
    xlab("Genomic position (Mbp)") + ylab("") +
    ggtitle(rname)
}


make_plots <- function(r_name) {
  sdf = copy(df[r == r_name & q == r_name])
  # set the colors
  if(!ONECOLORSCALE){
    sdf$discrete = get_colors(sdf)
  }
  
  # make the plots
  if(TRI){
    p_lone = make_tri(sdf, rname=r_name)
    scale = 2/3
  }else{
    scale = 1
    p_lone = make_dot(sdf, rname=r_name)
  }
  p_hist = make_hist(sdf)
  
  
  # setup the space
  dir.create(glue("{OUT}/pdfs/{r_name}/"))
  dir.create(glue("{OUT}/pngs/{r_name}/"))
  
  # save the plots
  plot = cowplot::plot_grid(
    p_lone, p_hist,
    ncol=1, 
    rel_heights = c(3*scale,1)
  )
  ggsave(plot=plot,
         file=glue("{OUT}/pdfs/{r_name}/{PRE}__{r_name}__tri.{TRI}__onecolorscale.{ONECOLORSCALE}.pdf"),
         height = 12*scale, width = 9)
  
  ggsave(plot=plot,
         file = glue("{OUT}/pngs/{r_name}/{PRE}__{r_name}__tri.{TRI}__onecolorscale.{ONECOLORSCALE}.png"), 
         height = 12*scale, width = 10, dpi = DPI)
  
  p_lone
}
########################################################################################################
#
# EDIT THIS SECION FOR YOUR INPUTS
#
parser <- ArgumentParser()
parser$add_argument("-b", "--bed",  help="bedfile with alignment information")
parser$add_argument("-p", "--prefix",  help="Prefix for the outputs")
args <- parser$parse_args()

PRE = args$prefix
GLOB = args$bed
OUT=glue("results/{PRE}_figures")
print(PRE) 
print(GLOB)

DPI=600
#
# STOP EDITING
#
########################################################################################################

dir.create(OUT)
dir.create(glue("{OUT}/pdfs"))
dir.create(glue("{OUT}/pngs"))
all.files = Sys.glob(GLOB)
df = read_bedpe(all.files)
#df=fread(GLOB)
Qs = unique(df$q)
N=length(Qs)
columns = ceiling(sqrt(N+1))
rows = ceiling( (N+1) / columns)

#
# big plot
#
if(T){
  facet_fig = cowplot::plot_grid(make_hist(df), make_dot(df), rel_heights = c(1,4), ncol=1)
  ggsave(plot=facet_fig, file=glue("{OUT}/pdfs/{PRE}.facet.all.pdf"), height = 20, width = 16)
  ggsave(plot=facet_fig, file=glue("{OUT}/pngs/{PRE}.facet.all.png"), height = 20, width = 16, dpi=DPI)
}
vals = c(TRUE, FALSE)
for(TRI in vals){
  if(TRI){
    scale = 2/3
  }else{
    scale = 1
  }
  for(ONECOLORSCALE in vals){
      plots = lapply(Qs, make_plots)
      plots[[N+1]] = make_hist(df) 
      p = cowplot::plot_grid(plotlist = plots, nrow=rows, ncol=columns, labels = "auto");
      ggsave(glue("{OUT}/pdfs/{PRE}.tri.{TRI}__onecolorscale.{ONECOLORSCALE}__all.pdf"), plot=p, height = 6*rows*scale, width = 6*columns)
      ggsave(glue("{OUT}/pngs/{PRE}.tri.{TRI}__onecolorscale.{ONECOLORSCALE}__all.png"), plot=p, height = 6*rows*scale, width = 6*columns, dpi=DPI)
  }
}

if(F){
  ONECOLORSCALE=FALSE
  TRI=TRUE
  z = make_plots(Qs[1])
  z
}

