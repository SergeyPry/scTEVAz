library(readr)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(paletteer)


# load the data 
data <- read.csv("temp.csv")

# specify value and label to use in the plot
y_value = "counts_sum"
y_label = "Normalized cell expression value"

gene <- "current_gene"

# sorting code


# the fonts were optimized during the initial testing but the user will
# need to adjust them further
ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col() +
          scale_y_continuous(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 14),
                strip.text.x = element_text(size = 11, margin = margin(0.12,0,0.12,0, "cm")),
                axis.text.y = element_text(size = 11),
                axis.title.y = element_text(size = 12),
                axis.title.x = element_text(size = 12),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 10, face="plain"),
                legend.title = element_text( size = 10, face = "bold"),
                legend.text = element_text( size = 10, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines")
          ) 


# optionally save the plot as desired including a new file name
# raster image format
ggsave("plot.png", dpi = 300)

# alternative vector formats (adjust the dimensions as needed)
ggsave("plot.svg", width = 12, height = 12, units = "in")
ggsave("plot.pdf", width = 12, height = 12, units = "in")

