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
y_value = "cell_counts"
y_label = "Normalized cell counts"

gene <- "current_gene"

# sorting code


# the fonts were optimized during the initial testing but the user will
# need to adjust them further
ggplot(data, aes(x = tissue, y = .data[[y_value]], fill = tissue)) +
          geom_col() +
          scale_y_sqrt(n.breaks = 5) +
          facet_wrap(~stage, scales = "free", ncol = 2) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_viridis_d(option = "turbo", direction = 1) +
            theme(plot.title = element_text(size = 14),
                  strip.text.x = element_text(size = 11, margin = margin(0.12,0,0.12,0, "cm")),
                  axis.text.y = element_text(size = 11),
                  axis.title.y = element_text(size = 12),
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, colour="grey20", 
                                             size= 8, face="plain"),
                  legend.title = element_text( size = 10, face = "bold"),
                  legend.text = element_text( size = 10, face = "plain"),
                  legend.key.size = unit(0.15, "cm"),
                  legend.position = "top", 
                  legend.direction = "vertical",
                  panel.spacing = unit(0.1, "lines")) +
            guides(fill=guide_legend(nrow = 5)) 



# optionally save the plot as desired
ggsave("plot.png", dpi = 300)