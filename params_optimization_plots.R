#!/usr/bin/Rscript
# r script to plot the metrics resulted from parameter optimization of stacks assembly 
# input file contains only one metric of all parameters. It should contain four columns, without header: parameter (e.g. m), param_value (e.g. m2), metric (e.g. assembled_loci) and metric_value (e.g. 2987)
# this script requires two arguments: an input file and an output name
#the output name should include the extention (e.g. .png)

#packages needed for plotting the results. It may be necessary to install them first
library(ggrepel)
library(ggplot2)
library(ggbreak)
library(lemon)
library(scales)
library(reprex)
library(ggtext)

args <- commandArgs(trailingOnly = TRUE)

#arguments
input_file <- args[1]
output_name <- args[2]

#import input file

input <- read.delim(file = input_file, sep = "\t", header = FALSE)

#create a header for the input file
colnames(input) <- c("parameter","param_value","metric","metric_value")
input$param <- as.factor(input$parameter)
input$metric <- as.factor(input$metric)


#plot the results as a bar plot

png(output_name)

ggplot(input, aes(parameter_value,metric_value)) + 
	facet_wrap(parameter, scales="free_x", strip.position = "left") +
	ylab(NULL) +
  	theme( strip.text.y = element_blank(), plot.title = element_text(size=30, hjust = .07, vjust = 3),
        	panel.background = element_rect(fill="grey96"),
        	legend.title = element_blank(), legend.text = element_blank(), legend.background = element_blank(),
        	axis.text = element_text(size=33), axis.text.x = element_text(size=30)) +
	geom_bar(stat = "identity", aes(fill=parameter))+
	#limits can be adjusted depending on the data
	expand_limits(y=c(31000))+
	geom_text(size=10,aes(label= value,y=value+600*sign( value)),vjust=.6, position=position_dodge(width = 1), angle=90)+
	xlab(NULL)+
	#scale breaks can be adjusted depending on the data, or even removed 
	scale_y_break(c(2100,24500), scales="fixed")+
	#title of the plot is optional
	ggtitle("a) Assembled loci")+
	scale_fill_brewer(palette = "Paired", labels=c("M","m","n"))
dev.off()
