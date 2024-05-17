setwd("~/")

source("~/AdlungLab_HelperFunctions.R")
library(ggplot2)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Profiles + Validation ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


r <- readRDS("Data/r_Frede_i6.Rdata")

result_df <- data.frame(column = character(0), min_value = numeric(0), max_value = numeric(0))

params <- c("kturnB", "kturnM", "kturnN","kturnT") #names(r$params)[5:16]

for (i in params) {
  pp <- as.data.frame(r$profiles[which(r$profiles$whichPar %in% i), ])
  
  subset_result <- data.frame(
    column = paste0(i),
    min_value = min(pp[,i]),
    max_value = max(pp[,i])
  )
  
  result_df <- rbind(result_df, subset_result)
  
}


result_df$best <- as.numeric(r$params[1,params])
result_ordered <- result_df[order(result_df$column), ]
#result_ordered$type <- c(rep("Death",6),rep("Proliferation",6))

p1 <- ggplot(result_ordered, aes(x=column, y=best)) +
  geom_crossbar(aes(ymin=min_value, ymax=max_value,
                    fill=column), width=0.66, size=0.25) +
  # scale_fill_manual(values=ALcols[rep(c(8,4,10,6,5,9),2)]) +
  scale_fill_manual(values=ALcols[c(8,10,6,9)]) +
  labs(title="Inferred from Frede et al.", x=NULL, 
       y="log10(estimated value)") +
  #  facet_wrap(~type, scales = "free", nrow = 1) +
  theme_adlunglab(base_size = 10) +
  theme(legend.position = "none", axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

p1


r <- readRDS("Data/r_Ho_fit_i6.Rdata")

result_df <- data.frame(column = character(0), min_value = numeric(0), max_value = numeric(0))

params <- c("kturnB", "kturnM", "kturnN","kturnT") #names(r$params)[5:16]

for (i in params) {
  pp <- as.data.frame(r$profiles[which(r$profiles$whichPar %in% i), ])
  
  subset_result <- data.frame(
    column = paste0(i),
    min_value = min(pp[,i]),
    max_value = max(pp[,i])
  )
  
  result_df <- rbind(result_df, subset_result)
  
}


result_df$best <- as.numeric(r$params[1,params])
result_ordered <- result_df[order(result_df$column), ]
#result_ordered$type <- c(rep("Death",6),rep("Proliferation",6))

p2 <- ggplot(result_ordered, aes(x=column, y=best)) +
  geom_crossbar(aes(ymin=min_value, ymax=max_value,
                    fill=column), width=0.66, size=0.25) +
  # scale_fill_manual(values=ALcols[rep(c(8,4,10,6,5,9),2)]) +
  scale_fill_manual(values=ALcols[c(8,10,6,9)]) +
  labs(title="Inferred from Ho et al.", x=NULL, 
       y="log10(estimated value)") +
  #  facet_wrap(~type, scales = "free", nrow = 1) +
  theme_adlunglab(base_size = 10) +
  theme(legend.position = "none", axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

p2


require(readxl)
ki <- as.data.frame(read_excel("Data/HU2023_23.xlsx"))

kiplot <- melt(ki, measure.vars = c("Median_Bcell","Median_Tcell",
                                    "Median_Mac","Median_Neutr"))

kiplot$celltype <- factor(kiplot$variable, levels=c("Median_Bcell",
                                                    "Median_Tcell",
                                                    "Median_Mac",
                                                    "Median_Neutr"), 
                          labels=c("Bcell","Tcell","Mac","Neutr"))

kiplot$celltype <- factor(kiplot$celltype, levels= c("Bcell","Mac","Neutr","Tcell"))

# Validation signature ####

require(ggbeeswarm)
p3 <- ggplot(kiplot, aes(x=celltype,y=log10(value), fill=celltype)) +
  #geom_boxplot(outlier.shape = NA, width=0.66) + 
  geom_violin(color="#000000", draw_quantiles = c(0.5),width=1,trim = F,size=0.5) +
  #geom_crossbar() +
  geom_quasirandom(color="#000000", size=1, width=0.25) +
  scale_fill_manual(values=ALcols[c(8,10,6,9)]) +
  labs(title="This work", x=NULL, 
       y="log10(KI-67::PeCy7 median)") +
  theme_adlunglab(base_size = 10) +
  theme(legend.position = "none", axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

p3


require(ggpubr)
figure3 <- ggarrange(
  p1, p2, p3, #r$plots$param, 
  ncol = 3,widths = c(1/3,1/3,1/3),labels=c("A", "B","C")) 

figure3

fw <- 200
presdir <- "~/Output/"
ggsave(filename = paste0(presdir,"Fig3_", Sys.Date()  ,".pdf"),
       width=fw, height=fw*9/16/2, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig3_", Sys.Date()  ,".png"),
       width=fw, height=fw*9/16/2, units="mm", dpi=300)
