require(readxl)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)

source("~/AdlungLab_HelperFunctions.R")
# load data ####

setwd("~/")

# load data from GSE171770 ####
meta.raw <- as.data.frame(read_excel("GSE171770_Clinical_Information_FUTURE.xlsx"))

# load processed counts ####
counts.ave <- read.csv("GSE171770_RawCounts_biopsy_Olamkicept_with_symbols.csv")
# # Convert to data frame
expr_matrix <- as.data.frame(counts.ave[,-1])
rownames(expr_matrix) = counts.ave[,1]
colnames(expr_matrix) = colnames(counts.ave)[-1]


# load deconvolution results ####
df <- readRDS("Data/results_decon_RNAseq_Ho.rds")

frac.out <- as.data.frame(t(df$out.all))
colnames(frac.out) <- colnames(expr_matrix)

rownames(meta.raw) <- meta.raw$Sample_Biopsy
meta.raw2 <- meta.raw[!rownames(meta.raw) %in% c("Sample2"),]
names(frac.out)[names(frac.out) == "H18624.L1_S55_L006_"] <- "H18624.L1_S55_L006"

frac.trans <- as.data.frame(t(frac.out))

meta.raw2$RowName <- rownames(meta.raw2)
frac.trans$RowName <- rownames(frac.trans)

meta.raw2 <- meta.raw2 %>%
  mutate(Time_seq = case_when(
    Time_seq == "0h"  ~ 0,
    Time_seq == "4h"  ~ 4,
    Time_seq == "24h" ~ 24,
    Time_seq == "2w"  ~ 336,   # 2 weeks * 7 days/week * 24 hours/day
    Time_seq == "6w"  ~ 1008,  # 6 weeks * 7 days/week * 24 hours/day
    Time_seq == "14w" ~ 2352,  # 14 weeks * 7 days/week * 24 hours/day
    TRUE ~ as.numeric(Time_seq) # Keeps the original value if none of the above conditions are met
  ))

merge.plot2 <- merge(meta.raw2, frac.trans, by = "RowName")

require(reshape2)

melt.plot2 <- melt(merge.plot2, measure.vars = c("Bcell","Mac","Neutr","Tcell","Epi","Other",
                                                 "pSTAT3_Epithel","pSTAT3_Lamina","Ki67_Epithel","Ki67_Lam",
                                                 "IL6","RIL6","STAT3_hIL6","CD68","CD3","MPO","Mayo","CDAI",
                                                 "Calprotectin"), 
                   id.vars = c("Patient","Time_seq","Remission","Response","Diagnose"))


celltypes <- read.csv("Output/CellTypeAnnotation_Ho_redone.csv")
color_assignment <- setNames(celltypes$celltypecolor[1:20], celltypes$CellType[1:20])

# Cross-Correlation ####

library(ggcorrplot)

cor.plot <- merge.plot2[merge.plot2$Remission=="R",c("Bcell","Mac","Neutr","Tcell",#"Epi", 
                                                     "pSTAT3_Epithel","pSTAT3_Lamina"#,"Ki67_Epithel","Ki67_Lam"
)]

corr <- round(cor(cor.plot,use = "pairwise.complete.obs"), 1)
p.mat <- cor_pmat(cor.plot)

#require(viridis)
fig4c <- ggcorrplot(corr, hc.order = F, #type = "lower",
                    p.mat = p.mat) +
  theme_adlunglab(base_size = 10) +
  labs(title = NULL, x= NULL, y=NULL) +
  scale_fill_gradient2(high = "red",
                       mid = "white",
                       low = "blue",
                       midpoint = 0,
                       name="Pearson's\nrho") +
  theme_adlunglab(base_size = 10) +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1),
        legend.position = "none")

fig4c

presdir <- "~/Dropbox/UKE/Files/Writing/Manuscripts/scMod_Manuscript/Figures/Revision/"
fw <- 130

ggsave(filename = paste0(presdir,"Fig4D_", Sys.Date()  ,".pdf"),
       width=fw, height=fw, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig4D_", Sys.Date()  ,".png"),
       width=fw, height=fw, units="mm", dpi=300)

# data extraction ####

plot <- melt.plot2[which(melt.plot2$variable %in% c("Bcell","Mac","Neutr","Tcell","Epi","Calprotectin") 
                         & melt.plot2$Remission %in% c("R")),] 
colnames(plot)[7] <- "freq"

require(tidyr)
require(dplyr)

schreiber.sum <- plot %>%
  group_by(variable, Time_seq) %>%
  summarize(value = mean(freq, na.rm = TRUE),
            sigma = sd(freq, na.rm = TRUE))

schreiber.filtered <- schreiber.sum %>% filter(sigma != 0)

saveRDS(schreiber.filtered, "Output/Schreiber_R_Data.Rds")

# dMod ####

# clear workspace otherwise profile function won't work

rm(list = ls())

library(ggplot2)

source("~/AdlungLab_HelperFunctions.R")
# load data ####

setwd("~/")

schreiber.filtered <- readRDS("Data/Schreiber_R_Data.Rds")

#plot <- readRDS("FrequencyTable.Rds")
check <- as.data.frame(schreiber.filtered)
#check$celltype <- as.character(check$celltype)

library(dplyr)
require(tidyr)
ag <- check %>% 
  rename(name=variable)  %>% rename(time=Time_seq)

library(dMod)
library(cOde)

require(reshape2)

# ag$valreshape2# ag$value <- ag$value/100
ag$sigma <- ag$sigma/sqrt(3)

#ag <- ag[-which(ag$time %in% c(4,24)),]
ag$time <- ag$time/24

ag$value[ag$name=="Calprotectin"] <- ag$value[ag$name=="Calprotectin"]/1000
ag$sigma[ag$name=="Calprotectin"] <- ag$sigma[ag$name=="Calprotectin"]/1000


ag <- ag[-which(ag$name %in% c("Calprotectin")),]

#ag <- ag[which(ag$time %in% c(0,14,42,98)),]

# ag <- ag[-which(ag$name == "Bcell" & ag$time==1),]

data <- datalist(
  Schreiber = as.data.frame(ag))

timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))


{
  f <- NULL
  
  # B cells ####                      
  f <- addReaction(f, from = "", to = "Bcell", 
                   rate = "Resolution*kturnB*Bcell",
                   description = "B_Proliferation")
  f <- addReaction(f, from = "Bcell", to = "", 
                   rate = "(1/kturnB)*Bcell",
                   description = "B_Death")
  
  # Macrophages ####                               
  f <- addReaction(f, from = "", to = "Mac", 
                   rate = "ActiveDisease*kturnM*Mac",
                   description = "Mac_Proliferaton")
  f <- addReaction(f, from = "Mac", to = "",
                   rate = "(1/kturnM)*Mac",
                   description = "Mac_Death")
  
  # Neutrophils ####                                 
  f <- addReaction(f, from = "", to = "Neutr", 
                   rate = "(DSS+1)*kturnN*Neutr",
                   description = "N_Activation")
  f <- addReaction(f, from = "Neutr", to = "", 
                   rate = "(1/kturnN)*Neutr",
                   description = "Neutr_Death")
  
  # T cells ####
  f <- addReaction(f, from = "", to = "Tcell", 
                   rate = "ActiveDisease*kturnT*Tcell",
                   description = "T_Proliferation")
  f <- addReaction(f, from = "Tcell", to = "",
                   rate = "(1/kturnT)*Tcell",
                   description = "T_Death")
  
  # # Locals ####                                  
  # f <- addReaction(f, from = "", to = "Cue",
  #                  rate = "kcue",
  #                  description = "Cue")
  # 
  # f <- addReaction(f, from = "Cue", to = "",
  #                  rate = "kcue*Cue",
  #                  description = "Uncue")
  
  f <- addReaction(f, from = "", to = "Resolution",
                   rate = "kres",
                   description = "Resolution")
  
  f <- addReaction(f, from = "Resolution", to = "",
                   rate = "1/(kres)*Resolution",
                   description = "Deresolution")
  
  
  f <- addReaction(f, from = "", to = "ActiveDisease",
                   rate = "kact",
                   description = "DiseaseActivation")
  
  f <- addReaction(f, from = "ActiveDisease", to = "",
                   rate = "1/(kact)*ActiveDisease",
                   description = "DiseaseDeactication")
  
  
  f <- addReaction(f, from = "Epi", to = "",
                   rate = "(DSS+1)*(1/kturnE)*Epi",
                   description = "Epi_Death")
  
  f <- addReaction(f, from = "", to = "Epi",
                   rate = "kturnE*Epi",
                   description = "Epi_Proliferation")
  
  f <- addReaction(f, from = "DSS", to = "2*DSS",
                   rate = "0",
                   description = "DSS") 
  
  print(f)             
  e <- NULL
  e <- eventlist() %>%
    addEvent(var = "DSS", time = 1, value = 0)
  
}

dmod_pipe <- function(f, e) { 
  model <- odemodel(f, modelname = "Best")
  x <- Xs(model, events = e)
  observables <- eqnvec(#DSS = "DSS", 
    ActiveDisease= "ActiveDisease", 
    Resolution = "Resolution", 
    DSS = "DSS", 
    Bcell = "Bcell", Epi = "Epi", #Stromal = "Stromal", 
    Mac = "Mac", Tcell = "Tcell", Neutr="Neutr")
  g <- Y(observables, x, compile = TRUE, modelname = "Best", attach.input = FALSE)
  innerpars <- getParameters(g*x)
  trafo <- repar("x~x", x = innerpars)
  #trafo <- repar(paste0("x~",DSSt0), x = c("DSS"), trafo)
  #trafo <- repar("x~0", x = c("healedMucosa"), trafo)
  #  trafo <- repar("x~0", x = c("damage"), trafo)
  trafo <- repar("x~exp(x)", x = innerpars, trafo)
  trafo1 <- trafo
  p <- P(trafo1, condition = "Schreiber")
  outerpars <- getParameters(p)
  pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)
  timesM <- unique(data$Schreiber$time) #seq(0,100,0.1)
  times <- unique(data$Schreiber$time) #0:100
  
  p_init <- plot((g*x*p)(times, pouter)) +
    scale_color_manual(values=c("#220589")) +
    labs(title= NULL, x= "Time (days", y="Relative frequency") +
    #  theme_adlunglab(base_size = 12) +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = "none")
  prior <- structure(rep(0, length(pouter)), names = names(pouter))
  obj <- normL2(data, g*x*p) + constraintL2(mu = prior, sigma = 10)
  myfit <- trust(obj, pouter, rinit = 1, rmax = 10)
  fitOutput <- mstrust(objfun=obj,center=pouter,studyname="fits", fits = 100, cores = 4)
  fitResults  <- as.parframe(fitOutput)
  
  p_waterfall <- plotValues(fitResults) +
    scale_color_gradient(low = "#220589", high = "#F0E634") +
    #  scale_color_gradient(guide = F, alpha=1) +
    scale_shape_discrete(guide=F) +
    labs(title= NULL, x= "run index (sorted by likelihood)", 
         y="objective function") +
    #  facet_wrap(nrow = 1) +
    theme_adlunglab(base_size = 10) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.22, "cm"),
          legend.key.width = unit(0.33, "cm"),
          legend.direction = "horizontal",
          legend.title.align = 1)
  
  simuplot <- (g*x*p)(timesM, myfit$argument)
  sp <- as.data.frame(simuplot$Schreiber)
  
  dss <- ggplot(sp[,c(1,2)]) +
    geom_line(aes(x=time, y=Calprotectin), color="#000000", linewidth=2) +
    labs(title= "Model input", x= "Time (days)", y="Calprotectin (g/g)") +
    theme_adlunglab(base_size = 10) +
    #scale_x_continuous(breaks=c(0,7,14)) +
    theme(legend.position = "none")
  
  # remove DSS
  # sp <- sp[,-c(2)]
  #sp <- sp[,-c(3)]
  
  sm <- melt(sp, id.vars = "time")
  sm$alpha <- "Model"
  dm <- as.data.frame(data$Schreiber)
  colnames(dm)[c(1)] <- c("variable")
  dm$alpha <- "Data"
  dm$variable <- factor(dm$variable, levels = sort(unique(dm$variable)))
  cols.model <- c(ALcols[c(8,10,6,9,4)],rep("#000000",3))
  
  # "#D90767" "#09A64A" "#F39F07" "#8D2668" "#EB5D12"
  
  p_fit <- ggplot() +
    #  geom_point(data=dm, aes(x=timepoint,y=value),color=met.brewer("Kandinsky",2)[2]) +
    geom_pointrange(data=dm, aes(x=time,ymin=(value-sigma), color=variable,
                                 ymax=(value+sigma), y=value, alpha=factor(alpha))) +
    geom_line(data=sm, aes(x=time, y=value, alpha=factor(alpha),
                           color=variable), linewidth=2) +
    labs(title= NULL, x= "Time (days)", y="Cell number (%)") +
    theme_adlunglab(base_size = 10) +
    scale_alpha_manual(name=NULL, values=c(1,0.6),
                       guide = guide_legend(override.aes = list(
                         linetype = c("blank", "solid"),
                         shape = c(16, NA)))) +
    scale_color_manual(values=cols.model, guide=F) +
    scale_x_continuous(breaks=unique(dm$time)) +
    facet_wrap(~variable, scales = "free", nrow = 1) +
    theme(legend.position = "none")
  
  k <- length(outerpars)
  n <- nrow(data$Schreiber)
  L <- fitResults[1]$value
  
  cAIC <- fitResults[1]$value + k * 2 + (2*k*(k+1))/(n-k-1) #\frac{2k(k+1)}{n-k-1}
  bestfitpars <- fitResults[1]
  
  
  bestfit <- myfit$argument
  profiles <- profile(obj, bestfit, names(bestfit), limits = c(-10, 10), cores = 4)
  
  
  simu <- as.data.frame((g*x*p)(timesD, pouter)$Schreiber)
  
  # Take a look at each parameter
  p_param <- plotProfile(profiles) + 
    scale_color_manual(values="#000000", guide=F) +
    scale_linetype_manual(values=c("dotted","longdash","solid"), name=NULL) +
    labs(title= NULL, x= "log10(parameter value)") +
    #  facet_wrap(nrow = 1) +
    theme_adlunglab(base_size = 6) +
    theme(legend.position = "right")
  plots <- list(init = p_init, fit = p_fit, param = p_param, waterfall = p_waterfall)
  return(list(dss = dss, aic = cAIC, plots = plots, 
              params = bestfitpars, profiles = profiles,
              k = k, n = n, L = L, simu = sp,f=f))  
}


for(i in 1:20){
  try({
    suppressMessages({suppressWarnings({r <- dmod_pipe(f, e)})})
    break #break/exit the for-loop
  }, silent = FALSE)
  message("Error! Re-Trying! #Try:", i)
}

r$plots$fit #+ theme_minimal()
r$plots$param #+ theme_minimal()
r$plots$waterfall
r$params
r$aic
r$L
r$n
r$k

rprint <- r

saveRDS(r,"Output/SchreiberFitBestModelFinal.Rdata")



# Frequency plots ####

r <- readRDS("Data/SchreiberFitBestModelFinal.Rdata")

dm <- data$Schreiber

A <- ggplot(data=dm[which(dm$name %in% c("Bcell","Mac","Neutr","Tcell","Epi")),],
            aes(x=factor(time*24),y=value)) +
  geom_col(aes(fill=name), position = "fill") +
  labs(x="Time",y="Relative abundance",title="Deconvoluted data") +
  scale_x_discrete(labels=c("0h","4h","24h","2w","6w","14w")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color_assignment,name=NULL) +
  theme_adlunglab(base_size = 10) +
  theme(axis.line.x = element_blank())

A

require(reshape2)
sm <- melt(r$simu, measure.vars = c("Bcell","Mac","Neutr","Tcell","Epi"), id.vars = "time")

B <- ggplot(data=sm,
            aes(x=factor(time*24),y=value)) +
  geom_col(aes(fill=variable), position = "fill") +
  labs(x="Time",y="Relative abundance",title="Model fit") +
  scale_x_discrete(labels=c("0h","4h","24h","2w","6w","14w")) +
  scale_fill_manual(values = color_assignment,name=NULL) +
  scale_y_continuous(expand = c(0,0)) +
  theme_adlunglab(base_size = 10) +
  theme(axis.line.x = element_blank())

B


require(ggpubr)
figure4 <- ggarrange(
  A, B, labels = c("C",""), 
  common.legend = T, legend = "right",
  ncol = 2,widths = c(1/2,1/2))

figure4

fw <- 200
presdir <- "~/Output/"
ggsave(filename = paste0(presdir,"Fig4_FreqBand_", Sys.Date(), ".pdf"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300, useDingbats=FALSE)

ggsave(filename = paste0(presdir,"Fig4_FreqBand_", Sys.Date(), ".png"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300)





r <- readRDS("Data/SchreiberFitBestModelFinal.Rdata")


pred <- as.data.frame(data$Schreiber)
simu <- r$simu
simu.m <- melt(simu, id.vars = "time")
#simu.m <- simu %>% gather(variable, value, -time)

pred$join <- paste0(pred$time, "_", pred$name)
simu.m$join <- paste0(simu.m$time, "_", simu.m$variable)

cormat <- full_join(pred, simu.m, by="join")

cormat.plot <- cormat#[which(cormat$time.x != 0),]

ct <- cor.test(cormat.plot$value.x,cormat.plot$value.y)



c <- ggplot(data=cormat.plot[-which(is.na(cormat.plot$name)),], aes(x=(value.x), y=(value.y))) +
  labs(title  = paste0("Pearson's rho=", round(ct$estimate,2), "; p=", signif(ct$p.value,1)), # plot title
       x      = "Model fit cell fraction (r.u.)" , # x-axis label
       y      = "Deconvoluted bulk RNA-seq data\nCell fraction (r.u.)") + # y-axis label
  geom_smooth(method = "lm", color="#000000", se = T, size=0.5) +
  # geom_abline(slope=1,intercept = 0, size=0.5, color="#000000", linetype="dashed") +
  geom_point(aes(fill=name), size=2,
             color="#000000", stroke=0.2, shape=21) +
  scale_fill_manual(values = color_assignment, name=NULL) +
  #coord_cartesian(expand = F) +
  # xlim(c(0,3)) + ylim(c(0,3)) +
  #  coord_fixed(expand = F) +
  #  expand_limits(y=0) +
  theme_adlunglab(base_size = 10) +
  theme(legend.key.size = unit(0.42, 'lines'), legend.position = "right") # +
# guides(fill = guide_legend(override.aes = list(size = 2)))                                                                                      color = c("#000000", "#000000"))))

c #+ theme_minimal()

fw <- 100

ggsave(filename = paste0(presdir,"FitvsData_Schreiber_", Sys.Date()  ,".pdf"), 
       width=fw, height=200*9/16/1.67, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FitvsData_Schreiber_", Sys.Date()  ,".png"), 
       width=fw, height=200*9/16/1.67, units="mm", dpi=300)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Neutrophil validation 




# library(GEOquery)
# gsm <- getGEO('GSE16879')


# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)


# load series and platform data from GEO

gset <- getGEO("GSE16879", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }


countmatrix <- gset@assayData$exprs
rownames(countmatrix) <- gset@featureData@data$Gene.Symbol

source("combMASlots.R")

meanmatrix <- combMASlots(countmatrix, func = mean)

saveRDS(meanmatrix,"meanmatrix_Arjis.Rds")

metadata <- gset@phenoData@data[,c(2,8,10:13)]
colnames(metadata) <- c("SampleID","Description","Tissue","Disease","Response","TimePoint")
metadata$Tissue <- sub(".*: ", "", metadata$Tissue)
metadata$Disease <- sub(".*: ", "", metadata$Disease)
metadata$Response <- sub(".*: ", "", metadata$Response)
metadata$TimePoint <- sub(".*: ", "", metadata$TimePoint)

saveRDS(metadata,"metadata_Arjis.Rds")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Vulcano abundance ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

g1 <- metadata$SampleID[metadata$TimePoint == "Before first infliximab treatment" & metadata$Response == "Yes"]
g3 <- metadata$SampleID[metadata$TimePoint == "After first infliximab treatment" & metadata$Response == "Yes"]

# Fold-change
logFC = as.vector(apply(meanmatrix, 1, function(x) log2((mean(2^x[g1]))/(mean(2^x[g3])))))
logFCx = as.vector(apply(meanmatrix, 1, function(x) (mean(x[g1]))))
logFCy = as.vector(apply(meanmatrix, 1, function(x) (mean(x[g3]))))

dim(meanmatrix)

# Wilcoxon-Mann-Whitney test
pval = as.vector(apply(meanmatrix, 1, function(x) wilcox.test(2^x[g1], 2^x[g3])$p.value#, exact=TRUE
))
wilc = p.adjust(pval, method = "fdr")


pb <- data.frame(rownames(meanmatrix), logFC, logFCx, logFCy, pval, wilc)
colnames(pb) <- c("gene", "fc", "x", "y","pval", "qval")

pb$qval[pb$qval=="NaN"] <- NA
pb <- na.omit(pb)

# write results to file
write.csv(pb,file=paste0("Output/DGE_Before_vs_After_Responders_mean.csv"))

# LAOK ####

#pb <- read.csv(file="Data/DGE_Before_vs_After_Responders_mean.csv")
#pb <- pb[-which(abs(pb$fc) == "Inf"),]
#pb <- pb[-which(pb$gene %in% c("")),]
duplicates <- duplicated(pb$gene)
pb <- pb[!duplicates, ]
# pb <- pb[-which(pb$pval > 0.1),]
rownames(pb) <- pb$gene #sub(" ///.*", "", pb$gene)

threshold <- 0.1#quantile(abs(pb$pval),1:20/20)[1]
pb[["mark"]] = ifelse(pb$qval < 0.1#pb$fc > 1.96*sd(pb$fc, na.rm=T)
                      , "#220589", "#E9EBED")
show <- c("S100A8","S100A9")

pb[["show"]] = ifelse(pb$gene %in% show, as.character(pb$gene), NA)
mark <- as.character((pb[,"mark"]))
names(mark) <- pb$mark
# pb <- pb[-which(pb$x < -3.950468 & pb$fc > 0 | pb$y < -3.950468 & pb$fc < 0),]
# pb <- pb[-which(pb$pval == 1),]
# pb$pval[pb$pval == 0] <- 6.891173e-281
pb$lp <- -log10(pb$qval)


require(ggrepel)

df <- data.frame(label=c("Before","After"), x = c(-2.5,2.5), y=c(7,7))

p <- ggplot(pb, aes(x=(-fc), y=(lp))) + 
  labs(title  = "IBD patients responding to infliximab", # plot title
       x      = paste0("log2(fold-change)") ,
       y      = paste0("-log10(p adj.)")) + 
  geom_point(aes(fill=mark), shape=21, size=1.5, color="#220589",stroke=0.15) +
  geom_vline(xintercept=0, colour="grey", na.rm = FALSE, show.legend = NA, linetype = "dotted", size=1) +
  geom_hline(yintercept=-log10(0.1), colour="grey", na.rm = FALSE, show.legend = NA, linetype = "dotted", size=1) +
  geom_label_repel(aes(label=show), nudge_x = -2 #, angle=0, size=2.5, vjust=0, hjust=-0.25
  ) +
  geom_text(data=df, aes(label=label, x=x,y=y), fontface="bold") +
  scale_y_continuous(limits = c(0,7.5), expand = c(0,0)) +
  xlim(-6,6) +
  scale_fill_manual(values=mark, guide="none") + 
  theme_adlunglab(base_size=12) +
  theme(plot.title = element_text(hjust = 0.5))

p

presdir <- "~/Output/"
fw <- 80

ggsave(filename = paste0(presdir,"Fig4_Vulcano_", Sys.Date()  ,".pdf"),
       width=fw, height=200*9/16/1.67, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig4_Vulcano_", Sys.Date()  ,".png"),
       width=fw, height=200*9/16/1.67, units="mm", dpi=300)


# LAOK ####
