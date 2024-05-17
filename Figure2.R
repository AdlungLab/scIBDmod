setwd("~/")

require(dMod)
require(readxl)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)

source("~/AdlungLab_HelperFunctions.R")

# load GSE data from Ho et al. ####
# tc_ibd
count_table <- read.table("raw_data/GSE148794_tc_ibd.count_table.tsv")
metadata <- read.table("raw_data/tc_ibd.metadata.tsv", fill = TRUE)

# Generate Seurat object ####
require(Seurat)
ho_colitis <- CreateSeuratObject(counts = count_table)

# edit metadata ----------------------------------------------------------------

# create a new column with the sample name
ho_colitis@meta.data$sample <- rownames(ho_colitis@meta.data)

# split sample column into different varibles
ho_colitis@meta.data <- separate(ho_colitis@meta.data, col ="sample", into = c("id", "id2", "id3", "id4", "id5", "id6", "barcode"))

# create column with time

# | 190426 = control & day 9 | 190417 = day 3 & day 6 | 190509 = day 12 & day 15 |
# | 1-10 control & 11-20 day 9 | 1-10 day 3 & 11-20 day 6 | 1-10 day 12 & 11-20 day 15 |

# using pattern I found i metadata to find time for probe... control = timepoint ??? 

ho_colitis@meta.data <- ho_colitis@meta.data %>% mutate(time = case_when(id == "R190426" & id6 <= 10 ~ "00",
                                                                         id == "R190426" & id6 > 10 ~ "09",
                                                                         id == "R190417" & id6 <= 10 ~ "03",
                                                                         id == "R190417" & id6 > 10 ~ "06",
                                                                         id == "R190509" & id6 <= 10 ~ "12",
                                                                         id == "R190509" & id6 > 10 ~ "15"))


saveRDS(ho_colitis, "Output/ho_colitis.rdata")


# load data ####
load("ho_colitis.rdata")

# # read annotations ####
require(readxl)
meta <- read_excel("~/Data/mmc1.xlsx",
                   sheet = "T1.IBD UMAP & Cluster", skip = 2)

ho_colitis@meta.data$name <- rownames(ho_colitis@meta.data)

ho_colitis@meta.data$celltype = meta$`cell name`[match(ho_colitis@meta.data$name,meta$cell_id)]
ho_colitis@meta.data$timepoint = meta$`time course`[match(ho_colitis@meta.data$name,meta$cell_id)]

Idents(ho_colitis) <- ho_colitis$seurat_clusters
Idents(ho_colitis) <- ho_colitis$celltype

ho_colitis[["percent.mt"]] <- PercentageFeatureSet(ho_colitis, pattern = "^MT-")

VlnPlot(ho_colitis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

DimPlot(ho_colitis, reduction = "umap")

FeaturePlot(ho_colitis, features = c("Trem2", "Cd9","Cd4","Cd8a","Cd68",
                                     "Tlr2","Ptprc","Epcam","Lgr5"), order = T)


ho.markers <- FindAllMarkers(ho_colitis, only.pos = TRUE)


ho_colitis$cluster_type <- paste0(ho_colitis$seurat_clusters, "_", ho_colitis$celltype)

unique(ho_colitis$cluster_type)


# redo clustering ####
FeatureScatter(ho_colitis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


ho.sub <- subset(ho_colitis, subset = nFeature_RNA > 250 & nFeature_RNA < 6000)

ho.sub <- NormalizeData(ho.sub, normalization.method = "LogNormalize", scale.factor = 10000)

ho.sub <- FindVariableFeatures(ho.sub, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ho.sub), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ho.sub)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(ho.sub)
ho.sub <- ScaleData(ho.sub, features = all.genes)

ho.sub <- RunPCA(ho.sub, features = VariableFeatures(object = ho.sub))

ElbowPlot(ho.sub)


ho.sub <- FindNeighbors(ho.sub, dims = 1:20)
ho.sub <- FindClusters(ho.sub, resolution = 0.5)


ho.sub <- RunUMAP(ho.sub, dims = 1:20)

DimPlot(ho.sub, reduction = "umap", label = T, label.size = 6) + NoLegend() + NoAxes() + 
  ggtitle(paste0(dim(ho.sub)[2], " cells"))


saveRDS(ho.sub,"Output/ho.sub.rdata")

# Annotation ####
readRDS("Output/ho.sub.rdata")


ho.sub <- readRDS("Output/ho_sub.rdata")

Idents(ho.sub) <- ho.sub$celltype

celltypes <- read.csv("Data/CellTypeAnnotation_Ho_redone.csv")
color_assignment <- setNames(celltypes$celltypecolor[1:20], celltypes$CellType[1:20])


fig2b <- DimPlot(ho.sub, cols = color_assignment, order = F, label=T, label.size = 3.5,label.box = T,
                 label.color = "#FFFFFF",repel = T) +
  NoAxes() + NoLegend()

fig2b


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data extraction ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hoplot <- as.data.frame(ho.sub@meta.data)

require(tidyr)
require(dplyr)
plot <- hoplot %>%
  group_by(id2,time,id6, celltype, .drop=FALSE) %>%
  summarise(n = n()) #%>%
# mutate(freq = n / sum(n))

plot$time <- as.numeric(plot$time)

plot$celltype <- factor(plot$celltype, levels = sort(levels(ho.sub)))

cols.cells <- c(cols.adlunglab[c(4,8)],rep("#000000",4),cols.adlunglab[6],
                "#000000", cols.adlunglab[10], rep("#000000",4),
                cols.adlunglab[c(5,9)])


# Data Plot ####
datasmooth <- ggplot(data= plot,color="#220589",fill="#220589", aes(x = time, y = n/100)) +
  labs(title  = NULL, # plot title
       x      = "Time (days)" , # x-axis label
       y      = "Cell number (1e2)") +
  geom_point(size=0.5, aes(fill=celltype), shape=21, color="#000000", 
             stroke=0.05) +
  geom_smooth(linewidth=0.5, aes(color=celltype, fill=celltype), 
              method="loess") +
  scale_fill_manual(values = color_assignment,name=NULL) +
  scale_color_manual(values = color_assignment,name=NULL) +
  # scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks=unique(plot$time)) +
  theme_adlunglab(base_size = 6) +
  facet_wrap(~celltype,nrow=3, scales = "free") +
  expand_limits(y=0) +
  theme(legend.position = "none", strip.text.x = element_text(size = 6)#,
        #axis.text.x = element_text(angle=45, hjust=1)#
        #axis.text = element_text(color="#000000", face="bold")
  )

datasmooth


#plot

ho.sum <- plot %>%
  group_by(celltype, time) %>%
  summarize(value = mean(n, na.rm = TRUE),
            sigma = sd(n, na.rm = TRUE))

# replace NAs in sigma with actual value
ho.sum$sigma <- ifelse(is.na(ho.sum$sigma), ho.sum$value, ho.sum$sigma)

saveRDS(ho.sum, "Output/Ho_DataNew_01.Rds")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fitting ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ho.sum <- readRDS("Output/Ho_DataNew_01.Rds")

#plot <- readRDS("FrequencyTable.Rds")
check <- as.data.frame(ho.sum[which(ho.sum$celltype %in% 
                                      c("Bcell", "Mac", "Neutr", 
                                        "Tcell", "Epi")),])
check$celltype <- as.character(check$celltype)

library(dplyr)
require(tidyr)
ag <- check %>% 
  rename(name=celltype) # %>% rename(value=freq)

library(dMod)
library(cOde)

ag$value <- ag$value/100
ag$sigma <- ag$sigma/100

data <- datalist(
  Ho = as.data.frame(ag))

timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))


{
  f <- NULL
  
  # B cells ####                      
  f <- addReaction(f, from = "", to = "Bcell", 
                   rate = "healedMucosa*kturnB*Bcell",
                   description = "B_Proliferation")
  f <- addReaction(f, from = "Bcell", to = "", 
                   rate = "(1/kturnB)*Bcell",
                   description = "B_Death")
  
  # Macrophages ####                               
  f <- addReaction(f, from = "", to = "Mac", 
                   rate = "(damage)*kturnM*Mac",
                   description = "Mac_Proliferaton")
  f <- addReaction(f, from = "Mac", to = "",
                   rate = "(1/kturnM)*Mac",
                   description = "Mac_Death")
  
  # Neutrophils ####                                 
  f <- addReaction(f, from = "", to = "Neutr", 
                   rate = "DSS*kturnN*Neutr",
                   description = "N_Activation")
  f <- addReaction(f, from = "Neutr", to = "", 
                   rate = "(1/kturnN)*Neutr",
                   description = "Neutr_Death")
  
  # T cells ####
  f <- addReaction(f, from = "", to = "Tcell", 
                   rate = "(damage)*kturnT*Tcell",
                   description = "T_Proliferation")
  f <- addReaction(f, from = "Tcell", to = "",
                   rate = "(1/kturnT)*Tcell",
                   description = "T_Death")
  
  # # Locals ####                                  
  f <- addReaction(f, from = "", to = "damage",
                   rate = "DSS",
                   description = "Delay_reaction")
  
  f <- addReaction(f, from = "damage", to = "healedMucosa",
                   rate = "kheal*damage",
                   description = "Healing")
  
  f <- addReaction(f, from = "Epi", to = "",
                   rate = "DSS*(1/kturnE)*Epi",
                   description = "Epi_Death")
  
  f <- addReaction(f, from = "", to = "Epi",
                   rate = "kturnE*Epi",
                   description = "Epi_Proliferation")
  
  
  # + DSS  ####                                
  f <- addReaction(f, from = "DSS", to = "2*DSS",
                   rate = "0",
                   description = "DSS") 
  
  print(f)             
  e <- NULL
  e <- eventlist() %>% 
    addEvent(var = "DSS", time = 6, value = 0)
  
}

dmod_pipe <- function(f, e, DSSt0=DSSt0) { 
  model <- odemodel(f, modelname = "Best")
  x <- Xs(model, events = e)
  observables <- eqnvec(DSS = "DSS", damage= "damage", healedMucosa = "healedMucosa" ,
                        Bcell = "Bcell", Epi = "Epi", #Stromal = "Stromal", 
                        Mac = "Mac", Tcell = "Tcell", Neutr="Neutr")
  g <- Y(observables, x, compile = TRUE, modelname = "Best", attach.input = FALSE)
  innerpars <- getParameters(g*x)
  trafo <- repar("x~x", x = innerpars)
  trafo <- repar(paste0("x~",DSSt0), x = c("DSS"), trafo)
  trafo <- repar("x~0", x = c("healedMucosa"), trafo)
  trafo <- repar("x~0", x = c("damage"), trafo)
  # r <- readRDS("r_Frede_i_2_2.Rdata")
  # 
  # for (i in c("kturnB","kturnM","kturnN","kturnT")) {
  #   trafo[i] <- exp(as.numeric(r$params[1,i]))
  # }
  trafo <- repar("x~exp(x)", x = innerpars, trafo)
  trafo1 <- trafo
  p <- P(trafo1, condition = "Ho")
  outerpars <- getParameters(p)
  pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)
  timesM <- seq(0,15,0.1)
  times <- 0:15
  
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
  sp <- as.data.frame(simuplot$Ho)
  
  dss <- ggplot(sp[,c(1,2)]) +
    geom_line(aes(x=time, y=DSS*100), color="#000000", linewidth=2) +
    labs(title= "Model input", x= "Time (days)", y="DSS (%)") +
    theme_adlunglab(base_size = 10) +
    scale_x_continuous(breaks=c(0,7,14)) +
    theme(legend.position = "none")
  
  # remove DSS
  # sp <- sp[,-c(2)]
  sp <- sp[,-c(2,3,4)]
  
  sm <- melt(sp, id.vars = "time")
  sm$alpha <- "Model"
  dm <- as.data.frame(data$Ho)
  colnames(dm)[c(1)] <- c("variable")
  dm$alpha <- "Data"
  dm$variable <- factor(dm$variable, levels = c("Bcell","Mac","Neutr","Tcell","Epi")) #sort(unique(dm$variable)))
  cols.model <- ALcols[c(8,10,6,9,4)]
  
  p_fit <- ggplot() +
    #  geom_point(data=dm, aes(x=timepoint,y=value),color=met.brewer("Kandinsky",2)[2]) +
    geom_pointrange(data=dm, aes(x=time,ymin=(value-sigma), color=variable,
                                 ymax=(value+sigma), y=value, alpha=factor(alpha))) +
    geom_line(data=sm, aes(x=time, y=value, alpha=factor(alpha),
                           color=variable), linewidth=2) +
    labs(title= NULL, x= "Time (days)", y="Cell number (1e2)") +
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
  n <- nrow(data$Ho)
  L <- fitResults[1]$value
  
  cAIC <- fitResults[1]$value + k * 2 + (2*k*(k+1))/(n-k-1) #\frac{2k(k+1)}{n-k-1}
  bestfitpars <- fitResults[1]
  
  
  bestfit <- myfit$argument
  profiles <- profile(obj, bestfit, names(bestfit), limits = c(-10, 10), cores = 4)
  
  
  simu <- as.data.frame((g*x*p)(timesD, pouter)$Ho)
  
  # Take a look at each parameter
  p_param <- plotProfile(profiles) + 
    scale_color_manual(values="#000000", guide=F) +
    scale_linetype_manual(values=c("dotted","longdash","solid"), name=NULL) +
    labs(title= NULL, x= "log10(parameter value)") +
    #  facet_wrap(nrow = 1) +
    theme_adlunglab(base_size = 10) +
    theme(legend.position = "right")
  plots <- list(init = p_init, fit = p_fit, param = p_param, waterfall = p_waterfall)
  return(list(dss = dss, aic = cAIC, plots = plots, 
              params = bestfitpars, profiles = profiles,
              k = k, n = n, L = L, simu = sp, f=f))  
}


for(i in 1:20){
  try({
    suppressMessages({suppressWarnings({r <- dmod_pipe(f, e, DSSt0=0.015)})})
    break #break/exit the for-loop
  }, silent = FALSE)
  message("Error! Re-Trying! #Try:", i)
}

r$plots$fit
r$plots$param
r$plots$waterfall
r$params
r$aic
r$L
r$n
r$k

rprint <- r

saveRDS(r,paste0("Output/r_Ho_fit_i6.Rdata"))



# Observed vs. Predicted plot ####

r <- readRDS("Output/r_Ho_fit_i6.Rdata")


pred <- as.data.frame(data$Ho)
simu <- r$simu
simu.m <- melt(simu[which(simu$time %in% c(0,3,6,9,12,15)),], id.vars = "time")
#simu.m <- simu %>% gather(variable, value, -time)

pred$join <- paste0(pred$time, "_", pred$name)
simu.m$join <- paste0(simu.m$time, "_", simu.m$variable)

cormat <- full_join(pred, simu.m, by="join")

cormat.plot <- cormat#[which(cormat$time.x != 0),]

ct <- cor.test(cormat.plot$value.x,cormat.plot$value.y)



c <- ggplot(data=cormat.plot, aes(x=(value.x), y=(value.y))) +
  labs(title  = paste0("Pearson's rho=", round(ct$estimate,2), "; p=", signif(ct$p.value,1)), # plot title
       y      = "Model fit cell number (1e2)" , # x-axis label
       x      = "scRNA-seq data cell number (1e2)") + # y-axis label
  geom_smooth(method = "lm", color="#000000", se = T, size=0) +
  geom_abline(slope=1,intercept = 0, size=0.5, color="#000000", linetype="dashed") +
  geom_point(aes(fill=name), size=2,
             color="#000000", stroke=0.2, shape=21) +
  scale_fill_manual(values = ALcols[c(8,4,10,6,9)], name=NULL) +
  #coord_cartesian(expand = F) +
  # xlim(c(0,3)) + ylim(c(0,3)) +
  #  coord_fixed(expand = F) +
  #  expand_limits(y=0) +
  theme_adlunglab(base_size = 10) +
  theme(legend.key.size = unit(0.42, 'lines'), legend.position = c(0.2,0.8)) # +
# guides(fill = guide_legend(override.aes = list(size = 2)))                                                                                      color = c("#000000", "#000000"))))

c

fw <- 70
presdir <- "~/Output/"
ggsave(filename = paste0(presdir,"FitvsData_Ho_", Sys.Date()  ,".pdf"), 
       width=fw, height=fw, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FitvsData_Ho_", Sys.Date()  ,".png"), 
       width=fw, height=fw, units="mm", dpi=300)


ples <- r$plots$param + theme_adlunglab(base_size = 10)

ples

fw <- 130

ggsave(filename = paste0(presdir,"PLEs_Ho_", Sys.Date()  ,".pdf"), 
       width=fw, height=fw, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"PLEs_Ho_", Sys.Date()  ,".png"), 
       width=fw, height=fw, units="mm", dpi=300)

r$plots$fit

fw <- 200
ggsave(filename = paste0(presdir,"Fig2D_", Sys.Date()  ,".pdf"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig2D_", Sys.Date()  ,".png"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Some more plots ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


r <- readRDS("Data/r_Ho_fit_i6.Rdata")

fig2c <- r$dss + scale_x_continuous(breaks = c(0,6,15))

fig2c  

fw <- 40
presdir <- "~/Output/"
ggsave(filename = paste0(presdir,"FigS2_DSS_", Sys.Date()  ,".pdf"), 
       width=fw, height=200*9/16/2.5, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FigS2_DSS_", Sys.Date()  ,".png"), 
       width=fw, height=200*9/16/2.5, units="mm", dpi=300)

ggarrange(fig2b, fig2c, ncol = 2, labels = c("B","C"),
          widths = c(0.68,0.32))

fw <- 200-75

ggsave(filename = paste0(presdir,"Fig2BC_", Sys.Date()  ,".pdf"), 
       width=fw, height=62, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig2BC_", Sys.Date()  ,".png"), 
       width=fw, height=62, units="mm", dpi=300)





show.markers <- c("Ptprc", "Cd79a", "H2-Ab1",
                  "Fcgr1", "Cd68", "Ly6c2","Itgax",# Itgax = "Cd11c", #Fcgr1="Cd64",
                  "Ly6g", "Itgam", "S100a8", "S100a9", # Itgam = "Cd11b",
                  "Trbc2", "Thy1","Cd3e", #Thy1 =Cd90
                  "Epcam")

levels(ho.sub) <- rev(c("Bcell","Mac","Neutr","Tcell","Epi","Other"))

s2c <- DotPlot(ho.sub, features = show.markers, split.by = "celltype",
               cols = color_assignment[rev(c(2,1,11,5,7,4))]) + RotatedAxis() +
  #scale_size_continuous(name="Percent\nExpressed")
  guides(size = guide_legend(title = "Percent\nExpressed")) +
  theme_adlunglab(base_size = 10) + theme(legend.position = "right") +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_discrete(labels=rev(c("Bcell","Mac","Neutr","Tcell","Epi","Other"))) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

s2c

fw <- 200

ggsave(filename = paste0(presdir,"FigS2_BubblePlot_", Sys.Date()  ,".pdf"),
       width=fw, height=fw*9/16/2.5, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FigS2_BubblePlot_", Sys.Date()  ,".png"),
       width=fw, height=fw*9/16/2.5, units="mm", dpi=300)


Idents(ho.sub) <- ho.sub$seurat_clusters


s2a <- DimPlot(ho.sub, reduction = "umap", label = T, label.size = 3.5) + NoLegend() + NoAxes() +
  ggtitle(paste0(dim(ho.sub)[2], " cells")) + theme(title = element_text(size=6))

Idents(ho.sub) <- ho.sub$celltype

p1 <- FeaturePlot(ho.sub, features = c("Ptprc"), order = F, cols = c("grey89","#220589")) + NoAxes() +
  theme(legend.text = element_text(size=10), plot.title = element_text(size=10),
        legend.key.width = unit(3,"mm"))

p2 <- FeaturePlot(ho.sub, features = c("Epcam"), order = F, cols = c("grey89","#220589")) + NoAxes() +
  theme(legend.text = element_text(size=10), plot.title = element_text(size=10),
        legend.key.width = unit(3,"mm"))

figures2 <- ggarrange(
  ggarrange(s2a, p1, p2, nrow = 1,widths = c(1/3,1/3,1/3), labels = c("A", "B", "")),
  s2c, nrow = 2, labels = c("","C"), heights = c(1/2,1/2))

figures2
fw <- 200

ggsave(filename = paste0(presdir,"FigS2ABC_", Sys.Date()  ,".pdf"),
       width=fw, height=fw*9/16, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FigS2ABC_", Sys.Date()  ,".png"),
       width=fw, height=fw*9/16, units="mm", dpi=300)

