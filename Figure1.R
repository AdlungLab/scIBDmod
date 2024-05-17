setwd("~/")

require(dMod)
require(readxl)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)

source("~/AdlungLab_HelperFunctions.R")

# load raw data ####
frede <- as.data.frame(read_excel("Frede_FACSexported.xlsx"))

# process data ####
rownames(frede) <- frede$cells_1e6
frede <- frede[,-c(1)]

# cells_sums <- colSums(frede, na.rm=T)
# normalized_data <- sweep(frede, 2, cells_sums, FUN = "/")
# frede.frac <- normalized_data

Bcells <- colSums(frede[c("B cells"), ],na.rm = T)
Tcells <- colSums(frede[c("IELs"#,"Teff","LP T cells"
), ],na.rm = T)
Mac <- colSums(frede[c(#"DCs (CD11b)",
  "Monocytes",
  "Macrophages"), ],na.rm = T)
Neutr <- colSums(frede[c(#"Neutrophils",
  "Granulocytes"), ],na.rm = T)
frede.combined <- rbind(frede, Bcells, Tcells, Mac, Neutr)
rownames(frede.combined)[c(14,15,16,17)] <- c("Bcell","Tcell","Mac","Neutr")

frede.combined$id <- row.names(frede.combined)

# melt all time points
frede.melt <- frede.combined %>% melt(id.vars = "id") 

# extract time points
frede.melt$time <- as.numeric(substr(frede.melt$variable, 2, 3))

# remove empty column
frede.melt <- na.omit(frede.melt)


frede.sel <- frede.melt[which(frede.melt$id %in% c("Bcell","Tcell","Mac","Neutr")),]


# calculate mean and standard deviation 
# across cells and time points
frede.sum <- frede.melt %>%
  group_by(id, time) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            sigma = sd(value, na.rm = TRUE)) %>% 
  rename(value=mean_value)

# # remove or replace all characters that could be problematic
# frede.sum$name <- gsub("\\s+", "", frede.sum$id)
# frede.sum$name <- gsub("/","", frede.sum$name)
# frede.sum$name <- gsub("\\(","_", frede.sum$name)
# frede.sum$name <- gsub("\\)","", frede.sum$name)

frede.fit <- frede.sum[which(frede.sum$id %in% c("Bcell","Mac","Neutr","Tcell")),]
# frede.nozero <- frede.fit[which(#frede.fit$time==0 | 
#   frede.fit$time==14),]

colnames(frede.fit)[1] <- "name"

# save processed data ####
#saveRDS(frede.fit,"Data/scMod_Frede_Fit_Data.rds")


# load processed data ####
frede.fit <- readRDS("Data/scMod_Frede_Fit_Data.rds")

data <- datalist(
  Frede = as.data.frame(frede.fit))

timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))



# dmod_pipe ####

dmod_pipe <- function(f, e, DSSt0=DSSt0) { 
  model <- odemodel(f, modelname = "Best")
  x <- Xs(model, events = e)
  observables <- eqnvec(DSS = "DSS", damage= "damage", healedMucosa = "healedMucosa",
                        Bcell = "Bcell", #Epi = "Epi", Stromal = "Stromal", 
                        Mac = "Mac", Tcell = "Tcell", Neutr="Neutr")
  g <- Y(observables, x, compile = TRUE, modelname = "Best", attach.input = FALSE)
  innerpars <- getParameters(g*x)
  trafo <- repar("x~x", x = innerpars)
  trafo <- repar(paste0("x~",DSSt0), x = c("DSS"), trafo)
  trafo <- repar("x~0", x = c("healedMucosa"), trafo)
    trafo <- repar("x~0", x = c("damage"), trafo)
  trafo <- repar("x~exp(x)", x = innerpars, trafo)
  trafo1 <- trafo
  p <- P(trafo1, condition = "Frede")
  outerpars <- getParameters(p)
  pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)
  timesM <- seq(0,14,0.1)
  times <- 0:14
  
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
  sp <- as.data.frame(simuplot$Frede)
  
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
  dm <- as.data.frame(data$Frede)
  colnames(dm)[c(1)] <- c("variable")
  dm$alpha <- "Data"
  dm$variable <- factor(dm$variable, levels = sort(unique(dm$variable)))
  cols.model <- ALcols[c(8,10,6,9)]
  
  p_fit <- ggplot() +
    #  geom_point(data=dm, aes(x=timepoint,y=value),color=met.brewer("Kandinsky",2)[2]) +
    geom_pointrange(data=dm, aes(x=time,ymin=(value-sigma), color=variable,
                                 ymax=(value+sigma), y=value, alpha=factor(alpha))) +
    geom_line(data=sm, aes(x=time, y=value, alpha=factor(alpha),
                           color=variable), linewidth=2) +
    labs(title= NULL, x= "Time (days)", y="Cell number (1e6)") +
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
  n <- nrow(data$Frede)
  L <- fitResults[1]$value
  
  cAIC <- fitResults[1]$value + k * 2 + (2*k*(k+1))/(n-k-1) #\frac{2k(k+1)}{n-k-1}
  bestfitpars <- fitResults[1]
  
  
  bestfit <- myfit$argument
  profiles <- profile(obj, bestfit, names(bestfit), limits = c(-10, 10), cores = 4)
  
  
  simu <- as.data.frame((g*x*p)(timesD, pouter)$Frede)
  
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


# Reactions ####
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
  # f <- addReaction(f, from = "Stromal", to = "",
  #                  rate = "(1/ktS)*Stromal",
  #                  description = "Stromal_Death")
  # 
  f <- addReaction(f, from = "", to = "damage",
                   rate = "DSS",
                   description = "Delay_reaction")
  # 
  # f <- addReaction(f, from = "", to = "Epi",
  #                  rate = "ktE*Epi",
  #                  description = "Epi_Proliferation")
  # 
  f <- addReaction(f, from = "damage", to = "healedMucosa",
                   rate = "kheal*damage",
                   description = "Healing")
  # 
  # f <- addReaction(f, from = "", to = "Stromal",
  #                  rate = "ktS*Stromal",
  #                  description = "Stromal_Proliferation")
  
  
  # + DSS  ####                                
  f <- addReaction(f, from = "DSS", to = "2*DSS",
                   rate = "0",
                   description = "DSS") 
  
  print(f)             
  e <- NULL
  e <- eventlist() %>% 
    addEvent(var = "DSS", time = 7, value = 0)
  
}


#Try til it works
for(i in 1:20){
  try({
    suppressMessages({suppressWarnings({r <- dmod_pipe(f, e, DSSt0=0.025)})})
    break #break/exit the for-loop
  }, silent = FALSE)
  message("Error! Re-Trying! #Try:", i)
}

r$plots$fit
r$plots$param
r$aic
r$L
r$plots$waterfall
r$params
r$n
r$k

rprint <- r

saveRDS(r,paste0("r_Frede_i6.Rdata"))


# Observed vs. Predicted plot ####

r <- readRDS("r_Frede_i6.Rdata")


r$plots$waterfall + theme_adlunglab(base_size = 10) + 
  labs(title=NULL, x= "Run index (sorted by likelihood", y="Objective function") +
  theme(legend.position = "bottom")

fw <- 66

presdir <- "~/Dropbox/UKE/Files/Writing/Manuscripts/scMod_Manuscript/Figures/Revision/"
ggsave(filename = paste0(presdir,"FigS1B_", Sys.Date()  ,".pdf"), 
       width=fw, height=60, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FigS1B_", Sys.Date()  ,".png"), 
       width=fw, height=60, units="mm", dpi=300)


pred <- as.data.frame(data$Frede)
simu <- r$simu
simu.m <- melt(simu[which(simu$time %in% c(0,3,7,10,14)),], id.vars = "time")
#simu.m <- simu %>% gather(variable, value, -time)

pred$join <- paste0(pred$time, "_", pred$name)
simu.m$join <- paste0(simu.m$time, "_", simu.m$variable)

cormat <- full_join(pred, simu.m, by="join")

cormat.plot <- cormat#[which(cormat$time.x != 0),]

ct <- cor.test(cormat.plot$value.x,cormat.plot$value.y)



c <- ggplot(data=cormat.plot, aes(x=(value.x), y=(value.y))) +
  labs(title  = paste0("Pearson's rho=", round(ct$estimate,2), "; p=", signif(ct$p.value,1)), # plot title
       y      = "Model fit cell number (1e6)" , # x-axis label
       x      = "FACS data cell number (1e6)") + # y-axis label
  geom_smooth(method = "lm", color="#000000", se = T, size=0) +
  geom_abline(slope=1,intercept = 0, size=0.5, color="#000000", linetype="dashed") +
  geom_point(aes(fill=name), size=2,
             color="#000000", stroke=0.2, shape=21) +
  scale_fill_manual(values = ALcols[c(8,10,6,9)], name=NULL) +
  #coord_cartesian(expand = F) +
    # xlim(c(0,3)) + ylim(c(0,3)) +
  #  coord_fixed(expand = F) +
  #  expand_limits(y=0) +
  theme_adlunglab(base_size = 10) +
  theme(legend.key.size = unit(0.42, 'lines'), legend.position = c(0.2,0.8)) # +
# guides(fill = guide_legend(override.aes = list(size = 2)))                                                                                      color = c("#000000", "#000000"))))

c


# Model selection plot ####

AIC <- read_excel("SupplementaryTable1.xlsx")

df <- AIC %>%
  filter(!is.na(cAIC)) %>%
  arrange((cAIC)) %>%
  mutate(i=factor(i, levels=i))

p3 <- ggplot(df, aes(x=i, y=(cAIC))) +
  geom_bar(stat = "identity", fill="#220589", width=0.89) +
  labs(title= NULL, x= "Model index", y="corrected\nAkaike Information Criterion") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  #  scale_fill_manual(name=NULL, values = c("#220589")) +
#  coord_flip(ylim = c(0,155), xlim=c(1,15)) +
  expand_limits(y=0, x=0) +
  theme_adlunglab(base_size = 10)  +
  theme(legend.position = "bottom")

p3

empty_plot <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "white") +
  theme_void()

ples <- r$plots$param + theme_adlunglab(base_size = 10)

require(ggpubr)
figureS1 <- ggarrange(
  ggarrange(p3, c, ncol = 2, widths = c(2/3,1/3), labels=c("A","B")), 
  ggarrange(ples,empty_plot, ncol = 2, widths = c(5/7,2/7), labels=c("C","")), 
  nrow = 2, heights = c(2/5,3/5), labels = c("",""))

figureS1


fw <- 200

presdir <- "~/Output"
ggsave(filename = paste0(presdir,"FigS1_ABC_", Sys.Date()  ,".pdf"), 
       width=fw, height=fw*3/4, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"FigS1_ABC_", Sys.Date()  ,".png"), 
       width=fw, height=fw*3/4, units="mm", dpi=300)


# Plot best fit ####

figure1 <- ggarrange(r$dss, r$plots$fit, ncol = 2, labels = c("C","D"), widths = c(1/5,4/5)) 

figure1

ggsave(filename = paste0(presdir,"Fig1CD_", Sys.Date()  ,".pdf"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300, useDingbats=FALSE)
ggsave(filename = paste0(presdir,"Fig1CD_", Sys.Date()  ,".png"), 
       width=fw, height=fw*9/16/2, units="mm", dpi=300)



