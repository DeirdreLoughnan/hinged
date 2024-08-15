# February 20, 2024
# Aim of this code is to test whether having a second slope better informs the synchrony model:

# This code was started June 22
# The aim is to load the phylogeny and change the species names to the class names
rm(list = ls()) 
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = FALSE)

#source('stan_utility_rstan.R', local=util)

library(colormap)
library(phytools)
library(ape)
require(rstan)
require(caper)
require(shinystan)
require(reshape2)
library(stringr)
library(ggplot2)
library(plyr)
library(dplyr)

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/Synchrony")
} else if(length(grep("Lizzie", getwd()) > 0)) {
  setwd("~/Documents/git/projects/others/deirdre/Synchrony")
} else{
  setwd("/home/deirdre/Synchrony") # for midge
}

#setwd("~/Documents/github/hinged")
##################################################################
# Taking the phylogney code from OSPREE and trying to adapt it for my synchrony project:
#get the synchrony data
#source("Rcode/combiningallphenodata.R")
dat.fin <- read.csv("Input/dat_fin_Dec2021_w_neg.csv")
most <- dat.fin[,c("year", "doy", "species", "phenophase", "studyid", "datasource", "yr1980","species.name")]
final.t <- most[complete.cases(most), ]  # there are five rows of the seabird data that doesn't 
final.t$species_fact <- as.numeric(as.factor(final.t$species.name))
final.t$study.fact <- as.numeric(as.factor(final.t$studyid))

#make a matching species name and spname pheno column
final.t$temp <- final.t$species.name
temp <- str_split_fixed(final.t$temp, " ", 2)
final.t$phylo.name <- paste(temp[,1], temp[,2], sep="_")

final.t$phenophase[final.t$phenophase == "egg laying"] <- "egg_laying"
final.t$phenophase[final.t$phenophase == "first cut"] <- "first_cut"
final.t$phenophase[final.t$phenophase == "juveniles first seen"] <- "juveniles_first_seen"
final.t$phenophase[final.t$phenophase == "population growth"] <- "population_growth"
final.t$phenophase[final.t$phenophase == "last cut"] <- "last_cut"
final.t$phenophase[final.t$phenophase == "switch date"] <- "switch_date"
final.t$phenophase[final.t$phenophase == "first appearance"] <- "first_appearance"
final.t$phenophase[final.t$phenophase == "first ripe fruit"] <- "first_ripe_fruit"
final.t$phenophase[final.t$phenophase == "gathering for departure"] <- "gathering_for_departure"
final.t$phenophase[final.t$phenophase == "last appearance"] <- "last_appearance"
final.t$phenophase[final.t$phenophase == "peak abundnace"] <- "abundance"
final.t$phenophase[final.t$phenophase == "reproduction" & final.t$species.name == "Ptychoramphus aleuticus"] <- "egg_laying"

final.t$phylo.name[final.t$phylo.name == "Hyacinthoides_x variabilis"] <- "Hyacinthoides_x_variabilis"
final.t$phylo.name[final.t$phylo.name == "Adelgidae_and Phylloxeridae spp"] <- "Adelgidae_and_Phylloxeridae_spp"
final.t$phylo.name[final.t$phylo.name == "Taraxacum_sect. Erythrosperma"] <- "Taraxacum_sect._Erythrosperma"
final.t$phylo.name[final.t$phylo.name == "Daucus_carota carota"] <- "Daucus_carota_carota"
final.t$phylo.name[final.t$phylo.name == "Polygonum_aviculare ss"] <- "Polygonum_aviculare_ss"
final.t$phylo.name[final.t$phylo.name == "Tilia_x europaea"] <- "Tilia_x_europaea"
final.t$phylo.name[final.t$phylo.name == "Populus_x canescens"] <- "Populus_x_canescens"
final.t$phylo.name[final.t$phylo.name == "Symphytum_x uplandicum"] <- "Symphytum_x_uplandicum"

final.t$sp.pheno <- paste(final.t$phylo.name, final.t$phenophase, sep = "_")
final.t$sp.pheno <- str_replace(final.t$sp.pheno, "__", "_")
final.t$sp.pheno <- str_replace(final.t$sp.pheno, " ", "")

# tree <- read.tree("Input/synchronyTree.tre")
# 
# #get spp and phylo matrix in the same order
# head(tree$tip.label)
# length(tree$tip.label) #1280
# 
# treeSpp <- tree$tip.label#1279 - bombus_first_appearance was duplicated
# datSpp <- unique(final.t$sp.pheno)#1279
# # Full phylogeny
# d <- final.t[match(tree$tip.label, final.t$phylo.name),]
# d$sp.pheno <- str_replace(d$sp.pheno, "__", "_")
# d$sp.pheno <- str_replace(d$sp.pheno, " ", "")
# 
# phymatch <- data.frame(sp.pheno = tree$tip.label, sppnum = c(1:length(tree$tip.label)))
# phymatch$sp.pheno <- str_replace(phymatch$sp.pheno, "__", "_")
# phymatch$sp.pheno <- str_replace(phymatch$sp.pheno, " ", "")
# 
# d <- merge(final.t, phymatch, by="sp.pheno")
# length(unique(final.t$sp.pheno))
# 
# temp <- phymatch[is.na(match(phymatch$sp.pheno,final.t$sp.pheno)),]; dim(temp)
# temp <- final.t[is.na(match(final.t$sp.pheno,tree$tip.label)),]; dim(temp)
# 
# d <- d[order(d$sppnum),]
# nspecies <- max(d$sppnum)
# #nspecies <- 1275
# cophen_tree <- cophenetic(tree)
# vcv_tree <- vcv(tree, cor = TRUE)
# 
# final.t$pheno.fact <- as.numeric(as.factor(final.t$sp.pheno))
# length(unique(final.t$sp.pheno)) #1275 w/o thackeray 789




datalist <- append(list(N = nrow(final.t),
                        N_grid = nspecies,
                        x0 = 1980,
                        # Nstudy = length(unique(final.t$studyid)),
                        y = final.t$doy,
                        species = final.t$pheno.fact,
                        # study = final.t$study.fact,
                        x = final.t$year), phypriors)


mdlHinge <- #stan("Stan/phylogeny_synchrony_cholesky_grandmean.stan",
  stan("Stan/stan_programs/fit_flex_hinge_dl.stan",
       data = datalist,
       init = simu_inits,
       iter = 4000,
       warmup = 3000,
       chains = 4,
       #seed = 62921,
       refresh = 20
  )

save( mdlHinge, file = "output/synchronyHinge.Rda")

### Adding partial pooling across species

datalistSp <- append(list(N = nrow(final.t),
  Nspp = nspecies,
  x0 = 1980,
  # Nstudy = length(unique(final.t$studyid)),
  y = final.t$doy,
  species = final.t$pheno.fact,
  # study = final.t$study.fact,
  x = final.t$year), phypriors)

mdlSp <- #stan("Stan/phylogeny_synchrony_cholesky_grandmean.stan",
  stan("Stan/stan_programs/fit_flex_hinge_sp.stan",
    data = datalistSp,
    iter = 4000,
    warmup = 3000,
    chains = 4,
    #seed = 62921,
    #refresh = 20
  )

save( mdlSp, file = "output/synchronyHingeSpecies.Rda")

sum <- summary(mdlSp)$summary

slopes <- sum[c("mu_beta1", "mu_beta2","mu_alpha","sigma_a", "sigma_b1", "sigma_b2", "sigma"), c("mean","2.5%", "97.5%", "n_eff", "Rhat")]

slopes

#######################################################
### Adding partial pooling across species and study

datalistSp <- append(list(N = nrow(final.t),
  Nspp = nspecies,
  x0 = 1980,
  # Nstudy = length(unique(final.t$studyid)),
  y = final.t$doy,
  species = final.t$pheno.fact,
  # study = final.t$study.fact,
  x = final.t$year), phypriors)

mdlSp <- #stan("Stan/phylogeny_synchrony_cholesky_grandmean.stan",
  stan("Stan/stan_programs/fit_flex_hinge_sp.stan",
    data = datalistSp,
    iter = 4000,
    warmup = 3000,
    chains = 4,
    #seed = 62921,
    #refresh = 20
  )

save( mdlSp, file = "output/synchronyHingeSpecies.Rda")


# plot slope with the data:
preSlope <- data.frame(sum[grep("beta1\\[", rownames(sum)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

postSlope <- data.frame(sum[grep("beta2\\[", rownames(sum)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

int <- data.frame(sum[grep("alpha\\[", rownames(sum)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

col.sp <- c(rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8), rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.5))

pdf("analyses/figures/post_pre1980_hist.pdf", width = 5, height = 5)
hist(postSlope$mean*10, col = col.sp[1], main = NA, breaks = 25, xlab = "Shift in phenology (days/decade)", ylim = c(0, 1000))
hist(preSlope$mean*10, col = col.sp[2], main = NA, add = T, breaks = 10)
dev.off()

legend("topright",legend = c("pre-1980", "post-1980"),
  col = c(col.sp[1], col.sp[2]),   bty = "n", pch = 19, cex =1.5)

mdlOut <- data.frame(speciesPheno = unique(sort(final.t$sp.pheno)), 
  pre1980 = preSlope$mean, 
  post1980 = postSlope$mean,
  alpha = int$mean)

pdf("analyses/figures/post_pre1980Diff.pdf", width = 5, height = 5)
hist((mdlOut$post1980-mdlOut$pre1980), main = NA)
dev.off()

# top three biggest changes pre-post climate change are all amphibians

mdlOut$diff <- mdlOut$post1980-mdlOut$pre1980
row.names(mdlOut) <- mdlOut$speciesPheno

#######################################################
# plot the raw data and whether the model output fits

dat <- subset(final.t, sp.pheno == "Acer_campestre_flowering")

pm.0 <- mdlOut["Acer_campestre_flowering", "pre1980"] * min(dat$yr1980) + mdlOut["Acer_campestre_flowering", "alpha"]
pm.1 <- mdlOut["Acer_campestre_flowering", "pre1980"] * 0 + mdlOut["Acer_campestre_flowering", "alpha"]

pm.2 <- mdlOut["Acer_campestre_flowering", "post1980"] * 0 + mdlOut["Acer_campestre_flowering", "alpha"]
pm.3 <- mdlOut["Acer_campestre_flowering", "post1980"] * max(dat$yr1980) + mdlOut["Acer_campestre_flowering", "alpha"]

pdf("analyses/figures/Acer_campestre_2hinge.pdf", width =3, height = 3)
plot(doy~year, data = dat, type="l", ylim = c(100,200), col = "darkslategrey", ylab = "Day of year", xlab = "Year", cex = 3,lwd =2, main = paste("Acer_campestre_flowering", round(mdlOut["Acer_campestre_flowering", "pre1980"],3),  round(mdlOut["Acer_campestre_flowering", "post1980"],3), sep = "_" ))
points(doy~year, data=dat, cex=0.6, col = "darkslategrey", pch =19)
abline((lm(doy~year, data=dat)), lty =2,lwd =2)
segments(x0 = min(dat$year), x1 = 1980, y0 = pm.0, y1 = pm.1 ,lwd =2)
segments(x0 = 1980, x1 = max(dat$year), y0 = pm.2, y1 = pm.3,lwd =2)
dev.off()


### 
spPheno <- "Bupalus_piniaria_abundance"
spPheno <- "Macoma_balthica_spawning"
spPheno <- "Pagodroma_nivea_egg_laying"
spPheno <- "Pseudacris_crucifer_first_appearance"
dat <- subset(final.t, sp.pheno == spPheno)

dat <- dat[order(dat$year),]
pm.0 <- mdlOut[spPheno, "pre1980"] * min(dat$yr1980) + mdlOut[spPheno, "alpha"]
pm.1 <- mdlOut[spPheno, "pre1980"] * 0 + mdlOut[spPheno, "alpha"]

pm.2 <- mdlOut[spPheno, "post1980"] * 0 + mdlOut[spPheno, "alpha"]
pm.3 <- mdlOut[spPheno, "post1980"] * max(dat$yr1980) + mdlOut[spPheno, "alpha"]

pdf(paste("analyses/figures/", spPheno, ".pdf", sep = ""), width = 5, height = 5)
plot(doy~year, data = dat, type="l", ylim = c(0,300), col = "darkslategrey", ylab = "Day of year", xlab = "Year", cex = 3,lwd =2, main = paste(spPheno, round(mdlOut[spPheno, "pre1980"],3),  round(mdlOut[spPheno, "post1980"],3), sep = "_" ))
points(doy~year, data=dat, cex=0.6, col = "darkslategrey", pch =19)
#abline((lm(doy~year, data=dat)), lty =2,lwd =2)
segments(x0 = min(dat$year), x1 = 1980, y0 = pm.0, y1 = pm.1 ,lwd =2)
segments(x0 = 1980, x1 = max(dat$year), y0 = pm.2, y1 = pm.3,lwd =2)
dev.off()
