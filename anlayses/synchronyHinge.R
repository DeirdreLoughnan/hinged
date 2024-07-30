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

tree <- read.tree("Input/synchronyTree.tre")

#get spp and phylo matrix in the same order
head(tree$tip.label)
length(tree$tip.label) #1280

treeSpp <- tree$tip.label#1279 - bombus_first_appearance was duplicated
datSpp <- unique(final.t$sp.pheno)#1279
# Full phylogeny
d <- final.t[match(tree$tip.label, final.t$phylo.name),]
d$sp.pheno <- str_replace(d$sp.pheno, "__", "_")
d$sp.pheno <- str_replace(d$sp.pheno, " ", "")

phymatch <- data.frame(sp.pheno = tree$tip.label, sppnum = c(1:length(tree$tip.label)))
phymatch$sp.pheno <- str_replace(phymatch$sp.pheno, "__", "_")
phymatch$sp.pheno <- str_replace(phymatch$sp.pheno, " ", "")

d <- merge(final.t, phymatch, by="sp.pheno")
length(unique(final.t$sp.pheno))

temp <- phymatch[is.na(match(phymatch$sp.pheno,final.t$sp.pheno)),]; dim(temp)
temp <- final.t[is.na(match(final.t$sp.pheno,tree$tip.label)),]; dim(temp)

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 1275
cophen_tree <- cophenetic(tree)
vcv_tree <- vcv(tree, cor = TRUE)

final.t$pheno.fact <- as.numeric(as.factor(final.t$sp.pheno))
length(unique(final.t$sp.pheno)) #1275 w/o thackeray 789

phypriors <- list(
  a_z_prior_mu = 150,
  a_z_prior_sigma = 50,
  astudy_prior_mu = 0,
  astudy_prior_sigma = 50,
  lam_interceptsa_prior_alpha = 1, #
  lam_interceptsa_prior_beta = 1, #
  sigma_interceptsa_prior_mu = 40,
  sigma_interceptsa_prior_sigma = 20,
  b_z_prior_mu = 0,
  b_z_prior_sigma = 10,
  lam_interceptsb_prior_alpha = 2, 
  lam_interceptsb_prior_beta = 5, #
  sigma_interceptsb_prior_mu = 5,
  sigma_interceptsb_prior_sigma = 5,
  sigma_y_prior_mu = 20,
  sigma_y_prior_sigma = 10
)


# Function for generating "good" initial values
simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  b_z.temp <- rnorm(n = nspecies, mean = phypriors[["b_z_prior_mu"]], sd = phypriors[["b_z_prior_sigma"]])
  return(append(list(
    a = a_z.temp,
    b_year = b_z.temp),
    phypriors))
}


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

# plot slope with the data:
preSlope <- data.frame(sum[grep("beta1\\[", rownames(sum)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

postSlope <- data.frame(sum[grep("beta2\\[", rownames(sum)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

col.sp <- c(rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8), rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.5))

hist(postSlope$mean*10, col = col.sp[1], main = NA, breaks = 25, xlab = "Shift in phenology (days/decade)", ylim = c(0, 1000))
hist(preSlope$mean*10, col = col.sp[2], main = NA, add = T, breaks = 10)

legend("topright",legend = c("pre-1980", "post-1980"),
  col = c(col.sp[1], col.sp[2]),   bty = "n", pch = 19, cex =1.5)

mdlOut <- data.frame(speciesPheno = unique(sort(final.t$sp.pheno)), 
  pre1980 = preSlope$mean, 
  post1980 = postSlope$mean)

hist((mdlOut$post1980-mdlOut$pre1980), main = NA)
# top three biggest changes pre-post climate change are all amphibians

mdlOut$diff <- mdlOut$post1980-mdlOut$pre1980
