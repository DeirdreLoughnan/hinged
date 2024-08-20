# February 20, 2024
# Aim of this code is to test whether having a second slope better informs the synchrony model:
# In addition, have since added partial pooling on species and study-level variation
rm(list = ls()) 
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = FALSE)

library(colormap)
# library(phytools)
# library(ape)
require(rstan)
# require(caper)
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

dat <- read.csv("synchronyData.csv")

datalist <- append(list(N = nrow(dat),
                        N_grid = length(unique(dat$sp.pheno)),
                        x0 = 1980,
                        # Nstudy = length(unique(dat$studyid)),
                        y = dat$doy,
                        species = dat$pheno.fact,
                        # study = dat$study.fact,
                        x = dat$year), phypriors)


mdlHinge <- 
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

datalistSp <- list(N = nrow(dat),
  Nspp = length(unique(dat$sp.pheno)),
  x0 = 1980,
  Nstudy = length(unique(dat$studyid)),
  y = dat$doy,
  species = dat$pheno.fact,
  # study = dat$study.fact,
  x = dat$year)

mdlSp <-
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

datalistStudy <- list(N = nrow(dat),
  Nspp = length(unique(dat$sp.pheno)),
  x0 = 1980,
  Nstudy = length(unique(dat$studyid)),
  y = dat$doy,
  species = dat$pheno.fact,
  study = dat$study.fact,
  x = dat$year)

mdlStudy <- 
  stan("Stan/stan_programs/fit_flex_hinge_sp_study.stan",
    data = datalistStudy,
    iter = 4000,
    warmup = 3000,
    chains = 4,
    #seed = 62921,
    refresh = 20
  )

save(mdlStudy, file = "..//hinged/analyses/output/hingeSpeciesStudyYpred.Rda")


##########################################################################################
# mikes utility plots and diagnostics
util <- new.env()
source('..//hinged/analyses/stan_utility_rstan.R', local=util)

diagnostics <- util$extract_hmc_diagnostics(mdlStudy)
util$check_all_hmc_diagnostics(diagnostics)
# All Hamiltonian Monte Carlo diagnostics are consistent with reliable Markov chain Monte Carlo.

samples <- util$extract_expectands(mdlStudy)
util$check_all_expectand_diagnostics(samples)
# All expectands checked appear to be behaving well enough for reliable Markov chain Monte Carlo estimation.

# Retrodictive checks
hist_retro <- function(obs, samples, pred_names,
  bin_min, bin_max, delta,
  xlab="", display_ylim=NULL, title="") {
  if (is.na(bin_min)) bin_min <- min(pred)
  if (is.na(bin_max)) bin_max <- max(pred)
  breaks <- seq(bin_min, bin_max, delta)
  B <- length(breaks) - 1
  idx <- rep(1:B, each=2)
  xs <- sapply(1:length(idx),
    function(b) if(b %% 2 == 0) breaks[idx[b] + 1]
    else                        breaks[idx[b]] )
  obs_counts <- hist(obs[bin_min < obs & obs < bin_max], breaks=breaks, plot=FALSE)$counts
  pad_obs_counts <- do.call(cbind,
    lapply(idx, function(n) obs_counts[n]))
  pred <- sapply(pred_names,
    function(name) c(t(samples[[name]]), recursive=TRUE))
  N <- dim(pred)[1]
  pred_counts <- sapply(1:N,
    function(n) hist(pred[n,][bin_min < pred[n,] & pred[n,] < bin_max],
      breaks=breaks,
      plot=FALSE)$counts)
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B,
    function(b) quantile(pred_counts[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
  
  if (is.null(display_ylim)) {
    display_ylim <- c(0, max(c(obs_counts, cred[9,])))
  }
  
  plot(1, type="n", main=title,
    xlim=c(bin_min, bin_max), xlab=xlab,
    ylim=display_ylim, ylab="Counts")
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
    col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
    col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
    col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
    col = c_mid_highlight, border = NA)
  lines(xs, pad_cred[5,], col=c_dark, lwd=2)
  lines(xs, pad_obs_counts, col="white", lty=1, lw=2.5)
  lines(xs, pad_obs_counts, col="black", lty=1, lw=2)
}

plot_cont_marginal_quantiles <- function(xs, preds, 
  display_xlim=NULL, display_ylim=NULL, 
  title="", x_name="", y_name="") {
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(preds, function(pred) quantile(c(t(pred), recursive=TRUE), probs=probs))
  
  if (is.null(display_xlim)) {
    display_xlim <- range(xs)
  }
  
  if (is.null(display_ylim)) {
    display_ylim <- c(min(cred[1,]), max(cred[9,]))
  }
  
  plot(1, type="n", main=title,
    xlim=display_xlim, xlab=x_name,
    ylim=display_ylim, ylab=y_name)
  
  polygon(c(xs, rev(xs)), c(cred[1,], rev(cred[9,])),
    col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(cred[2,], rev(cred[8,])),
    col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(cred[3,], rev(cred[7,])),
    col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(cred[4,], rev(cred[6,])),
    col = c_mid_highlight, border = NA)
  lines(xs, cred[5,], col=c_dark, lwd=2)
}

par(mfrow=c(2, 1), mar = c(5, 4, 2, 1))

pred_names <- grep('y_pred', names(samples), value=TRUE)
hist_retro(dat$doy, samples, pred_names, -3.5, 1, 0.25, "y")

pred_names <- grep('y_grid_pred', names(samples), value=TRUE)
preds <- samples[pred_names]
plot_cont_marginal_quantiles(data$x_grid, preds, 
  title="Conditional Check", 
  x_name="x", y_name="y",
  display_ylim=c(-3, 1))
points(data$x, data$y, col="white", pch=16, cex=1.2)
points(data$x, data$y, col="black", pch=16, cex=0.8)

# Marginal posterior distributions
par(mfrow=c(1, 2), mar = c(5, 4, 2, 1)) 
util$plot_expectand_pushforward(samples[["alpha"]], 25,
  display_name="alpha")
util$plot_expectand_pushforward(samples[["sigma"]], 25,
  display_name="sigma")

##########################################################################################
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

mdlOut <- data.frame(speciesPheno = unique(sort(dat$sp.pheno)), 
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

dat <- subset(dat, sp.pheno == "Acer_campestre_flowering")

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
dat <- subset(dat, sp.pheno == spPheno)

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
