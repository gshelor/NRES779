##### NRES 779 Bayesian Modeling Final Project #####
## script by Griffin Shelor

### This script will be used for building and evaluating models of Snow Water Equivalent (SWE) at SNOTEL sites
## data to be used was cleaned in SNOTELDaymetClean.R
## script by Griffin Shelor

##### Loading Packages, setting up script to use parallel package, reading in data #####
library(pacman)
p_load(snotelr, here, tidyverse, rstan, daymetr, parallel, isotracer, LaplacesDemon, rjags, loo, gt, gtExtras, RColorBrewer, ggpubr, webshot2, pals)

## setting script to maximise use of cores
options(mc.cores = detectCores())

## reading in csv and filtering to only years post 2000 for a specific site because any other model would 
SnotelDaymet <- read_csv(here("Project", "Data", "SNOTEL", "SnotelDaymetClean.csv")) |>
  filter(Year >= 2010) |>
  filter(state == "CA") |>
  filter(elev == min(elev) | elev == max(elev)) |>
  mutate(swe_plus = snow_water_equivalent + 1) |>
  arrange(date, site_id)

## moving data down a row so that it is easier to work into a model statement and Stan data list
SnotelDaymet$prevday_swe <- -999
SnotelDaymet$prevday_maxtemp <- -999
SnotelDaymet$prevday_mintemp <- -999
SnotelDaymet$prevday_meantemp <- -999
SnotelDaymet$prevday_precip <- -999
SnotelDaymet$prevday_cumulative_precip <- -999

SnotelDaymet574 <- SnotelDaymet |>
  filter(site_id == 574) |>
  arrange(date)
SnotelDaymet778 <- SnotelDaymet |>
  filter(site_id == 778) |>
  arrange(date)

for (i in 2:nrow(SnotelDaymet574)){
  SnotelDaymet574$prevday_swe[i] = SnotelDaymet574$snow_water_equivalent[i-1]
  SnotelDaymet574$prevday_maxtemp[i] <- SnotelDaymet574$daymet_maxtemp[i-1]
  SnotelDaymet574$prevday_mintemp[i] <- SnotelDaymet574$daymet_mintemp[i-1]
  SnotelDaymet574$prevday_meantemp[i] <- SnotelDaymet574$daymet_meantemp[i-1]
  SnotelDaymet574$prevday_precip[i] <- SnotelDaymet574$daymet_precip[i-1]
  SnotelDaymet574$prevday_cumulative_precip[i] <- SnotelDaymet574$precipitation_cumulative[i-1]
}
for (i in 2:nrow(SnotelDaymet778)){
  SnotelDaymet778$prevday_swe[i] = SnotelDaymet778$snow_water_equivalent[i-1]
  SnotelDaymet778$prevday_maxtemp[i] <- SnotelDaymet778$daymet_maxtemp[i-1]
  SnotelDaymet778$prevday_mintemp[i] <- SnotelDaymet778$daymet_mintemp[i-1]
  SnotelDaymet778$prevday_meantemp[i] <- SnotelDaymet778$daymet_meantemp[i-1]
  SnotelDaymet778$prevday_precip[i] <- SnotelDaymet778$daymet_precip[i-1]
  SnotelDaymet778$prevday_cumulative_precip[i] <- SnotelDaymet778$precipitation_cumulative[i-1]
}

## cutting off top row of each site subset since it does not contain true previous day information
SnotelDaymet574 <- SnotelDaymet574[2:nrow(SnotelDaymet574),]
SnotelDaymet778 <- SnotelDaymet778[2:nrow(SnotelDaymet778),]

## merging them back together
SnotelDaymet <- rbind(SnotelDaymet574, SnotelDaymet778) |>
  arrange(date, site_id)




##### Stan Models #####
##### Stan Model 1 #####
## making list of data to declare what goes into stan model
model1_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$elev, x2 = SnotelDaymet$prevday_precip)

## fitting stan model
set.seed(802)
model1_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel1.stan"),data = model1_datalist, chains = 3, iter = 25000, warmup = 5000)
model1_fit


## Extracting Parameters
model1_pars <- rstan::extract(model1_fit, c("b0", "b1", "b2", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model1_loglik <- extract_log_lik(model1_fit)

## performing leave-one-out CV
model1_loo <- loo(model1_loglik)
model1_loo$estimates

## looking at waic
model1WAIC <- waic(model1_loglik)
model1WAIC


##### Stan Model 2 #####
## making list of data to declare what goes into stan model
model2_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$prevday_cumulative_precip)

## fitting stan model
set.seed(802)
model2_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel2.stan"),data = model2_datalist, chains = 3, iter = 25000, warmup = 5000)
model2_fit


## Extracting Parameters
model2_pars <- rstan::extract(model2_fit, c("b0", "b1", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model2_loglik <- extract_log_lik(model2_fit)

## performing leave-one-out CV
model2_loo <- loo(model2_loglik)
model2_loo$estimates

## looking at waic
model2WAIC <- waic(model2_loglik)
model2WAIC


##### Stan Model 3 #####
## making list of data to declare what goes into stan model
model3_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$prevday_maxtemp)

## fitting stan model
set.seed(802)
model3_fit <- stan(file=here("Project", "Scripts", "Stan","ProjectModel3.stan"), data = model3_datalist, chains = 3, iter = 25000, warmup = 5000)
model3_fit

rstan::traceplot(model3_fit, pars = c("b0", "b1", "sigma"))


## Extracting Parameters
model3_pars <- rstan::extract(model3_fit, c("b0", "b1", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model3_loglik <- extract_log_lik(model3_fit)

## performing leave-one-out CV
model3_loo <- loo(model3_loglik)
model3_loo$estimates

## looking at waic
model3WAIC <- waic(model3_loglik)
model3WAIC

## making plots of posterior samples
par(mfrow = c(3,1))
hist(model3_pars$b0)
hist(model3_pars$b1)
hist(model3_pars$sigma)


##### Stan Model 4 #####
## making list of data to declare what goes into stan model
model4_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$prevday_mintemp)

## fitting stan model
set.seed(802)
model4_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel4.stan"),data = model4_datalist, chains = 3, iter = 25000, warmup = 5000)
model4_fit


## Extracting Parameters
model4_pars <- rstan::extract(model4_fit, c("b0", "b1", "sigma", "meanlog_gt", "sdlog_gt"))


## extracting log likelihood
model4_loglik <- extract_log_lik(model4_fit)

## performing leave-one-out CV
# model4_loo <- loo(model4_loglik)
# model4_loo$estimates

## looking at waic
model4WAIC <- waic(model4_loglik)
model4WAIC



##### Stan Model 5 #####
## making list of data to declare what goes into stan model
model5_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$prevday_meantemp, x2 = SnotelDaymet$prevday_precip)

## fitting stan model
set.seed(802)
model5_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel5.stan"),data = model5_datalist, chains = 3, iter = 25000, warmup = 5000)
model5_fit

## Extracting Parameters
model5_pars <- rstan::extract(model5_fit, c("b0", "b1", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model5_loglik <- extract_log_lik(model5_fit)

## performing leave-one-out CV
# model5_loo <- loo(model5_loglik)
# model5_loo$estimates

## looking at waic
model5WAIC <- waic(model5_loglik)
model5WAIC


##### Stan Model 6 #####
## making list of data to declare what goes into stan model
model6_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$prevday_precip)

## fitting stan model
set.seed(802)
model6_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel6.stan"),data = model6_datalist, chains = 3, iter = 25000, warmup = 5000)
model6_fit

## Extracting Parameters
model6_pars <- rstan::extract(model6_fit, c("b0", "b1", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model6_loglik <- extract_log_lik(model6_fit)

## performing leave-one-out CV
# model6_loo <- loo(model6_loglik)
# model6_loo$estimates

## looking at waic
model6WAIC <- waic(model6_loglik)
model6WAIC


##### Stan Model 7 #####
## making list of data to declare what goes into stan model
model7_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$elev, x2 = SnotelDaymet$prevday_cumulative_precip, x3 = SnotelDaymet$prevday_mintemp)

## fitting stan model
set.seed(802)
model7_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel7.stan"),data = model7_datalist, chains = 3, iter = 25000, warmup = 5000)
model7_fit

## Extracting Parameters
model7_pars <- rstan::extract(model7_fit, c("b0", "b1", "b2", "b3", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model7_loglik <- extract_log_lik(model7_fit)

## performing leave-one-out CV
# model7loo <- loo(model7_loglik)
# model7_loo$estimates

## looking at waic
model7WAIC <- waic(model7_loglik)
model7WAIC

##### Stan Model 8 #####
## making list of data to declare what goes into stan model
model8_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$elev, x2 = SnotelDaymet$prevday_maxtemp, x3 = SnotelDaymet$prevday_mintemp, x4 = SnotelDaymet$prevday_precip)

## fitting stan model
set.seed(802)
model8_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel8.stan"),data = model8_datalist, chains = 3, iter = 25000, warmup = 5000)
model8_fit

## Extracting Parameters
model8_pars <- rstan::extract(model8_fit, c("b0", "b1", "b2", "b3", "b4", "sigma", "meanlog_gt", "sdlog_gt"))

## extracting log likelihood
model8_loglik <- extract_log_lik(model8_fit)

## performing leave-one-out CV
# model8_loo <- loo(model8_loglik)
# model8_loo$estimates

## looking at waic
model8WAIC <- waic(model8_loglik)
model8WAIC



##### Stan Model 9, gamma dist #####
## making list of data to declare what goes into stan model
model9_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$swe_plus, x1 = SnotelDaymet$elev, x2 = SnotelDaymet$prevday_precip, x3 = SnotelDaymet$prevday_mintemp, x4 = SnotelDaymet$prevday_meantemp)

## fitting stan model
set.seed(802)
model9_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel9Gamma.stan"),data = model9_datalist, chains = 3, iter = 40000, warmup = 10000)
model9_fit

## Extracting Parameters
model9_pars <- rstan::extract(model9_fit, c("b0", "b1", "b2", "b3", "b4", "mean_gt", "sd_gt"))

## extracting log likelihood
model9_loglik <- extract_log_lik(model9_fit)

## performing leave-one-out CV
# model10_loo <- loo(model10_loglik)

## looking at waic
model9WAIC <- waic(model9_loglik)


##### Model 9 but in JAGS #####
## initializing model
# inits = list(
#   list(b0 = runif(1,0,5), b1 = runif(1,0,5), b2 = runif(1,0,5), b3 = runif(1,0,5), b4 = runif(1,0,5), sigma = runif(1,0,5)),
#   list(b0 = runif(1,0,5), b1 = runif(1,0,5), b2 = runif(1,0,5), b3 = runif(1,0,5), b4 = runif(1,0,5), sigma = runif(1,0,5)),
#   list(b0 = runif(1,0,5), b1 = runif(1,0,5), b2 = runif(1,0,5), b3 = runif(1,0,5), b4 = runif(1,0,5), sigma = runif(1,0,5)))
# data_swe = list(
#   n = nrow(SnotelDaymet),
#   x1 = as.double(SnotelDaymet$elev),
#   x2 = as.double(SnotelDaymet$prevday_cumulative_precip),
#   x3 = as.double(SnotelDaymet$prevday_mintemp),
#   x4 = as.double(SnotelDaymet$prevday_meantemp),
#   y = as.double(SnotelDaymet$snow_water_equivalent))
# n.adapt = 5000
# n.update = 20000
# n.iter = 20000
# # Call to JAGS
# set.seed(802)
# ## pop growth model
# jm_swe = jags.model(here("Project", "Scripts", "JAGS", "Model9GammaJAGS.R"), data = data_swe, inits = inits,
#                     n.chains = length(inits), n.adapt = n.adapt)
# update(jm_swe, n.iter = n.update)
# zm_swe = coda.samples(jm_swe, variable.names = c("b0", "b1", "b2", "b3", "b4", "mu", "sigma", "shape", "rate"),
#                       n.iter = n.iter, n.thin = 1)
# plot_swe = coda.samples(jm_swe, variable.names = c("shape", "rate"),
#                         n.iter = n.iter, n.thin = 1)
# zj_swe = jags.samples(jm_swe, variable.names = c("shape", "rate", "sigma", "mu"), n.iter = n.iter, n.thin = 1)
# ## estimating posterior dist
# post_df_swe = as.data.frame(rbind(zm_swe[[1]], zm_swe[[2]], zm_swe[[3]]))
# 
# ## plotting densities of params
# # par(mfrow = c(2,2))
# hist(post_df_swe$K, prob = TRUE)
# hist(post_df_swe$r, prob = TRUE)
# hist(post_df_swe$sigma, prob = TRUE)
# hist(post_df_swe$tau, prob = TRUE)
# ## plotting trace plots and density plots
# # par(mfrow = c(2,4), mar = c(2, 2, 2, 2))
# plot(plot_swe)


##### Stan Model 10 #####
## making list of data to declare what goes into stan model
model10_datalist <- list(N = nrow(SnotelDaymet), y = SnotelDaymet$snow_water_equivalent, x1 = SnotelDaymet$elev, x2 = SnotelDaymet$prevday_cumulative_precip, x3 = SnotelDaymet$prevday_mintemp)

## fitting stan model
set.seed(802)
model10_fit <- stan(file=here("Project", "Scripts", "Stan", "ProjectModel10.stan"),data = model10_datalist, chains = 3, iter = 30000, warmup = 10000)
model10_fit

## Extracting Parameters
model10_pars <- rstan::extract(model10_fit, c("b0", "b1", "b2", "b3", "sigma", "mean_gt", "sd_gt"))

## extracting log likelihood
model10_loglik <- extract_log_lik(model10_fit)

## performing leave-one-out CV
# model10_loo <- loo(model10_loglik)

## looking at waic
model10WAIC <- waic(model10_loglik)



##### Examining Fits #####
mean(model1_pars$meanlog_gt)
mean(model2_pars$meanlog_gt)
mean(model3_pars$meanlog_gt)
mean(model4_pars$meanlog_gt)
mean(model5_pars$meanlog_gt)
mean(model6_pars$meanlog_gt)
mean(model7_pars$meanlog_gt)
mean(model8_pars$meanlog_gt)
mean(model9_pars$meanlog_gt)
mean(model10_pars$meanlog_gt)


##### Comparing WAICs #####
WAICComp <- as.data.frame(loo_compare(model1WAIC, model2WAIC, model3WAIC, model4WAIC)) |>
  mutate(ModelID = c("Model 3", "Model 4", "Model 1", "Model 2"), .before = 1) |>
  mutate(DeltaWAIC = waic - min(waic))


WAIC_gtdf <- WAICComp |>
  select(ModelID, waic, DeltaWAIC, elpd_diff)


## Top 25 Table
# adding title and subtitle
WAICgt <- WAIC_gtdf |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_espn() |>
  tab_header(
    title = "WAIC for Top Models", # ...with this title
    subtitle = "calculated with Loo package")  |>  # and this subtitle
  ## tab_style(style = cell_fill("bisque"),
  ##           locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(waic), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 3 # With four decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(elpd_diff), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 3 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(waic),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = TRUE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(elpd_diff),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(ModelID = "Model ID", waic = "WAIC", DeltaWAIC = "Delta WAIC")

WAICgt
WAICgt |>
  gtsave(
    "WAICTable.png", expand = 5,
    path = here("Project", "Outputs")
  )

