library(rsample) # for bootstrapping
library(purrr) # for map_dbl function
library(infer)
library(Metrics)
library(tidyverse)
options(stringsAsFactors=T)

# Limitation Functions
LNO3 <- function(NO3, KNO3 = 0.5, KNH4 = 0.5, NH4) {
  return((NO3/(NO3 + KNO3)) / (1 + (NH4 / KNH4)))
}

LNH4 <- function(NH4, KNH4 = 0.5) {
  return(NH4/(NH4 + KNH4))
}

# Phytoplankton Primary Production Model 
# other models can be created, see equation 5 in the paper
# The log-transformed model. We use the property that the log of a product is the sum of the log's.
PP.model.szot.log <- function (Temp, NO3, NH4, Chla, psi_pmax = 0.08065, KNO3=0.5, KNH4=0.5, 
                               CChla = 38.2, gamma_P = 0.07, m_P = 0.05,
                               mu_max_min = 4, mu_max_0 = 0.55) {
  return((log(pmax(mu_max_min, mu_max_0 * exp(psi_pmax * Temp))) + 
            log((NO3 / (NO3 + KNO3) / (1 + NH4 / KNH4)) + NH4 / (NH4 + KNH4)) + 
            log((1 - gamma_P) - m_P)) + 
           log(Chla) + log(CChla / 24))  
}

# Different Regions of the Chesapeake Bay
levels <- c('Fresh Water', 'Oligohaline', 'Mesohaline', 
            'Polyhaline & Mixoeuhaline')

# Different Seasons of the year
seasons <- c('winter', 'spring', 'summer', 'autumn')

# creates model skill data frames w/out error column
subsets <- data.frame(Season = rep(NA, length(seasons) * length(levels)),
                      Salinity = rep(NA, length(seasons) * length(levels)),
                      Value1 = rep(NA, length(seasons) * length(levels)),
                      Error1 = rep(NA, length(seasons) * length(levels)),
                      Value2 = rep(NA, length(seasons) * length(levels)),
                      Error2 = rep(NA, length(seasons) * length(levels)),
                      Value3 = rep(NA, length(seasons) * length(levels)),
                      Error3 = rep(NA, length(seasons) * length(levels)),
                      St.Value1 = rep(NA, length(seasons) * length(levels)),
                      St.Value2 = rep(NA, length(seasons) * length(levels)),
                      St.Value3 = rep(NA, length(seasons) * length(levels)),
                      Conv = rep(NA, length(seasons) * length(levels)))


get_start_nls_values <- function (params_calibrate_) {
  params_init_ <- {}
  for (para_ in params_calibrate_) {
    if (para_ == "KNO3") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 0, max = 5)
    } else if (para_ == "KNH4") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 0, max = 5)
    } else if (para_ == "CChla") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 1, max = 150)
    } else if (para_ == "mu00") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 0.2, max = 5)
    } else if (para_ == "mu01") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 0.2, max = 4)
    } else if (para_ == "KPO4") {
      params_init_[[paste0(para_, "_")]] <- runif(1, min = 0, max = 5)
    } 
  }
  return(params_init_)
}

# Calibrating the model using the Nonlinear least squares function
# lines 74 - 147

params_calibrate <- c('CChla', 'KNO3', 'KNH4') # Parameters to calibrate 
number.of.start.values <- 200 # Max number of try start values in the nls function

# re-run this loop if there are any non-zero values under the conv. collumn 
# in the dataframe "subsets"
i <- 1
for (season_ in seasons) {
  for (level_ in levels) {
    
    sum.of.residuals <- Inf # Initial value of the sum of residuals (Inf is the worst model)
    
    if (!dir.exists("data")) {
      dir.create("data")
      
      dataset_name <- "train_data_VA_AP_MD_all_06102024.RData"
      
      # If dataset is in current working directory, move it into "data"
      if (file.exists(dataset_name)) {
        file.rename(dataset_name, file.path("data", dataset_name))
      } else {
        warning(paste("Dataset", dataset_name, "not found in current directory."))
      }
    }
    
    load("./data/train_data_VA_AP_MD_all_06102024.RData")
    
    train_data.alt <- train_data %>%
      filter(CARBFIX > 0 ) %>%
      filter(CHLA > 0 ) %>%
      mutate(SALINITY.lvl = case_when(SALINITY < 0.5 ~ "Fresh Water",
                                      SALINITY >= 0.5 & SALINITY < 5 ~ "Oligohaline",
                                      SALINITY >= 5 & SALINITY < 18 ~ "Mesohaline",
                                      SALINITY >= 18 ~ "Polyhaline & Mixoeuhaline"),
             lnCARBFIX = log(CARBFIX))
    train_data.alt <- train_data.alt[(train_data.alt$SALINITY.lvl == level_) & (train_data.alt$seasons == season_),]
    for (k in 1:number.of.start.values) {
      start_values.tmp <- get_start_nls_values(params_calibrate_ = params_calibrate)
      nls_fit_tmp <- nls(train_data.alt$lnCARBFIX ~  PP.model.szot.log(Temp = train_data.alt$WTEMP, 
                                                                 NO3 = train_data.alt$NO3F,
                                                                 NH4 = train_data.alt$NH4F,
                                                                 Chla = train_data.alt$CHLA,
                                                                 CChla = CChla_, 
                                                                 KNH4 = KNH4_,
                                                                 KNO3 = KNO3_,
                                                                 m_P = 0.0),
                         data = train_data.alt,
                         start = start_values.tmp,
                         algorithm = "port",
                         control = list(maxiter = 50000, minFactor = 1/2000, warnOnly=T),
                         lower = c(1, 0.01, 0.01),
                         upper = c(250, 5, 5))
      # Overwrite value if sum of squares is lower
      if (nls_fit_tmp$m$deviance() < sum.of.residuals) {
        # Update the minimum sum of residuals and optimal model
        sum.of.residuals <- nls_fit_tmp$m$deviance()
        nls_fit <- nls_fit_tmp
        start_values <- start_values.tmp
      }
    }
    subsets$Season[i] <- season_
    subsets$Salinity[i] <- level_
    subsets$Value1[i] <- coef(nls_fit)[1] 
    subsets$Value2[i] <- coef(nls_fit)[2]
    subsets$Value3[i] <- coef(nls_fit)[3]
    subsets$Error1[i] <- (summary(nls_fit)$parameters[, "Std. Error"][1]) 
    subsets$Error2[i] <- (summary(nls_fit)$parameters[, "Std. Error"][2]) 
    subsets$Error3[i] <- (summary(nls_fit)$parameters[, 'Std. Error'][3])
    subsets$St.Value1[i] <- start_values[1]
    subsets$St.Value2[i] <- start_values[2]
    subsets$St.Value3[i] <- start_values[3]
    subsets$Conv[i] <- nls_fit$convergence # 0 = converged, 1 = didn't converge
    i <- i + 1
  }
}
formula_used <- formula(nls_fit)
model_used <- as.character(formula_used[[3]])[1]
cat("Calibration complete. Model used:", model_used)

if (all(subsets$Conv == 0)) {
  cat("\nAll subsets converged successfully.")
} else {
  stop("\nOne or more of the subsets did not converge. Code ended.")
}

# Plotting all of the calibrated variables. You may need to comment out the code 
# for plotting certain variables if you alter which ones need calibrated 

library(ggplot2)
Palette.Cb.friendly <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Plotting Optimum C:Chla values
p <- subsets %>%
  mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline', 
                                'Polyhaline & Mixoeuhaline'),
         Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
  ggplot( aes(x = Season, y = Value1,
              fill = Salinity)) +
  geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
  theme_gray(base_size = 22)+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(size = 14, hjust = 0.5)     
  )+
  ylim(0, 50)+
  geom_errorbar(aes(x = Season, ymin = Value1 - as.numeric(Error1),
                    ymax = Value1 + as.numeric(Error1)), width = 0.5, 
                position = position_dodge(0.9)) +
  scale_fill_manual(values = Palette.Cb.friendly) +
  labs(title="Szot Carbon (C) to Chloraphyll-a (Chla) ratio", 
       x = "Seasons",
       y = "C:Chla") 
ggsave("figures/Szot_C:Chla.pdf", plot = p, width = 8, height =6)

# Plotting Optimum KNO3 Values
p <- subsets %>%
  mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline',
                                'Polyhaline & Mixoeuhaline'),
         Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
  ggplot( aes(x = Season, y = Value2,
              fill = Salinity)) +
  geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
  theme_gray(base_size = 22)+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(size = 14, hjust = 0.5)   
  )+
  ylim(0, 5) +
  # geom_errorbar(aes(x = Season, ymin = pmax(0, Value2 - as.numeric(Error2)),
  #                   ymax = pmin(5, Value2 + as.numeric(Error2))), width = 0.5,
  #               position = position_dodge(0.9)) +
  scale_fill_manual(values = Palette.Cb.friendly) +
  labs(title = expression(paste("Szot K"[NO3])),
       x = "Seasons",
       y = expression(paste("K"[NO3] ~ "(mmol-N"~ m^-3, ")", sep = ""))) 
ggsave("figures/Szot_KN03.pdf", plot = p, width = 8, height =6)

# Plotting Optimum KNH4 Values
p <-subsets %>%
  mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline',
                                'Polyhaline & Mixoeuhaline'),
         Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
  ggplot( aes(x = Season, y = Value3,
              fill = Salinity)) +
  geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
  theme_gray(base_size = 22)+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(size = 14, hjust = 0.5) 
  )+
  ylim(0, 5) +
  # geom_errorbar(aes(x = Season, ymin = pmax(0, Value3 - as.numeric(Error3)),
  #                   ymax = pmin(5, Value3 + as.numeric(Error3))), width = 0.5,
  #               position = position_dodge(0.9)) +
  scale_fill_manual(values = Palette.Cb.friendly) +
  labs(title = expression(paste("Szot K"[NH4])),
       x = "Seasons",
       y = expression(paste("K"[NH4] ~ "(mmol-N"~ m^-3, ")", sep = ""))) 
ggsave("figures/Szot_KNH4.pdf", plot = p, width = 8, height =6)

# Now that the Model is calibrated, the second half of this code is designed to 
# test and plot model skill, providing insight into the model's performance

subsets.skill <- data.frame(Season = rep(NA, length(seasons) * length(levels)),
                            Salinity = rep(NA, length(seasons) * length(levels)),
                            Value = rep(NA, length(seasons) * length(levels)),
                            Error = rep(NA, length(seasons) * length(levels)),
                            Lim = rep(NA, length(seasons) * length(levels)),
                            Err = rep(NA, length(seasons) * length(levels)),
                            Conv = rep(NA, length(seasons) * length(levels)))

# This loop only runs one statistical test at a time. Choose stat tests
# when you change the statistical test and re-run the code, make sure to comment 
# out the plots for the other tests that were not run (lines 388 - 460)

stat_test <- cor # functions: cor, bias, rrse, rmse etc 

if (identical(stat_test, cor)) {
  test_name <- "Correlation"
} else if (identical(stat_test, bias)) {
  test_name <- "Bias"
} else if (identical(stat_test, rrse)) {
  test_name <- "Root Relative Squared Error"
} else if (identical(stat_test, rmse)) {
  test_name <- "Root Mean Squared Error"
} else {
  test_name <- "unknown"
}

i <- 1
for (season_ in seasons) {
  for (level_ in levels) {
    load("./data/train_data_VA_AP_MD_all_06102024.RData")
    train_data.alt <- train_data %>%
      filter(CARBFIX > 0 ) %>%
      filter(CHLA > 0 ) %>%
      mutate(SALINITY.lvl = case_when(SALINITY < 0.5 ~ "Fresh Water",
                                      SALINITY >= 0.5 & SALINITY < 5 ~ "Oligohaline",
                                      SALINITY >= 5 & SALINITY < 18 ~ "Mesohaline",
                                      SALINITY >= 18 ~ "Polyhaline & Mixoeuhaline"),
             lnCARBFIX = log(CARBFIX))
    train_data.alt <- train_data.alt[(train_data.alt$SALINITY.lvl == level_) & (train_data.alt$seasons == season_),]
    
    resample_ <- bootstraps(train_data.alt, 
                            times = 2e2,
                            apparent = TRUE)
    
    values <- map_dbl(
      resample_$splits,
      function(x) {
        dat <- as.data.frame(x)
        nls_fit <- nls(dat$lnCARBFIX ~  PP.model.szot.log(Temp = dat$WTEMP, 
                                                    NO3 = dat$NO3F,
                                                    NH4 = dat$NH4F,
                                                    Chla = dat$CHLA,
                                                    CChla = CChla_,
                                                    KNO3 = KNO3_,
                                                    KNH4 = KNH4_,
                                                    m_P = 0.0), 
                       data = dat,
                       start = list(CChla_ = subsets$St.Value1[i],
                                    KNO3_ = subsets$St.Value2[i],
                                    KNH4_ = subsets$St.Value3[i]), 
                       algorithm = "port", 
                       control = list(warnOnly = TRUE),
                       lower = c(1, 0.01, 0.01),
                       upper = c(250, 20, 20))
        stat_test((dat$CARBFIX), exp(predict(nls_fit))) 
      } 
    )
    
    sampling_dist <- data.frame(values) %>% 
      specify(response = values) %>% 
      assume("t") 
    
    # Compute the mean of the values distribution
    sample_mean <-  data.frame(values) %>% 
      specify(response = values) %>%
      calculate(stat = "mean")
    
    subsets.skill$Season[i] <- season_
    subsets.skill$Salinity[i] <- level_
    
    subsets.skill$Lim[i] <- mean(LNO3(NO3 = train_data.alt$NO3F,
                                      NH4 = train_data.alt$NH4F,
                                      KNO3 = subsets$Value2[i],
                                      KNH4 = subsets$Value3[i]) + 
                                   LNH4(NH4 = train_data.alt$NH4F,
                                        KNH4 = subsets$Value3[i]))
    subsets.skill$Err[i] <- sd(LNO3(NO3 = train_data.alt$NO3F,
                                         NH4 = train_data.alt$NH4F,
                                         KNO3 = subsets$Value2[i],
                                         KNH4 = subsets$Value3[i]) + 
                                      LNH4(NH4 = train_data.alt$NH4F,
                                           KNH4 = subsets$Value3[i]))
    
    #error = 1 std.dev
    subsets.skill$Error[i] <- as.numeric(sd(values))
    #Value = Sample Mean
    subsets.skill$Value[i] <- as.numeric(sample_mean)
    # 0 = converged, 1 = didn't converge
    subsets.skill$Conv[i] <- nls_fit$convergence 
    i <- i + 1
  }
}

cat("\nModel skill test complete. Statistical test used:", test_name)

if (all(subsets.skill$Conv == 0)) {
  cat("\nAll model skill subsets converged successfully.")
} else {
  stop("\nOne or more model skill subsets did not converge. Code ended.")
}

# Plotting Mean Limitations (do not uncomment)
p <- subsets.skill %>%
  mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline', 
                                'Polyhaline & Mixoeuhaline'),
         Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
  ggplot( aes(x = Season, y = Lim, fill = Salinity)) +
  geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
  theme_gray(base_size = 22)+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(size = 14, hjust = 0.5) 
  )+
  ylim(0, 1) +
  geom_errorbar(aes(x = Season, ymin = Lim - Err, ymax = pmin(1, Lim + Err)), 
                width = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = Palette.Cb.friendly) +
  labs(title="Szot Mean Limitation",
       x = "Seasons",
       y = expression(paste("Mean Limitation (mmol-N"~ m^-3 , ")", sep = ""))) 
ggsave("figures/Szot_Mean_Limitation.pdf", plot = p, width = 8, height =6)

# Plotting Correlation (do not uncomment)
p <- subsets.skill %>%
  mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline', 
                                'Polyhaline & Mixoeuhaline'),
         Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
  ggplot( aes(x = Season, y = Value, fill = Salinity)) +
  geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
  theme_gray(base_size = 22)+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(size = 14, hjust = 0.5) 
  )+
  ylim(0, 1) + 
  geom_errorbar(aes(x = Season, ymin = Value - Error, ymax = Value + Error),
                width = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = Palette.Cb.friendly) +
  labs(title="Szot Correlation",
       x = "Seasons",
       y = expression(paste("Correlation"))) 
ggsave("figures/Szot_Correlation.pdf", plot = p, width = 8, height =6)

# Plotting Bias (do not uncomment)
# p <- subsets.skill %>%
#   mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline',
#               'Polyhaline & Mixoeuhaline'),
#          Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
#   ggplot( aes(x = Season, y = Value, fill = Salinity)) +
#   geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
#   theme_gray(base_size = 22)+
#   theme(
#     legend.text = element_text(size = 14),
#     legend.title = element_text(size = 16),
#     axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
#     plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
#     plot.subtitle = element_text(size = 14, hjust = 0.5) 
#   )+
#   ylim(-40, 10) + 
#   geom_errorbar(aes(x = Season, ymin = Value - Error, ymax = Value + Error),
#                 width = 0.5, position = position_dodge(0.9)) +
#   scale_fill_manual(values = Palette.Cb.friendly) +
#   labs(title="Szot Bias",
#        x = "Seasons",
#        y = expression(paste("Bias (mgC"~ m^-3 ~ h^-1, ")", sep = ""))) 
# ggsave("figures/Szot_Bias.pdf", plot = p, width = 8, height =6)

# Plotting RRSD (do not uncomment)
# p <- subsets.skill %>%
#   mutate(Salinity = fct_relevel(Salinity, 'Fresh Water', 'Oligohaline', 'Mesohaline',
#                                 'Polyhaline & Mixoeuhaline'),
#          Season = fct_relevel(Season, 'winter', 'spring', 'summer', 'autumn')) %>%
#   ggplot( aes(x = Season, y = Value, fill = Salinity)) +
#   geom_bar(stat="identity",width = 0.8, position = position_dodge(width = .9)) +
#   theme_gray(base_size = 22)+
#   theme(
#     legend.text = element_text(size = 14),
#     legend.title = element_text(size = 16),
#     axis.text.x = element_text(size = 14, angle = 30, hjust = 1),  
#     plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
#     plot.subtitle = element_text(size = 14, hjust = 0.5) 
#   )+
#   ylim(0, 1.5) + 
#   geom_errorbar(aes(x = Season, ymin = Value - Error, ymax = Value + Error),
#                 width = 0.5, position = position_dodge(0.9)) +
#   scale_fill_manual(values = Palette.Cb.friendly) +
#   labs(title="Szot Root Relative Squared Difference (RRSD)",
#        x = "Seasons",
#        y = expression(paste("RRSD"))) 
# ggsave("figures/Szot_RRSD.pdf", plot = p, width = 8, height =6)
dev.off()
