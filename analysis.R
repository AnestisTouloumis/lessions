## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(
  size = 'tiny',
  eval = TRUE,                          # TRUE; logical
  echo = TRUE,                          # TRUE; logical
  results = 'markup',                   # 'markup', 'asis', 'hold', 'hide' 
  collapse = TRUE,                      # FALSE; logical
  error = TRUE,                        # TRUE; logical
  message = TRUE,                      # TRUE; logical
  include = TRUE,                       # TRUE; logical
  warning = TRUE,                      # TRUE; logical
  strip.white = TRUE,                   # TRUE; logical
  tidy = TRUE,                          # FALSE; logical
  tidy.opts = list(width.cutoff = 70),  
  prompt = TRUE,                        # FALSE; logical
  comment = NA,                         # '##'; character or NA
  highlight = TRUE,                     # TRUE; logical
  fig.keep = 'high',                    # 'high', 'none', 'all', 'first', 'last'
  fig.show = 'asis',                    # 'asis', 'hold', 'animate', 'hide'
  #fig.width = 6.5,                     # 7; numeric
  #fig.height = 4.5,                    # 7; numeric
  out.width = '0.90\\textwidth',        # character
  out.height = '0.80\\textheight',      # character
  fig.align = 'center',                 # 'default', 'right', 'left', 'center' 
  purl = TRUE,                          # TRUE; logical
  R.options = list(width = 80)
)


## ---------------------------------------------------------------------------------------------
library("readr")
options(tidyverse.quiet = TRUE)
library("tidyverse")
## package testthat needed to check whether two R objects are equal/identical
library("testthat", warn.conflicts = FALSE) 
library("multgee", quietly = TRUE)
library("VGAM", quietly = TRUE)


## ---------------------------------------------------------------------------------------------
lesions <- 
  read_csv("xray_spa_dc_synd_Osteo_vert3s.csv", show_col_types = FALSE)
lesions <- 
  lesions |>
  group_by(id, vertebra) |>
  mutate(id_vertebra = cur_group_id()) |>
  ungroup()


## ---------------------------------------------------------------------------------------------
lesions |> 
  _$Osteo_synm_tot_corr_upp_low_2 |> 
  table(useNA = "always")

lesions |> 
  _$Osteo_synm_tot_corr_upp_low_2 |> 
  table() |> 
  prop.table()

lesions |> 
  filter(synm_tot_lag == 1) |> 
  _$Osteo_synm_tot_corr_upp_low_2 |> 
  table(useNA = "always")

lesions |> 
  filter(synm_tot_lag == 1) |> 
  _$Osteo_synm_tot_corr_upp_low_2 |> 
  table() |> 
  prop.table()

lesions |>
  with(table(synm_tot_lag, Osteo_synm_tot_corr_upp_low_2)) |>
  print()

lesions |>
  filter(!is.na(synm_tot_lag)) |>
  _$id |>
  n_distinct()

## This is different from the reported number of 3696
lesions |> 
  _$id_vertebra |> 
  n_distinct()


## ---------------------------------------------------------------------------------------------
lesions <- 
  lesions |>
  select(-id, -bdmard, -synm_tot, -vertebra)


## ----tidy=FALSE-------------------------------------------------------------------------------
lesions <-
  lesions |>
  mutate(Osteo_synm_tot_corr_upp_low_2 =
           factor(Osteo_synm_tot_corr_upp_low_2, levels = c("1", "2", "0")),
         Osteo_synm_tot_corr_upp_low_2_lag =
           factor(Osteo_synm_tot_corr_upp_low_2_lag),
         synm_tot_lag =
           factor(synm_tot_lag),
         sexe =
           factor(sexe),
         hla =
           factor(hla),
         tabac_10y =
           factor(tabac_10y),
         profession =
           factor(profession),
         bdmard_lag =
           factor(bdmard_lag))


## ---------------------------------------------------------------------------------------------
fit_null <- nomLORgee(formula = Osteo_synm_tot_corr_upp_low_2 ~ 1,
                      data = lesions,
                      id = id_vertebra,
                      repeated = t_new,
                      LORstr = "time.exch")
  
fit_null |> 
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ synm_tot_lag) |>
  summary()

fit_null |> 
  update(Osteo_synm_tot_corr_upp_low_2 ~ Osteo_synm_tot_corr_upp_low_2_lag) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ sexe) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ age_m0) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ bmi) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ hla) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ tabac_10y) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ profession) |>
  summary()

fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ bdmard_lag) |>
  summary()


## ----tidy=FALSE-------------------------------------------------------------------------------
model_gee_multinomial <- 
  fit_null |>
  update(formula = Osteo_synm_tot_corr_upp_low_2 ~ synm_tot_lag + 
           Osteo_synm_tot_corr_upp_low_2_lag + sexe + age_m0 + bmi + hla + 
           tabac_10y + profession + bdmard_lag) 
model_gee_multinomial |> summary()


## ---------------------------------------------------------------------------------------------
model_gee_multinomial |> 
  coefficients() |> 
  exp()


## ---------------------------------------------------------------------------------------------
model_gee_independence <- 
  model_gee_multinomial |>
  update(LORstr = "independence") 


## ---------------------------------------------------------------------------------------------
model_ml <-
  vglm(formula = Osteo_synm_tot_corr_upp_low_2 ~ synm_tot_lag +
         Osteo_synm_tot_corr_upp_low_2_lag + sexe + age_m0 + bmi + hla +
         tabac_10y + profession + bdmard_lag,
       data = lesions,
       family = multinomial(ref = "0"))


## ---------------------------------------------------------------------------------------------
model_gee_independence_coef <- 
  model_gee_independence |>
  coefficients() |> 
  matrix(12, 2)
model_ml_coef <-
  model_ml |> 
  coefficients() |> 
  matrix(12, 2, TRUE) 
expect_equal(model_gee_independence_coef, model_ml_coef)


## ---------------------------------------------------------------------------------------------
fitted_values_gee <- 
  model_gee_independence |> 
  fitted.values() |>
  as.numeric()
fitted_values_ml <-
  model_ml |> 
  fitted.values() |>
  as.numeric()
expect_equal(fitted_values_gee, fitted_values_ml)


## ----tidy=FALSE-------------------------------------------------------------------------------
lesions2 <-
  lesions |>
  mutate(Osteo_synm_tot_corr_upp_low_2 =
           factor(Osteo_synm_tot_corr_upp_low_2, levels = c("0", "1", "2")))
model_gee_multinomial2 <- 
  nomLORgee(formula = Osteo_synm_tot_corr_upp_low_2 ~ synm_tot_lag +
              Osteo_synm_tot_corr_upp_low_2_lag + sexe + age_m0 + bmi + hla +
              tabac_10y + profession + bdmard_lag,
            data = lesions2,
            id = id_vertebra,
            repeated = t_new,
            LORstr = "independence") 
model_gee_multinomial2 |> summary()


## ---------------------------------------------------------------------------------------------
fitted_values_gee2 <- 
  model_gee_multinomial2 |>
  fitted.values() |>
  _[, c(2, 3, 1)] |>
  as.numeric()
expect_equal(fitted_values_gee, fitted_values_gee2)


## ---------------------------------------------------------------------------------------------
model_gee_independence |> summary()

