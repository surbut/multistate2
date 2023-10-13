library(ggsci)
library("dplyr")
library("tidyverse")
library("survival")
library("survminer")
library("eulerr")
library("ggplot2")
library("reshape")
library("ggfortify")
library("gridExtra")
library("survMisc")
library("rsq")
library("pROC")
library("data.table")
library("RColorBrewer")
library("gt")
library("DT")
library(plotly)
library(dplyr)
library(carData)
library(gapminder)
library(babynames)
col_pce <- rev(brewer.pal(10, "Paired"))
col_p = c(col_pce[1:3],
          col_pce[5],
          col_pce[4],
          col_pce[c(6, 10)],
          "black",
          col_pce[9],
          "turquoise",
          col_pce[7])
tems = theme(
  axis.text.x = element_text(
    color = "grey20",
    size = 20,
    angle = 90,
    hjust = .5,
    vjust = .5,
    face = "plain"
  ),
  axis.text.y = element_text(
    color = "grey20",
    size = 20,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ),
  axis.title.x = element_text(
    color = "grey20",
    size = 20,
    angle = 0,
    hjust = .5,
    vjust = 0,
    face = "plain"
  ),
  axis.title.y = element_text(
    color = "grey20",
    size = 20,
    angle = 90,
    hjust = .5,
    vjust = .5,
    face = "plain"
  )
)

compute_CVrisk2=function (df, scores = c("ascvd_10y_accaha", "as2", 
                                         "ascvd_10y_frs_simple", "chd_10y_mesa", "chd_10y_mesa_cac"), 
                          age, gender, race, sbp = NULL, bmi = NULL, hdl = NULL, totchol = NULL, 
                          bp_med = NULL, smoker = NULL, diabetes = NULL, lipid_med = NULL, 
                          fh_heartattack = NULL, cac = NULL) 
{
  all_args <- as.list(environment())
  valid_pred <- c("age", "gender", "race", "sbp", "bmi", "hdl", 
                  "totchol", "bp_med", "smoker", "diabetes", "lipid_med", 
                  "fh_heartattack", "cac")
  pred_args <- list()
  for (var in valid_pred) {
    if (!is.null(all_args[[var]])) 
      pred_args[[var]] <- df[[all_args[[var]]]]
  }
  results <- sapply(scores, function(x) do.call(x, pred_args))
  row.names(results) <- NULL
  return(cbind(df, results))
}

as2=function (race = "white", gender = c("male", "female"), age, 
              totchol, hdl, sbp, bp_med, smoker, diabetes, ...) 
{
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  ascvd_pooled_coef <- NULL
  utils::data(ascvd_pooled_coef, envir = environment())
  race <- ifelse(race %in% c("white", "aa"), race, "white")
  race_sex <- data.frame(race, gender)
  race_sex$id <- as.numeric(row.names(race_sex))
  pooled_coef <- merge(race_sex, ascvd_pooled_coef)
  pooled_coef <- pooled_coef[order(pooled_coef$id), ]
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  indv_sum <- log(age) * pooled_coef$ln_age + log(age)^2 * 
    pooled_coef$ln_age_squared + log(totchol) * pooled_coef$ln_totchol + 
    log(age) * log(totchol) * pooled_coef$ln_age_totchol + 
    log(hdl) * pooled_coef$ln_hdl + log(age) * log(hdl) * 
    pooled_coef$ln_age_hdl + log(sbp_treated) * pooled_coef$ln_treated_sbp + 
    log(sbp_treated) * log(age) * pooled_coef$ln_age_treated_sbp + 
    log(sbp_untreated) * pooled_coef$ln_untreated_sbp + log(sbp_untreated) * 
    log(age) * pooled_coef$ln_age_untreated_sbp + smoker * 
    pooled_coef$smoker + smoker * log(age) * pooled_coef$ln_age_smoker + 
    diabetes * pooled_coef$diabetes
  risk_score <- round((1 - (pooled_coef$baseline_survival^exp(indv_sum - 
                                                                pooled_coef$group_mean))) * 100, 2)
  ifelse(risk_score < 1, 1, risk_score)
}


