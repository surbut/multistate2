

#
output_path = "~/Library/CloudStorage/Dropbox-Personal/ukbb-ehr-data/data/"
data_period <- readRDS(paste0(output_path, "data_period.rds"))
biomarkers <- readRDS(paste0(output_path, "biomarkers2.rds"))
prescriptions <- readRDS(paste0(output_path, "prescriptions.rds"))
records <- readRDS(paste0(output_path, "all_records.rds"))

library(data.table)
library(lubridate)
library(plyr)
library(ggplot2)
library(cowplot)

data_plot =
  function (period_data,
            record_data,
            labels,
            x_limits,
            y_title,
            x_axis = TRUE,
            legend = FALSE,
            top_margin = 5,
            bottom_margin = 5) {
    # Constants
    col_palette <- c(
      "#0077bb",
      # blue
      "#33bbee",
      # cyan
      "#009988",
      # teal
      "#ee7733",
      # orange
      "#cc3311",
      # red
      "#ee3377",
      # magenta
      "#dddddd",
      # light grey
      "#555555"
    )  # off-black
    # Format legend and x-axis
    theme_legend <- if (legend) {
      theme(legend.position = "bottom",
            legend.margin = margin(
              t = 0,
              l = 1,
              r = 1,
              b = 0
            ))
    } else {
      theme(legend.position = "none")
    }
    theme_x_axis <- if (x_axis) {
      NULL
    } else {
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    }
    # Return plot
    ggplot() +
      geom_rect(
        data = period_data,
        aes(
          xmin = from,
          xmax = to,
          fill = "Data collection"
        ),
        ymin = -0.3,
        ymax = 3.3
      ) +
      geom_linerange(data = record_data,
                     aes(
                       x = date,
                       y = 0,
                       ymin = ymin,
                       ymax = ymax,
                       colour = type
                     )) +
      scale_x_date(limits = x_limits) +
      scale_y_continuous(name = y_title, limits = c(-0.3, 3.3)) +
      scale_colour_manual(
        name = NULL,
        limits = c("event", "test", "presc"),
        labels = c("Diagnosis/event", "Test/observation", "Prescription"),
        values = c(col_palette[1], col_palette[5], col_palette[4])
      ) +
      scale_fill_manual(name = NULL,
                        limits = "Data collection",
                        values = col_palette[7]) +
      labels +
      theme_minimal() +
      theme(
        text = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(
          t = top_margin,
          l = 5,
          r = 5,
          b = bottom_margin
        )
      ) +
      theme_x_axis +
      theme_legend
  }

pheno_plot <- function (id, save = NULL) {
  # Constants
  
  drug_types <- c(
    "insulin" = "Insulin",
    "metformin" = "Metformin",
    "non_insulin" = "Other anti-diabetic",
    "anti_hypertensives" = "Anti-hypertensive",
    "antipsychotics" = "Atypical anti-psychotic",
    "statins" = "Statin",
    "steroids" = "Steroid"
  )
  col_palette <- c(
    "#0077bb",
    # blue
    "#33bbee",
    # cyan
    "#009988",
    # teal
    "#ee7733",
    # orange
    "#cc3311",
    # red
    "#ee3377",
    # magenta
    "#dddddd",
    # light grey
    "#555555"
  )  # off-black
  # Participant data
  data_period <- data_period[eid == id]
  biomarkers <- biomarkers[eid == id]
  prescriptions <- prescriptions[eid == id]
  prescriptions <-
    prescriptions[!(category == "diabetes" & type == "any")]
  records <- records[eid == id]
  records[type == "event", c("ymin", "ymax") := .(2, 3)]
  records[type == "test", c("ymin", "ymax") := .(1, 2)]
  records[type == "presc", c("ymin", "ymax") := .(0, 1)]
  # Plot limits
  x_limits <-
    c(min(data_period$from) - 365, max(data_period$to) + 365)
  # Plot theme
  theme_set(theme_minimal())
  theme_update(
    text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(
      t = 0,
      l = 1,
      r = 1,
      b = 0
    ),
    plot.margin = margin(
      t = 0,
      l = 5,
      r = 5,
      b = 0
    )
  )
  # EHR data plot
  ehr_plot <- data_plot(
    data_period,
    records,
    labs(
      title = paste("Participant summary: ID", id),
      subtitle = "Inferred period of data collection"
    ),
    x_limits,
    y_title = NULL,
    legend = TRUE,
    bottom_margin = 0
  )
  # Prescriptions plot
  prescriptions[, name := category]
  prescriptions[category == "diabetes", name := type]
  prescriptions[, name := revalue(name, drug_types, warn_missing = FALSE)]
  presc_categories <- prescriptions[, unique(name)]
  lapply(1:length(presc_categories), function (i) {
    prescriptions[name == presc_categories[i], c("ymin", "ymax") := .(i, i + 1)]
    NULL
  })
  drug_labels <- unique(prescriptions[, .(name, ymin)])
  drug_plot <- ggplot() +
    geom_rect(data = prescriptions,
              aes(
                xmin = from,
                xmax = to,
                ymin = ymin,
                ymax = ymax,
                fill = category
              )) +
    geom_text(
      data = drug_labels,
      aes(y = ymin + 0.5, label = name),
      x = x_limits[1],
      hjust = 0
    ) +
    scale_x_date(name = NULL,
                 limits = x_limits,
                 date_labels = "%Y") +
    scale_y_continuous(name = NULL) +
    scale_fill_discrete(type = col_palette) +
    labs(subtitle = "Prescription history") +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  # Biomarker plots
  ldl_plot <- ggplot() +
    geom_point(data = biomarkers[variable == "ldl"], aes(x = date, y = value), colour = col_palette[6]) +
    geom_line(data = biomarkers[variable == "ldl"], aes(x = date, y = value), colour = col_palette[6]) +
    scale_x_date(name = NULL,
                 limits = x_limits,
                 date_labels = "%Y") +
    scale_y_continuous(name = bquote(LDL ~ (ng / dL))) +
    labs(subtitle = "Biomarkers") +
    geom_vline(
      xintercept = min(cad$all_event_dt.Include_in_cases[identifier %in% id, eventdate]),
      color = "red",
      size = 5
    ) +
    theme(plot.margin = margin(
      t = 10,
      l = 5,
      r = 5,
      b = 5
    ))
  
  
  sbp_plot <- ggplot() +
    geom_point(data = biomarkers[variable == "sbp"], aes(x = date, y = value), colour = col_palette[5]) +
    geom_line(data = biomarkers[variable == "sbp"], aes(x = date, y = value), colour = col_palette[5]) +
    scale_x_date(name = NULL,
                 limits = x_limits,
                 date_labels = "%Y") +
    scale_y_continuous(name = bquote(SBP ~ (mmHg))) +
    labs(subtitle = "Biomarkers") +
    theme(plot.margin = margin(
      t = 10,
      l = 5,
      r = 5,
      b = 5
    ))
  # Composite plot
  all_plots <- list(ehr_plot, drug_plot, ldl_plot, sbp_plot)
  out <-
    plot_grid(
      plotlist = all_plots,
      ncol = 1,
      align = "hv",
      axis = "lr",
      rel_heights = c(1.45, 1.25, 1, 1.1)
    )
  print(out)
}




tenlb = function(start,
                 stop,
                 modelfit,
                 agesmooth,
                 agesint,
                 atrisk,
                 window_width = 20,
                 span = 0.75,
                 degree = 2,
                 prs_quants) {
  mat = coefplotsmooth2(
    agesmooth,
    start,
    stop,
    modelfit,
    window_width = window_width,
    span = span,
    degree = degree
  )$custom_smooth
  ten.year = matrix(NA,
                    nrow = length(agesint),
                    ncol = length(prs_quants) * 2)
  ten.year.ben = matrix(NA,
                        nrow = length(agesint),
                        ncol = length(prs_quants) * 2)
  lifetime = matrix(NA,
                    nrow = length(agesint),
                    ncol = length(prs_quants) * 2)
  lifetime.ben = matrix(NA,
                        nrow = length(agesint),
                        ncol = length(prs_quants) * 2)
  
  for (i in 1:length(agesint)) {
    age = agesint[i]
    ten.year[i,] = compute_prediction_product_matrix(
      atrisk = atrisk,
      agepredinterval = c(age:min(age + 10, 80)),
      coefmat = mat
    )$PredictedIntervalrisk
    ten.year.ben[i,] = compute_prediction_product_matrix(
      atrisk = atrisk,
      agepredinterval = c(age:min(age + 10, 80)),
      coefmat = mat
    )$Hazard_treated
    
    lifetime[i,] = compute_prediction_product_matrix(
      atrisk = atrisk,
      agepredinterval = c(age:(80)),
      coefmat = mat
    )$PredictedIntervalrisk
    lifetime.ben[i,] = compute_prediction_product_matrix(
      atrisk = atrisk,
      agepredinterval = c(age:(80)),
      coefmat = mat
    )$Hazard_treated
  }
  
  ten.year = data.frame(ten.year) * 100
  lifetime = data.frame(lifetime) * 100
  lifetime.ben = data.frame(lifetime.ben) * 100
  ten.year.ben = data.frame(ten.year.ben) * 100
  
  e = length(prs_quants)
  rownames(ten.year) = rownames(ten.year.ben) = agesint
  colnames(ten.year) = colnames(ten.year.ben) = paste0(rep(round(prs_quants, 1), 2), ":", rep(c("female", "male"), each =
                                                                                                e))
  ten.year$age = ten.year.ben$age = agesint
  
  rownames(lifetime) = rownames(lifetime.ben) = agesint
  colnames(lifetime) = colnames(lifetime.ben) = paste0(rep(round(prs_quants, 1), 2), ":", rep(c("female", "male"), each =
                                                                                                e))
  lifetime$age = lifetime.ben$age = agesint
  
  
  lookup_table <- data.frame(melt(ten.year, id.vars = c("age")))
  lookup_table_treat <-
    data.frame(melt(ten.year.ben, id.vars = c("age")))
  
  df <-
    inner_join(lookup_table, lookup_table_treat, by = c("age", "variable"))
  # colors <- c("darkred", "red", "lightcoral", "darkblue", "blue", "lightblue")
  #
  #
  
  
  d = df[, c(1, 2, 3, 4)]
  colnames(d) = c("year", "strata", "untreated_risk", "treated_risk")
  
  
  gtentreat = ggplot(d, aes(x = year)) +
    geom_ribbon(aes(
      ymin = treated_risk,
      ymax = untreated_risk,
      fill = strata
    ),
    alpha = 0.3) +
    geom_line(aes(
      y = treated_risk,
      color = strata,
      linetype = "treated"
    )) +
    geom_line(aes(
      y = untreated_risk,
      color = strata,
      linetype = "untreated"
    )) +
    labs(
      x = "Age",
      y = "Predicted Ten year risk (%)",
      fill = "Strata",
      color = "Strata",
      linetype = "Strategy"
    ) +
    theme_classic()
  
  
  
  
  
  names(lookup_table)[3] = "ten.year"
  gten = ggplot(lookup_table, aes(age, y = ten.year, color = as.factor(variable))) +
    stat_smooth() + labs(
      x = "Age",
      y = paste0("Ten Year Risk from ", start, " to ", stop),
      col = "PRS:Sex"
    )
  
  
  
  lookup_table2 <- data.frame(melt(lifetime, id.vars = c("age")))
  lookup_table2$strategy = rep("untreated_risk", nrow(lookup_table))
  lookup_table_treat2 <-
    data.frame(melt(lifetime.ben, id.vars = c("age")))
  lookup_table_treat2$strategy = rep("treated_risk", nrow(lookup_table_treat2))
  
  df_long_me = rbind(lookup_table2, lookup_table_treat2)
  df <-
    inner_join(lookup_table2, lookup_table_treat2, by = c("age", "variable"))
  
  
  d = df[, c(1, 2, 3, 5)]
  colnames(d) = c("year", "strata", "untreated_risk", "treated_risk")
  # gltreat = ggplot(d, aes(x = year)) +
  #   geom_ribbon(aes(
  #     ymin = treated_risk,
  #     ymax = untreated_risk,
  #     fill = strata
  #   ),
  #   alpha = 0.3) +
  #   geom_line(aes(
  #     y = treated_risk,
  #     color = strata,
  #     linetype = "treated"
  #   )) +
  #   geom_line(aes(
  #     y = untreated_risk,
  #     color = strata,
  #     linetype = "untreated"
  #   )) +
  #   scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(10)) +
  #   labs(
  #     title = "Predicted Lifetime Risk Curves",
  #     x = "Age",
  #     y = "Risk",
  #     fill = "Strata",
  #     color = "Strata",
  #     linetype = "Strategy"
  #   ) +
  #   theme_classic()
  
  
  gltreat <- ggplot(d, aes(x = year)) +
    geom_ribbon(aes(
      ymin = treated_risk,
      ymax = untreated_risk,
      fill = strata
    ),
    alpha = 0.3) +
    geom_line(aes(
      y = treated_risk,
      color = strata,
      linetype = "treated"
    )) +
    geom_line(aes(
      y = untreated_risk,
      color = strata,
      linetype = "untreated"
    )) +
    labs(
      x = "Age",
      y = "Predicted Lifetime Risk (%)",
      fill = "Strata",
      color = "Strata",
      linetype = "Strategy"
    ) +
    theme_classic()
  
  
  
  names(lookup_table2)[3] = "lifetime"
  glife = ggplot(lookup_table2, aes(age, y = lifetime, color = as.factor(variable))) +
    stat_smooth() + labs(
      x = "Age",
      y = paste0("Lifetme Risk from ", start, " to ", stop),
      col = "PRS:Sex"
    )
  
  return(
    list(
      "tenplot" = gten,
      "lifeplot" = glife,
      "gtentreat" = gtentreat,
      "lifetreat" = gltreat
    )
  )
}
