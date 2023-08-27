library(shiny)
library(ggplot2)

library(shiny)
library(shinyjs)
library(ggplot2)


model_names <- c("Health", "Hypertension", "Hyperlipidemia", "Diabetes", "Ht&HyperLip", "HyperLip&DM", "Ht&Dm", "Ht&HyperLip&Dm")

data_array = as.array(readRDS("~/risk/data_array.rds"))
predicted_risk_treated = as.array(readRDS("~/risk/predicted_risks_treated.rds"))

get_risk_profile=function(data_array, treated_array, 
                          start_cov_profile, change_age_cov = NULL, new_cov_profile = NULL, 
                          start_model, change_ages_model = NULL, new_models = NULL) {
  
  risk_profile <- numeric(length(40:70))
  risk_profile_treated <- numeric(length(40:70))
  
  age_range <- 40:70
  for(i in 1:length(age_range)) {
    age <- age_range[i]
    current_cov <- start_cov_profile
    current_model <- start_model
    
    if(!is.null(change_age_cov) && age >= change_age_cov) {
      current_cov <- new_cov_profile
    }
    
    if(!is.null(change_ages_model)) {
      change_idxs <- which(age >= change_ages_model)
      if(length(change_idxs) > 0) {
        current_model <- new_models[max(change_idxs)]
      }
    }
    
    risk_profile[i] <- data_array[current_cov, i, current_model]
    risk_profile_treated[i] <- treated_array[current_cov, i, current_model]
  }
  
  return(list(untreated = risk_profile, treated = risk_profile_treated))
}



ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  div(class = "title-banner", 
      h2("MSGene"),
      h4("Sarah Urbut, MD PhD")
  ),
  
  useShinyjs(),
  titlePanel("Risk Prediction App"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cadPrs", "CAD PRS:", choices = c(-2, -1, 0, 1, 2)),
      selectInput("sex", "Sex:", choices = c("0", "1")),
      selectInput("smoke", "Smoke:", choices = c("0", "1")),
      selectInput("antihtn", "Anti-HTN:", choices = c("0", "1")),
      selectInput("statin", "Statin:", choices = c("0", "1")),
      checkboxInput("changeProfile", "Change Profile?", FALSE),
      conditionalPanel(
        condition = "input.changeProfile == true",
        numericInput("changeAge", "Age of Covariate Change:", 55, min = 40, max = 70),
        selectInput("newCadPrs", "New CAD PRS:", choices = c(-2, -1, 0, 1, 2)),
        selectInput("newSex", "New Sex:", choices = c("0", "1")),
        selectInput("newSmoke", "New Smoke:", choices = c("0", "1")),
        selectInput("newAntihtn", "New Anti-HTN:", choices = c("0", "1")),
        selectInput("newStatin", "New Statin:", choices = c("0", "1"))
      ),
      selectInput("startModel", "Starting Model:", 
                  choices = c("Health" = 1, 
                              "Hypertension" = 2,
                              "Hyperlipidemia" = 3,
                              "Diabetes" = 4,
                              "Ht&HyperLip" = 5,
                              "HyperLip&DM" = 6,
                              "Ht&Dm" = 7,
                              "Ht&HyperLip&Dm" = 8)),
      actionButton("addModelChange", "Add another model change"),
      div(id = "modelChanges"),
      checkboxInput("showTreated", "Show treated risk profile", TRUE) # Add this line
    ),
    
    mainPanel(
      plotOutput("riskPlot")
    )
  )
)


server <- function(input, output, session) {
  observeEvent(input$addModelChange, {
    insertUI(
      selector = "#modelChanges",
      ui = tags$div(
        selectInput(paste0("newModel", input$addModelChange), "New Model", 
                    choices = c("Health" = 1, 
                                "Hypertension" = 2,
                                "Hyperlipidemia" = 3,
                                "Diabetes" = 4,
                                "Ht&HyperLip" = 5,
                                "HyperLip&DM" = 6,
                                "Ht&Dm" = 7,
                                "Ht&HyperLip&Dm" = 8)),
        numericInput(paste0("modelChangeAge", input$addModelChange), "Age of Model Change", value = 55, min = 40, max = 70)
      )
    )
  })
  
  
  
  output$riskPlot <- renderPlot({
    

    if(input$addModelChange > 0) {
      change_ages_model_values <- sapply(1:input$addModelChange, function(i) as.integer(input[[paste0("modelChangeAge", i)]]))
      new_models_values <- sapply(1:input$addModelChange, function(i) as.integer(input[[paste0("newModel", i)]]))
    } else {
      change_ages_model_values <- NULL
      new_models_values <- NULL
    }
    
    start_cov_row <- get_cov_row(input$cadPrs, input$sex, input$smoke, input$antihtn, input$statin)
    
    if (input$changeProfile) {
      new_cov_row <- get_cov_row(input$newCadPrs, input$newSex, input$newSmoke, input$newAntihtn, input$newStatin)
    } else {
      new_cov_row <- NULL
    }
    
    risk_values <- get_risk_profile(
      data_array,
      start_cov_profile = start_cov_row,
      change_age_cov = if(input$changeProfile) as.integer(input$changeAge) else NULL,
      new_cov_profile = new_cov_row,
      start_model = as.integer(input$startModel),
      change_ages_model = change_ages_model_values,
      new_models = new_models_values
    )
    
    
    df <- data.frame(Age = 40:70, Risk = risk_values)
    
    # Adding a group column for different models
    df$Group <- as.integer(input$startModel) # initialize with the start model
    
    # If there are model changes, adjust the group column accordingly
    if(!is.null(change_ages_model_values)) {
      for(i in 1:length(change_ages_model_values)) {
        age_of_change <- change_ages_model_values[i]
        new_model <- new_models_values[i]
        df$Group[df$Age >= age_of_change] <- as.integer(new_model)
      }
    }
    
    ggplot(df, aes(x = Age, y = Risk, group = Group, color = as.factor(Group))) + 
      geom_line(size = 1.5) +  
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18),  
        axis.title.x = element_text(size = 14, face = "bold"),  
        axis.title.y = element_text(size = 14, face = "bold"),  
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),  
        legend.text = element_text(size = 12)
      ) + 
      labs(
        title = "Predicted Lifetime Risk Profile over Age",
        y = "Predicted Lifetime Risk",
        x = "Age"
      ) +
      scale_color_brewer(palette = "Set1", name = "Model", labels = model_names) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))
  })
  
  
  df_untreated <- data.frame(Age = 40:70, Risk = risk_values, Type = "Untreated")
  
  # Add this section to calculate the treated risks
  risk_values_treated <- get_risk_profile(
    predicted_risk_treated, # Note that we're using the treated array here
    start_cov_profile = start_cov_row,
    change_age_cov = if(input$changeProfile) as.integer(input$changeAge) else NULL,
    new_cov_profile = new_cov_row,
    start_model = as.integer(input$startModel),
    change_ages_model = change_ages_model_values,
    new_models = new_models_values
  )
  
  df_treated <- data.frame(Age = 40:70, Risk = risk_values_treated, Type = "Treated")
  
  # Check if the "showTreated" checkbox is checked
  if(input$showTreated) {
    df <- rbind(df_untreated, df_treated)
  } else {
    df <- df_untreated
  }
  
  ggplot(df, aes(x = Age, y = Risk, color = Type, group = interaction(Group, Type))) + 
    geom_line(size = 1.5) +  
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),  
      axis.title.x = element_text(size = 14, face = "bold"),  
      axis.title.y = element_text(size = 14, face = "bold"),  
      axis.text.x = element_text(size = 12),  
      axis.text.y = element_text(size = 12),
      legend.title = element_blank(),  
      legend.text = element_text(size = 12)
    ) + 
    labs(
      title = "Predicted Lifetime Risk Profile over Age",
      y = "Predicted Lifetime Risk",
      x = "Age"
    ) +
    scale_color_brewer(palette = "Set1", name = "Model", labels = model_names
}

