# Load ggplot2
library(ggplot2)

# Data
timeline_data <- data.frame(
  year = c(2005, 2009, 2011, 2013, 2013, 2018, 2019, 2021),
  event = c("White Sox Win",
            "Stats with John Storey",
            "Chalk Thesis",
            "GTEx project",
            "PhD start",
            "PRS published in Nature Genetics",
            "Start residency",
            "Come to Boston!")
)

# Adjust y-values for overlapping events
timeline_data$y <- c(1, 2, 3, 4, 5, 6, 7, 8) # You can adjust these values as per your visual preference

# Plot
timeline_plot <- ggplot(timeline_data, aes(x = year, y = y)) + 
  geom_point(color = "blue", size = 4) + 
  geom_text(aes(label = event), hjust = 1.2) + 
  labs(title = "Timeline", x = "Year", y = "") + 
  theme_minimal() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlim(2003, 2022)  # extend x-axis limits

print(timeline_plot)