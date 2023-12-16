library(shiny)
library(ggplot2)

# Assuming 'seurat_obj' is your Seurat object
# Ensure you have already defined 'seurat_obj' before running this Shiny app code

# UI
ui <- fluidPage(
  titlePanel("Flow vs CITE-seq Correlation Plot Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("threshold", "Threshold:", min = 0, max = 100, value = 0),
      radioButtons("DSB", "dsb normalised:", choices = c("Y", "N"), selected = "N"),
   #   radioButtons("lane", "Select Lane:", choices = c(1, 2), selected = 1),
   radioButtons("lane", "Select Lane:", choices = c(1, 2), selected = 1),
      checkboxInput("logCorrelation", "Log Correlation", value = FALSE),
      actionButton("submit", "Submit"),
      radioButtons("correlationMethod", "Select Correlation Method:", choices = c("Pearson", "Spearman"), selected = "Pearson")
    ),
    mainPanel(
      plotOutput("correlationPlot")
    )
  )
)

# Server
server <- function(input, output) {
  observe({
    threshold <- input$threshold
    lane <- input$lane
    
    # Assuming 'cd19_positive_cells' is already defined
   ## cd19_positive_cells <- subset(seurat_obj, subset = CD19.1 > threshold & Lane == lane)
   
    ### DBS
    
    ###ADD READ FUNCTION TO READ IN THE RDS FILES
    if (input$DSB == "N"){
      DefaultAssay(seurat_obj)<-"CITE"
      cd19_positive_cells <- subset(seurat_obj, subset = CD19.1 > threshold & Lane == lane) 
    }else{
      DefaultAssay(seur)<-"CITE"
      cd19_positive_cells <- subset(seur, subset = CD19.1 > threshold & Lane == lane)    
    }
      
  
    # Define the markers
    markers <- c("CD19.1", "IgD", "CD27.1", "CD21", "CD11c", "CD85j", "CD32", "CD72.1", "CD307d", "CD307e", "Ig-light-chain-kappa", "Ig-light-chain-lambda")
    
    # Create a matrix to store counts
    marker_counts_table <- matrix(NA, nrow = length(markers), ncol = 2, dimnames = list(markers, c("Count", "Percentage")))
    
    # Loop through each marker
    for (marker in markers) {
      # Count cells expressing the marker
      marker_counts <- ncol(cd19_positive_cells[, GetAssayData(cd19_positive_cells)[marker, ] > threshold])
      
      # Calculate percentage
      marker_percentage <- (marker_counts / ncol(cd19_positive_cells)) * 100
      
      # Store counts and percentage in the table
      marker_counts_table[marker, "Count"] <- marker_counts
      marker_counts_table[marker, "Percentage"] <- marker_percentage
    }
    
    # Provided table
    flow_cytometry <- data.frame(
      id = c("CD19.1", "IgD", "CD27.1", "CD21", "CD11c", "CD85j", "CD32", "CD72.1", "CD307d", "CD307e", "Ig-light-chain-kappa", "Ig-light-chain-lambda"),
      percent = c(100.0, 93.0, 27.0, 70.0, 7.0, 96.0, 90.0, 70.0, 2.0, 0.1, 53.0, 44.0)
    )
    
    # Extract percentages from the marker_counts_table
    calculated_percentages <- marker_counts_table[, "Percentage"]
    
    # Combine both tables into a single data frame
    if (input$logCorrelation) {
      correlation_data <- data.frame(log10(flow_cytometry$percent), log10(calculated_percentages))
    } else {
      correlation_data <- data.frame(flow_cytometry$percent, calculated_percentages)
    }
    
    # Calculate correlation
    correlation <- reactive({
      cc<<- correlation_data
      if (input$logCorrelation) {
        if (input$correlationMethod == "Pearson") {
          cor(correlation_data, method = "pearson")
        } else {
          cor(correlation_data, method = "spearman")
        }
      } else {
        if (input$correlationMethod == "Pearson") {
          cor(correlation_data, method = "pearson")
        } else {
          cor(correlation_data, method = "spearman")
        }
      }
    })
    
    # Create the correlation plot
    output$correlationPlot <- renderPlot({
      req(input$submit)
      
      if (input$logCorrelation) {
        ggplot(correlation_data, aes(x = log10(flow_cytometry$percent), y = log10(calculated_percentages))) +
          geom_point() +
          geom_text(aes(label = rownames(correlation_data))) +  # Adjust positioning as needed
          geom_smooth(method = "lm", se = FALSE, color = "blue") +
          labs(title = paste("Correlation Plot", " Sample ", input$lane, sep=""),
               x = "Provided Percentages\n By Flow",
               y = "Calculated Percentages\n By Single Cell") +
          annotate("text", x = -1, y = 0.1,
                   label = 
                     paste("Correlation =", round(correlation(), 2)[2], "P-value =", cor.test(correlation_data[,1], correlation_data[,2])$p.value),
                   hjust = 0, vjust = 1, size = 6, color = "red") + theme_grey(base_size = 25)
      } else {
        ggplot(correlation_data, aes(x = flow_cytometry$percent, y = calculated_percentages)) +
          geom_point() +
          geom_text(aes(label = rownames(correlation_data))) +  # Adjust positioning as needed
          geom_smooth(method = "lm", se = FALSE, color = "blue") +
          labs(title = paste("Correlation Plot", " Sample ", input$lane, sep=""),
               x = "Provided Percentages\n By Flow",
               y = "Calculated Percentages\n By Single Cell") +
          annotate("text", x = 0, y = 10,
                   label = 
                     paste("Correlation =", round(correlation(), 2)[2], "P-value =", cor.test(correlation_data[,1], correlation_data[,2])$p.value),
                   hjust = 0, vjust = -3, size = 6, color = "red")+ theme_grey(base_size = 25)
      }
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)