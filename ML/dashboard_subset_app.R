<<<<<<< HEAD

# dashboard_subset_app.R
# Minimal Shiny dashboard for ML results (subset, enhanced) for Genomic Variants project
# Uses outputs from PART_AB_Subset_100k_Enhanced.R

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(DT)
  library(plotly)
  library(bslib)
})

# ---- Load Data ----
combined_candidates <- c(
  "cleanVD2_GVPS_Pathogenicity_Combined_subset.rds",
  "cleanVD2_GVPS_Pathogenicity_Combined_subset.csv"
)
combined_path <- combined_candidates[file.exists(combined_candidates)][1]
if (is.na(combined_path)) stop("No combined subset data file found")

if (grepl(".rds$", combined_path)) {
  ml_data <- readRDS(combined_path)
} else {
  ml_data <- read.csv(combined_path, stringsAsFactors = FALSE)
}
# If Pathogenicity_Level is missing, create it from .pred_class
if (!"Pathogenicity_Level" %in% names(ml_data) && ".pred_class" %in% names(ml_data)) {
  ml_data$Pathogenicity_Level <- ml_data$.pred_class
}
ml_data <- ml_data %>% mutate(
  Pathogenicity_Level = factor(Pathogenicity_Level,
                               levels = c("Benign", "Possibly_Deleterious", "Probably_Deleterious", "Highly_Deleterious")),
  GVPS_risk = factor(GVPS_risk, levels = c("Low Risk", "Medium Risk", "High Risk"))
)

# ---- UI ----
theme <- bs_theme(version = 5, bootswatch = "flatly")

ui <- fluidPage(
  theme = theme,
  titlePanel("Genomic Variant ML Results Dashboard (Subset, Enhanced)"),
  tabsetPanel(id = "tabs",
              tabPanel("Overview",
                       fluidRow(
                         column(8,
                                h4("ML Results for Genomic Variants (Subset)"),
                                p("This section summarizes machine learning predictions for variant functional impact (GVPS) and pathogenicity class, using a 100,000-row subset with enhanced modeling.")
                         ),
                         column(4,
                                wellPanel(
                                  h5("Dataset"),
                                  verbatimTextOutput("data_summary")
                                )
                         )
                       )
              ),
              tabPanel("ML Explorer",
                       fluidRow(
                         column(3,
                                selectInput("risk_filter", "GVPS Risk", choices = c("All", levels(ml_data$GVPS_risk)), selected = "All"),
                                selectInput("path_filter", "Pathogenicity Class", choices = c("All", levels(ml_data$Pathogenicity_Level)), selected = "All"),
                                sliderInput("gvps_range", "GVPS range", min(ml_data$GVPS_pred, na.rm = TRUE), max(ml_data$GVPS_pred, na.rm = TRUE), value = range(ml_data$GVPS_pred, na.rm = TRUE)),
                                sliderInput("percentile_range", "GVPS Percentile", min(ml_data$GVPS_percentile, na.rm = TRUE), max(ml_data$GVPS_percentile, na.rm = TRUE), value = range(ml_data$GVPS_percentile, na.rm = TRUE)),
                                hr(),
                                selectInput("scatter_x", "Scatterplot X", choices = names(ml_data), selected = "GVPS_pred"),
                                selectInput("scatter_y", "Scatterplot Y", choices = names(ml_data), selected = "GVPS_percentile"),
                                selectInput("scatter_col", "Color by", choices = c("None", names(ml_data)), selected = "GVPS_risk"),
                                downloadButton("download_filtered", "Download filtered")
                         ),
                         column(9,
                                tabsetPanel(
                                  tabPanel("Table", DTOutput("variant_table")),
                                  tabPanel("Scatter Plot", plotlyOutput("explorer_scatter"))
                                )
                         )
                       )
              ),
              tabPanel("Ranking",
                       fluidRow(
                         column(3,
                                selectInput("rank_metric", "Rank by", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                                numericInput("rank_n", "Top N", value = 25, min = 1, max = 1000),
                                radioButtons("rank_order", "Order", choices = c("Descending" = "desc", "Ascending" = "asc"), selected = "desc")
                         ),
                         column(9,
                                DTOutput("ranking_table")
                         )
                       )
              ),

              tabPanel("Distributions",
                       fluidRow(
                         column(6, plotlyOutput("gvps_hist")),
                         column(6, plotlyOutput("path_bar"))
                       ),
                       fluidRow(
                         column(6, plotlyOutput("risk_bar")),
                         column(6, plotlyOutput("percentile_hist"))
                       ),
                       hr(),
                       fluidRow(
                         column(6,
                                h5("GVPS ECDF Plot"),
                                plotlyOutput("gvps_ecdf_plot")
                         ),
                         column(6,
                                h5("Feature Distribution by Class"),
                                selectInput("feature_select", "Feature (numeric)", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                                plotlyOutput("feature_by_class_plot")
                         )
                       ),
                       fluidRow(
                         column(6,
                                h5("Top Gene/Consequence Barplot"),
                                selectInput("bar_type", "Barplot Type", choices = c("Gene" = "SYMBOL", "Consequence" = "Consequence")),
                                plotlyOutput("top_barplot")
                         )
                       )
              ),
              tabPanel("Performance",
                       fluidRow(
                         column(12,
                                h4("Confusion Matrix (Subset)"),
                                uiOutput("conf_matrix")
                         )
                       )
              ),
              tabPanel("Variant Detail",
                       fluidRow(
                         column(4,
                                selectizeInput("variant_select", "Select rsID", choices = NULL, selected = NULL, options = list(placeholder = "Type to search"))
                         ),
                         column(8,
                                h5("Variant summary"),
                                tableOutput("variant_summary"),
                                h5("Raw fields"),
                                tableOutput("variant_detail")
                         )
                         
                       )
              )
  )
)
# ---- Server ----
server <- function(input, output, session) {
  # ROC curve image display logic
  output$roc_curve_img <- renderUI({
    class_map <- c(
      "Benign" = "roc_curve_Benign_subset.png",
      "Possibly_Deleterious" = "roc_curve_Possibly_Deleterious_subset.png",
      "Probably_Deleterious" = "roc_curve_Probably_Deleterious_subset.png",
      "Highly_Deleterious" = "roc_curve_Highly_Deleterious_subset.png"
    )
    req(input$roc_class)
    img_file <- class_map[[input$roc_class]]
    if (!is.null(img_file) && file.exists(img_file)) {
      tags$img(src = img_file, style = "max-width:100%;height:auto;border:1px solid #ccc;")
    } else {
      tags$div(style = "color:red;", paste("ROC curve image not found for class:", input$roc_class))
    }
  })
  output$data_summary <- renderPrint({
    list(rows = nrow(ml_data), cols = ncol(ml_data), risk_levels = levels(ml_data$GVPS_risk), path_classes = levels(ml_data$Pathogenicity_Level))
  })
  
  filtered_data <- reactive({
    df <- ml_data
    if (!is.null(input$risk_filter) && input$risk_filter != "All") {
      df <- df %>% filter(GVPS_risk == input$risk_filter)
    }
    if (!is.null(input$path_filter) && input$path_filter != "All") {
      df <- df %>% filter(Pathogenicity_Level == input$path_filter)
    }
    if (!is.null(input$gvps_range)) {
      df <- df %>% filter(GVPS_pred >= input$gvps_range[1], GVPS_pred <= input$gvps_range[2])
    }
    if (!is.null(input$percentile_range)) {
      df <- df %>% filter(GVPS_percentile >= input$percentile_range[1], GVPS_percentile <= input$percentile_range[2])
    }
    df
  })
  
  output$variant_table <- renderDT({
    filtered_data() %>%
      select(rsid, GVPS_pred, GVPS_risk, GVPS_percentile, Pathogenicity_Level, starts_with(".pred_"))
  }, options = list(pageLength = 25), filter = "top")
  
  output$explorer_scatter <- renderPlotly({
    req(input$scatter_x, input$scatter_y)
    df <- filtered_data()
    x <- input$scatter_x
    y <- input$scatter_y
    col <- if (input$scatter_col == "None") NULL else input$scatter_col
    p <- ggplot(df, aes_string(x = x, y = y, color = col)) +
      geom_point(alpha = 0.6) +
      labs(title = paste("Scatterplot:", x, "vs", y), x = x, y = y, color = col) +
      theme_minimal()
    ggplotly(p)
  })
  # Ranking Tab
  output$ranking_table <- renderDT({
    req(input$rank_metric)
    df <- filtered_data()
    metric <- input$rank_metric
    n <- input$rank_n
    ord <- if (input$rank_order == "desc") desc else identity
    df_ranked <- df %>% arrange(ord(.data[[metric]])) %>% head(n)
    datatable(df_ranked, options = list(pageLength = n))
  })
  
  output$download_filtered <- downloadHandler(
      filename = function() paste0("ml_variants_filtered_subset_", Sys.Date(), ".csv"),
      content = function(file) write.csv(filtered_data(), file, row.names = FALSE)
    )
  
  output$gvps_hist <- renderPlotly({
      if (!"GVPS_pred" %in% names(ml_data) || all(is.na(ml_data$GVPS_pred))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_pred missing or empty") + theme_void()))
      }
      p <- ggplot(ml_data, aes(x = GVPS_pred)) +
        geom_histogram(bins = 40, fill = "#2c3e50", color = "white", alpha = 0.8) +
        labs(title = "GVPS distribution", x = "GVPS", y = "Count") +
        theme_minimal()
      ggplotly(p)
  })
  
  output$path_bar <- renderPlotly({
      if (!"Pathogenicity_Level" %in% names(ml_data) || all(is.na(ml_data$Pathogenicity_Level))) {
        return(ggplotly(ggplot() + labs(title = "Pathogenicity_Level missing or empty") + theme_void()))
      }
      p <- ml_data %>%
        count(Pathogenicity_Level) %>%
        ggplot(aes(x = Pathogenicity_Level, y = n, fill = Pathogenicity_Level)) +
        geom_col() +
        labs(title = "Pathogenicity class distribution", x = "Class", y = "Count") +
        theme_minimal() + theme(legend.position = "none")
      ggplotly(p)
  })
  
  output$risk_bar <- renderPlotly({
      if (!"GVPS_risk" %in% names(ml_data) || all(is.na(ml_data$GVPS_risk))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_risk missing or empty") + theme_void()))
      }
      p <- ml_data %>%
        count(GVPS_risk) %>%
        ggplot(aes(x = GVPS_risk, y = n, fill = GVPS_risk)) +
        geom_col() +
        labs(title = "GVPS risk stratification", x = "Risk Level", y = "Count") +
        theme_minimal() + theme(legend.position = "none")
      ggplotly(p)
  })
  
  output$percentile_hist <- renderPlotly({
      if (!"GVPS_percentile" %in% names(ml_data) || all(is.na(ml_data$GVPS_percentile))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_percentile missing or empty") + theme_void()))
      }
      p <- ggplot(ml_data, aes(x = GVPS_percentile)) +
        geom_histogram(bins = 40, fill = "#8e44ad", color = "white", alpha = 0.8) +
        labs(title = "GVPS percentile distribution", x = "Percentile", y = "Count") +
        theme_minimal()
      ggplotly(p)
  })
  
  output$outlier_density <- renderPlotly({
      if (!"GVPS_pred" %in% names(ml_data) || !"GVPS_outlier" %in% names(ml_data)) {
        return(ggplotly(ggplot() + labs(title = "GVPS_pred or GVPS_outlier missing") + theme_void()))
      }
      df <- ml_data %>% filter(!is.na(GVPS_pred), !is.na(GVPS_outlier))
      if (nrow(df) == 0) {
        return(ggplotly(ggplot() + labs(title = "No data for outlier density plot") + theme_void()))
      }
      p <- ggplot(df, aes(x = GVPS_pred, color = GVPS_outlier)) +
        geom_density() +
        labs(title = "GVPS outlier density", x = "GVPS", y = "Density") +
        theme_minimal()
      ggplotly(p)
})
  
  output$conf_matrix <- renderUI({
    cm_file <- "pathogenicity_confusion_matrix_subset.png"
    if (file.exists(cm_file)) {
      tags$img(src = cm_file, style = "max-width:100%;height:auto;border:1px solid #ccc;")
    } else {
      tags$div(style = "color:red;", "Confusion matrix image not found")
    }
  })

  # --- Distributions tab: ECDF, feature by class, top barplot ---
  output$gvps_ecdf_plot <- renderPlotly({
    if (!"GVPS_pred" %in% names(ml_data) || all(is.na(ml_data$GVPS_pred))) {
      return(ggplotly(ggplot() + labs(title = "GVPS_pred missing or empty") + theme_void()))
    }
    p <- ggplot(ml_data, aes(x = GVPS_pred)) +
      stat_ecdf(geom = "step", color = "#2c3e50") +
      labs(title = "GVPS ECDF", x = "GVPS", y = "F(x)") +
      theme_minimal()
    ggplotly(p)
  })

  output$feature_by_class_plot <- renderPlotly({
    feat <- input$feature_select
    if (is.null(feat) || !feat %in% names(ml_data) || !"Pathogenicity_Level" %in% names(ml_data)) {
      return(ggplotly(ggplot() + labs(title = "Feature or class missing") + theme_void()))
    }
    df <- ml_data %>% filter(!is.na(.data[[feat]]))
    if (nrow(df) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data for selected feature") + theme_void()))
    }
    p <- ggplot(df, aes(x = Pathogenicity_Level, y = .data[[feat]], fill = Pathogenicity_Level)) +
      geom_violin(alpha = 0.4, color = NA) +
      geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
      labs(title = paste(feat, "by Pathogenicity Class"), x = "Class", y = feat) +
      theme_minimal() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$top_barplot <- renderPlotly({
    bar_type <- input$bar_type
    if (is.null(bar_type) || !bar_type %in% names(ml_data)) {
      return(ggplotly(ggplot() + labs(title = "Barplot type missing or invalid") + theme_void()))
    }
    df <- ml_data %>% count(.data[[bar_type]]) %>% arrange(desc(n)) %>% head(20)
    if (nrow(df) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data for barplot") + theme_void()))
    }
    p <- ggplot(df, aes(x = reorder(.data[[bar_type]], n), y = n)) +
      geom_col(fill = if (bar_type == "SYMBOL") "#8e44ad" else "#16a085") +
      coord_flip() +
      labs(title = paste("Top", bar_type, "(by count)"), x = bar_type, y = "Count") +
      theme_minimal()
    ggplotly(p)
  })

  # --- Variant Detail tab ---
  observe({
    if ("rsid" %in% names(ml_data)) {
      updateSelectizeInput(session, "variant_select", choices = unique(ml_data$rsid), server = TRUE)
    }
  })

  output$variant_summary <- renderTable({
    req(input$variant_select)
    if (!"rsid" %in% names(ml_data)) return(data.frame(Message = "No rsID column in data"))
    row <- ml_data[ml_data$rsid == input$variant_select, , drop = FALSE]
    if (nrow(row) == 0) return(data.frame(Message = "No data for selected rsID"))
    # Show summary fields (customize as needed)
    summary_fields <- intersect(c("rsid", "GVPS_pred", "GVPS_risk", "GVPS_percentile", "Pathogenicity_Level", "SYMBOL", "Consequence"), names(row))
    row[1, summary_fields, drop = FALSE]
  })

  output$variant_detail <- renderTable({
    req(input$variant_select)
    if (!"rsid" %in% names(ml_data)) return(data.frame(Message = "No rsID column in data"))
    row <- ml_data[ml_data$rsid == input$variant_select, , drop = FALSE]
    if (nrow(row) == 0) return(data.frame(Message = "No data for selected rsID"))
    as.data.frame(t(row[1, , drop = FALSE]))
  }, rownames = TRUE)
}

shinyApp(ui, server)
=======

# dashboard_subset_app.R
# Minimal Shiny dashboard for ML results (subset, enhanced) for Genomic Variants project
# Uses outputs from PART_AB_Subset_100k_Enhanced.R

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(DT)
  library(plotly)
  library(bslib)
})

# ---- Load Data ----
combined_candidates <- c(
  "cleanVD2_GVPS_Pathogenicity_Combined_subset.rds",
  "cleanVD2_GVPS_Pathogenicity_Combined_subset.csv"
)
combined_path <- combined_candidates[file.exists(combined_candidates)][1]
if (is.na(combined_path)) stop("No combined subset data file found")

if (grepl(".rds$", combined_path)) {
  ml_data <- readRDS(combined_path)
} else {
  ml_data <- read.csv(combined_path, stringsAsFactors = FALSE)
}
# If Pathogenicity_Level is missing, create it from .pred_class
if (!"Pathogenicity_Level" %in% names(ml_data) && ".pred_class" %in% names(ml_data)) {
  ml_data$Pathogenicity_Level <- ml_data$.pred_class
}
ml_data <- ml_data %>% mutate(
  Pathogenicity_Level = factor(Pathogenicity_Level,
                               levels = c("Benign", "Possibly_Deleterious", "Probably_Deleterious", "Highly_Deleterious")),
  GVPS_risk = factor(GVPS_risk, levels = c("Low Risk", "Medium Risk", "High Risk"))
)

# ---- UI ----
theme <- bs_theme(version = 5, bootswatch = "flatly")

ui <- fluidPage(
  theme = theme,
  titlePanel("Genomic Variant ML Results Dashboard (Subset, Enhanced)"),
  tabsetPanel(id = "tabs",
              tabPanel("Overview",
                       fluidRow(
                         column(8,
                                h4("ML Results for Genomic Variants (Subset)"),
                                p("This section summarizes machine learning predictions for variant functional impact (GVPS) and pathogenicity class, using a 100,000-row subset with enhanced modeling.")
                         ),
                         column(4,
                                wellPanel(
                                  h5("Dataset"),
                                  verbatimTextOutput("data_summary")
                                )
                         )
                       )
              ),
              tabPanel("ML Explorer",
                       fluidRow(
                         column(3,
                                selectInput("risk_filter", "GVPS Risk", choices = c("All", levels(ml_data$GVPS_risk)), selected = "All"),
                                selectInput("path_filter", "Pathogenicity Class", choices = c("All", levels(ml_data$Pathogenicity_Level)), selected = "All"),
                                sliderInput("gvps_range", "GVPS range", min(ml_data$GVPS_pred, na.rm = TRUE), max(ml_data$GVPS_pred, na.rm = TRUE), value = range(ml_data$GVPS_pred, na.rm = TRUE)),
                                sliderInput("percentile_range", "GVPS Percentile", min(ml_data$GVPS_percentile, na.rm = TRUE), max(ml_data$GVPS_percentile, na.rm = TRUE), value = range(ml_data$GVPS_percentile, na.rm = TRUE)),
                                hr(),
                                selectInput("scatter_x", "Scatterplot X", choices = names(ml_data), selected = "GVPS_pred"),
                                selectInput("scatter_y", "Scatterplot Y", choices = names(ml_data), selected = "GVPS_percentile"),
                                selectInput("scatter_col", "Color by", choices = c("None", names(ml_data)), selected = "GVPS_risk"),
                                downloadButton("download_filtered", "Download filtered")
                         ),
                         column(9,
                                tabsetPanel(
                                  tabPanel("Table", DTOutput("variant_table")),
                                  tabPanel("Scatter Plot", plotlyOutput("explorer_scatter"))
                                )
                         )
                       )
              ),
              tabPanel("Ranking",
                       fluidRow(
                         column(3,
                                selectInput("rank_metric", "Rank by", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                                numericInput("rank_n", "Top N", value = 25, min = 1, max = 1000),
                                radioButtons("rank_order", "Order", choices = c("Descending" = "desc", "Ascending" = "asc"), selected = "desc")
                         ),
                         column(9,
                                DTOutput("ranking_table")
                         )
                       )
              ),

              tabPanel("Distributions",
                       fluidRow(
                         column(6, plotlyOutput("gvps_hist")),
                         column(6, plotlyOutput("path_bar"))
                       ),
                       fluidRow(
                         column(6, plotlyOutput("risk_bar")),
                         column(6, plotlyOutput("percentile_hist"))
                       ),
                       hr(),
                       fluidRow(
                         column(6,
                                h5("GVPS ECDF Plot"),
                                plotlyOutput("gvps_ecdf_plot")
                         ),
                         column(6,
                                h5("Feature Distribution by Class"),
                                selectInput("feature_select", "Feature (numeric)", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                                plotlyOutput("feature_by_class_plot")
                         )
                       ),
                       fluidRow(
                         column(6,
                                h5("Top Gene/Consequence Barplot"),
                                selectInput("bar_type", "Barplot Type", choices = c("Gene" = "SYMBOL", "Consequence" = "Consequence")),
                                plotlyOutput("top_barplot")
                         )
                       )
              ),
              tabPanel("Performance",
                       fluidRow(
                         column(12,
                                h4("Confusion Matrix (Subset)"),
                                uiOutput("conf_matrix")
                         )
                       )
              ),
              tabPanel("Variant Detail",
                       fluidRow(
                         column(4,
                                selectizeInput("variant_select", "Select rsID", choices = NULL, selected = NULL, options = list(placeholder = "Type to search"))
                         ),
                         column(8,
                                h5("Variant summary"),
                                tableOutput("variant_summary"),
                                h5("Raw fields"),
                                tableOutput("variant_detail")
                         )
                         
                       )
              )
  )
)
# ---- Server ----
server <- function(input, output, session) {
  # ROC curve image display logic
  output$roc_curve_img <- renderUI({
    class_map <- c(
      "Benign" = "roc_curve_Benign_subset.png",
      "Possibly_Deleterious" = "roc_curve_Possibly_Deleterious_subset.png",
      "Probably_Deleterious" = "roc_curve_Probably_Deleterious_subset.png",
      "Highly_Deleterious" = "roc_curve_Highly_Deleterious_subset.png"
    )
    req(input$roc_class)
    img_file <- class_map[[input$roc_class]]
    if (!is.null(img_file) && file.exists(img_file)) {
      tags$img(src = img_file, style = "max-width:100%;height:auto;border:1px solid #ccc;")
    } else {
      tags$div(style = "color:red;", paste("ROC curve image not found for class:", input$roc_class))
    }
  })
  output$data_summary <- renderPrint({
    list(rows = nrow(ml_data), cols = ncol(ml_data), risk_levels = levels(ml_data$GVPS_risk), path_classes = levels(ml_data$Pathogenicity_Level))
  })
  
  filtered_data <- reactive({
    df <- ml_data
    if (!is.null(input$risk_filter) && input$risk_filter != "All") {
      df <- df %>% filter(GVPS_risk == input$risk_filter)
    }
    if (!is.null(input$path_filter) && input$path_filter != "All") {
      df <- df %>% filter(Pathogenicity_Level == input$path_filter)
    }
    if (!is.null(input$gvps_range)) {
      df <- df %>% filter(GVPS_pred >= input$gvps_range[1], GVPS_pred <= input$gvps_range[2])
    }
    if (!is.null(input$percentile_range)) {
      df <- df %>% filter(GVPS_percentile >= input$percentile_range[1], GVPS_percentile <= input$percentile_range[2])
    }
    df
  })
  
  output$variant_table <- renderDT({
    filtered_data() %>%
      select(rsid, GVPS_pred, GVPS_risk, GVPS_percentile, Pathogenicity_Level, starts_with(".pred_"))
  }, options = list(pageLength = 25), filter = "top")
  
  output$explorer_scatter <- renderPlotly({
    req(input$scatter_x, input$scatter_y)
    df <- filtered_data()
    x <- input$scatter_x
    y <- input$scatter_y
    col <- if (input$scatter_col == "None") NULL else input$scatter_col
    p <- ggplot(df, aes_string(x = x, y = y, color = col)) +
      geom_point(alpha = 0.6) +
      labs(title = paste("Scatterplot:", x, "vs", y), x = x, y = y, color = col) +
      theme_minimal()
    ggplotly(p)
  })
  # Ranking Tab
  output$ranking_table <- renderDT({
    req(input$rank_metric)
    df <- filtered_data()
    metric <- input$rank_metric
    n <- input$rank_n
    ord <- if (input$rank_order == "desc") desc else identity
    df_ranked <- df %>% arrange(ord(.data[[metric]])) %>% head(n)
    datatable(df_ranked, options = list(pageLength = n))
  })
  
  output$download_filtered <- downloadHandler(
      filename = function() paste0("ml_variants_filtered_subset_", Sys.Date(), ".csv"),
      content = function(file) write.csv(filtered_data(), file, row.names = FALSE)
    )
  
  output$gvps_hist <- renderPlotly({
      if (!"GVPS_pred" %in% names(ml_data) || all(is.na(ml_data$GVPS_pred))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_pred missing or empty") + theme_void()))
      }
      p <- ggplot(ml_data, aes(x = GVPS_pred)) +
        geom_histogram(bins = 40, fill = "#2c3e50", color = "white", alpha = 0.8) +
        labs(title = "GVPS distribution", x = "GVPS", y = "Count") +
        theme_minimal()
      ggplotly(p)
  })
  
  output$path_bar <- renderPlotly({
      if (!"Pathogenicity_Level" %in% names(ml_data) || all(is.na(ml_data$Pathogenicity_Level))) {
        return(ggplotly(ggplot() + labs(title = "Pathogenicity_Level missing or empty") + theme_void()))
      }
      p <- ml_data %>%
        count(Pathogenicity_Level) %>%
        ggplot(aes(x = Pathogenicity_Level, y = n, fill = Pathogenicity_Level)) +
        geom_col() +
        labs(title = "Pathogenicity class distribution", x = "Class", y = "Count") +
        theme_minimal() + theme(legend.position = "none")
      ggplotly(p)
  })
  
  output$risk_bar <- renderPlotly({
      if (!"GVPS_risk" %in% names(ml_data) || all(is.na(ml_data$GVPS_risk))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_risk missing or empty") + theme_void()))
      }
      p <- ml_data %>%
        count(GVPS_risk) %>%
        ggplot(aes(x = GVPS_risk, y = n, fill = GVPS_risk)) +
        geom_col() +
        labs(title = "GVPS risk stratification", x = "Risk Level", y = "Count") +
        theme_minimal() + theme(legend.position = "none")
      ggplotly(p)
  })
  
  output$percentile_hist <- renderPlotly({
      if (!"GVPS_percentile" %in% names(ml_data) || all(is.na(ml_data$GVPS_percentile))) {
        return(ggplotly(ggplot() + labs(title = "GVPS_percentile missing or empty") + theme_void()))
      }
      p <- ggplot(ml_data, aes(x = GVPS_percentile)) +
        geom_histogram(bins = 40, fill = "#8e44ad", color = "white", alpha = 0.8) +
        labs(title = "GVPS percentile distribution", x = "Percentile", y = "Count") +
        theme_minimal()
      ggplotly(p)
  })
  
  output$outlier_density <- renderPlotly({
      if (!"GVPS_pred" %in% names(ml_data) || !"GVPS_outlier" %in% names(ml_data)) {
        return(ggplotly(ggplot() + labs(title = "GVPS_pred or GVPS_outlier missing") + theme_void()))
      }
      df <- ml_data %>% filter(!is.na(GVPS_pred), !is.na(GVPS_outlier))
      if (nrow(df) == 0) {
        return(ggplotly(ggplot() + labs(title = "No data for outlier density plot") + theme_void()))
      }
      p <- ggplot(df, aes(x = GVPS_pred, color = GVPS_outlier)) +
        geom_density() +
        labs(title = "GVPS outlier density", x = "GVPS", y = "Density") +
        theme_minimal()
      ggplotly(p)
})
  
  output$conf_matrix <- renderUI({
    cm_file <- "pathogenicity_confusion_matrix_subset.png"
    if (file.exists(cm_file)) {
      tags$img(src = cm_file, style = "max-width:100%;height:auto;border:1px solid #ccc;")
    } else {
      tags$div(style = "color:red;", "Confusion matrix image not found")
    }
  })

  # --- Distributions tab: ECDF, feature by class, top barplot ---
  output$gvps_ecdf_plot <- renderPlotly({
    if (!"GVPS_pred" %in% names(ml_data) || all(is.na(ml_data$GVPS_pred))) {
      return(ggplotly(ggplot() + labs(title = "GVPS_pred missing or empty") + theme_void()))
    }
    p <- ggplot(ml_data, aes(x = GVPS_pred)) +
      stat_ecdf(geom = "step", color = "#2c3e50") +
      labs(title = "GVPS ECDF", x = "GVPS", y = "F(x)") +
      theme_minimal()
    ggplotly(p)
  })

  output$feature_by_class_plot <- renderPlotly({
    feat <- input$feature_select
    if (is.null(feat) || !feat %in% names(ml_data) || !"Pathogenicity_Level" %in% names(ml_data)) {
      return(ggplotly(ggplot() + labs(title = "Feature or class missing") + theme_void()))
    }
    df <- ml_data %>% filter(!is.na(.data[[feat]]))
    if (nrow(df) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data for selected feature") + theme_void()))
    }
    p <- ggplot(df, aes(x = Pathogenicity_Level, y = .data[[feat]], fill = Pathogenicity_Level)) +
      geom_violin(alpha = 0.4, color = NA) +
      geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
      labs(title = paste(feat, "by Pathogenicity Class"), x = "Class", y = feat) +
      theme_minimal() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$top_barplot <- renderPlotly({
    bar_type <- input$bar_type
    if (is.null(bar_type) || !bar_type %in% names(ml_data)) {
      return(ggplotly(ggplot() + labs(title = "Barplot type missing or invalid") + theme_void()))
    }
    df <- ml_data %>% count(.data[[bar_type]]) %>% arrange(desc(n)) %>% head(20)
    if (nrow(df) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data for barplot") + theme_void()))
    }
    p <- ggplot(df, aes(x = reorder(.data[[bar_type]], n), y = n)) +
      geom_col(fill = if (bar_type == "SYMBOL") "#8e44ad" else "#16a085") +
      coord_flip() +
      labs(title = paste("Top", bar_type, "(by count)"), x = bar_type, y = "Count") +
      theme_minimal()
    ggplotly(p)
  })

  # --- Variant Detail tab ---
  observe({
    if ("rsid" %in% names(ml_data)) {
      updateSelectizeInput(session, "variant_select", choices = unique(ml_data$rsid), server = TRUE)
    }
  })

  output$variant_summary <- renderTable({
    req(input$variant_select)
    if (!"rsid" %in% names(ml_data)) return(data.frame(Message = "No rsID column in data"))
    row <- ml_data[ml_data$rsid == input$variant_select, , drop = FALSE]
    if (nrow(row) == 0) return(data.frame(Message = "No data for selected rsID"))
    # Show summary fields (customize as needed)
    summary_fields <- intersect(c("rsid", "GVPS_pred", "GVPS_risk", "GVPS_percentile", "Pathogenicity_Level", "SYMBOL", "Consequence"), names(row))
    row[1, summary_fields, drop = FALSE]
  })

  output$variant_detail <- renderTable({
    req(input$variant_select)
    if (!"rsid" %in% names(ml_data)) return(data.frame(Message = "No rsID column in data"))
    row <- ml_data[ml_data$rsid == input$variant_select, , drop = FALSE]
    if (nrow(row) == 0) return(data.frame(Message = "No data for selected rsID"))
    as.data.frame(t(row[1, , drop = FALSE]))
  }, rownames = TRUE)
}

shinyApp(ui, server)
>>>>>>> e2a8a6051a5efdd3d6763ccdda1f2b6f35629b77
