library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(DT)
library(plotly)
library(bslib)
library(tidyverse)
library(tibble)
library(stringr)
library(plotly)
library(htmltools)
library(rentrez)
library(htmltools)





## --------- Helper for AF donut plot ----------
af_donut_plot <- function(af_value) {
  af_value <- as.numeric(af_value)
  if (is.na(af_value) || af_value <= 0 || af_value >= 1) {
    af_value <- NA_real_
  }
  
  df <- data.frame(
    cat   = c("Alt AF", "Other"),
    value = c(af_value, 1 - af_value)
  )
  
  ggplot(df, aes(x = 2, y = value, fill = cat)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(fill = NULL)
}

pubmed_cache_path <- "pubmed_cache.rds"

read_pubmed_cache <- function() {
  if (file.exists(pubmed_cache_path)) readRDS(pubmed_cache_path) else list()
}

write_pubmed_cache <- function(cache) {
  saveRDS(cache, pubmed_cache_path)
}

pubmed_fetch <- function(query, retmax = 20) {
  s <- rentrez::entrez_search(db = "pubmed", term = query, retmax = retmax)
  if (length(s$ids) == 0) return(tibble::tibble())
  
  sums <- rentrez::entrez_summary(db = "pubmed", id = s$ids)
  
  get_name <- function(y) {
    # y might be a list with $name, or a named list, or a character
    if (is.list(y) && !is.null(y[["name"]])) return(y[["name"]] %||% "")
    if (is.character(y)) return(y[1] %||% "")
    ""
  }
  
  df <- tibble::tibble(
    pmid   = vapply(sums, function(x) x$uid %||% "", character(1)),
    title  = vapply(sums, function(x) x$title %||% "", character(1)),
    source = vapply(sums, function(x) x$source %||% "", character(1)),
    year   = vapply(sums, function(x) x$pubdate %||% "", character(1)),
    authors = vapply(sums, function(x) {
      a <- x$authors
      if (is.null(a) || length(a) == 0) return("")
      
      # If authors comes back as a character vector
      if (is.character(a)) return(paste(a, collapse = ", "))
      
      # If authors is a list of objects
      if (is.list(a)) {
        nm <- vapply(a, get_name, character(1))
        nm <- nm[nzchar(nm)]
        return(paste(nm, collapse = ", "))
      }
      
      ""
    }, character(1))
  )
  
  df$url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", df$pmid, "/")
  df
}

`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || is.na(a)) b else a


make_pubmed_query <- function(vr) {
  # If vr is already a gene symbol like "COL11A1"
  if (is.atomic(vr) && length(vr) >= 1) {
    gene <- as.character(vr)[1]
  } else {
    # vr is a row-like object (tibble/data.frame/list)
    gene <- ""
    if (!is.null(vr[["SYMBOL"]])) gene <- as.character(vr[["SYMBOL"]])[1]
    if (identical(gene, "") && !is.null(vr[["hgnc_symbol"]])) 
      gene <- as.character(vr[["hgnc_symbol"]])[1]
  }
  
  gene <- trimws(gene)
  if (is.na(gene) || !nzchar(gene)) return("")
  
  paste0("(", gene, "[Title/Abstract])")
}


## --------- Helper to build a variant “info card” from one-row df ----------

variant_risk_panel <- function(v) {
  v <- v[1, , drop = FALSE]
  
  tagList(
    risk_gauge_ui("phyloP (clipped)", v$phylop_clip,
                  min = -10, max = 10, low = -2, high = 2),
    
    risk_gauge_ui("REVEL", v$revel_max,
                  min = 0, max = 1, high = 0.5),
    
    risk_gauge_ui("SpliceAI", v$spliceai_ds_max,
                  min = 0, max = 1, high = 0.2)
  )
}


risk_gauge_ui <- function(label, value, min, max, fmt = "%.2f",
                          low = NULL, high = NULL) {
  
  # handle missing
  if (is.null(value) || length(value) == 0 || is.na(value)) {
    value <- NA_real_
  } else {
    value <- as.numeric(value)
  }
  
  # clamp + scale to 0..100
  if (!is.na(value)) value_clamp <- max(min(value, max), min) else value_clamp <- NA_real_
  pct <- if (!is.na(value_clamp)) 100 * (value_clamp - min) / (max - min) else 0
  
  # optional marker positions (e.g., "high risk" threshold)
  marker <- function(x) {
    if (is.null(x) || is.na(x)) return(NULL)
    mp <- 100 * (x - min) / (max - min)
    tags$div(style = paste0(
      "position:absolute; left:", mp, "%; top:-3px; bottom:-3px; width:2px; background:#111; opacity:0.8;"
    ))
  }
  
  tags$div(
    style="margin:10px 0;",
    tags$div(
      style="display:flex; justify-content:space-between; align-items:baseline;",
      tags$span(style="font-weight:600;", label),
      tags$span(style="color:#555;",
                if (!is.na(value)) sprintf(fmt, value) else "NA")
    ),
    tags$div(
      style="position:relative; height:12px; border-radius:8px; background:#eee; overflow:hidden;",
      tags$div(style=paste0(
        "height:100%; width:", pct, "%; background:#4C72B0; border-radius:8px;"
      )),
      marker(low),
      marker(high)
    ),
    tags$div(
      style="display:flex; justify-content:space-between; margin-top:4px; font-size:0.75em; color:#777;",
      tags$span(min),
      tags$span(max)
    )
  )
}



variant_card <- function(v) {
  v <- v[1, , drop = FALSE]
  
  # graceful fallbacks if some cols are missing
  variant_id <- if ("Existing_variation" %in% names(v)) {
    v$Existing_variation
  } else if (all(c("chrom", "pos") %in% names(v))) {
    paste0("chr", v$chrom, ":", v$pos)
  } else {
    ""
  }
  
  aa_change <- if ("Amino_acids" %in% names(v)) v$Amino_acids else NA
  codons    <- if ("Codons" %in% names(v)) v$Codons else NA
  
  cdna_pos  <- if ("cDNA_position" %in% names(v)) v$cDNA_position else NA
  cds_pos   <- if ("CDS_position" %in% names(v)) v$CDS_position else NA
  prot_pos  <- if ("Protein_position" %in% names(v)) v$Protein_position else NA
  biotype <- if ("BIOTYPE" %in% names(v)) v$BIOTYPE else NA
  
  
  # GO tags as chips
  go_tags <- c(
    if ("GO_BP_count" %in% names(v) && !is.na(v$GO_BP_count) && v$GO_BP_count > 0)
      paste0("BP: ", v$GO_BP_count),
    if ("GO_MF_count" %in% names(v) && !is.na(v$GO_MF_count) && v$GO_MF_count > 0)
      paste0("MF: ", v$GO_MF_count),
    if ("GO_CC_count" %in% names(v) && !is.na(v$GO_CC_count) && v$GO_CC_count > 0)
      paste0("CC: ", v$GO_CC_count)
  )
  
  htmltools::tags$div(
    style = "border: 1px solid #ddd; border-radius: 8px; padding: 12px; margin-bottom: 12px;",
    
    # Header: gene + variant ID
    htmltools::tags$h4(
      style = "margin-top: 0; margin-bottom: 4px;",
      paste0(v$SYMBOL, if (!is.na(variant_id) && variant_id != "") paste0("  —  ", variant_id) else "")
    ),
    htmltools::tags$p(
      style = "color: #666; margin-top: 0;",
      v$Consequence
    ),
    clinsig_val <- if ("CLIN_SIG" %in% names(v)) v$CLIN_SIG else if ("clnsig_raw" %in% names(v)) v$clnsig_raw else NA,
    clinsig_traffic_ui(clinsig_val),
    htmltools::tags$hr(),
    
    htmltools::tags$h4("Functional Risk Scores"),
    variant_risk_panel(v),
    
    
    htmltools::tags$div(
      style = "display: flex; gap: 16px;",
      
      # LEFT: transcript-level change info
      htmltools::tags$div(
        style = "flex: 2;",
        htmltools::tags$dl(
          style = "margin-bottom: 6px;",
          htmltools::tags$dt(htmltools::tags$b("Amino acid change:")), htmltools::tags$dd(aa_change),
          htmltools::tags$dt(htmltools::tags$b("Codons:")),              htmltools::tags$dd(codons),
          htmltools::tags$dt(htmltools::tags$b("cDNA position:")),       htmltools::tags$dd(cdna_pos),
          htmltools::tags$dt(htmltools::tags$b("CDS position:")),        htmltools::tags$dd(cds_pos),
          htmltools::tags$dt(htmltools::tags$b("Biotype:")),    htmltools::tags$dd(biotype),
          htmltools::tags$dt(htmltools::tags$b("Protein position:")),    htmltools::tags$dd(prot_pos)
        ),
        # GO chips
        if (length(go_tags) > 0) {
          htmltools::tags$div(
            style = "margin-top: 6px;",
            htmltools::tags$span(style = "font-weight: 600; font-size: 0.9em;", "GO tags: "),
            lapply(go_tags, function(tg) {
              htmltools::tags$span(
                style = paste(
                  "display:inline-block; margin:2px; padding:2px 6px;",
                  "border-radius:10px; background-color:#eef; font-size:0.8em;"
                ),
                tg
              )
            })
          )
        }
      ),
      
      # RIGHT: AF donut
      htmltools::tags$div(
        style = "flex: 1; text-align:center;",
        htmltools::tags$p(style = "margin-bottom: 4px; font-weight: 600;", "Allele Frequency"),
        plotOutput("variant_af_donut", height = "140px")
      )
    )
  )
}


### Define S3 class: Variant
Variant <- function(df_row) {
  if (!is.data.frame(df_row)) {
    stop("Input must be a data frame")
  }
  if (nrow(df_row) != 1) {
    stop("Data frame must contain exactly one row for Variant")
  }
  
  obj <- list(
    genomic = list(
      chromosome = df_row$chrom, #
      position   = df_row$pos, #
      ref        = df_row$REF_ALLELE, # df_row$ref
      alt        = df_row$Allele # df_row$alt #
    ),
    annotation = list(
      gene        = df_row$Gene, # df_row$gene, #
      symbol      = df_row$SYMBOL, # df_row$symbol, #
      transcript  = df_row$STRAND, # df_row$transcript, #
      consequence = df_row$Consequence
    ),
    clinical = list(
      phenotype = df_row$CLIN_SIG # df_row$phenotype #
      #  disease   = df_row$disease
    ),
    pathogenicity = list(
      SIFT      = df_row$SIFT,
      PolyPhen  = df_row$PolyPhen,
      FATHMM    = df_row$fathmm, # df_row$fathmm,
      CADD      = df_row$CADD_PHRED,
      GVPS      = df_row$GVPS_pred# df_row$CADD_PHRED
    )
  )
  
  class(obj) <- "Variant"
  return(obj)
}

## Print method
print.Variant <- function(x, ...) {
  cat("Variant Object\n")
  cat("----------------------------\n")
  cat("Genomic:\n")
  cat(sprintf(
    "  %s:%s %s>%s\n",
    x$genomic$chromosome,
    x$genomic$position,
    x$genomic$ref,
    x$genomic$alt
  ))
  
  cat("\nAnnotation:\n")
  cat(sprintf(
    "  Gene: %s\n  Symbol: %s\n  Consequence: %s\n",
    x$annotation$gene,
    x$annotation$symbol,
    x$annotation$consequence
  ))
  
  cat("\nClinical:\n")
  cat(sprintf(
    "  Phenotype: %s\n",
    x$clinical$phenotype
  ))
  
  invisible(x)
}

## Summary method
summary.Variant <- function(object, ...) {
  out <- list(
    genomic = object$genomic,
    annotation = object$annotation,
    clinical = object$clinical
  )
  class(out) <- "summary.Variant"
  out
}

print.summary.Variant <- function(x, ...) {
  cat("Summary of Variant\n")
  cat("----------------------------\n")
  str(x)
  invisible(x)
}

## Accessor methods
gene <- function(x) UseMethod("gene")
gene.Variant <- function(x) x$annotation$gene

clinical_significance <- function(x) UseMethod("clinical_significance")
clinical_significance.Variant <- function(x) x$clinical$phenotype


### Read variant dataset

variant_df <- readRDS('Version4_scores_goterm_data_470K.rds')
# variant_df <- data
variant_df <- variant_df %>%
  mutate(
    phylop_raw  = phylop,
    phylop_clip = pmax(pmin(phylop, 10), -10)
  )


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


### Create Shinydashboard


ui <- dashboardPage(
  dashboardHeader(title = "Variant Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("HomePage", tabName = "about", icon = icon("info-circle")), 
      menuItem("Variant Information", tabName = "variant_info", icon = icon("dna")),
      menuItem("Variant Effect Prediction", tabName = "pathogenicity", icon = icon("exclamation-triangle")),
      menuItem("ML Explorer", tabName = "ml_explorer", icon = icon("search")),
      menuItem("Ranking", tabName = "ranking", icon = icon("sort-numeric-desc")),
      menuItem("Distributions", tabName = "distributions", icon = icon("bar-chart")),
      menuItem("Variant Detail", tabName = "variant_detail", icon = icon("info-circle")),
      menuItem("Predict on Uploaded Data", tabName = "predict_upload", icon = icon("cloud-upload-alt"))
    )
  )
  ,
  dashboardBody(
    tabItems(
      # Page: Predict on Uploaded Data
      tabItem(tabName = "predict_upload",
              h2("Predict on Uploaded Data"),
              fluidRow(
                box(width = 4,
                    title = "Upload Data",
                    status = "primary",
                    solidHeader = TRUE,
                    helpText(
                      "To generate predictions and receive full quality control (QC) and interpretation, your CSV must contain the following columns:",
                      tags$ul(
                        tags$li("chrom"),
                        tags$li("pos"),
                        tags$li("REF_ALLELE"),
                        tags$li("Allele"),
                        tags$li("Gene"),
                        tags$li("SYMBOL"),
                        tags$li("STRAND"),
                        tags$li("Consequence"),
                        tags$li("SIFT"),
                        tags$li("PolyPhen"),
                        tags$li("fathmm"),
                        tags$li("CADD_PHRED"),
                        tags$li("SIFT_score"),
                        tags$li("FATHMM_score"),
                        tags$li("PolyPhen_score"),
                        tags$li("polyphen_max"),
                        tags$li("revel_max"),
                        tags$li("spliceai_ds_max"),
                        tags$li("phylop")
                      ),
                      "All columns above are required for the prediction and QC functions to run without error. For best results, ensure your data includes these features. See documentation for details."
                    ),
                    downloadButton("download_template_csv", "Download CSV Template", icon = icon("file-csv")),
                    br(),
                    fileInput("user_file", "Upload CSV File", accept = ".csv"),
                    actionButton("run_prediction", "Run Prediction", icon = icon("play")),
                    br(),
                    downloadButton("download_pred_results", "Download Results")
                ),
                box(width = 8,
                    title = "Prediction Results",
                    status = "info",
                    solidHeader = TRUE,
                    tabsetPanel(
                      tabPanel("Raw Results", DTOutput("user_pred_table")),
                      tabPanel("Pathogenicity QC", DTOutput("user_patho_summary")),
                      tabPanel("GVPS QC", DTOutput("user_gvps_summary"))
                    )
                )
              )
      ),
      
      # Page: About
      tabItem(
        tabName = "about",
        
        h2("Welcome to the Variant Explorer Dashboard"),
        p("This dashboard is designed to allow researchers or clinicians to explore genetic variants, retrieve annotations, and evaluate predicted effect and pathogenicity of variants."),
        tags$hr(),
        
        h2("About This Tool"),
        
        box(
          title = "Purpose",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          p(
            "This tool is designed to make it easier and quicker for clinicians, researchers, and students ",
            "to access and aggregate information about genomic variants. It provides integrated genomic, ",
            "functional, and clinical annotations to support interpretation of potential clinical significance."
          )
        ),
        
        box(
          title = "How to Use",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          p(
            "Users can look up variants by chromosome and position using the Variant Information page. ",
            "The current version of this tool is restricted to the first 470,000 known single nucleotide variants on chromosome 1."
          )
        ),
        
        box(
          title = "Data Sources",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          
          p("This dashboard integrates genomic, annotation, and clinical data from several widely used databases:"),
          
          tags$ul(
            tags$li(
              tags$strong("Ensembl: "),
              a("ensembl.org", href = "https://www.ensembl.org/index.html"),
              " — genomic context, variant annotations, and functional prediction scores ",
              "(SIFT, PolyPhen, CADD)."
            ),
            tags$li(
              tags$strong("gnomAD: "),
              a("gnomad.broadinstitute.org", href = "https://gnomad.broadinstitute.org/"),
              " — population allele frequencies and quality control metrics."
            ),
            tags$li(
              tags$strong("Gene Ontology: "),
              a("geneontology.org", href = "https://geneontology.org/"),
              " — structured descriptions of gene functions, biological processes, and cellular components."
            ),
            tags$li(
              tags$strong("ClinVar: "),
              a("ncbi.nlm.nih.gov/clinvar", href = "https://www.ncbi.nlm.nih.gov/clinvar/"),
              " — curated clinical interpretations (Pathogenic, Benign, Uncertain significance)."
            )
          )
        ),
        
        box(
          title = "ML Predictor",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          p(
            "Existing pathogenicity scores provide useful information but often lack consensus. ",
            "The machine learning predictor implemented in this tool integrates multiple established ",
            "scores to provide a more robust assessment of variant pathogenicity."
          )
        )
      ),
      
      
      # Page: ML Explorer
      tabItem(tabName = "ml_explorer",
              fluidRow(
                column(3,
                       box(
                         title = "Filters",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
                         selectInput("risk_filter", "GVPS Risk", choices = c("All", levels(ml_data$GVPS_risk)), selected = "All"),
                         selectInput("path_filter", "Pathogenicity Class", choices = c("All", levels(ml_data$Pathogenicity_Level)), selected = "All"),
                         sliderInput("gvps_range", "GVPS range", min(ml_data$GVPS_pred, na.rm = TRUE), max(ml_data$GVPS_pred, na.rm = TRUE), value = range(ml_data$GVPS_pred, na.rm = TRUE)),
                         sliderInput("percentile_range", "GVPS Percentile", min(ml_data$GVPS_percentile, na.rm = TRUE), max(ml_data$GVPS_percentile, na.rm = TRUE), value = range(ml_data$GVPS_percentile, na.rm = TRUE)),
                         hr(),
                         selectInput("scatter_x", "Scatterplot X", choices = names(ml_data), selected = "GVPS_pred"),
                         selectInput("scatter_y", "Scatterplot Y", choices = names(ml_data), selected = "GVPS_percentile"),
                         selectInput("scatter_col", "Color by", choices = c("None", names(ml_data)), selected = "GVPS_risk"),
                         downloadButton("download_filtered", "Download filtered")
                       )
                ),
                column(9,
                       tabBox(
                         width = 12,
                         title = NULL,
                         tabPanel(
                           title = "Table",
                           box(
                             width = 12,
                             style = "overflow-x: auto;",
                             DTOutput("variant_table")
                           )
                         ),
                         tabPanel(
                           title = "Scatter Plot",
                           box(
                             width = 12,
                             style = "overflow-x: auto;",
                             plotlyOutput("explorer_scatter")
                           )
                         )
                       )
                )
              )
      ),
      
      # Page: Ranking
      tabItem(tabName = "ranking",
              fluidRow(
                column(3,
                       box(
                         title = "Ranking Controls",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
                         selectInput("rank_metric", "Rank by", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                         numericInput("rank_n", "Top N", value = 25, min = 1, max = 1000),
                         radioButtons("rank_order", "Order", choices = c("Descending" = "desc", "Ascending" = "asc"), selected = "desc")
                       )
                ),
                column(9,
                       box(
                         title = "Ranking Table",
                         status = "info",
                         solidHeader = TRUE,
                         width = 12,
                         style = "overflow-x: auto;",
                         DTOutput("ranking_table")
                       )
                )
              )
      ),
      
      # Page: Distributions
      tabItem(tabName = "distributions",
              fluidRow(
                box(width = 12, style = "overflow-x: auto;", plotlyOutput("gvps_hist"))
              ),
              fluidRow(
                box(width = 12, style = "overflow-x: auto;", plotlyOutput("path_bar"))
              ),
              fluidRow(
                box(width = 12, style = "overflow-x: auto;", plotlyOutput("risk_bar"))
              ),
              fluidRow(
                box(width = 12, style = "overflow-x: auto;", plotlyOutput("percentile_hist"))
              ),
              hr(),
              fluidRow(
                box(width = 12,
                    h5("GVPS ECDF Plot"),
                    style = "overflow-x: auto;",
                    plotlyOutput("gvps_ecdf_plot")
                )
              ),
              fluidRow(
                box(width = 12,
                    h5("Feature Distribution by Class"),
                    selectInput("feature_select", "Feature (numeric)", choices = names(ml_data)[sapply(ml_data, is.numeric)], selected = "GVPS_pred"),
                    style = "overflow-x: auto;",
                    plotlyOutput("feature_by_class_plot")
                )
              ),
              fluidRow(
                box(width = 12,
                    h5("Top Gene/Consequence Barplot"),
                    selectInput("bar_type", "Barplot Type", choices = c("Gene" = "SYMBOL", "Consequence" = "Consequence")),
                    style = "overflow-x: auto;",
                    plotlyOutput("top_barplot")
                )
              )
      ),
      
      
      # Page: Variant Detail
      tabItem(tabName = "variant_detail",
              fluidRow(
                box(width = 4,
                    selectizeInput("variant_select", "Select rsID", choices = NULL, selected = NULL, options = list(placeholder = "Type to search"))
                ),
                box(
                  title = "Variant Summary",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 8,
                  #h5("Variant summary"),
                  style = "overflow-x: auto;",
                  tableOutput("variant_summary"),
                  h5("Raw fields"),
                  tableOutput("variant_detail")
                )
              )
      ),
      
      # Page: Variant Information 
      tabItem(
        tabName = "variant_info",
        h2("Variant Information"),
        fluidRow(
          box(
            width = 12,
            title = "Input Variant",
            status = "primary",
            solidHeader = TRUE,
            textInput("chrom_input", "Chromosome:", value = "1"),
            numericInput("pos_input", "Position:", value = "12345"),
            actionButton("lookup_btn", "Lookup Variant")
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Variant Details",
            status = "info",
            solidHeader = TRUE,
            uiOutput("variant_card_ui"),
            uiOutput("genomic_out"),
            uiOutput("annotation_out"),
            uiOutput("clinical_out")
          )
        ),
        
        fluidRow(
          box(
            width = 8,
            title = "Genomic Location",
            status = "success",
            solidHeader = TRUE,
            plotOutput("genomic_mini", height = "130px")
          ),
          box(
            width = 4,
            title = "Explore in IGB",
            status = "warning",
            solidHeader = TRUE,
            p("Copy this region into IGB to explore the local genomic context:"),
            verbatimTextOutput("igb_region_text"),
            tags$small(em("Region: ±50kb around the selected variant"))
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Literature (PubMed)",
            status = "primary",
            solidHeader = TRUE,
            actionButton("fetch_pubmed", "Get PubMed literature for selected variant", icon = icon("book")),
            tags$small(style="margin-left:8px; color:#666;",
                       "Fetches on demand and caches results locally."),
            tags$div(style="margin-top:10px;", uiOutput("pubmed_ui"))
          )
          
        ),
        fluidRow(
          box(
            width = 12,
            title = "Gene Ontology Overview",
            status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6, plotlyOutput("go_treemap", height = "300px")),   # Left: IMAGE
              column(6, uiOutput("go_term_list"))                        # Right: TERMS
            )
          )
        )
      ),
      
      # Page: Variant Effect Prediction
      tabItem(
        tabName = "pathogenicity",
        h2("Variant Effect Prediction"),
        fluidRow(
          box(
            width = 4,
            title = "Prediction Settings",
            status = "primary",
            solidHeader = TRUE,
            selectInput(
              "score_type",
              "Select Pathogenicity Score:",
              choices = c("SIFT", "PolyPhen", "FATHMM", "CADD", "GVPS")
            ),
            actionButton("predict_btn", "Show Score")
          ),
          box(
            width = 8,
            title = "Prediction Result",
            status = "warning",
            solidHeader = TRUE,
            verbatimTextOutput("score_out")
          )
        ),
        box(
          width = 8,
          title = "Consensus Interpretation",
          status = "info",
          solidHeader = TRUE,
          uiOutput("consensus_out")
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "All Pathogenicity Scores",
            status = "info",
            solidHeader = TRUE,
            plotOutput("score_plot", height = "300px")
          )
        )
      )
    )
  )
)




split_terms <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  x |>
    strsplit(";") |>
    unlist() |>
    trimws()
}

pubmed_results_ui <- function(df, query) {
  if (nrow(df) == 0) {
    return(tags$div(
      tags$b("PubMed results"),
      tags$p(style="margin:6px 0;", "No results found for:"),
      tags$code(query)
    ))
  }
  
  tags$div(
    tags$div(
      style="display:flex; justify-content:space-between; align-items:center;",
      tags$b("PubMed results"),
      tags$a("Open full search in PubMed",
             href = paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", URLencode(query, reserved = TRUE)),
             target = "_blank")
    ),
    tags$hr(style="margin:8px 0;"),
    lapply(seq_len(nrow(df)), function(i) {
      tags$div(
        style="margin-bottom:10px; padding-bottom:10px; border-bottom:1px solid #eee;",
        tags$a(df$title[i], href = df$url[i], target = "_blank", style="font-weight:600;"),
        tags$div(style="color:#666; font-size:0.9em;",
                 paste0(df$source[i], " • ", df$year[i])),
        if (nzchar(df$authors[i])) tags$div(style="color:#444; font-size:0.85em;", df$authors[i])
      )
    })
  )
}



go_term_list_ui <- function(bp, mf, cc) {
  
  htmltools::tags$div(
    style = "max-height: 300px; overflow-y: auto; padding: 6px;",
    
    # BP
    htmltools::tags$h4("Biological Process"),
    htmltools::tags$ul(
      lapply(bp, function(t) htmltools::tags$li(t))
    ),
    
    # MF
    htmltools::tags$h4("Molecular Function"),
    htmltools::tags$ul(
      lapply(mf, function(t) htmltools::tags$li(t))
    ),
    
    # CC
    htmltools::tags$h4("Cellular Component"),
    htmltools::tags$ul(
      lapply(cc, function(t) htmltools::tags$li(t))
    )
  )
}


info_card_ui <- function(title, items) {
  # items = named list: list("Label" = value, ...)
  safe <- function(x) {
    if (length(x) == 0 || is.null(x) || is.na(x) || x == "") "—" else as.character(x)
  }
  
  htmltools::tags$div(
    style = "border:1px solid #eee; border-radius:10px; padding:10px; margin-bottom:10px; background:#fff;",
    htmltools::tags$h4(style="margin:0 0 8px 0;", title),
    htmltools::tags$dl(
      style="margin:0; display:grid; grid-template-columns: 160px 1fr; row-gap:6px; column-gap:10px;",
      lapply(names(items), function(nm) {
        htmltools::tagList(
          htmltools::tags$dt(style="font-weight:600; color:#555;", nm),
          htmltools::tags$dd(style="margin:0; color:#111;", safe(items[[nm]]))
        )
      })
    )
  )
}



clinsig_bucket <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (is.na(x) || x == "" || x == "—") return("Unknown")
  
  if (grepl("pathogenic", x) && grepl("likely", x)) return("Likely Pathogenic")
  if (grepl("\\bpathogenic\\b", x)) return("Pathogenic")
  if (grepl("benign", x) && grepl("likely", x)) return("Likely Benign")
  if (grepl("\\bbenign\\b", x)) return("Benign")
  if (grepl("uncertain", x) || grepl("vus", x)) return("VUS")
  
  "Other"
}

clinsig_traffic_ui <- function(clnsig) {
  lvl <- clinsig_bucket(clnsig)
  
  levels_order <- c("Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic")
  colors <- c(
    "Benign"            = "#2ca02c",
    "Likely Benign"     = "#98df8a",
    "VUS"               = "#ffdd57",
    "Likely Pathogenic" = "#ff8c42",
    "Pathogenic"        = "#d62728"
  )
  
  htmltools::tags$div(
    style = "margin-top:10px;",
    htmltools::tags$div(
      style = "display:flex; gap:6px; align-items:center;",
      lapply(levels_order, function(L) {
        active <- identical(lvl, L)
        htmltools::tags$div(
          title = L,
          style = paste0(
            "flex:1; height:14px; border-radius:8px; background:", colors[[L]], ";",
            if (active) "outline:3px solid #111; outline-offset:1px;" else "opacity:0.35;",
            "transition:0.2s;"
          )
        )
      })
    ),
    htmltools::tags$div(
      style = "display:flex; justify-content:space-between; margin-top:6px; font-size:0.78em; color:#555;",
      htmltools::tags$span("Benign"),
      htmltools::tags$span("Pathogenic")
    ),
    htmltools::tags$div(
      style = "margin-top:6px; font-size:0.9em;",
      htmltools::tags$b("ClinVar: "),
      htmltools::tags$span(ifelse(lvl %in% c("Unknown","Other"), as.character(clnsig), lvl))
    )
  )
}




### Server
server <- function(input, output, session) {
  # --- Downloadable CSV template for user upload ---
  output$download_template_csv <- downloadHandler(
    filename = function() "prediction_template.csv",
    content = function(file) {
      template <- data.frame(
        chrom = "1",
        pos = 12345,
        REF_ALLELE = "A",
        Allele = "G",
        Gene = "GENE1",
        SYMBOL = "GENE1",
        STRAND = "+",
        Consequence = "missense_variant",
        SIFT = 0.01,
        PolyPhen = 0.8,
        fathmm = -2.5,
        CADD_PHRED = 25.1,
        SIFT_score = 0.99,
        FATHMM_score = 0.85,
        PolyPhen_score = 0.92,
        polyphen_max = 0.95,
        revel_max = 0.88,
        spliceai_ds_max = 0.12,
        phylop = 2.1
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
  
  # --- Predict on Uploaded Data ---
  user_data <- reactive({
    req(input$user_file)
    tryCatch({
      read.csv(input$user_file$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      showNotification("Error reading uploaded file", type = "error"); NULL
    })
  })
  
  pubmed_payload <- eventReactive(input$fetch_pubmed, {
    req(variant_row())
    vr <- variant_row()
    
    query <- make_pubmed_query(vr)
    
    cache <- read_pubmed_cache()
    key <- query
    
    if (!is.null(cache[[key]])) {
      return(list(query = query, results = cache[[key]]))
    }
    
    res <- pubmed_fetch(query, retmax = 20)
    
    cache[[key]] <- res
    write_pubmed_cache(cache)
    
    list(query = query, results = res)
  })
  
  output$pubmed_ui <- renderUI({
    p <- pubmed_payload()
    req(p)
    
    pubmed_results_ui(p$results, p$query)
  })
  
  
  user_pred_results <- eventReactive(input$run_prediction, {
    df <- user_data()
    if (is.null(df)) return(NULL)
    results <- df
    # --- Pathogenicity Dashboard Report ---
    patho_report <- NULL
    if (exists("generate_dashboard_predictions")) {
      try({
        patho_report <- generate_dashboard_predictions(df)
      }, silent = TRUE)
    }
    # --- GVPS Dashboard Report ---
    gvps_report <- NULL
    if (exists("generate_gvps_dashboard_report")) {
      try({
        gvps_report <- generate_gvps_dashboard_report(df)
      }, silent = TRUE)
    }
    # Attach to results for downstream use
    attr(results, "patho_report") <- patho_report
    attr(results, "gvps_report") <- gvps_report
    results
  })
  
  # --- Pathogenicity Dashboard Summary Table ---
  output$user_patho_summary <- renderDT({
    req(user_pred_results())
    patho_report <- attr(user_pred_results(), "patho_report")
    if (!is.null(patho_report) && !is.null(patho_report$summary)) {
      datatable(patho_report$summary, options = list(scrollX = TRUE, pageLength = 10))
    }
  })
  
  # --- GVPS Dashboard Summary Table ---
  output$user_gvps_summary <- renderDT({
    req(user_pred_results())
    gvps_report <- attr(user_pred_results(), "gvps_report")
    if (!is.null(gvps_report) && !is.null(gvps_report$summary)) {
      datatable(gvps_report$summary, options = list(scrollX = TRUE, pageLength = 10))
    }
  })
  
  output$user_pred_table <- renderDT({
    req(user_pred_results())
    datatable(user_pred_results(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$download_pred_results <- downloadHandler(
    filename = function() paste0("predictions_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(user_pred_results(), file, row.names = FALSE)
    }
  )
  # Reactive: Create Variant object
  variant_obj <- eventReactive(input$lookup_btn, {
    row <- variant_df %>%
      filter(chrom == input$chrom_input, pos == input$pos_input)
    
    if (nrow(row) == 0) {
      showNotification("Variant not found!", type = "error")
      return(NULL)
    }
    
    # Variant(row[1, ])
    # Variant(row[1,c(chrom, pos, REF_ALLELE, Allele, Gene, SYMBOL, STRAND, Consequence, CLIN_SIG, SIFT, PolyPhen, fathmm, CADD_PHRED)])
    Variant(row[1, c("chrom", "pos", "REF_ALLELE", "Allele", "Gene", "SYMBOL", "STRAND", "Consequence", "CLIN_SIG", "SIFT", "PolyPhen", "fathmm", "CADD_PHRED","GVPS_pred")])
  })
  
  
  
  variant_row <- eventReactive(input$lookup_btn, {
    row <- variant_df %>%
      filter(chrom == input$chrom_input, pos == input$pos_input)
    
    if (nrow(row) == 0) {
      return(NULL)
    }
    
    row[1, ]
  })
  
  
  output$variant_card_ui <- renderUI({
    req(variant_row())
    variant_card(variant_row())
  })
  
  ## AF donut plot used inside the card
  output$variant_af_donut <- renderPlot({
    req(variant_row())
    # change AF to whatever your AF column is called
    af_value <- variant_row()$AF.y
    af_donut_plot(af_value)
  })
  
  ### Display Output: Genomic + Annotation + Clinical
  output$genomic_out <- renderUI({
    req(variant_obj())
    g <- variant_obj()$genomic
    
    info_card_ui("Genomic Information", list(
      "Chromosome" = g$chromosome,
      "Position"   = g$position,
      "Ref"        = g$ref,
      "Alt"        = g$alt
    ))
  })
  
  output$annotation_out <- renderUI({
    req(variant_obj())
    a <- variant_obj()$annotation
    
    info_card_ui("Annotation Information", list(
      "Gene (Ensembl)" = a$gene,
      "Symbol"         = a$symbol,
      "Transcript"     = a$transcript,
      "Consequence"    = a$consequence
    ))
  })
  
  output$clinical_out <- renderUI({
    req(variant_obj())
    c <- variant_obj()$clinical
    
    info_card_ui("Clinical Information", list(
      "Clinical significance" = c$phenotype
    ))
  })
  
  
  ### Plot: Genomic Position
  output$genomic_mini <- renderPlot({
    
    req(variant_row())
    vr <- variant_row()
    
    # 🔧 FIX: shrink margins for Shiny box
    par(mar = c(3, 2, 2, 1))  # bottom, left, top, right
    
    chrom <- as.character(vr$chrom)
    pos   <- as.numeric(vr$pos)
    
    window <- 5e4
    start  <- max(pos - window, 1)
    end    <- pos + window
    
    plot(c(start, end), c(0, 0),
         type = "n",
         axes = FALSE,
         xlab = "",
         ylab = "",
         ylim = c(-0.5, 0.6))
    
    # backbone
    segments(start, 0, end, 0, lwd = 6, lend = "round", col = "#4C72B0")
    
    # ticks
    ticks <- pretty(c(start, end), n = 4)
    ticks <- ticks[ticks >= start & ticks <= end]
    segments(ticks, -0.08, ticks, 0.08)
    text(ticks, -0.2, format(ticks, big.mark=","), cex = 0.7)
    
    # variant marker
    segments(pos, -0.35, pos, 0.45, lwd = 2, col = "#DD8452")
    points(pos, 0, pch = 19, cex = 1.2, col = "#DD8452")
    
    # label
    text(start, 0.5,
         paste0("chr", chrom, ": ", format(start, big.mark=","), "–", format(end, big.mark=",")),
         adj = c(0, 0),
         cex = 0.8)
  },
  height = 130)  
  
  
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  
  parse_label_value <- function(x) {
    # numeric (e.g., CADD = 26.6, FATHMM = 0.91)
    if (is.numeric(x)) {
      return(list(label = NA_character_, value = as.numeric(x), raw = as.character(x)))
    }
    
    x <- as.character(x)
    if (length(x) == 0 || is.na(x) || x == "") {
      return(list(label = NA_character_, value = NA_real_, raw = x))
    }
    
    # label(value) e.g., "probably_damaging(1)"
    m <- regexec("^\\s*([^\\(]+)\\s*\\(([-0-9\\.eE]+)\\)\\s*$", x)
    reg <- regmatches(x, m)[[1]]
    
    if (length(reg) == 3) {
      return(list(
        label = trimws(reg[2]),
        value = suppressWarnings(as.numeric(reg[3])),
        raw   = x
      ))
    }
    
    # fallback: plain number in a string
    v <- suppressWarnings(as.numeric(x))
    if (!is.na(v)) {
      return(list(label = NA_character_, value = v, raw = x))
    }
    
    # final fallback
    list(label = trimws(x), value = NA_real_, raw = x)
  }
  
  interpret_tool_vote <- function(tool, label, value) {
    tool <- tolower(tool %||% "")
    lab  <- tolower(label %||% "")
    
    # SIFT strings: deleterious / deleterious_low_confidence / tolerated
    if (tool == "sift") {
      if (grepl("deleterious", lab)) return("Damaging")
      if (grepl("tolerated", lab))   return("Benign/Low")
      return("Uncertain")
    }
    
    # PolyPhen strings: probably_damaging / possibly_damaging / benign
    if (tool == "polyphen") {
      if (grepl("probably_damaging|possibly_damaging", lab)) return("Damaging")
      if (grepl("^benign", lab)) return("Benign/Low")
      return("Uncertain")
    }
    
    # CADD numeric (PHRED-like): higher = more deleterious
    if (tool == "cadd") {
      if (is.na(value)) return("Missing")
      if (value >= 20)  return("Damaging")
      if (value < 10)   return("Benign/Low")
      return("Uncertain")
    }
    
    # FATHMM numeric: depends on your scale/direction.
    # Common: higher => more damaging. If your results look flipped, tell me and we invert.
    if (tool == "fathmm") {
      if (is.na(value)) return("Missing")
      if (value >= 0.5) return("Damaging")
      if (value < 0.5)  return("Benign/Low")
      return("Uncertain")
    }
    
    if (is.na(label) && is.na(value)) return("Missing")
    "Uncertain"
  }
  
  
  
  ### Variant Effect Prediction
  selected_score <- eventReactive(input$predict_btn, {
    req(variant_obj(), input$score_type)
    
    # pull the selected tool from the pathogenicity list
    variant_obj()$pathogenicity[[input$score_type]]
  })
  
  output$score_out <- renderPrint({
    req(selected_score())
    
    # If GVPS is numeric, this prints just the number.
    # If you stored it as "label(value)", it prints that string.
    selected_score()
  })
  
  output$score_plot <- renderPlot({
    req(variant_obj())
    
    raw_all <- variant_obj()$pathogenicity
    req(!is.null(raw_all), length(raw_all) > 0)
    
    # remove GVPS from plot
    raw_list <- raw_all[setdiff(names(raw_all), "GVPS")]
    
    parsed <- lapply(raw_list, parse_label_value)
    
    df <- data.frame(
      Score = names(raw_list),
      Value = vapply(parsed, `[[`, numeric(1), "value"),
      stringsAsFactors = FALSE
    )
    
    ggplot(df, aes(x = Score, y = Value, fill = Score)) +
      geom_col(na.rm = TRUE) +
      theme_minimal() +
      labs(title = "Pathogenicity Score Comparison", y = "Score Value")
  })
  
  
  # Consensus interpretation panel (Feature 1)
  consensus <- eventReactive(input$predict_btn, {
    req(variant_obj())
    
    raw_all <- variant_obj()$pathogenicity
    req(!is.null(raw_all), length(raw_all) > 0)
    
    # remove GVPS from consensus
    raw_list <- raw_all[setdiff(names(raw_all), "GVPS")]
    
    parsed <- lapply(raw_list, parse_label_value)
    
    df <- data.frame(
      Tool  = names(raw_list),
      Raw   = vapply(parsed, `[[`, character(1), "raw"),
      Label = vapply(parsed, `[[`, character(1), "label"),
      Value = vapply(parsed, `[[`, numeric(1),  "value"),
      stringsAsFactors = FALSE
    )
    
    df$Vote <- mapply(
      interpret_tool_vote,
      tool  = df$Tool,
      label = df$Label,
      value = df$Value,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
    
    n_available <- sum(df$Vote != "Missing")
    n_damaging  <- sum(df$Vote == "Damaging")
    n_benign    <- sum(df$Vote == "Benign/Low")
    n_uncertain <- sum(df$Vote == "Uncertain")
    
    consensus_label <- dplyr::case_when(
      n_available == 0 ~ "No consensus (no scores available)",
      n_damaging >= 3  ~ "Likely Pathogenic (majority damaging)",
      n_benign >= 3    ~ "Likely Benign/Low impact (majority benign)",
      TRUE             ~ "Uncertain (mixed evidence)"
    )
    
    agree_prop <- if (n_available > 0) max(n_damaging, n_benign) / n_available else NA_real_
    confidence <- dplyr::case_when(
      is.na(agree_prop)  ~ "NA",
      agree_prop >= 0.75 ~ "High",
      agree_prop >= 0.60 ~ "Moderate",
      TRUE               ~ "Low"
    )
    
    list(
      label = consensus_label,
      agreement = paste0(n_damaging, " damaging / ", n_available, " available (", n_uncertain, " uncertain)"),
      confidence = confidence,
      df = df
    )
  })
  
  output$consensus_out <- renderUI({
    req(consensus())
    cns <- consensus()
    
    tagList(
      tags$h4(style = "margin-top:0; margin-bottom:6px;", cns$label),
      tags$p(style = "margin:0;", tags$strong("Agreement: "), cns$agreement),
      tags$p(style = "margin:0;", tags$strong("Confidence: "), cns$confidence),
      tags$hr(style = "margin:10px 0;"),
      tags$small("Per-tool votes:"),
      tableOutput("vote_table")
    )
  })
  
  output$vote_table <- renderTable({
    req(consensus())
    consensus()$df[, c("Tool", "Raw", "Vote")]
  }, striped = TRUE, spacing = "xs", width = "100%")
  
  
  output$igb_region_text <- renderText({
    req(variant_obj())
    chrom <- variant_obj()$genomic$chromosome
    pos   <- as.numeric(variant_obj()$genomic$position)
    
    # window around the variant for IGB
    window <- 5e4  # 50,000 bp on each side
    start  <- max(pos - window, 1L)
    end    <- pos + window
    sprintf("chr%s:%d-%d", chrom, start, end)
  })
  
  
  go_terms_raw <- reactive({
    req(variant_row())        # use the full merged row
    vr <- variant_row()
    
    
    gene_symbol <- vr$SYMBOL
    
    go_row <- vr |>
      filter(hgnc_symbol == gene_symbol) |>
      slice_head(n = 1)
    
    if (nrow(go_row) == 0) {
      return(NULL)
    }
    
    list(
      bp = split_terms(go_row$biological_process),
      mf = split_terms(go_row$molecular_function),
      cc = split_terms(go_row$cellular_component)
    )
  })
  
  
  
  output$go_term_list <- renderUI({
    terms <- go_terms_raw()
    req(terms)
    
    go_term_list_ui(terms$bp, terms$mf, terms$cc)
  })
  
  go_terms_for_variant <- reactive({
    req(variant_row())
    vr <- variant_row()
    
    bp <- if ("GO_BP_count" %in% names(vr)) vr$GO_BP_count else NA
    mf <- if ("GO_MF_count" %in% names(vr)) vr$GO_MF_count else NA
    cc <- if ("GO_CC_count" %in% names(vr)) vr$GO_CC_count else NA
    
    # if everything is missing or zero, bail
    if (all(is.na(c(bp, mf, cc))) || all(c(bp, mf, cc) == 0, na.rm = TRUE)) {
      return(NULL)
    }
    
    tibble::tibble(
      label  = c("GO",              "Molecular Function", "Biological Process", "Cellular Component"),
      parent = c("",                "GO",                 "GO",                 "GO"),
      value  = c(NA,                as.numeric(mf),       as.numeric(bp),       as.numeric(cc))
    )
  })
  
  output$go_treemap <- renderPlotly({
    df <- go_terms_for_variant()
    
    if (is.null(df)) {
      return(
        plotly::plot_ly(
          type = "scatter",
          x = 0, y = 0, mode = "text",
          text = "No GO summary available for this variant",
          textposition = "middle center"
        ) %>%
          plotly::layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE)
          )
      )
    }
    
    plot_ly(
      data = df,
      type = "sunburst",
      labels = ~label,
      parents = ~parent,
      values = ~value,
      branchvalues = "total"
    )
  })
  
  # Appended server dashboard
  
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


### run

shinyApp(ui, server)
