#!/usr/bin/env Rscript
################################################################################
# BIOPAC ACQ Batch Converter
################################################################################
# Converts AcqKnowledge (.acq) files to text format
# Supports versions 3.x - 5.x (uncompressed)
# 
# ACQ parser ported from bioread (Python) by Nate Vack, UW-Madison
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(shiny, shinyjs, shinyFiles, plotly)

# =============================================================================
# ACQ PARSER
# =============================================================================

V_370  <- 38L
V_400B <- 61L
V_430  <- 124L

VERSION_MAP <- c(

"30"="2.0a","34"="3.0","38"="3.7.0","41"="3.8.1","45"="3.9.0",
"61"="4.0.0B","68"="4.0.0","76"="4.0.1","83"="4.1.0","108"="4.2.0",
"124"="4.3.0","128"="4.4.0","132"="5.0.1"
)

read_int <- function(con, size=2, signed=TRUE, endian="little") {
raw_data <- readBin(con, "raw", n=size)
if (length(raw_data) < size) return(NA)
if (endian == "big") raw_data <- rev(raw_data)
val <- sum(as.integer(raw_data) * (256^(0:(size-1))))
if (signed && size > 1 && val >= 256^size/2) val <- val - 256^size
val
}

read_double <- function(con, endian="little") readBin(con, "double", n=1, size=8, endian=endian)

read_string <- function(con, n) {
raw_data <- readBin(con, "raw", n=n)
if (length(raw_data) < n) return("")
null_pos <- which(raw_data == 0x00)[1]
if (!is.na(null_pos) && null_pos > 1) raw_data <- raw_data[1:(null_pos-1)]
else if (!is.na(null_pos) && null_pos == 1) return("")
tryCatch(rawToChar(raw_data), error=function(e) "")
}

skip_bytes <- function(con, n) if (n > 0) readBin(con, "raw", n=n)

read_acq <- function(filename) {
if (!file.exists(filename)) stop(paste("File not found:", filename))
con <- file(filename, "rb")
on.exit(close(con))

seek(con, 2)
version_le <- read_int(con, 4, signed=FALSE, endian="little")
seek(con, 2)
version_be <- read_int(con, 4, signed=FALSE, endian="big")
byte_order <- if (version_le <= version_be) "little" else "big"
revision <- if (version_le <= version_be) version_le else version_be

seek(con, 0)
skip_bytes(con, 6)
lExtItemHeaderLen <- read_int(con, 4, signed=TRUE, endian=byte_order)
nChannels <- read_int(con, 2, signed=FALSE, endian=byte_order)
skip_bytes(con, 4)
dSampleTime <- read_double(con, endian=byte_order)

if (revision >= V_400B) {
  seek(con, 972)
  bCompressed <- read_int(con, 4, signed=TRUE, endian=byte_order)
  if (bCompressed != 0) stop("Compressed files not supported")
  seek(con, 2398)
  hExpectedPaddings <- if (revision >= V_430) read_int(con, 2, signed=TRUE, endian=byte_order) else 0
} else {
  hExpectedPaddings <- 0
}

samples_per_sec <- 1000 / dSampleTime
padding_offset <- lExtItemHeaderLen

if (revision >= V_430 && hExpectedPaddings > 0) {
  seek(con, lExtItemHeaderLen)
  for (p in seq_len(hExpectedPaddings)) {
    pad_len <- read_int(con, 4, signed=TRUE, endian=byte_order)
    padding_offset <- padding_offset + pad_len
  }
}

channels <- list()
pos <- padding_offset
for (i in seq_len(nChannels)) {
  seek(con, pos)
  ch <- list()
  ch$lChanHeaderLen <- read_int(con, 4, signed=TRUE, endian=byte_order)
  skip_bytes(con, 2)
  ch$name <- trimws(read_string(con, 40))
  skip_bytes(con, 22)
  ch$units <- trimws(read_string(con, 20))
  ch$point_count <- read_int(con, 4, signed=TRUE, endian=byte_order)
  ch$raw_scale <- read_double(con, endian=byte_order)
  ch$raw_offset <- read_double(con, endian=byte_order)
  skip_bytes(con, 4)
  
  if (revision >= V_400B) {
    skip_bytes(con, 40)
    ch$frequency_divider <- read_int(con, 2, signed=TRUE, endian=byte_order)
  } else if (revision >= V_370) {
    skip_bytes(con, 138)
    ch$frequency_divider <- read_int(con, 2, signed=TRUE, endian=byte_order)
  } else {
    ch$frequency_divider <- 1
  }
  if (is.na(ch$frequency_divider) || ch$frequency_divider < 1) ch$frequency_divider <- 1
  
  channels[[i]] <- ch
  pos <- pos + ch$lChanHeaderLen
}

seek(con, pos)
foreign_len <- if (revision >= V_400B) read_int(con, 4, TRUE, byte_order) else read_int(con, 2, TRUE, byte_order)

start_offset <- pos + foreign_len
dtype_headers <- NULL
data_start <- NULL

for (i in 0:4095) {
  seek(con, start_offset + i)
  valid <- TRUE
  headers <- list()
  for (ch in seq_len(nChannels)) {
    nSize <- read_int(con, 2, TRUE, byte_order)
    nType <- read_int(con, 2, TRUE, byte_order)
    expected <- if (nType %in% c(0,1)) 8 else if (nType == 2) 2 else NA
    if (is.na(expected) || nSize != expected) { valid <- FALSE; break }
    headers[[ch]] <- list(sample_size=nSize, dtype=nType)
  }
  if (valid && length(headers) == nChannels) {
    dtype_headers <- headers
    data_start <- start_offset + i + (nChannels * 4)
    break
  }
}
if (is.null(data_start)) stop("Could not locate data")

for (i in seq_len(nChannels)) {
  channels[[i]]$sample_size <- dtype_headers[[i]]$sample_size
  channels[[i]]$samples_per_sec <- samples_per_sec / channels[[i]]$frequency_divider
}

seek(con, data_start)
dividers <- sapply(channels, function(ch) ch$frequency_divider)
sizes <- sapply(channels, function(ch) ch$sample_size)
gcd <- function(a,b) { while(b!=0) {t<-b; b<-a%%b; a<-t}; a }
lcm <- function(a,b) (a*b) / gcd(a,b)
base_len <- Reduce(lcm, dividers)

sample_pat <- c()
for (slot in 0:(base_len-1)) {
  for (ch in seq_along(dividers)) {
    if ((slot %% dividers[ch]) == 0) sample_pat <- c(sample_pat, ch-1)
  }
}
byte_pattern <- rep(sample_pat, times=sizes[sample_pat+1])
pattern_len <- length(byte_pattern)
total_bytes <- sum(sapply(channels, function(ch) ch$point_count * ch$sample_size))

for (i in seq_len(nChannels)) {
  ch <- channels[[i]]
  channels[[i]]$raw_data <- if (ch$sample_size == 2) integer(ch$point_count) else numeric(ch$point_count)
  channels[[i]]$data_idx <- 1
}

chunk_size <- pattern_len * 256
bytes_read <- 0
while (bytes_read < total_bytes) {
  raw_bytes <- readBin(con, "raw", n=min(chunk_size, total_bytes - bytes_read))
  if (length(raw_bytes) == 0) break
  chunk_pat <- rep(byte_pattern, ceiling(length(raw_bytes)/pattern_len))[1:length(raw_bytes)]
  
  for (ch_idx in 0:(nChannels-1)) {
    ch <- channels[[ch_idx+1]]
    ch_bytes <- raw_bytes[chunk_pat == ch_idx]
    if (length(ch_bytes) == 0) next
    
    if (ch$sample_size == 2) {
      n_samp <- length(ch_bytes) %/% 2
      if (n_samp > 0) {
        vals <- readBin(ch_bytes, "integer", n=n_samp, size=2, signed=TRUE, endian=byte_order)
        start <- channels[[ch_idx+1]]$data_idx
        end <- min(start + length(vals) - 1, ch$point_count)
        channels[[ch_idx+1]]$raw_data[start:end] <- vals[1:(end-start+1)]
        channels[[ch_idx+1]]$data_idx <- end + 1
      }
    } else if (ch$sample_size == 8) {
      n_samp <- length(ch_bytes) %/% 8
      if (n_samp > 0) {
        vals <- readBin(ch_bytes, "double", n=n_samp, size=8, endian=byte_order)
        start <- channels[[ch_idx+1]]$data_idx
        end <- min(start + length(vals) - 1, ch$point_count)
        channels[[ch_idx+1]]$raw_data[start:end] <- vals[1:(end-start+1)]
        channels[[ch_idx+1]]$data_idx <- end + 1
      }
    }
  }
  bytes_read <- bytes_read + length(raw_bytes)
}

for (i in seq_len(nChannels)) {
  ch <- channels[[i]]
  channels[[i]]$data <- if (ch$sample_size == 2) (ch$raw_data * ch$raw_scale) + ch$raw_offset else ch$raw_data
}

version_str <- VERSION_MAP[as.character(revision)]
if (is.null(version_str) || is.na(version_str)) version_str <- paste("rev", revision)

list(
  filename = basename(filename),
  version = version_str,
  samples_per_sec = samples_per_sec,
  num_channels = nChannels,
  channels = channels
)
}

acq_to_df <- function(acq) {
max_points <- max(sapply(acq$channels, function(ch) ch$point_count))
df <- data.frame(time_sec = seq(0, by=1/acq$samples_per_sec, length.out=max_points))
for (i in seq_len(acq$num_channels)) {
  ch <- acq$channels[[i]]
  col_name <- gsub("[^[:alnum:]]", "_", ch$name)
  col_name <- gsub("_+", "_", gsub("^_|_$", "", col_name))
  if (nchar(col_name) == 0) col_name <- paste0("Ch", i)
  df[[col_name]] <- if (ch$point_count < max_points) {
    ch$data[floor(seq(1, ch$point_count, length.out=max_points))]
  } else ch$data[1:max_points]
}
df
}

build_metadata <- function(acq) {
lines <- c(
  paste0("# Source: ", acq$filename),
  paste0("# Version: ", acq$version),
  paste0("# Sample Rate: ", acq$samples_per_sec, " Hz"),
  paste0("# Channels: ", acq$num_channels),
  paste0("# Samples: ", acq$channels[[1]]$point_count),
  "#"
)
for (j in seq_len(acq$num_channels)) {
  ch <- acq$channels[[j]]
  lines <- c(lines, paste0("# Ch", j, ": ", ch$name, " (", ch$units, ")"))
}
c(lines, "#")
}

# =============================================================================
# STYLES
# =============================================================================

app_css <- '
@import url("https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap");
* { box-sizing: border-box; }
body {
font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
background: #f8f9fa; color: #212529; font-size: 13px; line-height: 1.5;
}
.container-fluid { padding: 16px 24px; }
.title-panel {
background: #ffffff; padding: 16px 24px; margin: -16px -24px 20px -24px;
border-bottom: 1px solid #dee2e6;
}
.title-panel h2 { font-size: 18px; font-weight: 600; color: #212529; margin: 0 0 4px 0; }
.title-panel .subtitle { font-size: 12px; color: #6c757d; }
h4 {
font-size: 11px; font-weight: 600; color: #212529; text-transform: uppercase;
letter-spacing: 0.75px; margin: 0 0 12px 0; padding-bottom: 8px;
border-bottom: 1px solid #dee2e6;
}
.well {
background: #ffffff; border: 1px solid #dee2e6; border-radius: 6px;
box-shadow: 0 1px 3px rgba(0,0,0,0.04); padding: 16px; margin-bottom: 12px;
}
.control-label { font-size: 11px; font-weight: 500; color: #495057; margin-bottom: 4px; display: block; }
.form-control {
font-size: 12px; padding: 8px 10px; border: 1px solid #ced4da;
border-radius: 4px; background: #ffffff;
}
.form-control:focus { border-color: #495057; box-shadow: 0 0 0 2px rgba(73,80,87,0.1); outline: none; }
.btn {
font-size: 11px; font-weight: 500; padding: 8px 14px; border-radius: 4px;
border: none; text-transform: uppercase; letter-spacing: 0.3px;
}
.btn-primary { background: #212529; color: #ffffff; }
.btn-primary:hover { background: #343a40; }
.btn-success { background: #495057; color: #ffffff; }
.btn-success:hover { background: #343a40; }
.btn-info { background: #6c757d; color: #ffffff; }
.btn-info:hover { background: #5a6268; }
.btn-block { display: block; width: 100%; margin-bottom: 8px; }
.plot-container {
background: #ffffff; border: 1px solid #dee2e6; border-radius: 6px;
padding: 16px; margin-bottom: 12px;
}
.plot-header {
font-size: 12px; font-weight: 600; color: #212529; margin-bottom: 12px;
padding-bottom: 8px; border-bottom: 1px solid #e9ecef;
}
.status-panel {
background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px;
padding: 12px; font-size: 11px; font-family: "SF Mono", "Consolas", monospace; color: #495057;
}
.info-panel {
background: #f8f9fa; border-left: 3px solid #6c757d; padding: 10px 12px;
font-size: 11px; color: #495057; margin-bottom: 12px; border-radius: 0 4px 4px 0;
}
.obs-table { max-height: 180px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 4px; }
.obs-row {
display: flex; padding: 8px 12px; border-bottom: 1px solid #f1f3f5;
cursor: pointer; font-size: 11px;
}
.obs-row:hover { background: #f8f9fa; }
.obs-row.selected { background: #495057; color: #ffffff; }
.obs-row .obs-name { flex: 2; font-weight: 500; }
.obs-row .obs-info { flex: 1; text-align: right; color: #6c757d; }
.obs-row.selected .obs-info { color: #adb5bd; }
.obs-header {
display: flex; padding: 8px 12px; background: #f8f9fa;
border-bottom: 1px solid #dee2e6; font-size: 10px; font-weight: 600;
text-transform: uppercase; letter-spacing: 0.5px; color: #495057;
position: sticky; top: 0; z-index: 1;
}
.obs-header .obs-name { flex: 2; }
.obs-header .obs-info { flex: 1; text-align: right; }
.data-scroll {
max-height: 350px; overflow-y: auto; border: 1px solid #dee2e6;
border-radius: 4px; font-family: "SF Mono", "Consolas", monospace; font-size: 11px;
}
.data-scroll table { width: 100%; border-collapse: collapse; }
.data-scroll th {
background: #f8f9fa; padding: 8px 12px; text-align: right; font-weight: 600;
font-size: 10px; text-transform: uppercase; letter-spacing: 0.5px;
border-bottom: 1px solid #dee2e6; position: sticky; top: 0; z-index: 1;
}
.data-scroll th:first-child { text-align: left; }
.data-scroll td { padding: 4px 12px; border-bottom: 1px solid #f1f3f5; text-align: right; }
.data-scroll td:first-child { text-align: left; color: #6c757d; }
.data-scroll tr:hover td { background: #fafafa; }
.shiny-file-input-progress { display: none; }
.form-group { margin-bottom: 12px; }
.selectize-input {
font-size: 12px !important; padding: 6px 10px !important;
min-height: 32px !important; border: 1px solid #ced4da !important; border-radius: 4px !important;
}
.selectize-dropdown { font-size: 12px; }
.selectize-dropdown-content .option { padding: 8px 10px; }
.plot-header .checkbox { display: inline-block; margin: 0; }
.plot-header .checkbox label { font-size: 11px; font-weight: 400; color: #6c757d; }
.plot-header .checkbox input { margin-right: 4px; }
'

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
useShinyjs(),
tags$head(tags$style(HTML(app_css))),

div(class = "title-panel",
    h2("BIOPAC ACQ Batch Converter"),
    div(class = "subtitle", "Convert AcqKnowledge files to text format")),

sidebarLayout(
  sidebarPanel(width = 3,
    wellPanel(
      h4("Data Input"),
      selectInput("input_mode", "Input Mode", choices = c("Select Files" = "files", "Select Folder" = "folder")),
      conditionalPanel("input.input_mode == 'files'",
        shinyFilesButton("in_files", "Browse Files", "Select .acq files", multiple = TRUE, class = "btn-info btn-block")),
      conditionalPanel("input.input_mode == 'folder'",
        shinyDirButton("in_folder", "Browse Folder", "Select folder", class = "btn-info btn-block")),
      uiOutput("input_status")
    ),
    wellPanel(
      h4("Output"),
      selectInput("format", "Format", choices = c("Tab-delimited (.txt)" = "txt", "CSV (.csv)" = "csv")),
      selectInput("meta", "Metadata", choices = c("Omit" = "omit", "Include" = "include")),
      shinyDirButton("out_dir", "Output Folder", "Select output folder", class = "btn-info btn-block"),
      uiOutput("out_dir_display")
    ),
    wellPanel(
      h4("Actions"),
      actionButton("convert", "Convert", class = "btn-success btn-block"),
      actionButton("save", "Save", class = "btn-primary btn-block"),
      actionButton("clear", "Clear", class = "btn-info btn-block")
    ),
    wellPanel(
      h4("Status"),
      div(class = "status-panel", uiOutput("status"))
    )
  ),
  
  mainPanel(width = 9,
    div(class = "plot-container",
        div(class = "plot-header", "Converted Files"),
        div(class = "obs-table", uiOutput("file_list"))),
    div(class = "plot-container",
        div(class = "plot-header", uiOutput("plot_header")),
        plotlyOutput("signal_plot", height = "350px")),
    div(class = "plot-container",
        div(class = "plot-header",
            span("Data Preview"),
            span(style = "float: right;", checkboxInput("show_all", "Show all rows", value = FALSE, width = "auto"))),
        div(class = "data-scroll", uiOutput("data_preview")))
  )
)
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

volumes <- c(Home = path.expand("~"), getVolumes()())
shinyFileChoose(input, "in_files", roots = volumes, filetypes = "acq")
shinyDirChoose(input, "in_folder", roots = volumes)
shinyDirChoose(input, "out_dir", roots = volumes)

rv <- reactiveValues(
  input_files = list(),
  input_dir = NULL,
  out_path = NULL,
  converted = list(),
  selected = NULL
)

# File selection
observeEvent(input$in_files, {
  if (is.integer(input$in_files)) return()
  files <- parseFilePaths(volumes, input$in_files)
  if (nrow(files) == 0) return()
  
  rv$input_files <- lapply(seq_len(nrow(files)), function(i) {
    list(name = files$name[i], path = as.character(files$datapath[i]))
  })
  names(rv$input_files) <- sapply(rv$input_files, function(x) x$name)
  rv$input_dir <- dirname(as.character(files$datapath[1]))
  rv$out_path <- rv$input_dir
  rv$converted <- list()
  rv$selected <- NULL
})

# Folder selection
observeEvent(input$in_folder, {
  if (is.integer(input$in_folder)) return()
  folder_path <- parseDirPath(volumes, input$in_folder)
  if (length(folder_path) == 0 || nchar(folder_path) == 0) return()
  
  acq_files <- list.files(folder_path, pattern = "\\.acq$", full.names = TRUE, ignore.case = TRUE)
  rv$input_files <- lapply(acq_files, function(f) list(name = basename(f), path = f))
  names(rv$input_files) <- sapply(rv$input_files, function(x) x$name)
  rv$input_dir <- folder_path
  rv$out_path <- folder_path
  rv$converted <- list()
  rv$selected <- NULL
})

# Output folder selection
observeEvent(input$out_dir, {
  if (is.integer(input$out_dir)) return()
  rv$out_path <- parseDirPath(volumes, input$out_dir)
})

# Input status
output$input_status <- renderUI({
  n <- length(rv$input_files)
  if (n == 0) return(NULL)
  div(class = "info-panel", style = "margin-top: 8px;",
      if (!is.null(rv$input_dir)) div(style = "word-break: break-all; font-weight: 500;", basename(rv$input_dir)),
      div(paste(n, ".acq file(s)")))
})

# Output folder display
output$out_dir_display <- renderUI({
  if (is.null(rv$out_path) || length(rv$out_path) == 0) return(NULL)
  div(class = "info-panel", style = "margin-top: 8px;",
      div(style = "word-break: break-all;", rv$out_path))
})

# Convert
observeEvent(input$convert, {
  req(length(rv$input_files) > 0)
  rv$converted <- list()
  
  for (fname in names(rv$input_files)) {
    f <- rv$input_files[[fname]]
    tryCatch({
      acq <- read_acq(f$path)
      df <- acq_to_df(acq)
      meta <- build_metadata(acq)
      rv$converted[[fname]] <- list(name = fname, acq = acq, df = df, meta = meta, error = NULL)
    }, error = function(e) {
      rv$converted[[fname]] <- list(name = fname, error = e$message)
    })
  }
  
  if (length(rv$converted) > 0) rv$selected <- names(rv$converted)[1]
})

# Clear
observeEvent(input$clear, {
  rv$converted <- list()
  rv$selected <- NULL
  rv$input_files <- list()
  rv$input_dir <- NULL
})

# File selection in list
observeEvent(input$select_file, { rv$selected <- input$select_file })

# File list
output$file_list <- renderUI({
  if (length(rv$converted) == 0) return(div(style = "padding: 20px; text-align: center; color: #6c757d;", "No files converted"))
  
  header <- div(class = "obs-header", span(class = "obs-name", "File"), span(class = "obs-info", "Channels"))
  rows <- lapply(names(rv$converted), function(fname) {
    item <- rv$converted[[fname]]
    is_sel <- identical(fname, rv$selected)
    n_ch <- if (is.null(item$error)) item$acq$num_channels else "Error"
    div(class = paste("obs-row", if (is_sel) "selected" else ""),
        onclick = sprintf("Shiny.setInputValue('select_file', '%s', {priority: 'event'})", fname),
        span(class = "obs-name", tools::file_path_sans_ext(fname)),
        span(class = "obs-info", n_ch))
  })
  tagList(header, rows)
})

# Plot header
output$plot_header <- renderUI({
  if (is.null(rv$selected)) return("Signal")
  item <- rv$converted[[rv$selected]]
  if (is.null(item) || !is.null(item$error)) return("Signal")
  HTML(paste0(tools::file_path_sans_ext(rv$selected),
              " &nbsp;<span style='color:#6c757d;font-weight:400;'>",
              item$acq$samples_per_sec, " Hz &middot; ",
              item$acq$num_channels, " ch &middot; ",
              format(item$acq$channels[[1]]$point_count, big.mark = ","), " samples</span>"))
})

# Signal plot
output$signal_plot <- renderPlotly({
  if (is.null(rv$selected) || !rv$selected %in% names(rv$converted)) {
    return(plot_ly(type = "scatter", mode = "lines") %>%
             layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                    annotations = list(list(text = "Select files and convert", x = 0.5, y = 0.5,
                                            xref = "paper", yref = "paper", showarrow = FALSE,
                                            font = list(size = 14, color = "#6c757d"))),
                    plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff") %>%
             config(displaylogo = FALSE))
  }
  
  item <- rv$converted[[rv$selected]]
  if (!is.null(item$error)) {
    return(plot_ly(type = "scatter", mode = "lines") %>%
             layout(annotations = list(list(text = paste("Error:", item$error), x = 0.5, y = 0.5,
                                            xref = "paper", yref = "paper", showarrow = FALSE,
                                            font = list(size = 12, color = "#dc3545"))),
                    plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff") %>%
             config(displaylogo = FALSE))
  }
  
  df <- item$df
  n <- nrow(df)
  if (n > 10000) df <- df[seq(1, n, length.out = 10000), ]
  
  colors <- c("#212529", "#0d6efd", "#198754", "#dc3545", "#ffc107",
              "#0dcaf0", "#6f42c1", "#fd7e14", "#20c997", "#6c757d")
  
  p <- plot_ly()
  ch_names <- names(df)[-1]
  for (i in seq_along(ch_names)) {
    ch <- ch_names[i]
    col <- colors[((i - 1) %% length(colors)) + 1]
    ch_units <- if (i <= length(item$acq$channels)) item$acq$channels[[i]]$units else ""
    legend_name <- if (nchar(ch_units) > 0) paste0(ch, " (", ch_units, ")") else ch
    p <- p %>% add_trace(x = df$time_sec, y = df[[ch]], type = "scatter", mode = "lines",
                         name = legend_name, line = list(color = col, width = 1))
  }
  
  p %>% layout(
    xaxis = list(title = list(text = "Time (s)", font = list(size = 11, color = "#495057")),
                 tickfont = list(size = 10, color = "#495057"), gridcolor = "#dee2e6", linecolor = "#adb5bd"),
    yaxis = list(title = list(text = "Amplitude", font = list(size = 11, color = "#495057")),
                 tickfont = list(size = 10, color = "#495057"), gridcolor = "#dee2e6", linecolor = "#adb5bd",
                 fixedrange = FALSE),
    legend = list(orientation = "h", x = 0, y = 1.1, font = list(size = 10)),
    plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff",
    margin = list(t = 40, b = 40, l = 60, r = 20), hovermode = "x unified"
  ) %>% config(displaylogo = FALSE, modeBarButtonsToRemove = list("lasso2d", "select2d"))
})

# Data preview
output$data_preview <- renderUI({
  if (is.null(rv$selected) || !rv$selected %in% names(rv$converted)) {
    return(div(style = "padding: 20px; text-align: center; color: #6c757d;", "No data"))
  }
  
  item <- rv$converted[[rv$selected]]
  if (!is.null(item$error)) return(div(style = "padding: 20px; text-align: center; color: #dc3545;", item$error))
  
  df <- item$df
  total_rows <- nrow(df)
  
  if (!isTRUE(input$show_all) && total_rows > 100) {
    df <- df[1:100, ]
    row_info <- paste0("Showing 1-100 of ", format(total_rows, big.mark = ","), " rows")
  } else {
    row_info <- paste0("Showing all ", format(total_rows, big.mark = ","), " rows")
  }
  
  header <- tags$tr(lapply(names(df), tags$th))
  rows <- lapply(seq_len(nrow(df)), function(i) {
    tags$tr(lapply(seq_along(df), function(j) {
      val <- df[i, j]
      tags$td(if (j == 1) sprintf("%.4f", val) else sprintf("%.6f", val))
    }))
  })
  
  tagList(
    div(style = "padding: 8px 12px; background: #f8f9fa; border-bottom: 1px solid #dee2e6; font-size: 11px; color: #6c757d;", row_info),
    tags$table(tags$thead(header), tags$tbody(rows))
  )
})

# Status
output$status <- renderUI({
  n_files <- length(rv$input_files)
  n_conv <- length(rv$converted)
  n_ok <- if (n_conv > 0) sum(sapply(rv$converted, function(x) is.null(x$error))) else 0
  
  lines <- c()
  if (n_files > 0) lines <- c(lines, paste("Input:", n_files, "file(s)"))
  if (n_conv > 0) lines <- c(lines, paste("Converted:", n_ok, "/", n_conv))
  if (!is.null(rv$selected)) lines <- c(lines, paste("Selected:", tools::file_path_sans_ext(rv$selected)))
  if (length(lines) == 0) lines <- "Ready"
  HTML(paste(lines, collapse = "<br>"))
})

# Save
observeEvent(input$save, {
  req(length(rv$converted) > 0, rv$out_path)
  
  ext <- if (input$format == "csv") ".csv" else ".txt"
  delim <- if (input$format == "csv") "," else "\t"
  saved <- 0
  
  for (fname in names(rv$converted)) {
    item <- rv$converted[[fname]]
    if (!is.null(item$error)) next
    
    out_path <- file.path(rv$out_path, paste0(tools::file_path_sans_ext(fname), ext))
    tryCatch({
      con <- file(out_path, "w")
      if (input$meta == "include") writeLines(item$meta, con)
      writeLines(paste(names(item$df), collapse = delim), con)
      write.table(item$df, con, sep = delim, row.names = FALSE, col.names = FALSE, quote = FALSE)
      close(con)
      saved <- saved + 1
    }, error = function(e) NULL)
  }
  
  showNotification(paste("Saved", saved, "file(s) to", rv$out_path), type = "message", duration = 4)
})
}

shinyApp(ui = ui, server = server)
