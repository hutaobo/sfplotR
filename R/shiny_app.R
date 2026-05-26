ui <- shiny::fluidPage(
  shiny::titlePanel("Cell-GPS heatmap app"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::fileInput("upload", "Upload CSV file", accept = c(".csv")),
      shiny::selectInput(
        "sortOrder",
        "Sort order",
        choices = c("Ascending" = "asc", "Descending" = "desc"),
        selected = "asc"
      ),
      shiny::numericInput("cellsize", "Cell size (pixels)", value = 20, min = 5, step = 1),
      shiny::numericInput("fontsize", "Font size", value = 10, min = 5, step = 1),
      shiny::numericInput("pdf_margin", "PDF extra margin (inches)", value = 2, min = 0.5, step = 0.5),
      shiny::downloadButton("downloadPDF", "Save as PDF"),
      shiny::hr(),
      shiny::div(
        style = "height:400px; overflow:auto; border:1px solid #ccc; padding:5px;",
        shiny::uiOutput("extraRowsUI")
      )
    ),
    shiny::mainPanel(
      shiny::div(
        style = "display: flex; justify-content: center; align-items: center; min-height: 600px;",
        shiny::plotOutput("heatmapPlot")
      )
    )
  )
)

server <- function(input, output, session) {
  mat_data <- shiny::reactive({
    shiny::req(input$upload)
    df <- utils::read.csv(input$upload$datapath, row.names = 1, check.names = FALSE)
    as.matrix(df)
  })

  overlap_rows <- shiny::reactive({
    intersect(rownames(mat_data()), colnames(mat_data()))
  })

  extra_rows <- shiny::reactive({
    setdiff(rownames(mat_data()), overlap_rows())
  })

  output$extraRowsUI <- shiny::renderUI({
    shiny::req(extra_rows())
    sorted_choices <- if (input$sortOrder == "asc") {
      sort(extra_rows())
    } else {
      sort(extra_rows(), decreasing = TRUE)
    }
    shiny::tagList(
      shiny::checkboxGroupInput("extraRows", "Additional rows", choices = sorted_choices),
      shiny::actionButton("selectAll", "Select all")
    )
  })

  shiny::observeEvent(input$selectAll, {
    shiny::req(extra_rows())
    sorted_choices <- if (input$sortOrder == "asc") {
      sort(extra_rows())
    } else {
      sort(extra_rows(), decreasing = TRUE)
    }
    shiny::updateCheckboxGroupInput(session, "extraRows", selected = sorted_choices)
  })

  reactive_data <- shiny::reactive({
    shiny::req(mat_data())
    rows <- intersect(rownames(mat_data()), c(overlap_rows(), input$extraRows))
    mat_data()[rows, , drop = FALSE]
  })

  output$heatmapPlot <- shiny::renderPlot(
    {
      shiny::req(reactive_data())
      pheatmap::pheatmap(
        reactive_data(),
        clustering_method = "average",
        fontsize = input$fontsize,
        cellwidth = input$cellsize,
        cellheight = input$cellsize,
        border_color = "white",
        color = grDevices::colorRampPalette(c("red", "white", "blue"))(100),
        breaks = seq(0, 1, length.out = 101)
      )
    },
    height = shiny::reactive({
      nrow(reactive_data()) * input$cellsize + 150
    }),
    width = shiny::reactive({
      ncol(reactive_data()) * input$cellsize + 150
    })
  )

  output$downloadPDF <- shiny::downloadHandler(
    filename = function() {
      paste("heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      mat <- reactive_data()
      pdf_width <- ncol(mat) * input$cellsize / 72 + input$pdf_margin
      pdf_height <- nrow(mat) * input$cellsize / 72 + input$pdf_margin
      grDevices::pdf(file, width = pdf_width, height = pdf_height)
      pheatmap::pheatmap(
        mat,
        clustering_method = "average",
        fontsize = input$fontsize,
        cellwidth = input$cellsize,
        cellheight = input$cellsize,
        border_color = "white",
        color = grDevices::colorRampPalette(c("red", "white", "blue"))(100),
        breaks = seq(0, 1, length.out = 101)
      )
      grDevices::dev.off()
    }
  )
}

run_cellgpsr_app <- function(...) {
  shiny::shinyApp(ui = ui, server = server, ...)
}
