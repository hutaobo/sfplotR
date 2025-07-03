library(shiny)
library(pheatmap)

ui <- fluidPage(
  titlePanel("热图和层次聚类的 Shiny App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("upload", "上传 CSV 文件",
                accept = c(".csv")),
      selectInput("sortOrder", "排序方式：",
                  choices = c("升序" = "asc", "降序" = "desc"), selected = "asc"),
      numericInput("cellsize", "单元格尺寸 (像素)", value = 20, min = 5, step = 1),
      numericInput("fontsize", "字体大小", value = 10, min = 5, step = 1),
      numericInput("pdf_margin", "PDF额外边距(英寸)", value = 2, min = 0.5, step = 0.5),
      downloadButton("downloadPDF", "保存为 PDF"),
      hr(),
      # 左侧额外行选择区域，带滚动条
      div(style = "height:400px; overflow:auto; border:1px solid #ccc; padding:5px;",
          uiOutput("extraRowsUI"))
    ),
    mainPanel(
      # 将热图放置在右侧空白区域的正中间
      div(style = "display: flex; justify-content: center; align-items: center; min-height: 600px;",
          plotOutput("heatmapPlot")
      )
    )
  )
)

server <- function(input, output, session) {

  # 读取 CSV 文件，保持列名不被修改
  mat_data <- reactive({
    req(input$upload)
    df <- read.csv(input$upload$datapath, row.names = 1, check.names = FALSE)
    as.matrix(df)
  })

  # 选取行名与列名均存在的行
  overlap_rows <- reactive({
    intersect(rownames(mat_data()), colnames(mat_data()))
  })

  # 计算上传数据中剩余的行（行名不在列名中）
  extra_rows <- reactive({
    setdiff(rownames(mat_data()), overlap_rows())
  })

  # 动态生成额外行选择的 UI，并添加“全选”按钮
  output$extraRowsUI <- renderUI({
    req(extra_rows())
    sorted_choices <- if (input$sortOrder == "asc") {
      sort(extra_rows())
    } else {
      sort(extra_rows(), decreasing = TRUE)
    }
    tagList(
      checkboxGroupInput("extraRows", "选择额外行：", choices = sorted_choices),
      actionButton("selectAll", "全选")
    )
  })

  # 点击“全选”按钮时，将所有额外行选中
  observeEvent(input$selectAll, {
    req(extra_rows())
    sorted_choices <- if (input$sortOrder == "asc") {
      sort(extra_rows())
    } else {
      sort(extra_rows(), decreasing = TRUE)
    }
    updateCheckboxGroupInput(session, "extraRows", selected = sorted_choices)
  })

  # 构建用于热图的数据集：初始使用 overlap_rows，
  # 如果用户选择额外行，则将其加入数据中
  reactive_data <- reactive({
    req(mat_data())
    rows <- intersect(rownames(mat_data()), c(overlap_rows(), input$extraRows))
    mat_data()[rows, , drop = FALSE]
  })

  # 绘制热图，颜色为红-白-蓝渐变，边框为白色
  output$heatmapPlot <- renderPlot({
    req(reactive_data())
    pheatmap(reactive_data(),
             clustering_method = "average",
             fontsize = input$fontsize,
             cellwidth = input$cellsize,
             cellheight = input$cellsize,
             border_color = "white",
             color = colorRampPalette(c("red", "white", "blue"))(100),
             breaks = seq(0, 1, length.out = 101))
  },
  height = reactive({
    nrow(reactive_data()) * input$cellsize + 150  # 为 dendrogram 和标签预留空间
  }),
  width = reactive({
    ncol(reactive_data()) * input$cellsize + 150
  }))

  # 下载 PDF：采用相同的热图参数输出 PDF 文件
  output$downloadPDF <- downloadHandler(
    filename = function() {
      paste("heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      mat <- reactive_data()
      pdf_width <- ncol(mat) * input$cellsize / 72 + input$pdf_margin
      pdf_height <- nrow(mat) * input$cellsize / 72 + input$pdf_margin
      pdf(file, width = pdf_width, height = pdf_height)
      pheatmap(mat,
               clustering_method = "average",
               fontsize = input$fontsize,
               cellwidth = input$cellsize,
               cellheight = input$cellsize,
               border_color = "white",
               color = colorRampPalette(c("red", "white", "blue"))(100),
               breaks = seq(0, 1, length.out = 101))
      dev.off()
    }
  )
}

shinyApp(ui, server)
