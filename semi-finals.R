library(shiny)
library(EBImage)  # readImage()
library(bslib)
library(bsicons) 
library(rentrez)
library(msa)
library(ape)
library(phangorn)
library(Biostrings)


# ---------- GUI ----------
GUI = fluidPage(
  theme = bs_theme(bg = "#2A3133", fg = "#FFFFFF", primary = "#007bff", secondary = "#6c757d", warning = "orange"),
  div("EvoBuild", style = "font-size: 42px;"),
  fluidRow(
    column(12, align = "center",
           actionButton("btn_phylo", "Phylogenetic Tree", class = "btn-primary btn-lg"),
           actionButton("btn_col", "Colony Counter", class = "btn-warning btn-lg"),
           actionButton("btn_genome", "De novo Genome Assembly", class = "btn-success btn-lg")
    )
  ),
  
  br(), br(),
  uiOutput("dynamic_ui"),
)


# ---------- server ----------
server = function(input, output, session) {
  
  
  # ---------- Phylogenetic tree ----------
  
  # Function to fetch sequences from NCBI
  SFID <- function(Sequence_id) {
    Sequence = entrez_fetch(db = "nuccore", id = Sequence_id, rettype = "fasta", sep = "\n")
    Sequence <- unlist(strsplit(Sequence, "\n"))
    Sequence <- Sequence[!grepl("^>", Sequence)]
    paste(Sequence, collapse = "")
  }
  
  sacc = function(Scientific_Name) {
    query = paste(Scientific_Name, "COX1")  # สร้างคำค้น
    result = entrez_search(db = "nucleotide", term = query, retmax = 1)  # ค้นหาในฐานข้อมูล
    
    # ตรวจสอบว่ามี ID ถูกคืนค่าหรือไม่
    if (length(result$id) == 0) {
      return(NA)  # ถ้าไม่พบ ID ให้คืนค่า NA
    }
    
    id = result$id  # ดึง ID
    return(id)  # คืนค่า ID
  }
  
  # Phylogenetic Tree UI
  observeEvent(input$btn_phylo, {
    output$dynamic_ui <- renderUI({
      fluidPage(
        fluidRow(
          column(6,
                 bslib::card(
                   style = "background-color: #353A3A; color: white;",
                   card_header("Phylogenetic Tree", style = "font-size: 24px;"),
                   card_body(
                     p("Enter the SCIENTIFIC NAMES for analysis", style = "color: white;"),
                     p("Comma to separated and underscore for space"),
                     textInput("input_name", "Example : Pinctada margaritifera, Octopus vulgaris, Lymnaea stagnalis, Cornu aspersum,Venerupis philippinarum", width = "100%"),
                     actionButton("run_phylo", "Run Analysis", class = "btn-primary"),
                   )
                 )
          ),
          column(6,
                 bslib::card(
                   style = "background-color: #353A3A; color: white;",
                   card_header("Sequence ID"),
                   card_body(
                     tableOutput("species_table"),
                     p("Please make sure that sequence ID match with your interested species !!",style = "color: red;")
                   ))
          ),
          plotOutput("tree_plot")  # เพิ่ม Output สำหรับ plot
        )
      )
    })
  })
  
  # Phylogenetic Tree Generation
  observeEvent(input$run_phylo, {
    
    species_names = strsplit(input$input_name, ",")[[1]]
    ids = unlist(lapply(species_names, sacc))
    
    # ดึงลำดับ DNA เฉพาะ ID ที่ไม่เป็น NA
    seq = sapply(ids, SFID)
    valid_sequences = seq[nzchar(seq)]  # กรองเฉพาะลำดับที่ไม่ว่างเปล่า
    names(valid_sequences) <- species_names[nzchar(seq)]
    
    print(seq)
    #สร้าง data frame ทำตาราง
    species_data <- data.frame(Species = names(valid_sequences),Sequence_ID = ids[nzchar(seq)],Link = paste("https://www.ncbi.nlm.nih.gov/nuccore/",ids[nzchar(seq)],sep = ""),stringsAsFactors = FALSE)
    output$species_table <- renderTable({
      species_data  # แสดง Data Frame ในตาราง
    })
    
    SeqStringset = DNAStringSet(valid_sequences)
    print(SeqStringset)
    alignment = msa(SeqStringset)
    print(alignment)
    DNA = as.phyDat(as.matrix(alignment), type = "DNA")  # สร้าง phyDat จาก alignment
    
    # --- NJ ---
    dm = dist.hamming(DNA)
    tree = NJ(dm)
    
    # --- Plotting ---
    output$tree_plot <- renderPlot({
      plot(tree, main = "Phylogenetic Tree (Neighbor Joining)")
    })
    
    set.seed(123)
    NJtrees = bootstrap.phyDat(DNA, FUN = function(x) NJ(dist.hamming(x)), bs = 100)
    plotBS(tree, NJtrees, "phylogram")
    
  })
  
  
  # ---------- Colony counter ----------
  observeEvent(input$btn_col, {
    output$dynamic_ui = renderUI({
      fluidPage(
        fluidRow(
          column(4,  # Left column for controls
                 bslib::card(
                   h2("Colony Counter"),
                   p("Colony Counter Function used for colony plating"),
                   
                   # เพิ่มปุ่มอัปโหลดภาพ
                   fileInput("file_image", "Choose your plate pic",
                             accept = c("image/jpeg", "image/png")),
                   
                   # สร้าง slider input เพื่อปรับภาพ
                   sliderInput("offset", "Threshold Offset:", min = 0, max = 0.1, value = 0.025, step = 0.005),
                   sliderInput("makebrush", "Make Brush:", min = 1, max = 15, value = 9, step = 2),
                   
                   # ปุ่มนับ colony
                   actionButton("count_colonies", "Count Colonies", class = "btn-warning btn-lg", w = "300", h = "300"),
                   style = "background-color: #353A3A; color: white;",
                   p("Please check your plate pic to make sure it is clear", style = "color: red;")
                 )
          ),
          column(6,  # Right column for image output
                 imageOutput("processed_nmaskt")  # เปลี่ยนขนาดที่นี่
          )
        )
      )
    })
  })
  
  
  # ฟังก์ชันสำหรับนับ colonies
  observeEvent(input$count_colonies, {
    req(input$file_image)  # ตรวจสอบว่ามีไฟล์ภาพแล้ว
    
    # เรียกใช้ฟังก์ชันนับโคโลนี
    colony_count = count_colonies_function(input$file_image$datapath,
                                           input$offset,
                                           input$makebrush)
    
    # แสดงผลลัพธ์
    showModal(modalDialog(
      title = "ผลลัพธ์การนับ colony",
      paste("จำนวน colony ที่นับได้:", colony_count),  # แสดงจำนวนโคโลนีที่นับได้
      easyClose = TRUE
    ))
  })
  
  # ฟังก์ชันเพื่อแสดงภาพที่ปรับแต่งแล้ว (nmaskf)
  output$processed_nmaskf = renderImage({
    req(input$file_image)
    
    # ปรับแต่งภาพ
    plate = readImage(input$file_image$datapath)
    nmaskf = thresh(plate, w = 15, h = 15, offset = input$offset)
    
    # บันทึกภาพ
    img_path = tempfile(fileext = ".png")
    writeImage(nmaskf, img_path)
    
    # ส่งคืนข้อมูลภาพ
    list(src = img_path, contentType = 'image/png', alt = "Processed nmaskf", weight = "400", height = "400")
  }, deleteFile = TRUE)
  
  # ฟังก์ชันเพื่อแสดงภาพที่ปรับแต่งแล้ว (nmaskt)
  output$processed_nmaskt = renderImage({
    req(input$file_image)
    
    # ปรับแต่งภาพ
    plate = readImage(input$file_image$datapath)
    nmaskf = thresh(plate, w = 15, h = 15, offset = input$offset)
    nmaskt = fillHull(opening(nmaskf, makeBrush(input$makebrush, shape = 'diamond')))
    
    # บันทึกภาพ
    img_path = tempfile(fileext = ".png")
    writeImage(nmaskt, img_path)
    
    # ส่งคืนข้อมูลภาพ
    list(src = img_path, contentType = 'image/png', alt = "Processed nmaskt", weight = "400", height = "400")
  }, deleteFile = TRUE)
  
  
  # ฟังก์ชันนับโคโลนี
  count_colonies_function = function(image_path, nmaskf, nmaskt) {
    # ดึงภาพมาเก็บไว้ในตัวแปร plate
    plate = readImage(image_path)
    
    # ปรับแต่งภาพ
    nmaskf = thresh(plate, w = 15, h = 15, offset = input$offset)  # กำหนด threshold
    nmaskt = fillHull(opening(nmaskf, makeBrush(input$makebrush, shape = 'diamond')))
    
    # นับ colony
    colony_labels = bwlabel(nmaskt)
    colony_count = max(colony_labels)
    
    return(colony_count)  # คืนค่าจำนวนโคโลนีที่นับได้
  }
  
  # ---------- Genome assembly ----------
  observeEvent(input$btn_genome, {
    output$dynamic_ui = renderUI({
      fluidPage(
        h2("ฟังก์ชัน: De novo Genome Assembly"),
        p("ใส่โค้ดการทำงานสำหรับสร้าง genome assembly ได้ที่นี่")
      )
    })
  })
  
}


# ---------- เรียกใช้แอป Shiny ----------
shinyApp(ui = GUI, server = server)
