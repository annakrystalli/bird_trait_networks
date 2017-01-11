shinyApp(

ui = fluidPage(
  h1("Configure data_log"),
  h3("current data_log"),
  dataTableOutput('table'),
  br(),
  actionButton("BUTnew", "Edit"),
  bsModal("modalnew", "Edit dataset entry", "BUTnew", size = "large",
          selectInput("csv.file", label = "Select csv.file to edit", choices = csv.files, selected = csv.files[1], multiple = FALSE,
                      selectize = TRUE, width = NULL, size = NULL),
          uiOutput("sliders"),
          actionButton("BUTsave", "Save changes")),
  actionButton("BUTreset", "Reset"),
  actionButton("BUTsaveD", "Save data_log"),
  br(),
  textOutput("text"),
  
  
  h1("Configure vnames"),
  actionButton("BUTupdate", "Update vnames"),
  h3("current vnames"),
  br(),
  actionButton("BUTnew2", "Edit"),
  bsModal("modalnew2", "Edit dataset entry", "BUTnew2", size = "large",
          selectInput("csv.file2", label = "Select vnames column to edit", choices = vname.cols, selected = vname.cols[1], 
                      multiple = FALSE,
                      selectize = TRUE, width = NULL, size = NULL),
          p("selecting an fcode creates an input box for each column name available in the corresponding file"),
          textOutput("error"),
          tableOutput("dup.df"),
          uiOutput("sliders2"),
          uiOutput("sliders3"),
          actionButton("BUTsave2", "Save changes")),
  actionButton("BUTreset2", "Reset"),
  actionButton("BUTsaveD2", "Save vnames"),
  br(),
  textOutput("text2"),
  dataTableOutput('vnames')
),



server = function(input, output, session) {
  
  v <- reactiveValues(data_log = data_log, msg = "", vnames = vnames, msg2 = "", error = "",
                      format = "wide", input.values = NULL, input.vars = NULL, inputs = NULL,
                      dup.df = NULL, dat = NULL,
                      vars = NULL, l.vars = NULL, tmp_vnames = NULL)

  var.exists <- reactive({is.na(vnames[,input$csv.file2][vnames$code == "var"])})


  # ---- data-log ---- 
  output$table <- renderDataTable(data_log)
    
  output$sliders <- renderUI({
    lapply(1:numCols, function(i) {
      column <- columns[i]
      names(column) <- names(meta.var_columns)[meta.var_columns == column]
      if(column %in% meta.var_columns){
        files <- list.files(paste(input.folder, "pre/",
                                  names(meta.var_columns)[meta.var_columns == column],
                                  "/", sep = "")) %>%
          grep(pattern = "Icon\r", inv=T, value=T)
        
        files <- files[!files %in% data_log[,column]]
        current <- data_log[data_log$dcode == input$csv.file, column]
        
        empty <- NA
        files <- c(current, empty, files)
        
        selectInput(paste("col", i, sep = ""), label = column, choices = files, selected = current, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      }else{
        if(column == "format"){
          selectInput(paste("col", i, sep = ""), label = column, choices = c("wide", "long"), selected = NULL, multiple = FALSE,
                      selectize = TRUE, width = NULL, size = NULL)
        }
        current <- data_log[data_log$dcode == input$csv.file, column]
        textInput(paste0("col", i), label = column, value = current)
      }
    })
  })

  observeEvent(input$BUTsave, {
    if(is.null(input$BUTsave)){data_log}else{
    tmp_data_log <- v$data_log
    for(i in 1:numCols){
      cat("Showing",paste0("col", i), input[[paste0("col", i)]], "\n")
      column <- columns[i]
      tmp_data_log[tmp_data_log$dcode == input$csv.file, column] <- input[[paste0("col", i)]]}
    
    v$data_log <- tmp_data_log
    v$data_log
    v$msg <- "unsaved changes in browser data_log"
    }
  })
  
  output$table <- renderDataTable({
    v$data_log
  }) 
  
  observeEvent(input$BUTreset,{
    v$data_log <- read.csv(file = paste(input.folder, "metadata/data_log.csv", sep = ""),
                         stringsAsFactors = F, fileEncoding = fileEncoding, 
                         na.strings=na.strings, strip.white = T, 
                         blank.lines.skip = T)
    v$msg <- "browser data_log reset to file"
  })
  
  observeEvent(input$BUTsaveD,{
    v$msg <- "file data_log synced to browser"
    write.csv(v$data_log, file = paste(input.folder, "metadata/data_log.csv", sep = ""),
              row.names = F)
  })

  output$text <- renderText({ 
    v$msg
  })
  
  
  # ---- vnames ----
  
  output$vnames <- renderDataTable(vnames)
  
  output$sliders2 <- renderUI({
    data_log <- v$data_log
    vnames <- v$vnames
    
    if(is.na(data_log$format[data_log$dcode == paste0("D", substring(input$csv.file2,2))])){
    }else{
      v$format <- data_log$format[data_log$dcode == paste0("D", substring(input$csv.file2,2))]
    }
    
    fcode <- fcodes[names(fcodes) == gsub('[0-9]+', '',input$csv.file2)]
    path <- paste0(input.folder, "pre/", fcode, "/",
                   data_log[data_log$dcode == paste0("D", substring(input$csv.file2,2)), 
                            paste0(fcode, "_file.name")])
    dat <- read.csv(path, 
             header = T, strip.white = T, stringsAsFactors = F,
             na.strings = na.strings, blank.lines.skip = TRUE,
             fileEncoding = fileEncoding, check.names = F)
    
    dat <-  dat[,!sapply(dat, function(x)all(is.na(x)))]
    vars <- names(dat)
    vars <- setNames(make.names(vars), vars)
    v$vars <- vars
    v$dat <- dat
    
    numCols2 <- as.integer(length(vars))
    
    if(v$format == "long"){
      choices2 <- c(NA, master.vars, taxo.vars)
    }
    if(v$format == "wide"){
      choices2 <- c(NA, vnames$code)
    }

    lapply(1:numCols2, function(i) {
      var <- vars[i]
        current <- if(!var %in% vnames[,input$csv.file2]){NA}else{
        vnames$code[vnames[,input$csv.file2] == names(var)]}

        selectInput(paste("i_", var, sep = ""), label = names(var), choices = choices2, 
                    selected = current, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
    })
    
  
  })
  
  output$sliders3 <- renderUI({
    
    vnames <- v$vnames
    data_log <- v$data_log
    dat <- v$dat
    
    
    if(is.na(v$format)){return()}
    if(v$format == "wide"){
      v$l.vars <- NULL
      return()}
     var.name <- vnames[,input$csv.file2][vnames$code == "var"]
     if(is.na(var.name)){return()}else{

      choices3 <- c(NA, vnames$code[!vnames$code %in% c(master.vars, taxo.vars)])
      
        l.vars <- na.omit(unique(dat[, var.name]))
        l.vars <-  setNames(make.names(l.vars), l.vars)
        
        v$l.vars <- l.vars
        
        lapply(1:length(l.vars), FUN = function(i){
          l.var <- l.vars[i]
          current <- if(!names(l.var) %in% vnames[,input$csv.file2]){NA}else{
            vnames$code[vnames[,input$csv.file2] == names(l.var)]
            }
          selectInput(paste("i_", l.var, sep = ""), label = names(l.var), choices = choices3, 
                      selected = current, multiple = FALSE,
                      selectize = TRUE, width = NULL, size = NULL)
        })
      }
  })
  
  get_inputs <- reactive({
    current.inputs <- paste0("i_", c(v$vars, v$l.vars))
    all.inputs <- names(input)[grep("i_",names(input))]
    all.inputs[all.inputs %in% current.inputs]
  })
  
  get_input.vars <- reactive({
    all.vars <- c(v$vars, v$l.vars)
    inputs <- v$inputs
    all.vars[match(gsub("i_", "",inputs), all.vars)]
  })
  
  get_input.values <- reactive({
    inputs <- v$inputs
    input.vars <- v$input.vars
    input.values <- NULL
    
    for(j in 1:length(input.vars)){
      if(is.null(input[[inputs[j]]])){}else{
      if(any(is.na(input[[inputs[j]]]), input[[inputs[j]]] == "NA")){
        add <- NA 
        names(add) <- names(input.vars[j])
        input.values <- c(input.values, add)
      }else{
        add <- input[[inputs[j]]]
        names(add) <- names(input.vars[j])
        input.values <- c(input.values, add)
      }
        }
    }
    input.values
  })
  
  get_tmp_vnames <- reactive({
    inputs <- v$inputs
    input.vars <- v$input.vars
    input.values <- v$input.values
    tmp_vnames <- v$vnames
    
    for(j in 1:length(input.vars)){
      if(is.null(input.values[j])){
      }
      if(is.na(input.values[j])){
        tmp_vnames[, input$csv.file2][tmp_vnames[, input$csv.file2] == names(input.vars[j])] <- NA
      }else{
        tmp_vnames[tmp_vnames$code == input.values[j], 
                   input$csv.file2] <- names(input.vars[j])
      }
    }
    tmp_vnames[tmp_vnames == "NA"] <- NA
    tmp_vnames
  })
  
  any_dups <- reactive({any(duplicated(na.omit(v$input.values)))})
  
  get_dups <- reactive({
    if(any_dups()){
      dups <- na.omit(v$input.values[duplicated(v$input.values)])
      dup.df <- data_frame(code = v$input.values[v$input.values %in% dups],
                           file = names(v$input.values)[v$input.values %in% dups])
      dup.df
    }else{NULL}
  })
  
  update_vnames <- reactive({
    v$vnames <- v$tmp_vnames
    v$msg2 <- "unsaved changes in browser vnames"
  })
  
  observeEvent(input$BUTsave2, {
    if(is.null(input$BUTsave2)){}else{
      v$inputs <- get_inputs()
      v$input.vars <- get_input.vars()
      v$input.values <- get_input.values()
    
      v$tmp_vnames <- get_tmp_vnames()
      if(any_dups()){
        v$error <- "ERROR: duplicate file vnames attempted to be assigned to single var code"
        }else{
          v$error <- "" }
    }
  })
  
  output$dup.df <- renderTable({
  get_dups()
  })
  
  output$error <-  renderText(v$error)
   
  
  output$text2 <- renderText({ 
    v$msg2
  })
  observeEvent(input$BUTsave2,{
    if(!any_dups()){
    v$vnames <- v$tmp_vnames}
    v$inputs <- get_inputs()
    v$input.vars <- get_input.vars()
    v$input.values <- get_input.values()})
  
  output$vnames <- renderDataTable({
    v$vnames
  }) 
  
  observeEvent(input$BUTreset2,{
    v$vnames <- read.csv(file = paste(input.folder, "metadata/vnames.csv", sep = ""),
                         stringsAsFactors = F, fileEncoding = fileEncoding, 
                         na.strings=na.strings, strip.white = T, 
                         blank.lines.skip = T)
    v$msg2 <- "browser vnames reset to file"
  })
  
  observeEvent(input$BUTsaveD2,{
    v$msg2 <- "file vnames synced to browser"
    write.csv(v$vnames, file = paste(input.folder, "metadata/vnames.csv", sep = ""),
              row.names = F)
  })
  

}
)
