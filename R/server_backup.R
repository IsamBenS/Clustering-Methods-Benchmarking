library("shiny")
library("shinyjs")
library(flowCore)
library(gplots)

source("../FilePreHandling.R")
source("../StatisticalAnalysis.R")
source("../TimeAnalysis.R")
source("../EfficiencyAnalysis.R")

server <- function(input, output, session)
{
    useShinyjs()
    #======================================================================================================================
    #======================REACTIVE VALUES=================================================================================
    #======================================================================================================================

    global.values <- reactiveValues(
        fcs.files.fg.ref = NULL,
        fcs.files.fg.proj.1 = NULL,
        fcs.files.fg.proj.2 = NULL,
        
        fcs.files.mfi.proj.1 = NULL,
        fcs.files.mfi.proj.2 = NULL,
        
        current.panel = "FG_panel",
        plot = TRUE
    )

    analysis.variables <- reactiveValues(
        clusters.relations.list.1 = NULL,
        clusters.relations.id.list.1 = NULL,
        theta.coef.list.1 = NULL,
        F.beta.coef.list.1 = NULL,
        G.coef.list.1 = NULL,
        FG.accuracy.coef.list.1 = NULL,
        raw.FG.accuracy.coef.list.1 = NULL,
        
        clusters.relations.list.2 = NULL,
        clusters.relations.id.list.2 = NULL,
        theta.coef.list.2 = NULL,
        F.beta.coef.list.2 = NULL,
        G.coef.list.2 = NULL,
        FG.accuracy.coef.list.2 = NULL,
        raw.FG.accuracy.coef.list.2 = NULL
    )
    
    mfi.variables <- reactiveValues(
        clusters.resemblance.list = NULL,
        s.id.list = NULL,
        sent.list = NULL,
        smfi.list = NULL,
        s.list = NULL
    )
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #=======================TABS CHANGING==================================================================================
    #======================================================================================================================
    
    observe(
    {
        if(!is.null(input$tabs))
        {
            if(input$tabs != global.values$current.panel)
            {
                removeUI("#left_menu_items")
                insertUI("#left_menu",
                         "beforeEnd",
                         div(id="left_menu_items")
                )
                
                removeUI("#plots_section_fg")
                insertUI("#FG_description",
                         "afterEnd",
                         div(id="plots_section_fg")
                )
                
                removeUI("#plots_section_mfi")
                insertUI("#MFI_description",
                         "afterEnd",
                         div(id="plots_section_mfi")
                )
                
                global.values$current.panel = input$tabs
                global.values$fcs.files.fg.ref <- NULL
                global.values$fcs.files.fg.proj.1 <- NULL
                global.values$fcs.files.fg.proj.2 <- NULL
                global.values$fcs.files.mfi.proj.1 <- NULL
                global.values$fcs.files.mfi.proj.2 <- NULL
                global.values$plot <- FALSE
                
                analysis.variables$clusters.relations.list.1 = NULL
                analysis.variables$clusters.relations.id.list.1 = NULL
                analysis.variables$theta.coef.list.1 = NULL
                analysis.variables$F.beta.coef.list.1 = NULL
                analysis.variables$G.coef.list.1 <- NULL
                analysis.variables$FG.accuracy.coef.list.1 <- NULL
                analysis.variables$raw.FG.accuracy.coef.list.1 <- NULL
                    
                analysis.variables$clusters.relations.list.2 <- NULL
                analysis.variables$clusters.relations.id.list.2 <- NULL
                analysis.variables$theta.coef.list.2 <- NULL
                analysis.variables$F.beta.coef.list.2 <- NULL
                analysis.variables$G.coef.list.2 <- NULL
                analysis.variables$FG.accuracy.coef.list.2 <- NULL
                analysis.variables$raw.FG.accuracy.coef.list.2 = NULL
                
                mfi.variables$clusters.resemblance.list <- NULL
                mfi.variables$s.id.list <- NULL
                mfi.variables$sent.list <- NULL
                mfi.variables$smfi.list <- NULL
                mfi.variables$s.list = NULL
            }
        }
    })
    
    
    
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #=======================COMPARE FILES==================================================================================
    #======================================================================================================================
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.ref))
        {
            ref.col <- 1:ncol(global.values$fcs.files.fg.ref[[1]]@exprs)
            names(ref.col) <- colnames(global.values$fcs.files.fg.ref[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_fg_r", "Clusters column", choices=ref.col, selected=length(ref.col))
            
            lapply(global.values$fcs.files.fg.ref, function(f)
            {
                f.info <- FPH.get.file.information(f,ref.col[[length(ref.col)]])
                cl.info <- "<p><BLOCKQUOTE>"
                size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                lapply(1:length(f.info),function(i)
                {
                    cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                })
                cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                insertUI("#ref_files_list",
                         "beforeEnd",
                         ui = div(id=paste0("ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                  h5(paste0("REFERENCE: ", basename(f@description$FILENAME)), " ", size, " events"),
                                  div(HTML(cl.info),class="menu_file_info"),
                                  class="menu_file"
                                )
                )
            })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_fg_r", "Clusters column", choices=list(" "=0))
        }
    })
    
    observeEvent(input$ref_selection_fg,
    {
         shinyjs::disable("ref_selection_fg")
         shinyjs::disable("compare")
         
         removeUI("#ref_files_list")
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "FlowFrames"
         m[1,2] = "*.csv;*.fcs"
         temp.files <- choose.files(filters = m,multi = FALSE)
         global.values$fcs.files.fg.ref <- NULL
         if(length(temp.files) != 0)
         {
             insertUI("#left_menu_items",
                      "afterBegin",
                      div(id="ref_files_list",
                          class="menu_set"))
             global.values$fcs.files.fg.ref <- NULL
             progress <- Progress$new()
             progress$set(message="LOADING FILE", value=0)
             
             lapply(temp.files, function(f)
             {
                 l <- length(f)
                 x <- NULL
                 if(grepl("csv",f))
                 {
                     x <- as.matrix(read.csv(f))
                     x <- flowFrame(x)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.ref <<- c(global.values$fcs.files.fg.ref, x)
                 }
                 else
                 {
                     x <- read.FCS(f,emptyValue = FALSE)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.ref <<- c(global.values$fcs.files.fg.ref, x)
                 }
                 names(global.values$fcs.files.fg.ref)[[length(global.values$fcs.files.fg.ref)]] <<- substr(f,1,nchar(f)-4)
             })
             progress$inc(1, detail="FILE LOADED")
             
             progress$close()
         }
         
         shinyjs::enable("ref_selection_fg")
         shinyjs::enable("compare")
    })
    
    observeEvent(input$ref_update_fg,
    {
        if(!is.null(global.values$fcs.files.fg.ref) && !is.null(input$clust_col_selection_fg_r) && input$clust_col_selection_fg_r != "")
        {
            lapply(global.values$fcs.files.fg.ref, function(f)
            {
                f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_fg_r))
                if(length(f.info)>0)
                {
                    cl.info <- "<p><BLOCKQUOTE>"
                    removeUI(paste0("#ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
                    size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                    lapply(1:length(f.info),function(i)
                    {
                        cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                    })
                    cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                    insertUI("#ref_files_list",
                             "beforeEnd",
                             ui = div(id=paste0("ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                      h5(paste0("REFERENCE: ", basename(f@description$FILENAME)), " ", size, " events"),
                                      div(HTML(cl.info),class="menu_file_info"),
                                      class="menu_file"
                             )
                    )
                }
            })
        }
    })
    
    
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.proj.1))
        {
            p1.col <- 1:ncol(global.values$fcs.files.fg.proj.1[[1]]@exprs)
            names(p1.col) <- colnames(global.values$fcs.files.fg.proj.1[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_fg_1", "Clusters column", choices=p1.col, selected=length(p1.col))
            lapply(global.values$fcs.files.fg.proj.1, function(f)
            {
                f.info <- FPH.get.file.information(f,p1.col[[length(p1.col)]])
                cl.info <- "<p><BLOCKQUOTE>"
                size <- sum(unlist(sapply(f.info, function(f.i){return(f.i[[1]])})))
                lapply(1:length(f.info),function(i)
                {
                    cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                })
                cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                insertUI("#p1_files_list",
                         "beforeEnd",
                         ui = div(id=paste0("p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                  h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
                                  div(HTML(cl.info),class="menu_file_info"),
                                  class="menu_file"
                         )
                )
            })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_fg_1", "Clusters column", choices=list(" "=0))
        }
    })
    
    observeEvent(input$proj1_selection_fg,
    {
        shinyjs::disable("proj1_selection_fg")
        shinyjs::disable("compare")
        
        removeUI("#p1_files_list")
        
        
        m <- matrix(nrow=1,ncol=2)
        m[1,1] = "FlowFrames"
        m[1,2] = "*.csv;*.fcs"
        temp.files <- choose.files(filters = m,multi = F)
        global.values$fcs.files.fg.proj.1 <- NULL
        if(length(temp.files) != 0)
        {
            insertUI("#left_menu_items",
                     "beforeEnd",
                     div(id="p1_files_list",
                         class="menu_set"))
             progress <- Progress$new()
             progress$set(message="LOADING FILE")
             
             lapply(temp.files, function(f)
             {
                 l <- length(f)
                 if(grepl("csv",f))
                 {
                     x <- as.matrix(read.csv(f))
                     x <- flowFrame(x)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.proj.1 <<- c(global.values$fcs.files.fg.proj.1, x)
                 }
                 else
                 {
                     x <- read.FCS(f,emptyValue = FALSE)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.proj.1 <<- c(global.values$fcs.files.fg.proj.1, x)
                 }
                 names(global.values$fcs.files.fg.proj.1)[[length(global.values$fcs.files.fg.proj.1)]] <<- substr(f,1,nchar(f)-4)
             })
            progress$inc(1, detail="FILE LOADED")
            progress$close()
        }
         
         shinyjs::enable("proj1_selection_fg")
         shinyjs::enable("compare")
         
    })
    
    observeEvent(input$proj1_update_fg,
    {
        if(!is.null(global.values$fcs.files.fg.proj.1) && !is.null(input$clust_col_selection_fg_1) && input$clust_col_selection_fg_1 != "")
        {
            lapply(global.values$fcs.files.fg.proj.1, function(f)
            {
                f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_fg_1))
                if(length(f.info)>0)
                {
                    cl.info <- "<p><BLOCKQUOTE>"
                    removeUI(paste0("#p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
                    size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                    lapply(1:length(f.info),function(i)
                    {
                        cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                    })
                    cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                    insertUI("#p1_files_list",
                             "beforeEnd",
                             ui = div(id=paste0("p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                      h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
                                      div(HTML(cl.info),class="menu_file_info"),
                                      class="menu_file"
                             )
                    )
                }
            })
        }
    })
    
    
    
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.proj.2))
        {
            p2.col <- 1:ncol(global.values$fcs.files.fg.proj.2[[1]]@exprs)
            names(p2.col) <- colnames(global.values$fcs.files.fg.proj.2[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_fg_2", "Clusters column", choices=p2.col, selected=length(p2.col))
            
            lapply(global.values$fcs.files.fg.proj.2, function(f)
            {
                f.info <- FPH.get.file.information(f,p2.col[[length(p2.col)]])
                cl.info <- "<p><BLOCKQUOTE>"
                size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                lapply(1:length(f.info),function(i)
                {
                    cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                })
                cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                insertUI("#p2_files_list",
                         "beforeEnd",
                         ui = div(id=paste0("p2_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                  h5(paste0("FILE 2: ", basename(f@description$FILENAME)), " ", size, " events"),
                                  div(HTML(cl.info),class="menu_file_info"),
                                  class="menu_file"
                         )
                )
            })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_fg_2", "Clusters column", choices=list(" "=0))
        }
    })
    
    observeEvent(input$proj2_selection_fg,
    {
        shinyjs::disable("proj2_selection_fg")
        shinyjs::disable("compare")
        
        removeUI("#p2_files_list")
        
        
        m <- matrix(nrow=1,ncol=2)
        m[1,1] = "FlowFrames"
        m[1,2] = "*.csv;*.fcs"
        temp.files <- choose.files(filters = m,multi = FALSE)
        global.values$fcs.files.fg.proj.2 <- NULL
        if(length(temp.files) > 0)
        {
            insertUI("#left_menu_items",
                     "beforeEnd",
                     div(id="p2_files_list",
                         class="menu_set"))
            progress <- Progress$new()
            progress$set(message="LOADING FILE")
            
            lapply(temp.files, function(f)
            {
                 if(grepl("csv",f))
                 {
                     x <- as.matrix(read.csv(f))
                     x <- flowFrame(x)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.proj.2 <<- c(global.values$fcs.files.fg.proj.2, x)
                 }
                 else
                 {
                     x <- read.FCS(f,emptyValue = FALSE)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.fg.proj.2 <<- c(global.values$fcs.files.fg.proj.2, x)
                 }
                names(global.values$fcs.files.fg.proj.2)[[length(global.values$fcs.files.fg.proj.2)]] <<- substr(f,1,nchar(f)-4)
            })
            progress$inc(1, detail="FILE LOADED")
            progress$close()
        }
        
        shinyjs::enable("proj2_selection_fg")
        shinyjs::enable("compare")
    })
    
    observeEvent(input$proj2_update_fg,
    {
         if(!is.null(global.values$fcs.files.fg.proj.2) && !is.null(input$clust_col_selection_fg_2) && input$clust_col_selection_fg_2 != "")
         {
             lapply(global.values$fcs.files.fg.proj.2, function(f)
             {
                 f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_fg_2))
                 if(length(f.info)>0)
                 {
                     cl.info <- "<p><BLOCKQUOTE>"
                     removeUI(paste0("#p2_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
                     size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                     lapply(1:length(f.info),function(i)
                     {
                         cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                     })
                     cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                     insertUI("#p2_files_list",
                              "beforeEnd",
                              ui = div(id=paste0("p2_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                       h5(paste0("FILE 2: ", basename(f@description$FILENAME)), " ", size, " events"),
                                       div(HTML(cl.info),class="menu_file_info"),
                                       class="menu_file"
                              )
                     )
                 }
             })
         }
    })
    
    
    
    observeEvent(input$compare,
    {
        analysis.variables$clusters.relations.list.2 = NULL
        analysis.variables$clusters.relations.id.list.2 = NULL
        analysis.variables$theta.coef.list.2 = NULL
        analysis.variables$F.beta.coef.list.2 = NULL
        analysis.variables$G.coef.list.2 = NULL
        analysis.variables$FG.accuracy.coef.list.2 = NULL
        analysis.variables$raw.FG.accuracy.coef.list.2 = NULL
        
        analysis.variables$clusters.relations.list.1 = NULL
        analysis.variables$clusters.relations.id.list.1 = NULL
        analysis.variables$theta.coef.list.1 = NULL
        analysis.variables$F.beta.coef.list.1 = NULL
        analysis.variables$G.coef.list.1 = NULL
        analysis.variables$FG.accuracy.coef.list.1 = NULL
        analysis.variables$raw.FG.accuracy.coef.list.1 = NULL
        
        global.values$plot <- TRUE
        
        
        shinyjs::disable("compare")
        
        
        #UPDATE PLOT SECTION===============================================================
        removeUI("#plots_section_fg")
        insertUI("#FG_description",
                 "afterEnd",
                 div(id="plots_section_fg")
        )
        
        fcs.files <- NULL
        used.proj <- 0
        if(!is.null(global.values$fcs.files.fg.ref))
        {
            fcs.files <- global.values$fcs.files.fg.ref
            used.proj <- 3
        }
        else
        {
            if(length(global.values$fcs.files.fg.proj.2) <= length(global.values$fcs.files.fg.proj.1) && !is.null(global.values$fcs.files.fg.proj.2))
            {
                fcs.files <- global.values$fcs.files.fg.proj.2
                used.proj <- 2
            }
            else
            {
                if(!is.null(global.values$fcs.files.fg.proj.1))
                {
                    fcs.files <- global.values$fcs.files.fg.proj.1
                    used.proj <- 1
                }
            }
        }
        
        if(used.proj>0)
        {
            lapply(c(1:length(names(fcs.files))), function(f.id)
            {
                
                f.name <- names(fcs.files)[[f.id]]
                fcs.ref <- fcs.files[[f.id]]
                fcs.ref.info <- 0
                
                
                insertUI(selector = "#plots_section_fg",
                         where = "beforeEnd",
                         ui = div(id=paste0("pl_",f.id,"_set"),
                                  h4(basename(f.name)),
                                  actionButton(paste0("pl_",f.id,"_set_h"), "Hide"),
                                  actionButton(paste0("pl_",f.id,"_set_s"), "Show"),
                                  class="pl_set"
                              )
                )
                
                observeEvent(input[[paste0("pl_",f.id,"_set_h")]],
                {
                    shinyjs::hide(paste0("pl1_",f.id))
                    shinyjs::hide(paste0("pl2_",f.id))
                })
                
                observeEvent(input[[paste0("pl_",f.id,"_set_s")]],
                {
                    shinyjs::show(paste0("pl1_",f.id))
                    shinyjs::show(paste0("pl2_",f.id))
                })
                
                rangeX <- c(max(min(0,fcs.ref@exprs),-100),min(max(10,fcs.ref@exprs),10))
                rangeY <- rangeX
                #VAL.1=======================================================================================================
                    if( (used.proj < 3) || (used.proj==3 && !is.null(global.values$fcs.files.fg.proj.1)) )
                    {
                        progress <- Progress$new()
                        progress$set(message="COMPARING FILES")
                        col.cl.1 <- 0
                        col.cl.2 <- 0
                        if(used.proj == 1)
                        {
                            col.cl.1 <- as.integer(input$clust_col_selection_fg_1)
                            col.cl.2 <- as.integer(input$clust_col_selection_fg_2)
                        }
                        else if(used.proj == 2)
                        {
                            col.cl.2 <- as.integer(input$clust_col_selection_fg_1)
                            col.cl.1 <- as.integer(input$clust_col_selection_fg_2)
                        }
                        else
                        {
                            col.cl.1 <- as.integer(input$clust_col_selection_fg_r)
                            col.cl.2 <- as.integer(input$clust_col_selection_fg_1)
                        }
                        
                        fcs.test.1 <- 0
                        fcs.test.1.info <- 0
                        val.1 <- 0
                        
                        if(used.proj == 1)
                        {
                            fcs.ref.info <- FPH.get.file.information(fcs.ref, as.integer(input$clust_col_selection_fg_1))
                            if(!is.null(global.values$fcs.files.fg.proj.2[[f.id]]))
                            {
                                fcs.test.1 <- global.values$fcs.files.fg.proj.2[[f.id]]
                                fcs.test.1 <- FPH.create.artificials.clusters(fcs.ref, fcs.test.1, col.cl.1, col.cl.2)
                                fcs.test.1.info <- FPH.get.file.information(fcs.test.1, col.cl.2)
                                
                                val.1 <- FPH.get.clusters.associations(fcs.ref.info, fcs.test.1.info)
                                analysis.variables$clusters.relations.list.1[[f.name]] <- val.1[[3]]
                                analysis.variables$clusters.relations.id.list.1[[f.name]] <- val.1[[2]]
                                analysis.variables$theta.coef.list.1[[f.name]] <- val.1[[1]]
                            }
                        }
                        else if(used.proj == 2)
                        {
                            fcs.ref.info <- FPH.get.file.information(fcs.ref, as.integer(input$clust_col_selection_fg_2))
                            if(!is.null(global.values$fcs.files.fg.proj.1[[f.id]]))
                            {
                                fcs.test.1 <- global.values$fcs.files.fg.proj.1[[f.id]]
                                fcs.test.1 <- FPH.create.artificials.clusters(fcs.ref, fcs.test.1, col.cl.1, col.cl.2)
                                fcs.test.1.info <- FPH.get.file.information(fcs.test.1, col.cl.2)
                                
                                val.1 <- FPH.get.clusters.associations(fcs.ref.info, fcs.test.1.info)
                                analysis.variables$clusters.relations.list.1[[f.name]] <- val.1[[3]]
                                analysis.variables$clusters.relations.id.list.1[[f.name]] <- val.1[[2]]
                                analysis.variables$theta.coef.list.1[[f.name]] <- val.1[[1]]
                            }
                            
                        }
                        else
                        {
                            fcs.ref.info <- FPH.get.file.information(fcs.ref, as.integer(input$clust_col_selection_fg_r))
                            if(!is.null(global.values$fcs.files.fg.proj.1[[f.id]]))
                            {
                                fcs.test.1 <- global.values$fcs.files.fg.proj.1[[f.id]]
                                fcs.test.1 <- FPH.create.artificials.clusters(fcs.ref, fcs.test.1, col.cl.1, col.cl.2)
                                fcs.test.1.info <- FPH.get.file.information(fcs.test.1, col.cl.2)
                                
                                val.1 <- FPH.get.clusters.associations(fcs.ref.info, fcs.test.1.info)
                                analysis.variables$clusters.relations.list.1[[f.name]] <- val.1[[3]]
                                analysis.variables$clusters.relations.id.list.1[[f.name]] <- val.1[[2]]
                                analysis.variables$theta.coef.list.1[[f.name]] <- val.1[[1]]
                            }
                        }
                        
                        cl.to.plot.test.1 <- 1:length(fcs.test.1.info)
                        names(cl.to.plot.test.1) <- sapply(1:length(cl.to.plot.test.1), function(i){return(paste0("C",i))})
                        
                        #========================================================================================================
                        prec.coef.1 <- 0
                        prec.coef.2 <- 0
                        rec.coef.1 <- 0
                        rec.coef.2 <- 0
                        
                        if(used.proj < 3)
                        {
                            prec.coef.1 <- mean(sapply(1:(length(val.1[[2]])-1), function(c)
                            {
                                return(val.1[[3]][[c+1]][[1]])
                            }))
                            rec.coef.1 <- mean(sapply(1:(length(val.1[[2]])-1), function(c)
                            {
                                return(val.1[[3]][[c+1]][[2]])
                            }))
                        }
                        else
                        {
                            if(!is.null(global.values$fcs.files.fg.proj.1))
                            {
                                prec.coef.1 <- mean(sapply(1:(length(val.1[[2]])-1), function(c)
                                {
                                    return(val.1[[3]][[c+1]][[1]])
                                }))
                                rec.coef.1 <- mean(sapply(1:(length(val.1[[2]])-1), function(c)
                                {
                                    return(val.1[[3]][[c+1]][[2]])
                                }))
                            }
                        }
                        
                        analysis.variables$F.beta.coef.list.1[[f.name]] <<- calculate.F.beta.coef(prec.coef.1, rec.coef.1)
                        analysis.variables$G.coef.list.1[[f.name]] <<- calculate.G.coef(prec.coef.1, rec.coef.1)
                        
                        raw.fg.coef <- calculate.FG.fiability.coef(analysis.variables$F.beta.coef.list.1[[f.name]],
                                                                   analysis.variables$G.coef.list.1[[f.name]])
                        
                        analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]] <<- raw.fg.coef
                        analysis.variables$FG.accuracy.coef.list.1[[f.name]] <<- min(raw.fg.coef, 
                                                                                   abs(raw.fg.coef - analysis.variables$theta.coef.list.1[[f.name]]))
                        if(!is.null(analysis.variables$clusters.relations.id.list.1) && !is.null(analysis.variables$clusters.relations.id.list.1[[f.name]]))
                        {
                            ref.markers <- 1:ncol(fcs.ref@exprs)
                            names(ref.markers) <- colnames(fcs.ref@exprs)
                            test.markers <- 1:ncol(fcs.test.1@exprs)
                            names(test.markers) <- colnames(fcs.test.1@exprs)
                            cl.to.plot.ref <- 1:length(fcs.ref.info)
                            names(cl.to.plot.ref) <- sapply(1:length(cl.to.plot.ref), function(i){return(paste0("C",i))})
                            
                            ##=================================================================================================
                            insertUI(selector = paste0("#pl_",f.id,"_set"),
                                     where = "beforeEnd",
                                     ui = div(id=paste0("pl1_",f.id),
                                              div(id=paste0("pl1_",f.id,"_ref_container"),
                                                  imageOutput(paste0("pl1_",f.id,"_ref")),
                                                  div(id=paste0("pl1_",f.id,"_test_container_m"),
                                                      div(
                                                          selectInput(paste0("pl1_",f.id,"_ref_cl"), "Reference populations/clusters", 
                                                                      choices = cl.to.plot.ref, 
                                                                      multiple = TRUE, 
                                                                      selected = cl.to.plot.ref),
                                                          class="pl_cl_highlight"
                                                      ),
                                                      
                                                      class="pl_m"
                                                  ),
                                                  class="plot_container"
                                              ),
                                              div(id=paste0("pl1_",f.id,"_test_container"),
                                                  imageOutput(paste0("pl1_",f.id,"_test")),
                                                  div(id=paste0("pl1_",f.id,"_test_container_m"),
                                                      div(
                                                          selectInput(paste0("pl1_",f.id,"_test_cl"), "Annotated populations", 
                                                                      choices = cl.to.plot.test.1, 
                                                                      multiple = TRUE, 
                                                                      selected = cl.to.plot.test.1),
                                                          class="pl_cl_highlight"
                                                      ),
                                                      class="pl_m"
                                                  ),
                                                  class="plot_container"
                                              ),
                                              div(id=paste0("pl1_",f.id,"_scores"),
                                                  selectInput(paste0("pl1_",f.id,"_m1"), "Marker 1", choices = ref.markers, selected = 1),
                                                  selectInput(paste0("pl1_",f.id,"_m2"), "Marker 2", choices = ref.markers, selected = 1),
                                                  class = "pl_score"),
                                              class="pl_div"
                                        )
                            )
                            ##=================================================================================================
                            
                            output[[paste0("pl1_",f.id,"_ref")]] <- renderImage(
                            {
                                if(global.values$plot)
                                {
                                    outfile <- tempfile(fileext = ".png")
                                    if(!is.null(paste0("pl1_",f.id,"_ref_cl")) && length(input[[paste0("pl1_",f.id,"_ref_cl")]]) > 0)
                                    {
                                        rangeX <- c(max(min(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                               as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),-100),
                                                    min(max(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                                as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),10))
                                        rangeY <- rangeX
                                        mat.ref <- fcs.ref@exprs
                                        nmb.el <- nrow(mat.ref)
                                        rows.to.p <- unlist(sapply(fcs.ref.info[which(cl.to.plot.ref%in%input[[paste0("pl1_",f.id,"_ref_cl")]])],
                                                            function(t)
                                        {
                                            return(t[[2]])
                                        }))
                                        pl.color <- rep(3,nrow(mat.ref))
                                        pl.color[c(rows.to.p)] <- 2
                                        
                                        png(outfile)
                                        plot(mat.ref[, c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                as.integer(input[[paste0("pl1_",f.id,"_m2")]]))], 
                                             pch=".", main = "REFERENCE",
                                             col=pl.color, xlim=rangeX,ylim=rangeY)
                                        lapply(as.integer(input[[paste0("pl1_",f.id,"_ref_cl")]]), function(i)
                                        {
                                            xco <- mean(fcs.ref@exprs[c(fcs.ref.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]))])
                                            yco <- mean(fcs.ref@exprs[c(fcs.ref.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m2")]]))])
                                            text(xco,yco, paste0("C",i))
                                        })
                                        dev.off()
                                    }
                                        
                                    list(src=outfile)
                                }
                               
                                
                            }, deleteFile = T)
                            if(used.proj < 3 || (used.proj==3 && !is.null(global.values$fcs.files.fg.proj.1)) )
                            {
                                progress$inc(0.45, detail="REFERENCE PLOTTED")
                            }
                            else
                            {
                                progress$inc(0.2, detail="REFERENCE PLOTTED")
                            }
                            
                            output[[paste0("pl1_",f.id,"_test")]] <- renderImage(
                            {
                                if(global.values$plot)
                                {
                                    outfile <- tempfile(fileext = ".png")
                                    if(!is.null(paste0("pl1_",f.id,"_test_cl")) && length(input[[paste0("pl1_",f.id,"_test_cl")]]) > 0)
                                    {
                                        mat.test.1 <- fcs.test.1@exprs
                                        nmb.el <- nrow(mat.test.1)
                                        rangeX <- c(max(min(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                               as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),-100),
                                                    min(max(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                                as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),10))
                                        rangeY <- rangeX
                                        rows.to.p <- unlist(sapply(fcs.test.1.info[which(cl.to.plot.test.1%in%input[[paste0("pl1_",f.id,"_test_cl")]])],
                                                           function(t)
                                       {
                                           return(t[[2]])
                                       }))
                                        pl.color <- rep(3,nrow(mat.test.1))
                                        pl.color[c(rows.to.p)] <- 2
                                        
                                        
                                        png(outfile)
                                        plot(mat.test.1[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                           as.integer(input[[paste0("pl1_",f.id,"_m2")]]))],
                                             pch=".", main = "1st FILE",
                                             col=pl.color, xlim=rangeX,ylim=rangeY)
                                        lapply(as.integer(input[[paste0("pl1_",f.id,"_test_cl")]]), function(i)
                                        {
                                            xco <- mean(fcs.test.1@exprs[c(fcs.test.1.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]))])
                                            yco <- mean(fcs.test.1@exprs[c(fcs.test.1.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m2")]]))])
                                            text(xco,yco, paste0("C",i))
                                        })
                                        dev.off()
                                    }
                                    
                                    list(src=outfile)
                                }
                            }, deleteFile = T)
                            if(used.proj < 3 || (used.proj==3 && !is.null(global.values$fcs.files.fg.proj.1)) )
                            {
                                progress$inc(0.45, detail="1st FILE PLOTTED")
                            }
                            else
                            {
                                progress$inc(0.2, detail="1st FILE PLOTTED")
                            }
                            
                            
                            if(!is.null(analysis.variables$FG.accuracy.coef.list.1[[f.name]]) && global.values$plot)
                            {
                                t <- HTML("<p><strong>Good</strong></p>")
                                if(analysis.variables$FG.accuracy.coef.list.1[[f.name]] < 0.6)
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Higly Innacurate</u></strong></p>")
                                }
                                else if(analysis.variables$FG.accuracy.coef.list.1[[f.name]] > 0.6 && analysis.variables$FG.accuracy.coef.list.1[[f.name]] < 0.9)
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Probably Inneficient</u></strong></p>")
                                }
                                else
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Highly Accurate</u></strong></p>")
                                }
                                insertUI(selector=paste0("#pl1_",f.id,"_scores"),
                                         where = "beforeEnd",
                                         ui = div(h4(paste0("ACCURACY (FG) = ",analysis.variables$FG.accuracy.coef.list.1[[f.name]])),
                                              t)
                                )
                                
                                if(!is.null(analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]]))
                                {
                                    insertUI(selector=paste0("#pl1_",f.id,"_scores"),
                                             where = "beforeEnd",
                                             ui = h5(paste0("Raw FG score = ",analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]]))
                                    )
                                }
                                if(!is.null(analysis.variables$theta.coef.list.1[[f.name]]))
                                {
                                    insertUI(selector=paste0("#pl1_",f.id,"_scores"),
                                             where = "beforeEnd",
                                             ui = h5(paste0("Theta correction = ",analysis.variables$theta.coef.list.1[[f.name]]))
                                    )
                                }
                            }
                            if(!is.null(analysis.variables$F.beta.coef.list.1[[f.name]]) && global.values$plot)
                            {
                                insertUI(selector=paste0("#pl1_",f.id,"_scores"),
                                         where = "beforeEnd",
                                         ui = h5(paste0("Fbeta score = ",analysis.variables$F.beta.coef.list.1[[f.name]]))
                                )
                            }
                            if(!is.null(analysis.variables$G.coef.list.1[[f.name]]) && global.values$plot)
                            {
                                insertUI(selector=paste0("#pl1_",f.id,"_scores"),
                                         where = "beforeEnd",
                                         ui = h5(paste0("G score = ",analysis.variables$G.coef.list.1[[f.name]]))
                                )
                            }
                            
                            progress$inc(0.1,detail="SCORES FILE 1 CALCULATED")
                            
                        }
                        if(used.proj == 3)
                        {
                            progress$close()
                        }
                    }
                
                #VAL.2=======================================================================================================
                    if(used.proj == 3 && !is.null(global.values$fcs.files.fg.proj.2))
                    {
                        progress <- Progress$new()
                        progress$set(message="COMPARING FILES")
                        col.cl.1 <- as.integer(input$clust_col_selection_fg_r)
                        col.cl.2 <- as.integer(input$clust_col_selection_fg_2)
                        
                        fcs.test.2 <- 0
                        fcs.test.2.info <- 0
                        val.2 <- 0
                        
                        if(used.proj == 3)
                        {
                            fcs.ref.info <- FPH.get.file.information(fcs.ref, as.integer(input$clust_col_selection_fg_r))
                            if(!is.null(global.values$fcs.files.fg.proj.2[[f.id]]))
                            {
                                fcs.test.2 <- global.values$fcs.files.fg.proj.2[[f.id]]
                                fcs.test.2 <- FPH.create.artificials.clusters(fcs.ref, fcs.test.2, col.cl.1, col.cl.2)
                                fcs.test.2.info <- FPH.get.file.information(fcs.test.2, col.cl.2)
                                
                                val.2 <- FPH.get.clusters.associations(fcs.ref.info, fcs.test.2.info)
                                analysis.variables$clusters.relations.list.2[[f.name]] <- val.2[[3]]
                                analysis.variables$clusters.relations.id.list.2[[f.name]] <- val.2[[2]]
                                analysis.variables$theta.coef.list.2[[f.name]] <- val.2[[1]]
                            }
                        }
                        
                        cl.to.plot.test.2 <- 1:length(fcs.test.2.info)
                        names(cl.to.plot.test.2) <- sapply(1:length(cl.to.plot.test.2), function(i){return(paste0("C",i))})
                        
                        
                        
                        #========================================================================================================
                        prec.coef.1 <- 0
                        prec.coef.2 <- 0
                        rec.coef.1 <- 0
                        rec.coef.2 <- 0
                        
                        if(!is.null(global.values$fcs.files.fg.proj.2))
                        {
                            prec.coef.2 <- mean(sapply(1:(length(val.2[[2]])-1), function(c)
                            {
                                return(val.2[[3]][[c+1]][[1]])
                            }))
                            rec.coef.2 <- mean(sapply(1:(length(val.2[[2]])-1), function(c)
                            {
                                return(val.2[[3]][[c+1]][[2]])
                            }))
                        }
                        
                        analysis.variables$F.beta.coef.list.2[[f.name]] <<- calculate.F.beta.coef(prec.coef.2, rec.coef.2)
                        analysis.variables$G.coef.list.2[[f.name]] <<- calculate.G.coef(prec.coef.2, rec.coef.2)
                        
                        raw.fg.coef <- calculate.FG.fiability.coef(analysis.variables$F.beta.coef.list.2[[f.name]],
                                                                   analysis.variables$G.coef.list.2[[f.name]])
                        
                        analysis.variables$raw.FG.accuracy.coef.list.2[[f.name]] <<- raw.fg.coef
                        analysis.variables$FG.accuracy.coef.list.2[[f.name]] <<- min(raw.fg.coef, 
                                                                                   abs(raw.fg.coef - analysis.variables$theta.coef.list.2[[f.name]]))
                        
                        if(!is.null(analysis.variables$clusters.relations.id.list.2) && !is.null(analysis.variables$clusters.relations.id.list.2[[f.name]]))
                        {
                            ref.markers <- 1:ncol(fcs.ref@exprs)
                            names(ref.markers) <- colnames(fcs.ref@exprs)
                            test.markers <- 1:ncol(fcs.test.2@exprs)
                            names(test.markers) <- colnames(fcs.test.2@exprs)
                            cl.to.plot.ref <- 1:length(fcs.ref.info)
                            names(cl.to.plot.ref) <- sapply(1:length(cl.to.plot.ref), function(i){return(paste0("C",i))})
                            
                            #==============================================================================================
                            insertUI(selector = paste0("#pl_",f.id,"_set"),
                                     where = "beforeEnd",
                                     ui = div(id=paste0("pl2_",f.id),
                                              div(id=paste0("pl2_",f.id,"_ref_container"),
                                                  imageOutput(paste0("pl2_",f.id,"_ref")),
                                                  div(id=paste0("pl2_",f.id,"_test_container_m"),
                                                      div(
                                                          selectInput(paste0("pl2_",f.id,"_ref_cl"), "Reference populations/clusters", 
                                                                      choices = cl.to.plot.ref, 
                                                                      multiple = TRUE, 
                                                                      selected = cl.to.plot.ref),
                                                          class="pl_cl_highlight"
                                                      ),
                                                      
                                                      class="pl_m"
                                                  ),
                                                  class="plot_container"
                                              ),
                                              div(id=paste0("pl2_",f.id,"_test_container"),
                                                  imageOutput(paste0("pl2_",f.id,"_test")),
                                                  div(id=paste0("pl2_",f.id,"_test_container_m"),
                                                      div(
                                                          selectInput(paste0("pl2_",f.id,"_test_cl"), "Annotated populations", 
                                                                      choices = cl.to.plot.test.2, 
                                                                      multiple = TRUE, 
                                                                      selected = cl.to.plot.test.2),
                                                          class="pl_cl_highlight"
                                                      ),
                                                      class="pl_m"
                                                  ),
                                                  class="plot_container"
                                              ),
                                              div(id=paste0("pl2_",f.id,"_scores"),
                                                  selectInput(paste0("pl2_",f.id,"_m1"), "Marker 1", choices = ref.markers, selected = 1),
                                                  selectInput(paste0("pl2_",f.id,"_m2"), "Marker 2", choices = ref.markers, selected = 1),
                                                  class = "pl_score"),
                                              class="pl_div"
                                     )
                            )
                            #==============================================================================================
                            
                            output[[paste0("pl2_",f.id,"_ref")]] <- renderImage(
                            {
                                if(global.values$plot)
                                {
                                    outfile <- tempfile(fileext = ".png")
                                    if(!is.null(paste0("pl2_",f.id,"_ref_cl")) && length(input[[paste0("pl2_",f.id,"_ref_cl")]]) > 0)
                                    {
                                        mat.ref <- fcs.ref@exprs
                                        nmb.el <- nrow(mat.ref)
                                        rows.to.p <- unlist(sapply(fcs.ref.info[which(cl.to.plot.ref%in%input[[paste0("pl2_",f.id,"_ref_cl")]])],
                                                                   function(t)
                                        {
                                            return(t[[2]])
                                        }))
                                        pl.color <- rep(3,nrow(mat.ref))
                                        pl.color[c(rows.to.p)] <- 2
                                        rangeX <- c(max(min(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                               as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),-100),
                                                    min(max(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                                as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),10))
                                        rangeY <- rangeX
                                        
                                        
                                        png(outfile)
                                        plot(mat.ref[,c(as.integer(input[[paste0("pl2_",f.id,"_m1")]]),
                                                        as.integer(input[[paste0("pl2_",f.id,"_m2")]]))], 
                                             pch=".", main = "REFERENCE",
                                             col=pl.color, xlim=rangeX,ylim=rangeY)
                                        lapply(as.integer(input[[paste0("pl2_",f.id,"_ref_cl")]]), function(i)
                                        {
                                            xco <- mean(fcs.ref@exprs[c(fcs.ref.info[[i]][[2]]), c(as.integer(input[[paste0("pl2_",f.id,"_m1")]]))])
                                            yco <- mean(fcs.ref@exprs[c(fcs.ref.info[[i]][[2]]), c(as.integer(input[[paste0("pl2_",f.id,"_m2")]]))])
                                            text(c(xco,yco), paste0("C",i))
                                        })
                                        dev.off()
                                    }
                                    
                                    list(src=outfile)
                                }
                                
                            }, deleteFile = T)
                            progress$inc(0.2, detail="REFERENCE PLOTTED")
                            
                            output[[paste0("pl2_",f.id,"_test")]] <- renderImage(
                            {
                                if(global.values$plot)
                                {
                                    outfile <- tempfile(fileext = ".png")
                                    if(!is.null(paste0("pl2_",f.id,"_test_cl")) && length(input[[paste0("pl2_",f.id,"_test_cl")]]) > 0)
                                    {
                                        mat.test.2 <- fcs.test.2@exprs
                                        nmb.el <- nrow(mat.test.2)
                                        rows.to.p <- unlist(sapply(fcs.test.2.info[which(cl.to.plot.test.2%in%input[[paste0("pl2_",f.id,"_test_cl")]])],
                                                                   function(t)
                                        {
                                             return(t[[2]])
                                        }))
                                        pl.color <- rep(3,nrow(mat.test.2))
                                        pl.color[c(rows.to.p)] <- 2
                                        rangeX <- c(max(min(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                               as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),-100),
                                                    min(max(0,fcs.ref@exprs[,c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]),
                                                                                as.integer(input[[paste0("pl1_",f.id,"_m2")]]))]),10))
                                        rangeY <- rangeX
                                        
                                        
                                        png(outfile)
                                        plot(mat.test.2[,c(as.integer(input[[paste0("pl2_",f.id,"_m1")]]),
                                                           as.integer(input[[paste0("pl2_",f.id,"_m2")]]))],
                                             pch=".", main = "2nd FILE",
                                             col=pl.color, xlim=rangeX,ylim=rangeY)
                                        lapply(as.integer(input[[paste0("pl2_",f.id,"_test_cl")]]), function(i)
                                        {
                                            xco <- mean(fcs.test.2@exprs[c(fcs.test.2.info[[i]][[2]]), c(as.integer(input[[paste0("pl2_",f.id,"_m1")]]))])
                                            yco <- mean(fcs.test.2@exprs[c(fcs.test.2.info[[i]][[2]]), c(as.integer(input[[paste0("pl2_",f.id,"_m2")]]))])
                                            text(c(xco,yco), paste0("C",i))
                                        })
                                        dev.off()
                                    }
                                    
                                    list(src=outfile)
                                }
                                
                            }, deleteFile = T)
                            progress$inc(0.2, detail="2nd FILE PLOTTED")
                            
                            if(!is.null(analysis.variables$FG.accuracy.coef.list.2[[f.name]]) && global.values$plot)
                            {
                                t <- HTML("<p><strong>Good</strong></p>")
                                if(analysis.variables$FG.accuracy.coef.list.2[[f.name]] < 0.6)
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Higly Innacurate</u></strong></p>")
                                }
                                else if(analysis.variables$FG.accuracy.coef.list.2[[f.name]] > 0.6 && analysis.variables$FG.accuracy.coef.list.2[[f.name]] < 0.9)
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Probably Inneficient</u></strong></p>")
                                }
                                else
                                {
                                    t <- HTML("<p><strong>Clustering Method <u>Highly Accurate</u></strong></p>")
                                }
                                insertUI(selector=paste0("#pl2_",f.id,"_scores"),
                                         where = "beforeEnd",
                                         ui = div(h4(paste0("ACCURACY (FG) = ",analysis.variables$FG.accuracy.coef.list.2[[f.name]])),
                                                  t)
                                )
                                if(!is.null(analysis.variables$raw.FG.accuracy.coef.list.2[[f.name]]))
                                {
                                    insertUI(selector=paste0("#pl2_",f.id,"_scores"),
                                             where = "beforeEnd",
                                             ui = h5(paste0("Raw FG score = ",analysis.variables$raw.FG.accuracy.coef.list.2[[f.name]]))
                                    )
                                }
                                if(!is.null(analysis.variables$theta.coef.list.2[[f.name]]))
                                {
                                    insertUI(selector=paste0("#pl2_",f.id,"_scores"),
                                             where = "beforeEnd",
                                             ui = h5(paste0("Theta correction = ",analysis.variables$theta.coef.list.2[[f.name]]))
                                    )
                                }
                            }
                            if(!is.null(analysis.variables$F.beta.coef.list.2[[f.name]]) && global.values$plot)
                            {
                                insertUI(selector=paste0("#pl2_",f.name,"_scores"),
                                         where = "beforeEnd",
                                         ui = h5(paste0("Fbeta score = ",analysis.variables$F.beta.coef.list.2[[f.name]]))
                                )
                            }
                            if(!is.null(analysis.variables$G.coef.list.2[[f.name]]) && global.values$plot)
                            {
                                insertUI(selector=paste0("#pl2_",f.name,"_scores"),
                                         where = "beforeEnd",
                                         ui = h5(paste0("G score = ",analysis.variables$G.coef.list.2[[f.name]]))
                                )
                            }
                            progress$inc(0.1,detail="SCORES FILE 2 CALCULATED")
                            progress$close()
                        }
                    }
            })
        }
        
        
        shinyjs::enable("compare")
    })
    
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #=======================MFI BASED SIMILARITY===========================================================================
    #======================================================================================================================
    
    observe(
    {
        if(!is.null(global.values$fcs.files.mfi.proj.1))
        {
            p1.col <- 1:ncol(global.values$fcs.files.mfi.proj.1[[1]]@exprs)
            names(p1.col) <- colnames(global.values$fcs.files.mfi.proj.1[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_mfi_1", "Clusters column", choices=p1.col, selected=length(p1.col)-1)
            updateSelectInput(session, "markers_selection_mfi_1", "Markers to use", choices=p1.col, selected=p1.col[1:(length(p1.col)-2)])
            updateSelectInput(session, "ev_col_selection_mfi_1", "Clusters size column", choices=p1.col, selected=length(p1.col))
            
            lapply(global.values$fcs.files.mfi.proj.1, function(f)
            {
                f.info <- FPH.get.file.information(f,p1.col[[length(p1.col)]]-1)
                cl.info <- "<p><BLOCKQUOTE>"
                size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                lapply(1:length(f.info),function(i)
                {
                    cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                })
                cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                insertUI("#p1_files_list",
                         "beforeEnd",
                         ui = div(id=paste0("p1_file_",f@description$FILENAME,"_info"),
                                  h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
                                  div(HTML(cl.info),class="menu_file_info"),
                                  class="menu_file"
                         )
                )
            })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_mfi_1", "Clusters column", choices=list(" "=0))
            updateSelectInput(session, "markers_selection_mfi_1", "Markers to use", choices=list(" "=0))
            updateSelectInput(session, "ev_col_selection_mfi_1", "Clusters size column", choices=list(" "=0))
        }
    })
    
    observeEvent(input$proj1_selection_mfi,
    {
         shinyjs::disable("proj1_selection_mfi")
         shinyjs::disable("calculate")
         
         removeUI("#p1_files_list")
         insertUI("#left_menu_items",
                  "beforeEnd",
                  div(id="p1_files_list",
                      class="menu_set"))
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "FlowFrames"
         m[1,2] = "*.csv;*.fcs"
         temp.files <- choose.files(filters = m,multi = FALSE)
         global.values$fcs.files.mfi.proj.1 <- NULL
         if(length(temp.files) != 0)
         {
             progress <- Progress$new()
             progress$set(message="LOADING FILE")
             
             lapply(temp.files, function(f)
             {
                 l <- length(f)
                 if(grepl("csv",f))
                 {
                     x <- as.matrix(read.csv(f))
                     x <- flowFrame(x)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.mfi.proj.1 <<- c(global.values$fcs.files.mfi.proj.1, x)
                 }
                 else
                 {
                     x <- read.FCS(f,emptyValue = FALSE)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.mfi.proj.1 <<- c(global.values$fcs.files.mfi.proj.1, x)
                 }
                 names(global.values$fcs.files.mfi.proj.1)[[length(global.values$fcs.files.mfi.proj.1)]] <<- substr(f,1,nchar(f)-4)
             })
             progress$inc(1,detail="FILE LOADED")
             progress$close()
         }
         
         shinyjs::enable("proj1_selection_mfi")
         shinyjs::enable("compare")
     })
    
    observeEvent(input$proj1_update_mfi,
    {
         if(!is.null(global.values$fcs.files.mfi.proj.1) && !is.null(input$clust_col_selection_mfi_1) && input$clust_col_selection_mfi_1 != "")
         {
             lapply(global.values$fcs.files.mfi.proj.1, function(f)
             {
                 f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_mfi_1))
                 if(length(f.info)>0)
                 {
                     cl.info <- "<p><BLOCKQUOTE>"
                     removeUI(paste0("#p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
                     size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                     lapply(1:length(f.info),function(i)
                     {
                         cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                     })
                     cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                     insertUI("#p1_files_list",
                              "beforeEnd",
                              ui = div(id=paste0("p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                       h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
                                       div(HTML(cl.info),class="menu_file_info"),
                                       class="menu_file"
                              )
                     )
                 }
             })
         }
    })
    
    
    
    observe(
    {
        if(!is.null(global.values$fcs.files.mfi.proj.2))
        {
            p2.col <- 1:ncol(global.values$fcs.files.mfi.proj.2[[1]]@exprs)
            names(p2.col) <- colnames(global.values$fcs.files.mfi.proj.2[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_mfi_2", "Clusters column", choices=p2.col, selected=length(p2.col)-1)
            updateSelectInput(session, "markers_selection_mfi_2", "Markers to use", choices=p2.col, selected=p2.col[1:(length(p2.col)-2)])
            updateSelectInput(session, "ev_col_selection_mfi_2", "Clusters size column", choices=p2.col, selected=length(p2.col))
            
            lapply(global.values$fcs.files.mfi.proj.2, function(f)
            {
                f.info <- FPH.get.file.information(f,p2.col[[length(p2.col)-1]])
                cl.info <- "<p><BLOCKQUOTE>"
                size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                lapply(1:length(f.info),function(i)
                {
                    cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                })
                cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                insertUI("#p2_files_list",
                         "beforeEnd",
                         ui = div(id=paste0("p2_file_",f@description$FILENAME,"_info"),
                                  h5(paste0("FILE 2: ", basename(f@description$FILENAME)), " ", size, " events"),
                                  div(HTML(cl.info),class="menu_file_info"),
                                  class="menu_file"
                         )
                )
            })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_mfi_2", "Clusters column", choices=list(" "=0))
            updateSelectInput(session, "markers_selection_mfi_2", "Markers to use", choices=list(" "=0))
            updateSelectInput(session, "ev_col_selection_mfi_2", "Clusters size column", choices=list(" "=0))
        }
    })

    observeEvent(input$proj2_selection_mfi,
    {
         shinyjs::disable("proj2_selection_mfi")
         shinyjs::disable("calculate")
         
         removeUI("#p2_files_list")
         insertUI("#left_menu_items",
                  "beforeEnd",
                  div(id="p2_files_list",
                      class="menu_set"))
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "FlowFrames"
         m[1,2] = "*.csv;*.fcs"
         temp.files <- choose.files(filters = m,multi = FALSE)
         global.values$fcs.files.mfi.proj.2 <- NULL
         if(length(temp.files) > 0)
         {
             progress <- Progress$new()
             progress$set(message="LOADING FILES")
             
             lapply(temp.files, function(f)
             {
                 if(grepl("csv",f))
                 {
                     x <- as.matrix(read.csv(f))
                     x <- flowFrame(x)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.mfi.proj.2 <<- c(global.values$fcs.files.mfi.proj.2, x)
                 }
                 else
                 {
                     x <- read.FCS(f,emptyValue = FALSE)
                     lapply(1:ncol(x@exprs), function(i)
                     {
                         nx <- x@description[[paste0("$P",i,"S")]]
                         if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                         {
                             colnames(x)[i] <<- nx
                         }
                     })
                     global.values$fcs.files.mfi.proj.2 <<- c(global.values$fcs.files.mfi.proj.2, x)
                 }
                 names(global.values$fcs.files.mfi.proj.2)[[length(global.values$fcs.files.mfi.proj.2)]] <<- substr(f,1,nchar(f)-4)
             })
             
             progress$inc(1,detail="FILE LOADED")
             progress$close()
         }
         
         shinyjs::enable("proj2_selection_mfi")
         shinyjs::enable("calculate")
     })
    
    observeEvent(input$proj2_update_mfi,
    {
         if(!is.null(global.values$fcs.files.mfi.proj.2) && !is.null(input$clust_col_selection_mfi_2) && input$clust_col_selection_mfi_2 != "")
         {
             lapply(global.values$fcs.files.mfi.proj.2, function(f)
             {
                 f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_mfi_2))
                 if(length(f.info)>0)
                 {
                     cl.info <- "<p><BLOCKQUOTE>"
                     removeUI(paste0("#p2_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
                     size <- sum(sapply(f.info, function(f){return(f[[1]])}))
                     lapply(1:length(f.info),function(i)
                     {
                         cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
                     })
                     cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
                     insertUI("#p2_files_list",
                              "beforeEnd",
                              ui = div(id=paste0("p2_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
                                       h5(paste0("FILE 2: ", basename(f@description$FILENAME)), " ", size, " events"),
                                       div(HTML(cl.info),class="menu_file_info"),
                                       class="menu_file"
                              )
                     )
                 }
             })
         }
    })
    
    
    
    
    observeEvent(input$calculate,
    {
         mfi.variables$sent.list <- NULL
         mfi.variables$smfi.list <- NULL
         mfi.variables$s.list <- NULL
         mfi.variables$s.id.list <- NULL
         global.values$plot <- TRUE
         
         shinyjs::disable("calculate")
         
         #UPDATE PLOT SECTION===============================================================
         removeUI("#plots_section_mfi")
         insertUI("#MFI_description",
                  "afterEnd",
                  div(id="plots_section_mfi")
         )
         
         if(!is.null(global.values$fcs.files.mfi.proj.1) && !is.null(global.values$fcs.files.mfi.proj.2))
         {
             fcs.files <- 0
             if(length(global.values$fcs.files.mfi.proj.1) <= length(global.values$fcs.files.mfi.proj.2))
             {
                 fcs.files <- global.values$fcs.files.mfi.proj.1
             }
             else
             {
                 fcs.files <- global.values$fcs.files.mfi.proj.2
             }
             
             progress <- Progress$new()
             progress$set(message="CALCULATING RESEMBLANCE")
             lapply(c(1:length(names(fcs.files))), function(f.id)
             {
                 #HIDE AND SEEK UI============================================================================================
                 f.name <- names(global.values$fcs.files.mfi.proj.1)[[f.id]]
                 insertUI(selector = "#plots_section_mfi",
                          where = "beforeEnd",
                          ui = div(id=paste0("pl_",f.id,"_set"),
                                   h4(basename(f.name)),
                                   actionButton(paste0("pl_",f.id,"_set_h"), "Hide"),
                                   actionButton(paste0("pl_",f.id,"_set_s"), "Show"),
                                   class="pl_set"
                          )
                 )
                 
                 observeEvent(input[[paste0("pl_",f.id,"_set_h")]],
                 {
                      shinyjs::hide(paste0("pl_",f.id))
                 })
                 
                 observeEvent(input[[paste0("pl_",f.id,"_set_s")]],
                 {
                      shinyjs::show(paste0("pl_",f.id))
                 })
                 
                 
                 #FILES INIT==================================================================================================
                 fcs.1 <- global.values$fcs.files.mfi.proj.1[[f.id]]
                 fcs.1.info <- FPH.get.file.information(fcs.1, as.integer(input$clust_col_selection_mfi_1))
                 fcs.2 <- global.values$fcs.files.mfi.proj.2[[f.id]]
                 fcs.2.info <- FPH.get.file.information(fcs.2, as.integer(input$clust_col_selection_mfi_1))
                 
                 #VALUES INIT=================================================================================================
                 col.cl.1 <- as.integer(input$clust_col_selection_mfi_1)
                 col.cl.2 <- as.integer(input$clust_col_selection_mfi_2)
                 markers.cols.1 <- as.integer(input$markers_selection_mfi_1)
                 names(markers.cols.1) <- colnames(fcs.1@exprs)[as.integer(markers.cols.1)]
                 markers.cols.2 <- as.integer(input$markers_selection_mfi_2)
                 names(markers.cols.2) <- colnames(fcs.2@exprs)[as.integer(markers.cols.2)]
                 nmb.ev.col.1 <- as.integer(input$ev_col_selection_mfi_1)
                 nmb.ev.col.2 <- as.integer(input$ev_col_selection_mfi_2)
                 
                 val <- FPH.get.clusters.resemblance(fcs.1,fcs.2,markers.cols.1,markers.cols.2,col.cl.1,col.cl.2,
                                                     nmb.ev.col.1,nmb.ev.col.2)
                 mfi.variables$s.id.list[[f.name]] <<- val[[1]]
                 mfi.variables$s.list[[f.name]] <<- val[[2]]
                 mfi.variables$s.mfi.list[[f.name]] <<- val[[3]]
                 mfi.variables$s.ent.list[[f.name]] <<- val[[4]]
                 
                 cl.to.plot.1 <- 1:length(fcs.1.info)
                 names(cl.to.plot.1) <- sapply(1:length(cl.to.plot.1), function(i){return(paste0("C",i))})
                 cl.to.plot.2 <- 1:length(fcs.2.info)
                 names(cl.to.plot.2) <- sapply(1:length(cl.to.plot.2), function(i){return(paste0("C",i))})
                 progress$inc(0.7,detail="FILES ANALYSED")
                 
                 
                 #PLOTS UI===================================================================================================
                 if(!is.null(mfi.variables$s.id.list) && !is.null(mfi.variables$s.id.list[[f.name]]))
                 {
                     ##GENERATE UI===========================================================================================
                     insertUI(selector = paste0("#pl_",f.id,"_set"),
                              where = "beforeEnd",
                              ui = div(id=paste0("pl_",f.id),
                                       div(id=paste0("pl_",f.id,"_mfi"),
                                           div(id=paste0("pl_",f.id,"_mfi_int_plot"),
                                               imageOutput(paste0("pl_",f.id,"_s")),
                                               div(id=paste0("pl_",f.id,"_mfi_m"),
                                                   sliderInput(paste0("pl_",f.id,"_2_ky"), "S_MFI Threshold", 
                                                               min = 0, max = 1, value = 0, step = 0.001),
                                                   sliderInput(paste0("pl_",f.id,"_2_kx"), "S_EV Threshold", 
                                                               min = 0, max = 1, value = 0, step = 0.001),
                                                   class="pl_m"),
                                               class="plot_inter_plot"
                                           ),
                                           div(id=paste0("pl_",f.id,"_mfi_cl")),
                                           class="plot_mfi"
                                       ),
                                       class="pl_div"
                              )
                     )
                     ##PLOT=============================================================================================
                     
                     output[[paste0("pl_",f.id,"_s")]] <- renderImage(
                     {
                         if(global.values$plot)
                         {
                             removeUI(paste0("#pl_",f.id,"_mfi_cl"))
                             
                             outfile <- tempfile(fileext = ".png")
                             s.id <- unlist(mfi.variables$s.id.list[[f.name]])
                             s.mfi <- -1/(log(unlist(mfi.variables$s.mfi.list[[f.name]]))-1)
                             s.ent <- -1/(log(unlist(mfi.variables$s.ent.list[[f.name]]))-1)
                             s <- mfi.variables$s.list[[f.name]]
                             
                             insertUI(paste0("#pl_",f.id,"_mfi_int_plot"),
                                  "afterEnd",
                                  div(id=paste0("pl_",f.id,"_mfi_cl"),
                                      h4(paste0("Clusters resemblance")),
                                      class="plot_mfi_res"
                                  )
                             )
                             
                             png(outfile)
                             plot(s.ent,s.mfi, pch="x", main = "MFI similarity VS Cluster size",xlim=c(0,1), ylim=c(0,1))
                             text(0.2,0.1,"No resemblance",col = "blue")
                             text(0.8,0.9,"High resemblance", col="blue")
                             rect(0,0,0.6,0.6,col="red")
                             rect(0,0.6,0.6,1,col="orange")
                             rect(0.6,0,1,0.6,col="orange")
                             rect(0.6,0.6,1,1,col="green")
                             
                             points(s.ent,s.mfi)
                             x0 <- as.numeric(input[[paste0("pl_",f.id,"_2_ky")]])
                             y0 <- as.numeric(input[[paste0("pl_",f.id,"_2_kx")]])
                             abline(v=x0)
                             abline(h=y0)
                             nmb.cl <- 0
                             lapply(1:length(s.mfi), function(i)
                             {
                                 if(s.mfi[i] >= y0 && s.ent[i] >= x0)
                                 {
                                     if(s.id[1] == 1)
                                     {
                                         insertUI(paste0("#pl_",f.id,"_mfi_cl"),
                                                  "beforeEnd",
                                                  div(paste0("C",i,"--C",s.id[i+1],": s=",s[i]))
                                         )
                                     }
                                     else
                                     {
                                         insertUI(paste0("#pl_",f.id,"_mfi_cl"),
                                                  "beforeEnd",
                                                  div(h5(paste0("C",s.id[i+1],"--C",i,":")),
                                                      h6("s=",s[i]))
                                         )
                                     }
                                     nmb.cl <<- nmb.cl + 1
                                 }
                                 
                             })
                             text(0.2,0.85,paste0(nmb.cl," clusters(",100*nmb.cl/length(s),"%)"))
                             dev.off()
                             
                             list(src=outfile)
                         }
                         
                     }, deleteFile = T)
                     progress$inc(0.3, detail="RESEMBLANCE COMPUTED")
                 }
             })
             progress$inc(0.3,detail="RESEMBLANCE COMPUTED")
             progress$close()
         }
         
         
         shinyjs::enable("calculate")
     })
    
    

}
