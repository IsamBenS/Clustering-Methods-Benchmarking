library(shiny)
library(shinydashboard)
library(flowCore)
library(gplots)
library(shinyjs)


server <- function(input, output, session)
{
    useShinyjs()
    #======================================================================================================================
    #======================REACTIVE VALUES=================================================================================
    #======================================================================================================================
    
    global.values <- reactiveValues(
        fcs.files.fg.ref = NULL,
        fcs.files.fg.proj.1 = NULL,
        fcs.files.fg.proj.1.annotated = NULL,
        fcs.files.fg.mapping = NULL,
        
        fcs.files.mfi.proj.1 = NULL,
        fcs.files.mfi.proj.2 = NULL,
        
        plot = TRUE,
        working.directory = getwd(),
        use.enriched.file = FALSE
        
    )
    
    analysis.variables <- reactiveValues(
        annotations.relations.list.1 = NULL,
        annotations.relations.id.list.1 = NULL,
        theta.coef.list.1 = NULL,
        F.beta.coef.list.1 = NULL,
        G.coef.list.1 = NULL,
        FG.accuracy.coef.list.1 = NULL,
        raw.FG.accuracy.coef.list.1 = NULL,
        purity.matrix = NULL,
        scores.table = NULL,
        params.table = NULL
    )
    
    mfi.variables <- reactiveValues(
        clusters.resemblance.list = NULL,
        s.id.list = NULL,
        sent.list = NULL,
        smfi.list = NULL,
        s.list = NULL
    )
    
    clustering.variables <- reactiveValues(
        available.methods.names = NULL,
        available.methods.parameters = NULL,
        added.keyword = NULL
    )
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #==========================================LOAD FILES==================================================================
    #======================================================================================================================
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.ref))
        {
            ref.col <- 1:ncol(global.values$fcs.files.fg.ref[[1]]@exprs)
            names(ref.col) <- colnames(global.values$fcs.files.fg.ref[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_fg_r", "Clusters column", choices=ref.col, selected=length(ref.col))
            updateSelectInput(session, "method_markers", "Markers to use", choices=ref.col, selected=NULL)
            
            # lapply(global.values$fcs.files.fg.ref, function(f)
            # {
            #     f.info <- FPH.get.file.information(f,ref.col[[length(ref.col)]])
            #     cl.info <- "<p><BLOCKQUOTE>"
            #     size <- sum(sapply(f.info, function(f){return(f[[1]])}))
            #     lapply(1:length(f.info),function(i)
            #     {
            #         cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
            #     })
            #     cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
            #     insertUI("#ref_files_list",
            #              "beforeEnd",
            #              ui = div(id=paste0("ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
            #                       h5(paste0("REFERENCE: ", basename(f@description$FILENAME)), " ", size, " events"),
            #                       div(HTML(cl.info),class="menu_file_info"),
            #                       class="menu_file"
            #              )
            #     )
            # })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_fg_r", "Clusters column", choices=list(" "=0))
        }
    })
    
    observe(
    {
        if(global.values$use.enriched.file)
        {
            shinyjs::disable("ref_use_same_fg")
            shinyjs::disable("proj1_selection_fg")
        }
        else
        {
            shinyjs::enable("ref_use_same_fg")
            shinyjs::enable("proj1_selection_fg")
        }
    })
    
    observeEvent(input$ref_use_same_fg,
    {
         removeUI("#download_enriched_file")
         clustering.variables$added.keyword <- NULL
         if(input$ref_use_same_fg)
         {
             disable("proj1_selection_fg")
             if(!is.null(global.values$fcs.files.fg.ref[[1]]))
             {
                 global.values$fcs.files.fg.proj.1 <- list()
                 global.values$fcs.files.fg.proj.1[[1]] <- global.values$fcs.files.fg.ref[[1]]
                 insertUI("#left_menu_items",
                          "beforeEnd",
                          div(id="p1_files_list",
                              class="menu_set"))
             }
             
         }
         else
         {
             enable("proj1_selection_fg")
             global.values$fcs.files.fg.proj.1 <- NULL
             removeUI("#p1_files_list")
         }
         
     })
    
    observeEvent(input$ref_selection_fg,
    {
         shinyjs::disable("ref_selection_fg")
         shinyjs::disable("compare")
         shinyjs::disable("save_results")
         
         removeUI("#ref_files_list")
         removeUI("#download_enriched_file")
         shinyjs::enable("run_clustering_button")
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "FlowFrames"
         m[1,2] = "*.csv;*.fcs"
         temp.files <- choose.files(filters = m,multi = FALSE)
         global.values$fcs.files.fg.ref <- NULL
         global.values$fcs.files.fg.mapping <- NULL
         clustering.variables$added.keyword <- NULL
         global.values$use.enriched.file <- F
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
         
         clustering.variables$available.methods.names <- NULL
         clustering.variables$available.methods.parameters <- NULL
         if(is.null(clustering.variables$available.methods))
         {
             temp.dir <- paste0(global.values$working.directory,"/methods_Folder/")
             methods.files <- list.files(temp.dir,pattern = ".R")
             lapply(methods.files, function(f)
             {
                 source(paste0(temp.dir,f))
                 clustering.variables$available.methods.names <- c(clustering.variables$available.methods.names, strsplit(f,".R"))
                 clustering.variables$available.methods.parameters <- c(clustering.variables$available.methods.parameters, list(fct.parameters))
                 names(clustering.variables$available.methods.parameters) <- clustering.variables$available.methods.names
             })
         }
         
         updateSelectInput(session, "method_selection_fg_r", "Select Methods", choices=clustering.variables$available.methods.names,
                           selected = clustering.variables$available.methods.names)
         
         
         shinyjs::enable("ref_selection_fg")
         shinyjs::enable("compare")
     })
    
    observeEvent(input$ref_mapping_fg,
    {
         shinyjs::disable("ref_mapping_fg")
         
         removeUI("#ref_files_list")
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "Mapping Matrix"
         m[1,2] = "*.csv"
         temp.file <- choose.files(filters = m,multi = FALSE)
         if(length(temp.file) != 0)
         {
             progress <- Progress$new()
             progress$set(message="LOADING FILE", value=0)
             global.values$fcs.files.fg.mapping <- as.matrix(read.csv(temp.file))
             names(global.values$fcs.files.fg.mapping) <<- substr(temp.file,1,nchar(temp.file)-4)
             progress$inc(1, detail="FILE LOADED")
             progress$close()
         }
         
         shinyjs::enable("ref_mapping_fg")
     })
    
    # observeEvent(input$ref_update_fg,
    # {
    #      if(!is.null(global.values$fcs.files.fg.ref) && !is.null(input$clust_col_selection_fg_r) && input$clust_col_selection_fg_r != "")
    #      {
    #          lapply(global.values$fcs.files.fg.ref, function(f)
    #          {
    #              f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_fg_r))
    #              if(length(f.info)>0)
    #              {
    #                  cl.info <- "<p><BLOCKQUOTE>"
    #                  removeUI(paste0("#ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
    #                  size <- sum(sapply(f.info, function(f){return(f[[1]])}))
    #                  lapply(1:length(f.info),function(i)
    #                  {
    #                      cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
    #                  })
    #                  cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
    #                  insertUI("#ref_files_list",
    #                           "beforeEnd",
    #                           ui = div(id=paste0("ref_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
    #                                    h5(paste0("REFERENCE: ", basename(f@description$FILENAME)), " ", size, " events"),
    #                                    div(HTML(cl.info),class="menu_file_info"),
    #                                    class="menu_file"
    #                           )
    #                  )
    #              }
    #          })
    #      }
    #  })
    
    observe(
    {
        if(!is.na(input$method_selection_fg_r) && !is.null(input$method_selection_fg_r) && input$method_selection_fg_r!="")
        {
            removeUI("#t_2_b_1")
            
            insertUI("#t_2_fr",
                     "beforeEnd",
                     ui = div(id="t_2_b_1", width="100%"))
            
            
            curr.method.name <- input$method_selection_fg_r
            insertUI("#t_2_b_1",
                     where = "afterBegin",
                     ui = h3(paste0("Parameters - ",curr.method.name), id=paste0("rmp_",curr.method.name,"_param"))
            )
            if(length(clustering.variables$available.methods.parameters[[curr.method.name]]) > 0)
            {
                lapply(1:length(clustering.variables$available.methods.parameters[[curr.method.name]]), function(curr.param.id)
                {
                    insertUI("#t_2_b_1",
                             where = "beforeEnd",
                             ui = box
                             (
                                 textInput(paste0("cm_",curr.method.name,"_param_",curr.param.id),
                                           paste0(names(clustering.variables$available.methods.parameters[[curr.method.name]])[curr.param.id]),
                                           value = clustering.variables$available.methods.parameters[[curr.method.name]][curr.param.id])
                             )
                    )
                })
            }
        }
    })
    
    
    
    
    
    
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.proj.1))
        {
            p1.col <- 1:ncol(global.values$fcs.files.fg.proj.1[[1]]@exprs)
            names(p1.col) <- colnames(global.values$fcs.files.fg.proj.1[[1]]@exprs)
            updateSelectInput(session, "clust_col_selection_fg_1", "Clusters column", choices=p1.col, selected=length(p1.col))
            # lapply(global.values$fcs.files.fg.proj.1, function(f)
            # {
            #     f.info <- FPH.get.file.information(f,p1.col[[length(p1.col)]])
            #     cl.info <- "<p><BLOCKQUOTE>"
            #     size <- sum(unlist(sapply(f.info, function(f.i){return(f.i[[1]])})))
            #     lapply(1:length(f.info),function(i)
            #     {
            #         cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
            #     })
            #     cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
            #     insertUI("#p1_files_list",
            #              "beforeEnd",
            #              ui = div(id=paste0("p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
            #                       h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
            #                       div(HTML(cl.info),class="menu_file_info"),
            #                       class="menu_file"
            #              )
            #     )
            # })
        }
        else
        {
            updateSelectInput(session, "clust_col_selection_fg_1", "Clusters column", choices=list(" "=""))
        }
    })
    
    observeEvent(input$proj1_selection_fg,
    {
         shinyjs::disable("proj1_selection_fg")
         shinyjs::disable("compare")
         shinyjs::disable("save_results")
         
         removeUI("#p1_files_list")
         removeUI("#download_enriched_file")
         shinyjs::enable("run_clustering_button")
         
         
         m <- matrix(nrow=1,ncol=2)
         m[1,1] = "FlowFrames"
         m[1,2] = "*.csv;*.fcs"
         temp.files <- choose.files(filters = m,multi = F)
         global.values$fcs.files.fg.proj.1 <- NULL
         clustering.variables$added.keyword <- NULL
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
    
    # observeEvent(input$proj1_update_fg,
    # {
    #      if(!is.null(global.values$fcs.files.fg.proj.1) && !is.null(input$clust_col_selection_fg_1) && input$clust_col_selection_fg_1 != "")
    #      {
    #          lapply(global.values$fcs.files.fg.proj.1, function(f)
    #          {
    #              f.info <- FPH.get.file.information(f,as.integer(input$clust_col_selection_fg_1))
    #              if(length(f.info)>0)
    #              {
    #                  cl.info <- "<p><BLOCKQUOTE>"
    #                  removeUI(paste0("#p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"))
    #                  size <- sum(sapply(f.info, function(f){return(f[[1]])}))
    #                  lapply(1:length(f.info),function(i)
    #                  {
    #                      cl.info <<- paste0(cl.info,"    Cluster ", i, " - ", f.info[[i]][[1]], " events (",as.integer(f.info[[i]][[1]]*10000000/size)/100000,"%) <br>")
    #                  })
    #                  cl.info <- paste0(cl.info,"</BLOCKQUOTE></p>")
    #                  insertUI("#p1_files_list",
    #                           "beforeEnd",
    #                           ui = div(id=paste0("p1_file_",substr(basename(f@description$FILENAME),1,nchar(basename(f@description$FILENAME))-4),"_info"),
    #                                    h5(paste0("FILE 1: ", basename(f@description$FILENAME)), " ", size, " events"),
    #                                    div(HTML(cl.info),class="menu_file_info"),
    #                                    class="menu_file"
    #                           )
    #                  )
    #              }
    #          })
    #      }
    #  })
    
    
    
    
    
    
    
    
    
    observeEvent(input$compare,
    {
         analysis.variables$annotations.relations.list.1 = NULL
         analysis.variables$annotations.relations.id.list.1 = NULL
         analysis.variables$theta.coef.list.1 = NULL
         analysis.variables$F.beta.coef.list.1 = NULL
         analysis.variables$G.coef.list.1 = NULL
         analysis.variables$FG.accuracy.coef.list.1 = NULL
         analysis.variables$raw.FG.accuracy.coef.list.1 = NULL
         
         global.values$plot <- TRUE
         
         
         shinyjs::disable("compare")
         shinyjs::disable("save_results")
         
         
         #UPDATE PLOT SECTION===============================================================
         removeUI("#t_3_2")
         insertUI("#t_3_1",
                  "afterEnd",
                  fluidRow(id="t_3_2")
         )
         removeUI("#t_3_3")
         insertUI("#t_3_2",
                  "afterEnd",
                  fluidRow(id="t_3_3")
         )
         
         fcs.files <- NULL
         used.proj <- 0
         if(!is.null(global.values$fcs.files.fg.ref))
         {
             fcs.files <- global.values$fcs.files.fg.ref
             lapply(c(1:length(names(fcs.files))), function(f.id)
             {
                 f.name <- names(fcs.files)[[f.id]]
                 fcs.ref <- fcs.files[[f.id]]
                 fcs.ref.info <- 0
                 
                 insertUI(selector = "#t_3_2",
                          where = "afterBegin",
                          ui = h3(paste("FILE:",basename(f.name)))
                 )
                 
                 rangeX <- c(max(min(0,fcs.ref@exprs),-100),min(max(10,fcs.ref@exprs),10))
                 rangeY <- rangeX
                 
                 if( !is.null(global.values$fcs.files.fg.proj.1) )
                 {
                     progress <- Progress$new()
                     progress$set(message="COMPARING FILES")
                     col.cl.1 <- as.integer(input$clust_col_selection_fg_r)
                     col.cl.2 <- as.integer(input$clust_col_selection_fg_1)
                     fcs.test.1 <- 0
                     fcs.test.1.info <- 0
                     val.1 <- 0
                     
                     fcs.ref.info <- FPH.get.file.information(fcs.ref, as.integer(input$clust_col_selection_fg_r))
                     if(!is.null(global.values$fcs.files.fg.proj.1[[f.id]]))
                     {
                         fcs.test.1.temp <- global.values$fcs.files.fg.proj.1[[f.id]]
                         p.matrix <- FPH.get.purity.matrix(fcs.ref, fcs.test.1.temp, col.cl.1, col.cl.2)
                         fcs.test.1 <- FPH.create.artificials.clusters(fcs.test.1.temp, p.matrix,col.cl.2) #NEW F1
                         #fcs.test.1 <- fcs.test.1.temp #Classic F1
                         fcs.test.1.info <- FPH.get.file.information(fcs.test.1, ncol(fcs.test.1@exprs))
                         fcs.test.1.temp.info <- FPH.get.file.information(fcs.test.1.temp, col.cl.2)
                         
                         global.values$fcs.files.fg.proj.1.annotated <<- fcs.test.1
                         
                         val.1 <- FPH.get.clusters.associations(fcs.ref.info, fcs.test.1.info)
                         analysis.variables$annotations.relations.list.1[[f.name]] <<- val.1[[3]]
                         analysis.variables$annotations.relations.id.list.1[[f.name]] <<- val.1[[2]]
                         analysis.variables$theta.coef.list.1[[f.name]] <<- val.1[[1]]
                         analysis.variables$purity.matrix <<- p.matrix
                     }
                     
                     
                     cl.to.plot.test.1 <- 1:length(fcs.test.1.info)
                     ord <- order(val.1[[2]][-1])
                     lapply(1:length(fcs.ref.info), function(i)
                     {
                         t <- which(val.1[[2]][-1]==i)
                         if(length(t)>0)
                         {
                             if(length(t)>1)
                             {
                                names(cl.to.plot.test.1)[t] <<- sapply(1:length(t), function(j){return(paste0("Annot-",i,"_",j))})
                             }
                             else
                             {
                                 names(cl.to.plot.test.1)[t] <<- paste0("Annot-",i)
                             }
                         }
                     })
                     cl.to.plot.test.1 <- cl.to.plot.test.1[ord]
                     
                     #========================================================================================================
                     prec.coef.1 <- 0
                     rec.coef.1 <- 0
                     
                     if(!is.null(global.values$fcs.files.fg.proj.1))
                     {
                         prec.coef.1 <- sapply(1:(length(val.1[[2]])-1), function(c)
                         {
                             return(val.1[[3]][[c+1]][[1]])
                         })
                         #prec.coef.1 <- mean(prec.coef.1)
                         rec.coef.1 <- sapply(1:(length(val.1[[2]])-1), function(c)
                         {
                             return(val.1[[3]][[c+1]][[2]])
                         })
                         #rec.coef.1 <- mean(rec.coef.1)
                     }
        
                     
                     threshold <- 0
                     if ( !is.na(input[[paste0("pl1_p_fct_ref",f.id,"_F1")]]) && !is.null(input[[paste0("pl1_p_fct_ref",f.id,"_F1")]]) )
                     {
                         threshold <- as.numeric(input[[paste0("pl1_p_fct_ref",f.id,"_F1")]])
                     }
                     analysis.variables$F.beta.coef.list.1[[f.name]] <<- calculate.F.beta.coef(prec.coef.1, rec.coef.1, threshold)
                     analysis.variables$G.coef.list.1[[f.name]] <<- calculate.G.coef(prec.coef.1, rec.coef.1, threshold)
                     
                     raw.fg.coef <- calculate.FG.fiability.coef(analysis.variables$F.beta.coef.list.1[[f.name]],
                                                                analysis.variables$G.coef.list.1[[f.name]])
                     
                     analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]] <<- raw.fg.coef
                     analysis.variables$FG.accuracy.coef.list.1[[f.name]] <<- min(raw.fg.coef, 
                                                                                  abs(raw.fg.coef - analysis.variables$theta.coef.list.1[[f.name]]))
                     
                     if(!is.null(analysis.variables$annotations.relations.id.list.1) && !is.null(analysis.variables$annotations.relations.id.list.1[[f.name]]))
                     {
                         ref.markers <- 1:ncol(fcs.ref@exprs)
                         names(ref.markers) <- colnames(fcs.ref@exprs)
                         test.markers <- 1:ncol(fcs.test.1@exprs)
                         names(test.markers) <- colnames(fcs.test.1@exprs)
                         cl.to.plot.ref <- 1:length(fcs.ref.info)
                         names(cl.to.plot.ref) <- sapply(1:length(cl.to.plot.ref), function(i){return(paste0("REF_POP-",i))})
                         
                         ##=================================================================================================
                         insertUI(selector = "#t_3_2",
                                  where = "beforeEnd",
                                  ui = div
                                  (
                                      # box
                                      # (
                                      #      id="t_3_2_1"
                                      #      imageOutput("pl_div_1_1_img")
                                      # ),
                                      box
                                      (
                                          id="t_3_2_2"
                                          #sliderInput("pl_div_1_2_scores_tp", "purity threshold", min=0, max=1,value=0.5, step=0.01),
                                          #sliderInput("pl_div_1_2_scores_fbeta", "Fbeta threshold", min=0, max=1,value=0.5, step=0.01),
                                      )
                                      
                                      # box
                                      # (
                                      #     id="t_3_2_3",
                                      #     imageOutput("pl_div_2_1_img")
                                      # ),
                                      # box
                                      # (
                                      #     id="t_3_2_4",
                                      #     tableOutput("pl_div_2_2_mat")
                                      # ),
                                      # 
                                      # box
                                      # (
                                      #     id="t_3_2_5",
                                      #     tableOutput("pl_div_3_1_mat")
                                      # ),
                                      # box
                                      # (
                                      #     id="t_3_2_6",
                                      #     tableOutput("pl_div_3_2_mat")
                                      # )
                                  )
                         )
                         insertUI(selector = "#t_3_3",
                                  where = "beforeEnd",
                                  ui = div
                                  (
                                      box
                                      (
                                          id="t_3_3_1",
                                          imageOutput("t_3_3_1_img"),
                                          div
                                          (
                                              selectInput(paste0("pl1_",f.id,"_ref_cl"), "Reference populations/clusters",
                                                          choices = cl.to.plot.ref,
                                                          multiple = TRUE,
                                                          selected = cl.to.plot.ref),
                                              class="img_cl"
                                          ),
                                          width=5
                                      ),
                                      box
                                      (
                                          id="t_3_3_2",
                                          imageOutput("t_3_3_2_img"),
                                          div
                                          (
                                              selectInput(paste0("pl1_",f.id,"_test_cl"), "Annotated populations",
                                                          choices = cl.to.plot.test.1,
                                                          multiple = TRUE,
                                                          selected = cl.to.plot.test.1),
                                              class="img_cl"
                                          ),
                                          width=5
                                      ),
                                      box
                                      (
                                          id="t_3_3_3",
                                          selectInput(paste0("pl1_",f.id,"_m1"), "Marker 1", choices = ref.markers, selected = 1),
                                          selectInput(paste0("pl1_",f.id,"_m2"), "Marker 2", choices = ref.markers, selected = 1),
                                          width=2
                                      )
                                  )
                         )
                         ##=================================================================================================
                         
                         output[["t_3_3_1_img"]] <- renderImage(
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
                                         text(xco,yco, paste0("REF_POP-",i))
                                     })
                                     dev.off()
                                 }
                                 
                                 list(src=outfile)
                             }
                             
                             
                         }, deleteFile = T)
                         progress$inc(0.45, detail="REFERENCE PLOTTED")
                         
                         output[["t_3_3_2_img"]] <- renderImage(
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
                                     rows.to.p <- unlist(sapply(fcs.test.1.info[as.numeric(input[[paste0("pl1_",f.id,"_test_cl")]])],
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
                                         j <- order(cl.to.plot.test.1)[i]
                                         xco <- mean(fcs.test.1@exprs[c(fcs.test.1.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m1")]]))])
                                         yco <- mean(fcs.test.1@exprs[c(fcs.test.1.info[[i]][[2]]), c(as.integer(input[[paste0("pl1_",f.id,"_m2")]]))])
                                         text(xco,yco, names(cl.to.plot.test.1[j]))
                                     })
                                     dev.off()
                                 }
                                 
                                 list(src=outfile)
                             }
                         }, deleteFile = T)
                         progress$inc(0.45, detail="TEST FILE PLOTTED")
                         
                         if(!is.null(analysis.variables$FG.accuracy.coef.list.1[[f.name]]) && global.values$plot)
                         {
                             removeUI("#FG_score")
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
                             insertUI(selector="#t_3_2_2",
                                      where = "beforeEnd",
                                      ui = box(id="FG_score",h4(paste0("ACCURACY (FG) = ",analysis.variables$FG.accuracy.coef.list.1[[f.name]])),
                                               t)
                             )
                             
                             if(!is.null(analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]]))
                             {
                                 removeUI("#RFG_score")
                                 insertUI(selector="#t_3_2_2",
                                          where = "beforeEnd",
                                          ui = h5(id="RFG_score", paste0("Raw FG score = ",analysis.variables$raw.FG.accuracy.coef.list.1[[f.name]]))
                                 )
                             }
                             if(!is.null(analysis.variables$theta.coef.list.1[[f.name]]))
                             {
                                 removeUI("#theta_corr")
                                 insertUI(selector="#t_3_2_2",
                                          where = "beforeEnd",
                                          ui = h5(id="theta_corr",paste0("Theta correction = ",analysis.variables$theta.coef.list.1[[f.name]]))
                                 )
                             }
                         }
                         if(!is.null(analysis.variables$F.beta.coef.list.1[[f.name]]) && global.values$plot)
                         {
                             removeUI("#Fbeta_score")
                             insertUI(selector="#t_3_2_2",
                                      where = "beforeEnd",
                                      ui = h5(id="Fbeta_score", paste0("Fbeta score = ",mean(analysis.variables$F.beta.coef.list.1[[f.name]])))
                             )
                         }
                         if(!is.null(analysis.variables$G.coef.list.1[[f.name]]) && global.values$plot)
                         {
                             removeUI("#G_score")
                             insertUI(selector="#t_3_2_2",
                                      where = "beforeEnd",
                                      ui = h5(id="G_score", paste0("G score = ",mean(analysis.variables$G.coef.list.1[[f.name]])))
                             )
                         }
                         
                         
                         
                     }
                     progress$close()
                 }
                 
             })
             shinyjs::enable("save_results")
         }
         shinyjs::enable("compare")
     })
    
    
    
    
    
    observeEvent(input$run_clustering_button,
    {
         shinyjs::disable("save_results")
         if( !is.null(global.values$fcs.files.fg.ref) )
         {
             if(!is.na(input$method_selection_fg_r) && !is.null(input$method_selection_fg_r) && input$method_selection_fg_r != "")
             {
                 global.values$fcs.files.fg.proj.1 <<- NULL
                 global.values$use.enriched.file <<- TRUE
                 curr.method.name <- input$method_selection_fg_r
                 
                 lapply(1:length(global.values$fcs.files.fg.ref), function(fcs.file.id)
                 {
                     temp.dir <- paste0(global.values$working.directory,"/GeneratedFCS/")
                     params <- list()
                     if(length(clustering.variables$available.methods.parameters[[curr.method.name]]) > 0)
                     {
                         params <- lapply(1:length(clustering.variables$available.methods.parameters[[curr.method.name]]), function(p)
                         {
                             x <- NULL
                             if(!is.null(input[[paste0("cm_",curr.method.name,"_param_",p)]]))
                             {
                                 x <- input[[paste0("cm_",curr.method.name,"_param_",p)]]
                             }
                             return(x)
                         })
                         names(params) <- names(clustering.variables$available.methods.parameters[[curr.method.name]])
                     }
                     markers_col <- 0
                     if( !is.null(input$method_markers))
                     {
                         if( !is.na(input$method_markers))
                         {
                             markers_col <- colnames(global.values$fcs.files.fg.ref[[1]]@exprs)[as.numeric(input$method_markers)]
                         }
                     }
                     method.output <- benchmark.method(curr.method.name, global.values$fcs.files.fg.ref[[1]], temp.dir, params, markers_col)
                     
                     global.values$fcs.files.fg.proj.1 <<- c(global.values$fcs.files.fg.proj.1, method.output[[1]])
                     names(global.values$fcs.files.fg.proj.1)[[length(global.values$fcs.files.fg.proj.1)]] <<- 
                         names(global.values$fcs.files.fg.ref)[[length(global.values$fcs.files.fg.ref)]]
                     
                     temp.file.path <- names(global.values$fcs.files.fg.proj.1)[[length(global.values$fcs.files.fg.proj.1)]]
                     temp.file.path <- paste0(temp.file.path, "-_-",curr.method.name)
                     clustering.variables$added.keyword <- paste0(ncol(global.values$fcs.files.fg.ref[[1]]@exprs)+1,"__")
                     if(length(params > 0))
                     {
                         lapply(1:length(params), function(i)
                         {
                             temp.file.path <<- paste0(temp.file.path,names(params)[i],"-", params[[i]],"__")
                             clustering.variables$added.keyword <<- paste0(clustering.variables$added.keyword,names(params)[i],"-", params[[i]],"__")
                         })
                     }
                     temp.file.path <- paste0(temp.file.path,"____")
                     names(clustering.variables$added.keyword) <- paste0("CLMETH_",ncol(global.values$fcs.files.fg.proj.1[[1]]@exprs),"_",curr.method.name)
                     global.values$fcs.files.fg.proj.1[[1]] <- add.keyword.to.fcs(global.values$fcs.files.fg.proj.1[[1]], clustering.variables$added.keyword, names(clustering.variables$added.keyword))
                     global.values$fcs.files.fg.ref[[1]] <- add.keyword.to.fcs(global.values$fcs.files.fg.ref[[1]], clustering.variables$added.keyword, names(clustering.variables$added.keyword))
                     names(global.values$fcs.files.fg.proj.1)[length(global.values$fcs.files.fg.proj.1)] <<- basename(temp.file.path)
                     
                     fcs.temp <- global.values$fcs.files.fg.proj.1[[1]]
                     
                     fcs.temp.name <- paste0(basename(names(global.values$fcs.files.fg.proj.1)[[length(global.values$fcs.files.fg.proj.1)]]),".fcs")
                     fcs.temp.dir <- tempdir()
                     write.FCS.CIPHE(fcs.temp,paste0(fcs.temp.dir,"/",fcs.temp.name))
                     
                     global.values$fcs.files.fg.proj.1[[1]] <- read.FCS(paste0(fcs.temp.dir,"/",fcs.temp.name), emptyValue = F)
                     
                     removeUI("#download_enriched_file")
                     insertUI("#t_1_b_1",
                              "beforeEnd",
                              downloadButton("download_enriched_file", "Download Enriched File"))
                     
                     
                     output$download_enriched_file <- downloadHandler(
                         filename = function()
                         {
                             paste0(names(global.values$fcs.files.fg.proj.1)[1],".fcs")
                         },
                         content = function(file)
                         {
                             write.FCS(global.values$fcs.files.fg.proj.1[[1]], file, delimiter="#")
                         }
                     )
                     
                     updateTabItems(session,'tabs',"t_1")
                 })
             }
             
         }
     })
    
    observe(
    {
        if(!is.null(global.values$fcs.files.fg.out))
        {
            enable("save_results")
        }
        else
        {
            disable("save_results")
        }
    })
    
    output$save_results <- downloadHandler(
        filename = function()
        {
            paste0("RESULTS_",names(global.values$fcs.files.fg.ref)[1],".zip")
        },
        content = function(file)
        {
            view(description(global.values$fcs.files.fg.ref[[1]]))
            shinyjs::disable("save_results")
            progress <- Progress$new()
            progress$set(message="EXPORTING RESULTS", value=0)
            temp.dir <- ""
            fnames <- c(1,2)
            
            out.mat <- NULL
            if(!is.null(global.values$fcs.files.fg.mapping))
            {
                out.mat <- global.values$fcs.files.fg.mapping
            }
            else
            {
                mat <- matrix(1:ncol(analysis.variables$purity.matrix), ncol = 1)
                out.mat <- cbind(mat,mat)
                global.values$fcs.files.fg.mapping <- out.mat
                names(global.values$fcs.files.fg.mapping) <- "generated_map"
            }
            out.fcs <- global.values$fcs.files.fg.ref[[1]]
            
            ann.fbeta.scores <- analysis.variables$F.beta.coef.list.1[[1]]
            ordered.annot.associated.pop <- analysis.variables$annotations.relations.id.list.1[[1]][-1]
            ref.info <- FPH.get.file.information(global.values$fcs.files.fg.ref[[1]],as.numeric(input$clust_col_selection_fg_r))
            ref.list <- rep(NA,length(ref.info))
            for (i in 1:length(ref.list)) 
            {
                t <- which(ordered.annot.associated.pop==i)
                if(length(t)>0)
                {
                    if(length(t)>1)
                    {
                        alert("MORE REFERENCES THAN ANNOTATIONS")
                    }
                    ref.list[i] <- ann.fbeta.scores[t[[1]]]
                }
                
            }
            ref.list.name <- colnames(global.values$fcs.files.fg.proj.1[[1]])[as.numeric(input$clust_col_selection_fg_1)]
            real.pop.order <- c()
            for(i in 1:length(ref.info))
            {
                real.pop.order <- c(real.pop.order,
                                    unique(global.values$fcs.files.fg.ref[[1]]@exprs[ref.info[[i]][[2]],as.numeric(input$clust_col_selection_fg_r)])[[1]])
            }
            ref.list <- ref.list[real.pop.order]
            
            
            
            meth.col <- ncol(global.values$fcs.files.fg.proj.1[[1]]@exprs)
            meth.name <- colnames(global.values$fcs.files.fg.proj.1[[1]])[meth.col] 
            if(!is.null(clustering.variables$added.keyword))
            {
                meth.name <- strsplit(names(clustering.variables$added.keyword),"_")[[1]]
                meth.name <- meth.name[length(meth.name)]
                out.fcs@exprs <- cbind(out.fcs@exprs,global.values$fcs.files.fg.proj.1[[1]]@exprs[,meth.col])
                colnames(out.fcs@exprs)[ncol(out.fcs@exprs)] <- paste0(meth.name,"-",ncol(out.fcs@exprs),"_clusters")
            }
            else
            {
                meth.col <- as.numeric(input$clust_col_selection_fg_1)
                clustering.variables$added.keyword <- paste0(ncol(out.fcs@exprs)+1,"__",meth.name)
                names(clustering.variables$added.keyword) <- paste0("CLMETH_",ncol(out.fcs@exprs)+1,"_",meth.name)
                out.fcs@exprs <- cbind(out.fcs@exprs,global.values$fcs.files.fg.proj.1[[1]]@exprs[,meth.col])
                colnames(out.fcs@exprs)[ncol(out.fcs@exprs)] <- paste0(meth.name,"-",ncol(out.fcs@exprs),"_clusters")
            }
            
            if(!is.null(out.mat))
            {
                out.mat <- cbind(out.mat,ref.list)
                colnames(out.mat)[ncol(out.mat)] <- paste0(meth.name,".",ncol(out.fcs@exprs),".fgcol")
            }
            
            out.fcs <- add.keyword.to.fcs(out.fcs, trunc(mean(analysis.variables$F.beta.coef.list.1[[1]])*10^4)/10^4, paste0("FSCLMETH_",meth.col,"_",meth.name))
            out.fcs@exprs <- cbind(out.fcs@exprs, NA)
            colnames(out.fcs@exprs)[ncol(out.fcs@exprs)] <- paste0(meth.name,"-",ncol(out.fcs@exprs)-1,"_ANNOTATIONS")
            fcs.annotations <- FPH.get.file.information(global.values$fcs.files.fg.proj.1.annotated,meth.col+1)
            
            lapply(1:length(fcs.annotations), function(i) 
            {
                if(!is.na(ordered.annot.associated.pop[i]))
                {
                    out.fcs@exprs[fcs.annotations[[i]][[2]],ncol(out.fcs@exprs)] <<- rep(ordered.annot.associated.pop[i],fcs.annotations[[i]][[1]])
                }
            })
            
            
            t <- basename(names(global.values$fcs.files.fg.ref)[1])
            fnames[1] <- paste0(temp.dir, t,"__",meth.name,"_.fcs")
            write.FCS.CIPHE(out.fcs, fnames[1])
            fnames[2] <- paste0(temp.dir, basename(names(global.values$fcs.files.fg.mapping)),"__",meth.name,"_.csv")
            write.csv(out.mat, fnames[2], col.names = T, row.names = F)
            
            zip(zipfile=file,files=fnames)
            file.remove(fnames)
            progress$inc(1, detail="FILES EXPORTED")
            progress$close()
            shinyjs::enable("save_results")
        },
        contentType = "application/zip"
    )
    
    
    
    
    
    observeEvent(input$t_4_1_refresh, 
    {
        shinyjs::disable("t_4_1_refresh")
        if( !is.null(global.values$fcs.files.fg.mapping) && ncol(global.values$fcs.files.fg.mapping)>2 )
        {
            analysis.variables$scores.table <- global.values$fcs.files.fg.mapping
            meth.names <- colnames(global.values$fcs.files.fg.mapping)[c(-1,-2)]
            if("popID.Scaffold"%in%colnames(global.values$fcs.files.fg.mapping) &&
               "pop"%in%colnames(global.values$fcs.files.fg.mapping))
            {
                cols.to.add <- c()
                for(i in 1:ncol(global.values$fcs.files.fg.mapping))
                {
                    if( length(grep(".fgcol",colnames(global.values$fcs.files.fg.mapping)[i]))>0 )
                    {
                        cols.to.add <- c(cols.to.add,i)
                    }
                }
                print(cols.to.add)
                analysis.variables$scores.table <- cbind(global.values$fcs.files.fg.mapping[,c("popID.Scaffold","pop")],
                                                         global.values$fcs.files.fg.mapping[,cols.to.add])
                meth.names <- colnames(global.values$fcs.files.fg.mapping)[cols.to.add]
            }
            colnames(analysis.variables$scores.table) <- c("ID","POP",meth.names)
            analysis.variables$scores.table <- rbind(analysis.variables$scores.table,0)
                
            param.ncol <- max(sapply(clustering.variables$available.methods.parameters, function(k){return(length(k))}))+1
            analysis.variables$params.table <- matrix("",ncol=param.ncol)
                
            lapply(2:ncol(analysis.variables$scores.table), function(i)
            {
                
                meth.id <- strsplit(meth.names[i-1],".")[[1]][2]
                key.list <- get.keywords.with.keypart.FCS(global.values$fcs.files.fg.ref[[1]],paste0("CLMETH_",meth.id))
                if(length(key.list)>0)
                {
                    lapply(1:length(key.list), function(k)
                    {
                        if(strsplit(names(key.list[k]),"_")[[1]][1] != "FSCLMETH")
                        {
                            t <- strsplit(key.list[k],"__")[[1]]
                            if(length(t)>2)
                            {
                                col.id <- as.numeric(strsplit(names(key.list[k]),"_")[[1]][2])
                                method.name <- paste(strsplit(names(key.list[k]),"_")[[1]][3], strsplit(names(key.list[k]),"_")[[1]][2], sep="_")
                                if(nrow(analysis.variables$params.table)>1)
                                {
                                    analysis.variables$params.table <<- rbind(analysis.variables$params.table,0)
                                }
                                else
                                {
                                    analysis.variables$params.table <- matrix("",ncol=param.ncol)
                                }
                                analysis.variables$params.table[nrow(analysis.variables$params.table),1] <<- method.name
                                lapply(2:length(t), function(l)
                                {
                                    pa.name <- strsplit(t[l],"-")[[1]][1]
                                    pa <- strsplit(t[l],"-")[[1]][2]
                                    analysis.variables$params.table[nrow(analysis.variables$params.table),l] <<- paste0(pa.name,": ",pa)
                                })
                            }
                        }
                    })
                }
                val <- analysis.variables$scores.table[-nrow(analysis.variables$scores.table),i]
                analysis.variables$scores.table[nrow(analysis.variables$scores.table),i] <<- trunc(mean(as.numeric(val[!is.na(val)]))*10^4)/10^4
            })
            analysis.variables$scores.table[nrow(analysis.variables$scores.table),1] <- "GLOBAL (mean)"
            colnames(analysis.variables$params.table) <- c("method","parameters",rep(" ",-2+param.ncol))
        }
        shinyjs::enable("t_4_1_refresh")
    })
    
    observeEvent(input$t_4_1_generate,
    {
        shinyjs::disable("t_4_1_generate")
        analysis.variables$scores.table <- NULL
        if( !is.null(global.values$fcs.files.fg.ref) && keyword.exists.FCS(global.values$fcs.files.fg.ref[[1]],"FSCLMETH_") )
        {
            key.list <- get.keywords.with.keypart.FCS(global.values$fcs.files.fg.ref[[1]],"FSCLMETH_")
            analysis.variables$scores.table <- c("GLOBAL (mean)")
            lapply(1:length(key.list), function(k)
            {
                t <- as.numeric(key.list[k])
                col.id <- as.numeric(strsplit(names(key.list[k]),"_")[[1]][2])+1
                method.name <- paste(strsplit(names(key.list[k]),"_")[[1]][3], strsplit(names(key.list[k]),"_")[[1]][2], sep="_")
                analysis.variables$scores.table <- cbind(analysis.variables$scores.table, t)
                colnames(analysis.variables$scores.table)[ncol(analysis.variables$scores.table)] <- method.name
            })
            colnames(analysis.variables$scores.table)[1] <- "POP"
            
            param.ncol <- max(sapply(clustering.variables$available.methods.parameters, function(k){return(length(k))}))+1
            analysis.variables$params.table <- matrix("",ncol=param.ncol)
            key.list <- get.keywords.with.keypart.FCS(global.values$fcs.files.fg.ref[[1]],"CLMETH_")
            if(length(key.list)>0)
            {
                lapply(1:length(key.list), function(k)
                {
                    if(strsplit(names(key.list[k]),"_")[[1]][1] != "FSCLMETH")
                    {
                        t <- strsplit(key.list[k],"__")[[1]]
                        if(length(t)>2)
                        {
                            col.id <- as.numeric(strsplit(names(key.list[k]),"_")[[1]][2])
                            method.name <- paste(strsplit(names(key.list[k]),"_")[[1]][3], strsplit(names(key.list[k]),"_")[[1]][2], sep="_")
                            if(nrow(analysis.variables$params.table)>1)
                            {
                                analysis.variables$params.table <<- rbind(analysis.variables$params.table,0)
                            }
                            else
                            {
                                analysis.variables$params.table <- matrix("",ncol=param.ncol)
                            }
                            analysis.variables$params.table[nrow(analysis.variables$params.table),1] <<- method.name
                            lapply(2:length(t), function(l)
                            {
                                pa.name <- strsplit(t[l],"-")[[1]][1]
                                pa <- strsplit(t[l],"-")[[1]][2]
                                analysis.variables$params.table[nrow(analysis.variables$params.table),l] <<-paste0(pa.name,": ",pa)
                            })
                        }
                    }
                })
            }
            colnames(analysis.variables$params.table) <- c("method","parameters",rep(" ",-2+param.ncol))
            
            
        }
        shinyjs::enable("t_4_1_generate")
    })
    
    
    output$t_4_1_table <- renderTable(analysis.variables$scores.table)
    output$t_4_1_params <- renderTable(analysis.variables$params.table)
}