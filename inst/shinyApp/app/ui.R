library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
    
    dashboardHeader
    (
        title="Clustering Methods Comparision"
    ),
    
    dashboardSidebar
    (
        sidebarMenu
        (
            id="tabs",
            menuItem("Files Selection", tabName="t_1"),
            menuItem("Clustering", tabName="t_2"),
            menuItem("Method Evaluation", tabName="t_3"),
            menuItem("Methods Comparision", tabName="t_4")
        )
    ),
    
    dashboardBody
    (
        useShinyjs(),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
        tabItems
        (
            tabItem
            (
                tabName="t_1",
                h2("Files Selection"),
                fluidRow
                (
                    box
                    (
                        title="Reference File", height = "28vh", id="t_1_b_1",
                        actionButton("ref_selection_fg", "Reference File"),
                        actionButton("ref_mapping_fg", "Populations Mapping File"),
                        #actionButton("ref_update_fg", "Update Clusters List"),
                        checkboxInput("ref_use_same_fg", "Use As Test",value = F),
                        selectInput("clust_col_selection_fg_r", "Clusters column", choices=NULL)
                    ),
                    
                    box
                    (
                        title="File 1", height = "28vh", id="t_1_b_2",
                        actionButton("proj1_selection_fg", "File 1"),
                        #actionButton("proj1_update_fg", "Update Clusters List"),
                        selectInput("clust_col_selection_fg_1", "Clusters column", choices=NULL)
                    )
                ),
                div(id="FG_description",
                    HTML("<p>
                                The <b>Reference file</b> is a usual FCS - or csv - file.<br />
                                <b>File 1</b> can be chosen in 3 different ways<br />
                                <BLOCKQUOTE>
                                    1. Selecting a file with the same event as the <b>Reference</b> to which the clusters labels from a clustering method were added.<br />
                                    2. Using the <b>Reference</b> file and choosing another column<br />
                                    3. Using the <b>Reference</b> file and running a clustering algorithm within the app.
                                </BLOCKQUOTE>
                                <b>The Mapping File</b> is a csv file containing the Ids attributed to each population (e.g NKt --> 1). It can be obtained from Scaffold.<br />
                           </p>"),
                    class="description_div"
                )
                
                
            ),
            tabItem
            (
                tabName="t_2",
                h2("Clustering methods"),
                fluidRow
                (
                    id="t_2_fr",
                    selectInput("method_selection_fg_r", "Clustering method", choices=NULL),
                    selectInput("method_markers", "Markers to use", choices=NULL, multiple = TRUE),
                    shinyjs::disabled
                    (
                        actionButton("run_clustering_button", "Run Clustering")
                    )
                )
            ),
            tabItem
            (
                tabName="t_3",
                h2("Efficiency scores"),
                fluidRow
                (
                    id="t_3_1",
                    div(id="FG_description",
                        HTML("<p>
                                Comparison steps:
                                <BLOCKQUOTE>
                                    1. Depending on the selected <u>Clusters column</u>, either the clusters or populations from the <b>Reference File</B> are extracted and used as landmarks.<br />
                                    2. Each cluster from <b>File 1</b> is annotated, using purity and the FG score to compare them to the landmarks from step 1.<br />
                                    3. The clusters from <b>File 1</b> are merged into <u>Annotated populations</u> based on their... annotations !<br />
                                    4. The <u>Reference clusters/populations</u> and the <u>Annotated populations</u> are plotted. The markers used to plots can be changed.<br />
                                    5. The FG score is printed and stored in the description of the file. <br />
                                    <strong>For methods using random numbers at some points (kmeans, clara, ...), it is recommanded to use <u>Reference populations</u> instead of <u>Reference clusters.</u></strong><br />
                                </BLOCKQUOTE>
                           </p>
                           <p>
                                FG Score: mean of F-beta score and G score, both giving the accuracy of the annotations attributed to the clusters from
                                the tested files, compared to the labels from the reference files.
                                <BLOCKQUOTE>
                                        -  F-Beta score = harmonic mean of precision and recall.<br />
                                        -  G score = geometric mean of precision and recall.<br />
                                        -  The <u>raw FG score</u> is the FG score computed as is - no correction applied to the score<br />
                                        <strong>A score of at least 0.9 indicates a good clustering method giving a proper annotation of most clusters from <b>File 1</b></strong><br />
                                </BLOCKQUOTE>
                           </p>"),
                        actionButton("compare", "Compare"),
                        downloadButton("save_results", "Export results"),
                        class="description_div"
                    )
                ),
                fluidRow
                (
                    id="t_3_2"
                ),
                fluidRow
                (
                    id="t_3_3"
                )
            ),
            tabItem
            (
                tabName="t_4",
                h2("F1 score based Methods Comparison"),
                box
                (
                    id="t_4_1",
                    actionButton("t_4_1_refresh", "Print scores from mapping file"),
                    actionButton("t_4_1_generate", "Print scores from previous analyses (keywords)"),
                    tableOutput("t_4_1_table"),
                    br(),br(),
                    tableOutput(outputId="t_4_1_params")
                )
            )
            
        )
        
    )
    
)