library(shiny)
library(shinyBS)
library(shinydashboard)
library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggpubr)
Coordinates_chum<- read.csv("Chumberas_Coords2.csv", header = TRUE, sep = ";",
                       dec = ".", col.names = c("Y", "X"))
Coordinates_alc<- read.csv("Coordenadas_alcornoques_def.csv", header = TRUE, sep = ";", dec = ".", 
                       col.names = c("Y", "X"))
Coordinates_vid<- read.csv("CoordenadasVid.csv", header = TRUE, sep = ";", dec = ".", 
                       col.names = c("Y", "X"))
Coordinates_oliv <- read.csv("Coordenadas_olivos.csv", header = TRUE, sep = ";",
                             dec = ".", col.names = c("Y", "X"))
world <- ne_countries(scale = "medium", returnclass = "sf")
# Define UI for application that draws a histogram
ui <- fluidPage(
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h1(strong('CONNECTIVITY BETWEEN NODES')),
            tags$hr(),
            selectInput(inputId = "data",
                        label = "Select which database to work with",
                        choices = c("Chumberas",
                                    "Alcornoques",
                                    "Vinedos",
                                    "Olivos")),
            sliderInput(inputId ="Ncoords", 
                        label="Choose number of Coordinates",
                        value=100, min=1, max=2000),
            tags$hr(),
            sliderInput(inputId ="Thres",
                        label ="Choose Distance Threshold (km)",
                        value=40, min=0, max=200),
            
            tags$hr(),
            plotOutput("plot1")
        ),  
        
        # Show a plot of the generated distribution
        mainPanel(
            tags$hr(),
            selectInput(inputId = "graph",
                        label = "Select which graph to plot",
                        choices = c("Degree Connection Distribution",
                                    "Bilogaritmic regresion of Degree",
                                    "Accumulated distribution of Degree",
                                    "Network Parameters v Threshold")),
            plotOutput("plotgraph"),
            tags$hr(),
            column(6,
                   h2(strong('Clustering GRAPH')),
                   tags$head(tags$style("#modal1 .modal-dialog{ width:1100px}")),
                   tags$head(tags$style("#modal1 .modal-body{ min-height:500px}")),
                   bsModal("modal1", "Clustering GRAPH", "display1", size = "large", 
                           column(6,
                                  h2("Clustering Graph"),
                                  plotOutput("plot2"), 
                                  verbatimTextOutput("text1")), 
                           column(6,
                                  h2("Modularity Graph"),
                                  plotOutput("plot3"),
                                  verbatimTextOutput("text2"))
                   ),
                   actionButton("display1", "Open Clustering Graphs")
            ),
            column(6,
                   h2(strong('Connection MAP')),
                   tags$head(tags$style("#modal2 .modal-dialog{ width:1100px}")),
                   tags$head(tags$style("#modal2 .modal-body{ min-height:1000px}")),
                   bsModal("modal2", "Connection MAP", "display", size = "large",
                           imageOutput("myimage")),
                   actionButton("display", "Open Connection MAP")
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    ##VARIABLES REACTIVAS.
    datasetInput <- reactive({
        Sample <- switch(input$data,
                        "Chumberas" = Coordinates_chum,
                        "Alcornoques" = Coordinates_alc,
                        "Vinedos" = Coordinates_vid,
                        "Olivos" = Coordinates_oliv)
        Workdata <- Sample[sample(nrow(Sample), input$Ncoords),]
        return(Workdata)
    })
    datasetInput3 <- reactive({
        Workdata <- datasetInput()
        FirstValues <- Workdata[which(Workdata$X<(-5)),]
        v1 <- FirstValues[order(-FirstValues$Y),]
        SecondValues <- Workdata[which(Workdata$X>(-5)),]
        v2 <- SecondValues[order(SecondValues$Y),]
        Workdata <- rbind(v1,v2)
        Workdata <- as.data.frame(Workdata)
        return(Workdata)
    })
    datasetInput5 <- reactive({
        Workdata <- datasetInput3()
        Distancias <- matrix(0, nrow(Workdata), nrow(Workdata))
        X=0
        Y=0
        for (i in 1:nrow(Workdata)){
            for (j in 1:nrow(Workdata)){
                X <- Workdata$X[i]-Workdata$X[j]
                Y <- Workdata$Y[i]-Workdata$Y[j]
                Distancias[i,j]<- sqrt((X)^2+(Y)^2)*111
            }
        }
        return(Distancias)
    })
    datasetInput6 <- reactive({
        Distancias <- datasetInput5()
        Conexiones <- (matrix(0, nrow(Distancias), nrow(Distancias)))
        for (i in 1:nrow(Distancias)){
            for (j in 1:nrow(Distancias)){
                if (Distancias[i,j]<(input$Thres)&&Distancias[i,j]!=0){
                    Conexiones[i,j]<- 1
                }
            }
        }
        return(Conexiones)
    })
    datasetInput7 <- reactive({
        Conexiones <- datasetInput6()
        Grado <- matrix(0, nrow(Conexiones), ncol=1)
        for (i in 1:nrow(Conexiones)){
            for (j in 1:nrow(Conexiones)){
                if  (Conexiones[i,j]==1){
                    Grado[i]<- Grado[i]+1
                }
            }
        }
        Distrib <- as.data.frame(table(Grado))
        return(Distrib)
    })
    datasetInput8 <- reactive({
        Workdata <- datasetInput3()
        nodes <- matrix(0, nrow(Workdata), 3)
        colnames(nodes) <- c("id", "X", "Y")
        for (i in 1:nrow(Workdata)){
            nodes[i,1] <- i
            nodes[i,2] <- Workdata$X[i]
            nodes[i,3] <- Workdata$Y[i]
        }
        nodes <- as.data.frame(nodes)
        return(nodes)
    })
    datasetInput9 <- reactive({
        Distancias <- datasetInput5()
        Conexiones <- datasetInput6()
            a=0
            for (i in 1:nrow(Conexiones)){
                for (j in 1:nrow(Conexiones)){
                    if (Conexiones[i,j]==1){
                        a<-a+1
                    }
                }
            }
            edges <- matrix(0, a, 3)
            colnames(edges) <- c("from", "to", "weight")
            means <- mean(Distancias)
            k=1
            for (i in 1:nrow(Distancias)){
                for (j in 1:nrow(Distancias)){
                    if (Conexiones[i,j]==1){
                        edges[k,1] <- i
                        edges[k,2] <- j
                        edges[k,3] <- means/Distancias[i,j]
                        k <- k + 1
                    }
                }
            }
            edges <- as.data.frame(edges)
        return(edges)
    })
    datasetInput10 <- reactive({
        Distrib <- datasetInput7()
        Log_dis <- matrix(0, nrow(Distrib), 2)
        for (i in 1:nrow(Distrib)){
            Log_dis[i,1] <- log(as.numeric(Distrib[i,1]))
            Log_dis[i,2] <- log(as.numeric(Distrib[i,2]))
        }
        colnames(Log_dis) <- c("Grado", "Freq")
        Log_dis <- as.data.frame(Log_dis)
        return(Log_dis)
    })
    datasetInput11 <- reactive({
        Log_dis <- datasetInput10()
        m <- lm(Freq~Grado, Log_dis)
        eq <- substitute(y == b %.% x+a*","~~italic(R)^2~"="~r2, 
                         list(a = format(unname(coef(m)[1]), digits = 2),
                              b = format(unname(coef(m)[2]), digits = 2),
                              r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq))
        return(as.character(as.expression(eq)))
    })
    datasetInput12 <- reactive({
        Workdata <- datasetInput3()
        Distrib <- datasetInput7()
        Norm_Distrib <- matrix(0, nrow(Distrib), ncol(Distrib))
        for (i in 1:nrow(Distrib)){
            Norm_Distrib[i,1] <- Distrib[i,1]
            Norm_Distrib[i,2] <- Distrib[i,2]/nrow(Workdata)
        }
        colnames(Norm_Distrib) <- c("Grado", "Freq")
        Norm_Distrib <- as.data.frame(Norm_Distrib)
        return(Norm_Distrib)
    })
    datasetInput13 <- reactive({
        Conexiones <- datasetInput6()
        g.graph <- graph_from_adjacency_matrix(Conexiones, "undirected")
        Cluster <- cluster_louvain(g.graph)$membership
        Clustering <- Conexiones
        for (i in 1:nrow(Conexiones)){
            for (j in 1:nrow(Conexiones)){
                if (Conexiones[i,j] == 1){
                    if (Cluster[i] == Cluster[j]){
                        Clustering[i,j] = Cluster[i]
                    }
                    else {
                        Clustering[i,j] = 0
                    }
                }
                else if (Conexiones[i,j] == 0){
                    Clustering[i,j] = NA
                }
            }
        }
        return(Clustering)
    })
    Value <- reactiveValues()
    Value$Con  <- data.frame(X = numeric(), 
                            Y = numeric())
    Value$Clus <- data.frame(X = numeric(), 
                             Y = numeric())
    Value$nClus <- data.frame(X = numeric(), 
                              Y = numeric())
    observeEvent(input$Thres, {
        ##Anadimos la ultima fila a la tabla de las conexiones.
        Clustering <- datasetInput13()
        Tabla <- table(Clustering)
        lcon <- length(Tabla[names(Tabla)==0])
        if (lcon == 0){
            ncon = 0
        } else {
            ncon = Tabla[names(Tabla)==0][[1]]
        }
        new_row <- data.frame(X = input$Thres,
                              Y = ncon)
        Value$Con <- rbind(Value$Con,new_row)
        
        ##Anadimos nueva fila a la tabla del coef.clustering.
        Conexiones <- datasetInput6()
        g.graph <- graph_from_adjacency_matrix(Conexiones, "undirected")
        com <- cluster_louvain(g.graph)
        mod <- round(modularity(g.graph, membership = com$membership),4)
        new_row2 <- data.frame(X = as.numeric(input$Thres),
                               Y = as.numeric(mod))
        Value$Clus <- rbind(Value$Clus, new_row2)
        
        ##Anadimos filas al numero de clusters.
        new_row3 <- data.frame(X = input$Thres,
                               Y = length(unique(com$membership)))
        Value$nClus <- rbind(Value$nClus, new_row3)
    })
    observeEvent(input$Ncoords, {
        Value$Con <- Value$Con[0,]
        Value$Clus <- Value$Clus[0,]
        Value$nClus <- Value$nClus[0,]
    })
    observeEvent(input$data, {
        Value$Con <- Value$Con[0,]
        Value$Clus <- Value$Clus[0,]
        Value$nClus <- Value$nClus[0,]
    })
    
    ##--------------------------------------------------------------------------
    ##                              PLOTS DEL PANEL.
    ##--------------------------------------------------------------------------
    output$plot1 <- renderPlot({
        Workdata <- datasetInput3()
        ggplot() + 
            geom_sf(data = world) +
            coord_sf(ylim = c(30, 50), xlim = c(-20, 10), expand = TRUE) +
            geom_point(data = Workdata, aes(x = X, y = Y), 
                       colour = "darkolivegreen",
                       position = "jitter")
    })
    output$plot2 <- renderPlot({
        Conexiones <- datasetInput6()
        net = graph_from_adjacency_matrix(Conexiones, "undirected")
        com <- cluster_louvain(net)
        group <- com$membership
        l = layout.fruchterman.reingold(net)
        
        plot(net, vertex.label = NA, vertex.color = group, vertex.size = 4,
             edge.color = "skyblue", layout = l)
    })
    output$text1 <- renderText({
        nodes <- datasetInput8()
        edges <- datasetInput9()
        routes_igraph <- graph_from_data_frame(d = edges,
                                               vertices = nodes,
                                               directed = TRUE)
        value <- round(transitivity(routes_igraph, type = "average"),3)
        diam <- round(diameter(routes_igraph, directed = TRUE), 3)
        first <- paste("The clustering coefficient is", value) 
        second <- paste("The diameter of the clusters is", diam, "km")
        paste(first, second, sep = "\n")
    })
    output$plot3 <- renderPlot({
        Workdata <- datasetInput3()
        Conexiones <- datasetInput6()
        g.graph <- graph_from_adjacency_matrix(Conexiones, "undirected")
        com <- cluster_louvain(g.graph)
        Nodos <- as.data.frame(cbind(Workdata[,2:3], 
                                     modulo = as.factor(as.vector(com$memberships))))
        ggplot() + 
            geom_sf(data = world) +
            coord_sf(ylim = c(30, 50), xlim = c(-20, 10), expand = TRUE) +
            geom_point(data = Nodos, aes(x = X, y = Y, 
                                         group = modulo,
                                         color = modulo), 
                       position = "jitter")
    })
    output$text2 <- renderText({
        Conexiones <- datasetInput6()
        g.graph <- graph_from_adjacency_matrix(Conexiones, "undirected")
        com <- cluster_louvain(g.graph)
        mod <- round(modularity(g.graph, membership = com$membership),4)
        paste("The value of modularity is", mod)
    })
    output$myimage <- renderImage({
        Clustering <- datasetInput13()
        colfunc <- colorRampPalette(rainbow(15, start = 0.5))
        png(filename="Clustering.png", width = 1000, height = 1000)
        image(t(apply(Clustering, 2, rev)), axes = FALSE, useRaster = TRUE,
              col = colfunc(length(unique(as.vector(Clustering)))))
        dev.off()
        return(list(
            src = "Clustering.png",
            filetype = "png",
            alt = "Clustering MAP"))
    }, deleteFile = TRUE)
    Dist_grado <- reactive({
        Distrib <- datasetInput7()
        ggplot() + geom_bar(data=Distrib,
                            aes(x = Grado, y = Freq),
                            stat = "identity")
    })
    Bilog_grado <- reactive({
        Log_dis <- datasetInput10()
        reg <- datasetInput11()
        ggplot(data = Log_dis, aes(x=Grado, y=Freq)) +
            geom_smooth(method = "lm", se=FALSE,
                        color="black",
                        formula = y ~ x) +
            geom_point(stat = "identity") +
            geom_text(x = mean(Log_dis$Grado)/2, y = mean(Log_dis$Freq)/2,
                      label = reg , parse = TRUE)
    })
    Acc_grado <- reactive({
        Norm_Distrib <- datasetInput12()
        ggplot() +
            geom_line(data = Norm_Distrib,
                      aes(x=Grado, y=cumsum(Freq), group = 1),
                      color = "navy") +
            geom_line(data = Norm_Distrib,
                      aes(x=Grado, y=(1-cumsum(Freq)), group = 1),
                      color = "firebrick") +
            xlab("Grado Conexiones") + ylab("Ocurrencias") +
            scale_fill_manual(name = "Distribution",
                              labels=c("Accumulated", "Inverse_Accumulated"))
    })
    Conx_Thres <- reactive({
        Conn <- ggplot(Value$Con, aes(x = X, y = Y)) + geom_point() + geom_line() +
                xlab("Distance Threshold") + ylab("Number of Connections")
        Clust <- ggplot(Value$Clus, aes(x = X, y = Y)) + geom_point() + geom_line() +
                xlab("Distance Threshold") + ylab("Clustering Coefficient")
        nClust <- ggplot(Value$nClus, aes(x = X, y = Y)) + geom_point() + geom_line() +
            xlab("Distance Threshold") + ylab("Number of Clusters")
        
        log_Conn <- ggplot(Value$Con, aes(x = log(X), y = log(Y))) + geom_point() + 
            geom_line() + xlab("Distance Threshold") + ylab("Number of Connections")
        log_Clust <- ggplot(Value$Clus, aes(x = X, y = log(Y))) + geom_point() + 
            geom_line() + xlab("Distance Threshold") + ylab("Clustering Coefficient")
        log_nClust <- ggplot(Value$nClus, aes(x = log(X), y = log(Y))) + geom_point() + 
            geom_line() + xlab("Distance Threshold") + ylab("Number of Clusters")
        ggarrange(Conn, Clust, nClust, log_Conn, log_Clust, log_nClust,
                  labels = c("Threshold - Connections",
                             "Threshold - Clustering",
                             "Threshold - N Clusters",
                             "log(Threshold) - log(Connections)",
                             "Threshold - log(Clustering)",
                             "log(Threshold) - log(N Clusters)"),
                  font.label = list(size = 10, color = "darkred", face = "bold"),
                  ncol = 3, nrow = 2)
    })
    output$plotgraph <- renderPlot({
        switch (input$graph,
                "Degree Connection Distribution" = Dist_grado(),
                "Bilogaritmic regresion of Degree" = Bilog_grado(),
                "Accumulated distribution of Degree" = Acc_grado(),
                "Network Parameters v Threshold" = Conx_Thres()
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
