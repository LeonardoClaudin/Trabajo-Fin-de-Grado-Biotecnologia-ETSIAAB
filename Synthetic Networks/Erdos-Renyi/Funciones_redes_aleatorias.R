##------------------------------------------------------------------------------
##                        FUNCIONES DE LA MODELIZACIÓN
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##                          LIBRERÍAS NECESARIAS
list.of.packages <- c("ggplot2", "igraph", "ggraph", "shinyBS", "shinydashboard",
                      "sf", "rnaturalearth", "rnaturalearthdata", "tidygraph",
                      "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only = TRUE)

##------------------------------------------------------------------------------
##                        REDES ALEATORIAS SINTÉTICAS
##------------------------------------------------------------------------------
##                      FUNCIÓN LCC ERDOS RENYI NETWORK

LCC_erdos_renyi <- function(size, prob, prom){
  lsize <- length(size)
  lprob <- length(prob)
  Result <- data.frame("Size" = c(),
                       "LCC" = c(),
                       "Prob" = c())
  for (k in 1:lprob){
    # EL tamaño de la red es de 10.000 nodos.
    net <- erdos.renyi.game(5000, prob[k])
    Resultados <- data.frame("Size" = size,
                             "LCC" = vector(length = lsize),
                             "Prob" = rep(prob[k], lsize))
    res <- matrix(0, nrow = lsize, ncol = prom)
    for (i in 1:lsize){
      for (j in 1:prom){
        #Tomamos muestras de tamaño definido por un vector.
        vertex_sel <- sample(5000, size[i])
        net_subgraph <- induced_subgraph(net, vertex_sel)
        res[i,j] <- max(clusters(net_subgraph)$csize)
      }
      Resultados$LCC <- apply(res, 1, mean, na.rm = TRUE)
    }
    Result <- rbind(Result, Resultados)
  }
  return(Result)
}

##------------------------------------------------------------------------------
##                FUNCIÓN MODULARIDAD ERDOS RENYI NETWORK

Mod_erdos_renyi <- function(size, prob, prom){
  lsize <- length(size)
  lprob <- length(prob)
  Result <- data.frame("Size" = c(),
                       "Mod" = c(),
                       "Prob" = c())
  for (k in 1:lprob){
    net <- erdos.renyi.game(5000, prob[k])
    Resultados <- data.frame("Size" = size,
                             "Mod" = vector(length = lsize),
                             "Prob" = rep(prob[k], lsize))
    res <- matrix(0, nrow = lsize, ncol = prom)
    for (i in 1:lsize){
      for (j in 1:prom){
        vertex_sel <- sample(5000, size[i])
        net_subgraph <- induced_subgraph(net, vertex_sel)
        res[i,j] <- modularity(net_subgraph, 
                               membership = clusters(net_subgraph)$membership)
      }
      Resultados$Mod <- apply(res, 1, mean, na.rm = TRUE)
    }
    Result <- rbind(Result, Resultados)
  }
  return(Result)
}

##------------------------------------------------------------------------------
##           FUNCIÓN COEFICIENTE CLUSTERING ERDOS RENYI NETWORK

CCoef_erdos_renyi <- function(size, prob, prom){
  lsize <- length(size)
  lprob <- length(prob)
  Result <- data.frame("Size" = c(),
                       "CCoef" = c(),
                       "Prob" = c())
  for (k in 1:lprob){
    net <- erdos.renyi.game(5000, prob[k])
    Resultados <- data.frame("Size" = size,
                             "CCoef" = vector(length = lsize),
                             "Prob" = rep(prob[k], lsize))
    res <- matrix(0, nrow = lsize, ncol = prom)
    for (i in 1:lsize){
      for (j in 1:prom){
        vertex_sel <- sample(5000, size[i])
        net_subgraph <- induced_subgraph(net, vertex_sel)
        res[i,j] <- transitivity(net_subgraph)
      }
      Resultados$CCoef <- apply(res, 1, mean, na.rm = TRUE)
    }
    Result <- rbind(Result, Resultados)
  }
  return(Result)
}

##------------------------------------------------------------------------------
##               FUNCIÓN DIÁMETRO ERDOS RENYI NETWORK

Dist_erdos_renyi <- function(size, prob, prom){
  lsize <- length(size)
  lprob <- length(prob)
  Result <- data.frame("Size" = c(),
                       "Dist" = c(),
                       "Prob" = c())
  for (k in 1:lprob){
    net <- erdos.renyi.game(5000, prob[k])
    Resultados <- data.frame("Size" = size,
                             "Dist" = vector(length = lsize),
                             "Prob" = rep(prob[k], lsize))
    res <- matrix(0, nrow = lsize, ncol = prom)
    for (i in 1:lsize){
      for (j in 1:prom){
        vertex_sel <- sample(5000, size[i])
        net_subgraph <- induced_subgraph(net, vertex_sel)
        res[i,j] <- diameter(net_subgraph)
      }
      Resultados$Dist <- apply(res, 1, mean, na.rm = TRUE)
    }
    Result <- rbind(Result, Resultados)
  }
  return(Result)
}

##------------------------------------------------------------------------------
##             FUNCIÓN EDGE DENSITY ERDOS RENYI NETWORK

EdgeD_erdos_renyi <- function(size, prob, prom){
  lsize <- length(size)
  lprob <- length(prob)
  Result <- data.frame("Size" = c(),
                       "EdgeD" = c(),
                       "Prob" = c())
  for (k in 1:lprob){
    net <- erdos.renyi.game(5000, prob[k])
    Resultados <- data.frame("Size" = size,
                             "EdgeD" = vector(length = lsize),
                             "Prob" = rep(prob[k], lsize))
    res <- matrix(0, nrow = lsize, ncol = prom)
    for (i in 1:lsize){
      for (j in 1:prom){
        vertex_sel <- sample(5000, size[i])
        net_subgraph <- induced_subgraph(net, vertex_sel)
        res[i,j] <- edge_density(net_subgraph)
      }
      Resultados$EdgeD <- apply(res, 1, mean, na.rm = TRUE)
    }
    Result <- rbind(Result, Resultados)
  }
  return(Result)
}
