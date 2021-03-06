##------------------------------------------------------------------------------
##                        FUNCIONES DE LA MODELIZACI�N
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##                          LIBRER�AS NECESARIAS
list.of.packages <- c("ggplot2", "igraph", "ggraph", "shinyBS", "shinydashboard",
                      "sf", "rnaturalearth", "rnaturalearthdata", "tidygraph",
                      "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only = TRUE)


##------------------------------------------------------------------------------
##                    OBTENCI�N DE LA MATRIZ DE DISTANCIAS

data_network <- function(data.csv, nCoord) {
  
  #Cargamos los datos
  Coordenadas <- read.csv(data.csv, header = TRUE, sep = ";",
                              dec = ".", col.names = c("Y", "X"))
  
  #Los datos de la pen�nsula deben ordenarse siguiendo la l�nea costera
  coord_orden <- Coordenadas[sample(nrow(Coordenadas), nCoord),]
  costa_izq <- coord_orden[which(coord_orden$X<(-5)),]
  costa_izq <- costa_izq[order(-costa_izq$Y),]
  
  costa_dcha <- coord_orden[which(coord_orden$X>=(-5)),]
  costa_dcha <- costa_dcha[order(costa_dcha$Y),]
  
  coord_orden <- as.data.frame(rbind(costa_izq, costa_dcha))
  return(coord_orden)
}

distance_matrix <- function(data_network){
  
  #Valores iniciales de la matriz de distancias
  nCoord <- nrow(data_network)
  dist_matrix <- matrix(0, nCoord, nCoord)
  X = 0
  Y = 0
  
  #distancias eucl�deas y equivalencia entre grados - kilometros (1gd = 111.1 Km)
  for (i in 1:nCoord){
    for (j in 1:nCoord){
      X <- data_network$X[i] - data_network$X[j]
      Y <- data_network$Y[i] - data_network$Y[j]
      dist_matrix[i,j] <- sqrt(X^2 + Y^2)*111.1
    }
  }
  return(dist_matrix)
}


##------------------------------------------------------------------------------
##                   OBTENCI�N DE LA MATRIZ DE ADYACENCIA

adjacency_matrix <- function(dist_matrix, dist_thres){
  
  #Dimensiones y valores iniciales de la matriz de adyacencia
  nCoord <- nrow(dist_matrix)
  adj_matrix <- matrix(0, nCoord, nCoord)
 
  #Asignaci�n de conexiones en funci�n de la distancia umbral
  for (i in 1:nCoord){
    for (j in 1:nCoord){
      if (dist_matrix[i,j] <= dist_thres && dist_matrix[i,j] != 0){
        adj_matrix[i,j] = 1
      }
    }
  }
  return(adj_matrix)
}


##------------------------------------------------------------------------------
##                 CALCULO DE LOS PESOS EN LAS CONEXIONES

weight_adj_matrix <- function(dist_matrix, max_disp, imp_disp){
  
  #Generamos una matriz de pesos (mismas dimensiones que distancias)
  nCoord <- nrow(dist_matrix)
  weight_matrix <- matrix(0, nCoord, nCoord)
  
  #Calculamos cada uno de los pesos con una funci�n exponenical.
  for (i in 1:nCoord){
    for (j in 1:nCoord){
      if (dist_matrix[i,j] < imp_disp){
        weight_matrix[i,j] = exp(-dist_matrix[i,j]/max_disp)
      }
      else {
        weight_matrix[i,j] = 0
      }
    }
  }
  
  return(weight_matrix)
}


##------------------------------------------------------------------------------
##               FUNCI�N DE MODULARIDAD (Algoritmo de Louvain)

mod_value <- function(adj_matrix){
  
  #Generamos la gr�fica con una matriz de adyacencia indirecta
  graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
  
  #Implementamos el algoritmo de Louvain
  louv_mod <- cluster_louvain(graph = graph_adj)
  
  #Obtenemos el valor de modularidad
  mod <- round(modularity(graph_adj, membership = louv_mod$membership),3)
  
  return(mod)
}


##------------------------------------------------------------------------------
##                     FUNCI�N DE MODULARIDAD MEMBERSHIP

mod_membership <- function(adj_matrix){
  
  #Generamos la gr�fica con una matriz de adyacencia indirecta
  graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
  
  #Implementamos el algoritmo de Louvain
  louv_mod <- cluster_louvain(graph = graph_adj)
  
  return(louv_mod$membership)
}


##------------------------------------------------------------------------------
##                  FUNCI�N MATRIZ DE MODULARIDAD

mod_matrix <- function(adj_matrix, mod_membership){
  
  #Generamos una matriz con las mismas dimensiones que la de adyacencia
  cluster_mat <- adj_matrix
  
  #Asignamos los m�dulos a los que pertenece cada nodo en la matriz
  for (i in 1:nrow(adj_matrix)){
    for (j in 1:ncol(adj_matrix)){
      if (adj_matrix[i,j] == 1){
        if (mod_membership[i] == mod_membership[j]){
          cluster_mat[i,j] = mod_membership[i]
        }
        else {
          cluster_mat[i,j] = 1
        }
      }
      else {
        cluster_mat[i,j] = NA
      }
    }
  }
  
  return(cluster_mat)
}

##------------------------------------------------------------------------------
##                            REDES REALES
##------------------------------------------------------------------------------
##                            FUNCI�N LCC
LCC_real_network <- function(doc.csv, nrep, sample_size, thres_vector){
  
  #Definimos la matriz en la que almacenar los valores
  LCC_tama�os <- matrix(0, 
                        nrow = length(thres_vector), 
                        ncol = length(sample_size))
  colnames(LCC_tama�os) <- c("100", "500", "1000", "1500", "2000")
  
  #Almacenamos en la matriz cada valor de LCC para cada tama�o y threshold
  for (i in 1:length(sample_size)){
    prom_val <- matrix(0, nrow = length(thres_vector), ncol = nrep)
    dist_matrix <- distance_matrix(doc.csv, sample_size[i])
    for (k in 1:nrep){
      for (j in 1:length(thres_vector)){
        
        #Calculamos la matriz de adyacencia y los pesos de la matriz
        adj_matrix <- adjacent_matrix(dist_matrix, thres_vector[j])
        
        #Calculamos el tama�o de los clusters de la red
        cluster_mem <- mod_membership(adj_matrix)
        prom_val[j,k] <- max(table(cluster_mem))[[1]]
      }
    }
    LCC_tama�os[,i] = apply(prom_val, 1, mean, na.rm = TRUE)
  }
  
  #Generamos una matriz con threshold, tama�os y LCC para cada combinaci�n
  Resultados <- cbind(thres_vector, LCC_tama�os)
  return(Resultados)
}

##------------------------------------------------------------------------------
##                            FUNCI�N MODULARIDAD
Mod_real_network <- function(doc.csv, nrep, sample_size, thres_vector){
  
  #Definimos la matriz en la que almacenar los valores
  Mod_values <- matrix(0, 
                       nrow = length(thres_vector), 
                       ncol = length(sample_size))
  colnames(Mod_values) <- c("100", "500", "1000", "1500", "2000")
  
  #Almacenamos en la matriz la Modularidad para cada tama�o y threshold
  for (i in 1:length(sample_size)){
    prom_val <- matrix(0, nrow = length(thres_vector), ncol = nrep)
    dist_matrix <- distance_matrix(doc.csv, sample_size[i])
    for (k in 1:nrep){
      for (j in 1:length(thres_vector)){
        
        #Calculamos la matriz de adyacencia y los pesos de la matriz
        adj_matrix <- adjacent_matrix(dist_matrix, thres_vector[j])
        
        #Generamos la gr�fica con una matriz de adyacencia indirecta
        graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
        
        #Implementamos el algoritmo de Louvain
        louv_mod <- cluster_louvain(graph = graph_adj)
        
        #Obtenemos el valor de modularidad
        prom_val[j,k] <- round(modularity(graph_adj, 
                                          membership = louv_mod$membership),
                               3)
      }
    }
    Mod_values[,i] = apply(prom_val, 1, mean, na.rm = TRUE)
  }
  
  #Generamos una matriz con threshold, tama�os la modularidad
  #para cada combinaci�n
  Resultados <- cbind(thres_vector, Mod_values)
  return(Resultados)
}

##------------------------------------------------------------------------------
##                      FUNCI�N COEFICIENTE DE CLUSTERING
CCoef_real_network <- function(doc.csv, nrep, sample_size, thres_vector){
  
  #Definimos la matriz en la que almacenar los valores
  CCoef_values <- matrix(0, 
                       nrow = length(thres_vector), 
                       ncol = length(sample_size))
  colnames(CCoef_values) <- c("100", "500", "1000", "1500", "2000")
  
  #Almacenamos en la matriz el coef. clustering para cada tama�o y threshold
  for (i in 1:length(sample_size)){
    prom_val <- matrix(0, nrow = length(thres_vector), ncol = nrep)
    dist_matrix <- distance_matrix(archivo.csv, sample_size[i])
    for (k in 1:nrep){
      for (j in 1:length(thres_vector)){
        
        #Calculamos la matriz de adyacencia y los pesos de la matriz
        adj_matrix <- adjacent_matrix(dist_matrix, thres_vector[j])
        
        #Generamos la gr�fica con una matriz de adyacencia indirecta
        graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
        
        #Calculamos el coeficiente de clustering para cada combinaci�n
        prom_val[j,k] <- round(transitivity(graph_adj, type = "average"),
                               3)
      }
    }
    CCoef_values[,i] = apply(prom_val, 1, mean, na.rm = TRUE)
  }
  
  #Generamos una matriz con threshold, tama�os el coef. clustering
  #para cada combinaci�n
  Resultados <- cbind(thres_vector, CCoef_values)
  return(Resultados)
}

##------------------------------------------------------------------------------
##                             FUNCI�N DI�METRO
Diameter_real_network <- function(archivo.csv, nrep, sample_size, thres_vector){
  
  #Definimos la matriz en la que almacenar los valores
  Diam_values <- matrix(0, 
                         nrow = length(thres_vector), 
                         ncol = length(sample_size))
  colnames(Diam_values) <- c("100", "500", "1000", "1500", "2000")
  
  #Almacenamos en la matriz el di�metro para cada tama�o y threshold
  for (i in 1:length(sample_size)){
    prom_val <- matrix(0, nrow = length(thres_vector), ncol = nrep)
    dist_matrix <- distance_matrix(archivo.csv, sample_size[i])
    for (k in 1:nrep){
      for (j in 1:length(thres_vector)){
        
        #Calculamos la matriz de adyacencia y los pesos de la matriz
        adj_matrix <- adjacent_matrix(dist_matrix, thres_vector[j])
        
        #Generamos la gr�fica con una matriz de adyacencia indirecta
        graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
        
        #Calculamos el di�metro para cada combinaci�n
        prom_val[j,k] <- diameter(graph_adj)
      }
    }
    Diam_values[,i] = apply(prom_val, 1, mean, na.rm = TRUE)
  }
  
  #Generamos una matriz con threshold, tama�os el di�metro para cada combinaci�n
  Resultados <- cbind(thres_vector, Diam_values)
  return(Resultados)
}

##------------------------------------------------------------------------------
##                             FUNCI�N EDGE DENSITY

EdgeDensity_real_network <- function(archivo.csv, nrep, sample_size, thres_vector){
  
  #Definimos la matriz en la que almacenar los valores
  Diam_values <- matrix(0, 
                        nrow = length(thres_vector), 
                        ncol = length(sample_size))
  colnames(Diam_values) <- c("100", "500", "1000", "1500", "2000")
  
  #Almacenamos en la matriz la densidad de ejes para cada tama�o y threshold
  for (i in 1:length(sample_size)){
    prom_val <- matrix(0, nrow = length(thres_vector), ncol = nrep)
    dist_matrix <- distance_matrix(archivo.csv, sample_size[i])
    for (k in 1:nrep){
      for (j in 1:length(thres_vector)){
        
        #Calculamos la matriz de adyacencia y los pesos de la matriz
        adj_matrix <- adjacent_matrix(dist_matrix, thres_vector[j])
        
        #Generamos la gr�fica con una matriz de adyacencia indirecta
        graph_adj <- graph_from_adjacency_matrix(adj_matrix, "Undirected")
        
        #Calculamos la densidad de ejes para cada combinaci�n
        prom_val[j,k] <- edge_density(graph_adj)
      }
    }
    Diam_values[,i] = apply(prom_val, 1, mean, na.rm = TRUE)
  }
  
  #Generamos una matriz con threshold, tama�os la densidad para cada combinaci�n
  Resultados <- cbind(thres_vector, Diam_values)
  return(Resultados)
}


