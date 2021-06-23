# Trabajo-Fin-de-Grado-Biotecnologia-ETSIAAB
MODELIZACION DE DATOS DE GEOLOCALIZACION DE ESPECIES VEGETALES MEDIANTE REDES COMPLEJAS

This github repository contains the programming files in R, corresponding to the Degree Thesis **Modelization of ocurrence data from plant species with complex networks**.
## Complex Networks.
The ocurrence data was obtained from [GBIF](https://www.gbif.org/). Once filtered, the data can be found in the folder **Ocurrence Data**. 
In **Funciones_redes_reales.R** there is a list of functions to obtain:
- Distance Matrix (*distance_matrix(data_network*)
- Adjancency Matrix (*adjacency_matrix(dis_matrix, dist_threshold)*)
- Weights matrix (*weight_adj_matrix(dist_matrix, max_disp, imp_disp)*)
- Modularity value (*mod_value(adj_matrix)*)
- Modularity matrix (membership matrix) (*mod_matrix(adj_matrix, mod_membership)*)

And extra functions to study properties of the network based on a distance vector and a sampling size vector, calculating mean values of *nrep* repetitions.
- Giant component (*LCC_real_network(doc.csv, nrep, sample_size, thres_vector)*)
- Modularity (*Mod_real_network(doc.csv, nrep, sample_size, thres_vector)*)
- Clustering coefficient (*CCoef_real_network(doc.csv, nrep, sample_size, thres_vector)*)
- Diameter (*Diameter_real_network(doc.csv, nrep, sample_size, thres_vector)*)
- Degree Connection (*EdgeD_real_network(doc.csv, nrep, sample_size, thres_vector)*)

The resulting data can be saved as a .csv document and plotted with [ggplot2](https://ggplot2.tidyverse.org/) using **Resultados_redes_reales.rmd** programme. Results will be saved in PDF format. There are two types of graphs:
- Property v Sample Size (with Size Vector in Legend)
- Property v Size Vector (with Sample Size in Legend)
