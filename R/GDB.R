
#### Create api for neo4j session
#' Create api for neo4j session
#'
#' @param user
#' @param password
#'
#' @return
#' @import neo4r
#'
#' @examples
get_neo_api <- function(user = "neo4j", password = "password"){
  gdb <- neo4j_api$new(url = "http://localhost:7474", user = user, password = password)
  return(gdb)
}

#### Checking connection to neo4j
#' Checking connection to neo4j
#'
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
neo_api_ping <- function(gdb){
  ping <- NULL
  if(class(try(gdb$ping(), silent = TRUE)) == "try-error"){
    return(paste0("There is a problem to connect to neo4j server.\n",
                  "Check that neo4j database is created and activated."))}
  ping <- gdb$ping()
  if(as.numeric(ping) == 200){return(NULL)}
  if(as.numeric(ping) != 200){return(paste0("There is a problem to connect to neo4j server.\n",
                                           "Ping answer: ", ping))}
}

#### Add samples as nodes to neo4j database
#' Add samples as nodes to neo4j database
#'
#' @param fcs_raw
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_sample_GDB <- function(fcs_raw, gdb){
  call_neo4j('CREATE CONSTRAINT ON (b:Sample) ASSERT b.name IS UNIQUE;', gdb)
  for(i in 1:length(fcs_raw)){
    call_neo4j(paste0('CREATE (:Sample { name: \'', sampleNames(fcs_raw)[i],
                      '\', cell_number: ', nrow(fcs_raw[[i]]), '})'), gdb)
  }
}

#### Add markers as nodes to neo4j database
#' Add markers as nodes to neo4j database
#'
#' @param use_markers
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_marker_GDB <- function(use_markers, gdb){
  call_neo4j('CREATE CONSTRAINT ON (b:Marker) ASSERT b.name IS UNIQUE;', gdb)
  for(i in 1:length(use_markers)){
    call_neo4j(paste0('CREATE (:Marker { name: \'', use_markers[i],
                      '\', antigen: \'', names(use_markers)[i], '\' })'), gdb)
  }
}

#### Add clusters as nodes to neo4j database
#' Add clusters as nodes to neo4j database
#'
#' @param cell_clustering
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_cluster_GDB <- function(cell_clustering, gdb){
  call_neo4j('CREATE CONSTRAINT ON (b:Cluster) ASSERT b.name IS UNIQUE;', gdb)
  cluster <- table(cell_clustering)
  for(i in 1:length(cluster)){
    call_neo4j(paste0('CREATE (:Cluster { name: \'', names(cluster)[i],
                      '\', cell_number: \'', cluster[[i]], '\' })'), gdb)
  }
}

#### Allocate cell populations within each sample based on clusterisation
#' Allocate cell populations within each sample based on clusterisation
#'
#' @param cell_clustering_list
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_populatio_GDB <- function(cell_clustering_list, gdb){
  call_neo4j('CREATE CONSTRAINT ON (b:Population) ASSERT b.name IS UNIQUE;', gdb)
  for(i in 1:length(cell_clustering_list)){
    cluster <- table(cell_clustering_list[[i]])
    for(j in 1:length(cluster)){
      call_neo4j(paste0('MATCH (s:Sample { name: \'',  names(cell_clustering_list)[[i]],
                        '\' }), (c:Cluster { name: \'', names(cluster)[j],
                        '\' }) CREATE (s)-[:Contains]->(p:Population { name: \'',
                        names(cluster)[j], '_', names(cell_clustering_list)[[i]],
                        '\', sample: \'', names(cell_clustering_list)[[i]],
                        '\', cluster: \'', names(cluster)[j],
                        '\', cell_number: \'', cluster[[i]], '\' }) ',
                        'CREATE (c)-[:Contains]->(p)'), gdb)
    }

    }
}

#### Allocate observations as nodes based on marker and cell population nodes
#' Allocate observations as nodes based on marker and cell population nodes
#'
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_observation_GDB <- function(gdb){
  call_neo4j('CREATE CONSTRAINT ON (o:Observation) ASSERT o.name IS UNIQUE;', gdb)
  p <- unnest_nodes(call_neo4j('MATCH (p:Population) RETURN p', gdb, type = 'graph')$nodes)
  m <- unnest_nodes(call_neo4j('MATCH (m:Marker) RETURN m', gdb, type = 'graph')$nodes)
  for (i in 1:nrow(m)) {
    for (j in 1:nrow(p)) {
      call_neo4j(paste0('MATCH (m:Marker { name: \'', m$name[i],
                        '\' }), (p:Population { name: \'', p$name[j],
                        '\' }) CREATE (m)-[:Contains]->(o:Observation { name: \'',
                        m$name[i], '_', p$name[j],
                        '\', marker: \'', m$name[i],
                        '\', antigen: \'', m$antigen[i],
                        '\', sample: \'', p$sample[j],
                        '\', cluster: \'', p$cluster[j],
                        '\', cell_number: \'', p$cell_number[j], '\' }) ',
                        'CREATE (p)-[:Contains]->(o)'), gdb)
    }
  }
}

#### Allocate phenounite (phenotype units) nodes based on markers and clusters in database
#' Allocate phenounite, phenotype units, nodes based on markers and clusters in database
#'
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_phenounite_GDB <- function(gdb){
  call_neo4j('CREATE CONSTRAINT ON (u:Phenounite) ASSERT u.name IS UNIQUE;', gdb)
  m <- as.data.frame(unnest_nodes(call_neo4j('MATCH (m:Marker) RETURN m', gdb, type = 'graph')$nodes))
  c <- as.data.frame(unnest_nodes(call_neo4j('MATCH (c:Cluster) RETURN c', gdb, type = 'graph')$nodes))
  for (i in 1:nrow(m)) {
    for (j in 1:nrow(c)) {
      call_neo4j(paste0('MATCH (m:Marker { name: \'', m$name[i],
                        '\' }), (c:Cluster { name: \'', c$name[j],
                        '\' }) CREATE (m)-[:Contains]->(u:Phenounite { name: \'',
                        c$name[j], '_', m$name[i],
                        '\', marker: \'', m$name[i],
                        '\', antigen: \'', m$antigen[i],
                        '\', cluster: \'', c$name[j], '\' }) ',
                        'CREATE (c)-[:Contains]->(u)'), gdb)
    }
  }
  o <- unnest_nodes(call_neo4j('MATCH (o:Observation) RETURN o', gdb, type = 'graph')$nodes)
  for (i in 1:nrow(o)) {
    call_neo4j(paste0('MATCH (u:Phenounite { name: \'',
                      o$cluster[i], '_', o$marker[i],
                      '\' }), (o:Observation { name: \'', o$name[i],
                      '\' }) CREATE (o)-[:Contains]->(u)'), gdb)
  }

}

#### Add correlation between clusters as relationships to the database
#' Add correlation between clusters as relationships to the database
#'
#' @param signals_between_clusters
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_signals_between_clusters_GDB <- function(signals_between_clusters, gdb){
  for (i in 1:nrow(signals_between_clusters)){
    call_neo4j(paste0('MATCH (c:Cluster { name: \'', as.character(signals_between_clusters$signaling_cluster[i]),
                      '\' }), (u:Phenounite { name: \'',
                      as.character(signals_between_clusters$targetet_cluster[i]),
                      '_', as.character(signals_between_clusters$gene_in_target_cluster[i]),
                      '\' }) CREATE (u)-[:Correlation { cor_significance: ',
                      signals_between_clusters$cor_significance[i],
                      ', cor_coefficient: ', signals_between_clusters$cor_coefficient[i],
                      ', adj_pValue: ', signals_between_clusters$adj_pValue[i], ' }]->(c)'), gdb)
  }
}

#### Add correlation within clusters as relationships to the database
#' Add correlation within clusters as relationships to the database
#'
#' @param signals_in_cluster
#' @param gdb
#'
#' @return
#' @import neo4r
#'
#' @examples
add_signals_in_cluster_GDB <- function(signals_in_cluster, gdb){
  for (i in 1:nrow(signals_in_cluster)){
    call_neo4j(paste0('MATCH (c:Cluster { name: \'', as.character(signals_in_cluster$cluster[i]),
                      '\' }), (u:Phenounite { name: \'',
                      as.character(signals_in_cluster$cluster[i]),
                      '_', as.character(signals_in_cluster$gene[i]),
                      '\' }) CREATE (u)-[:Correlation { cor_significance: ',
                      signals_in_cluster$cor_significance[i],
                      ', cor_coefficient: ', signals_in_cluster$cor_coefficient[i],
                      ', adj_pValue: ', signals_in_cluster$adj_pValue[i], ' }]->(c)'), gdb)
  }
}
