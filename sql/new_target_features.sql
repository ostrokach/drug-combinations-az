# PPI

drop table if exists ppi_full_vertex_gbd;
create table ppi_full_vertex_gbd as (
  SELECT
  d2g.drug d,

  max(f.degree) ppi_vertex_degree_gbd_max,
  min(f.degree) ppi_vertex_degree_gbd_min,
  avg(f.degree) ppi_vertex_degree_gbd_avg,
  std(f.degree) ppi_vertex_degree_gbd_std,

  max(f.closeness) ppi_vertex_closeness_gbd_max,
  min(f.closeness) ppi_vertex_closeness_gbd_min,
  avg(f.closeness) ppi_vertex_closeness_gbd_avg,
  std(f.closeness) ppi_vertex_closeness_gbd_std,

  max(f.betweenness) ppi_vertex_betweenness_gbd_max,
  min(f.betweenness) ppi_vertex_betweenness_gbd_min,
  avg(f.betweenness) ppi_vertex_betweenness_gbd_avg,
  std(f.betweenness) ppi_vertex_betweenness_gbd_std,

  max(f.clustering_coef) ppi_vertex_clustering_coef_gbd_max,
  min(f.clustering_coef) ppi_vertex_clustering_coef_gbd_min,
  avg(f.clustering_coef) ppi_vertex_clustering_coef_gbd_avg,
  std(f.clustering_coef) ppi_vertex_clustering_coef_gbd_std,

  count(*) ppi_vertex_gbd_count

  FROM az_dream_2015.drug_to_hgnc_target d2g
  JOIN az_dream_2015_features.ppi_full_vertex_gbg f ON (d2g.hgnc_name = f.g)
  GROUP BY d2g.drug
);
show warnings;
create index a on ppi_full_vertex_gbd (d);
show warnings;


drop table if exists ppi_full_edge_gbdd;
create table ppi_full_edge_gbdd as (
  SELECT
  dd2gg.d_1,
  dd2gg.d_2,

  max(f.edge_betweenness) ppi_full_edge_betweenness_gbdd_max,
  min(f.edge_betweenness) ppi_full_edge_betweenness_gbdd_min,
  avg(f.edge_betweenness) ppi_full_edge_betweenness_gbdd_avg,
  std(f.edge_betweenness) ppi_full_edge_betweenness_gbdd_std,

  count(*) ppi_full_edge_betweenness_gbdd_count

  FROM az_dream_2015.drug_pair_to_gene_pair dd2gg
  JOIN az_dream_2015_features.ppi_full_edge_gbgg f USING (g_1, g_2)
  GROUP BY dd2gg.d_1, dd2gg.d_2
);
show warnings;
create index a on ppi_full_edge_gbdd (d_1, d_2);
show warnings;


drop table if exists ppi_full_all_edge_gbdd;
create table ppi_full_all_edge_gbdd as (
  SELECT
  dd2gg.d_1,
  dd2gg.d_2,

  min(f.distance_min) ppi_all_full_distance_min_gbdd_min,
  avg(f.distance_min) ppi_all_full_distance_min_gbdd_avg,

  max(f.similarity_jaccard) ppi_full_all_similarity_jaccard_gbdd_max,
  min(f.similarity_jaccard) ppi_full_all_similarity_jaccard_gbdd_min,
  avg(f.similarity_jaccard) ppi_full_all_similarity_jaccard_gbdd_avg,
  std(f.similarity_jaccard) ppi_full_all_similarity_jaccard_gbdd_std,

  max(f.similarity_inverse_log_weighted) ppi_full_all_similarity_inverse_log_weighted_gbdd_max,
  min(f.similarity_inverse_log_weighted) ppi_full_all_similarity_inverse_log_weighted_gbdd_min,
  avg(f.similarity_inverse_log_weighted) ppi_full_all_similarity_inverse_log_weighted_gbdd_avg,
  std(f.similarity_inverse_log_weighted) ppi_full_all_similarity_inverse_log_weighted_gbdd_std,

  count(*) ppi_full_all_gbdd_count

  FROM az_dream_2015.drug_pair_to_gene_pair dd2gg
  JOIN az_dream_2015_features.ppi_full_all_edge_gbgg f USING (g_1, g_2)
  GROUP BY dd2gg.d_1, dd2gg.d_2
);
show warnings;
create index a on ppi_full_all_edge_gbdd (d_1, d_2);
show warnings;
