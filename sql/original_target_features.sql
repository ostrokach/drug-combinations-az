drop table if exists az_dream_2015_features.clare_gbdd_stitch;
show warnings;
create table az_dream_2015_features.clare_gbdd_stitch as (

-- drop table if exists az_dream_2015_features.clare_gbdd;
-- show warnings;
-- create table az_dream_2015_features.clare_gbdd as (
select
p.unique_id,
p.d_1 d_1,
p.d_2 d_2,

-- Biogrid
avg(b.bg_degree) bg_degree_mean_stitch,
max(b.bg_degree) bg_degree_max_stitch,
min(b.bg_degree) bg_degree_min_stitch,
std(b.bg_degree) bg_degree_std_stitch,

avg(b.bg_clustering_coef) bg_clustering_coef_mean_stitch,
max(b.bg_clustering_coef) bg_clustering_coef_max_stitch,
min(b.bg_clustering_coef) bg_clustering_coef_min_stitch,
std(b.bg_clustering_coef) bg_clustering_coef_std_stitch,

avg(b.bg_betweenness) bg_betweenness_mean_stitch,
max(b.bg_betweenness) bg_betweenness_max_stitch,
min(b.bg_betweenness) bg_betweenness_min_stitch,
std(b.bg_betweenness) bg_betweenness_std_stitch,

avg(b.bg_closeness) bg_closeness_mean_stitch,
max(b.bg_closeness) bg_closeness_max_stitch,
min(b.bg_closeness) bg_closeness_min_stitch,
std(b.bg_closeness) bg_closeness_std_stitch,

avg(b.bg_neighbor_sharing) bg_neighbor_sharing_mean_stitch,
max(b.bg_neighbor_sharing) bg_neighbor_sharing_max_stitch,
min(b.bg_neighbor_sharing) bg_neighbor_sharing_min_stitch,
std(b.bg_neighbor_sharing) bg_neighbor_sharing_std_stitch,

avg(b.bg_shortest_path_length) bg_shortest_path_length_mean_stitch,
min(b.bg_shortest_path_length) bg_shortest_path_length_min_stitch,

max(b_eb.bg_eb_max) bg_eb_max_max_stitch,
min(b_eb.bg_eb_min) bg_eb_min_min_stitch,
avg(b_eb.bg_eb_mean) bg_eb_mean_mean_stitch,
avg(b_eb.bg_eb_fraction) bg_eb_fraction_mean_stitch,

avg(b_nsp.bg_number_of_shortest_paths) bg_number_of_shortest_paths_mean_stitch,
max(b_nsp.bg_number_of_shortest_paths) bg_number_of_shortest_paths_max_stitch,
min(b_nsp.bg_number_of_shortest_paths) bg_number_of_shortest_paths_min_stitch,
std(b_nsp.bg_number_of_shortest_paths) bg_number_of_shortest_paths_std_stitch,


-- Co-expression and essentiality
avg(gc.coexpression) coexpression_mean_stitch,
max(gc.coexpression) coexpression_max_stitch,
min(gc.coexpression) coexpression_min_stitch,
std(gc.coexpression) coexpression_std_stitch,

chemical_interactions_v2.gene_pair_essentiality(ge1.gene_essentiality, ge2.gene_essentiality) gene_essentiality_max_stitch,


-- Genetic interaction
avg(gi.gi_degree) gi_degree_mean_stitch,
max(gi.gi_degree) gi_degree_max_stitch,
min(gi.gi_degree) gi_degree_min_stitch,
std(gi.gi_degree) gi_degree_std_stitch,

avg(gi.gi_clustering_coef) gi_clustering_coef_mean_stitch,
max(gi.gi_clustering_coef) gi_clustering_coef_max_stitch,
min(gi.gi_clustering_coef) gi_clustering_coef_min_stitch,
std(gi.gi_clustering_coef) gi_clustering_coef_std_stitch,

avg(gi.gi_betweenness) gi_betweenness_mean_stitch,
max(gi.gi_betweenness) gi_betweenness_max_stitch,
min(gi.gi_betweenness) gi_betweenness_min_stitch,
std(gi.gi_betweenness) gi_betweenness_std_stitch,

avg(gi.gi_closeness) gi_closeness_mean_stitch,
max(gi.gi_closeness) gi_closeness_max_stitch,
min(gi.gi_closeness) gi_closeness_min_stitch,
std(gi.gi_closeness) gi_closeness_std_stitch,

avg(gi.gi_neighbor_sharing) gi_neighbor_sharing_mean_stitch,
max(gi.gi_neighbor_sharing) gi_neighbor_sharing_max_stitch,
min(gi.gi_neighbor_sharing) gi_neighbor_sharing_min_stitch,
std(gi.gi_neighbor_sharing) gi_neighbor_sharing_std_stitch,

avg(gi.gi_shortest_path_length) gi_shortest_path_length_mean_stitch,
min(gi.gi_shortest_path_length) gi_shortest_path_length_min_stitch,

max(gi_eb.gi_eb_max) gi_eb_max_max_stitch,
min(gi_eb.gi_eb_min) gi_eb_min_min_stitch,
avg(gi_eb.gi_eb_mean) gi_eb_mean_mean_stitch,
avg(gi_eb.gi_eb_fraction) gi_eb_fraction_mean_stitch,

avg(gi_snp.gi_number_of_shortest_paths) gi_number_of_shortest_paths_mean_stitch,
max(gi_snp.gi_number_of_shortest_paths) gi_number_of_shortest_paths_max_stitch,
min(gi_snp.gi_number_of_shortest_paths) gi_number_of_shortest_paths_min_stitch,
std(gi_snp.gi_number_of_shortest_paths) gi_number_of_shortest_paths_std_stitch,


-- Gene ontology
avg(go_all.go_all_sem_sim) go_all_sem_sim_mean_stitch,
max(go_all.go_all_sem_sim) go_all_sem_sim_max_stitch,
min(go_all.go_all_sem_sim) go_all_sem_sim_min_stitch,
std(go_all.go_all_sem_sim) go_all_sem_sim_std_stitch,

avg(go_bp.go_bp_sem_sim) go_bp_sem_sim_mean_stitch,
max(go_bp.go_bp_sem_sim) go_bp_sem_sim_max_stitch,
min(go_bp.go_bp_sem_sim) go_bp_sem_sim_min_stitch,
std(go_bp.go_bp_sem_sim) go_bp_sem_sim_std_stitch,

avg(go_cc.go_cc_sem_sim) go_cc_sem_sim_mean_stitch,
max(go_cc.go_cc_sem_sim) go_cc_sem_sim_max_stitch,
min(go_cc.go_cc_sem_sim) go_cc_sem_sim_min_stitch,
std(go_cc.go_cc_sem_sim) go_cc_sem_sim_std_stitch,

avg(go_mf.go_mf_sem_sim) go_mf_sem_sim_mean_stitch,
max(go_mf.go_mf_sem_sim) go_mf_sem_sim_max_stitch,
min(go_mf.go_mf_sem_sim) go_mf_sem_sim_min_stitch,
std(go_mf.go_mf_sem_sim) go_mf_sem_sim_std_stitch,

avg(go_slim.go_slim_sem_sim) go_slim_sem_sim_mean_stitch,
max(go_slim.go_slim_sem_sim) go_slim_sem_sim_max_stitch,
min(go_slim.go_slim_sem_sim) go_slim_sem_sim_min_stitch,
std(go_slim.go_slim_sem_sim) go_slim_sem_sim_std_stitch,


-- Phylogenetic similarity
avg(phylo.phylogenic_similarity) phylogenic_similarity_mean_stitch,
max(phylo.phylogenic_similarity) phylogenic_similarity_max_stitch,
min(phylo.phylogenic_similarity) phylogenic_similarity_min_stitch,
std(phylo.phylogenic_similarity) phylogenic_similarity_std_stitch,


-- String
avg(s.s_degree) s_degree_mean_mean_stitch,
max(s.s_degree) s_degree_mean_max_stitch,
min(s.s_degree) s_degree_mean_min_stitch,
std(s.s_degree) s_degree_mean_std_stitch,

avg(s.s_clustering_coef) s_clustering_coef_mean_mean_stitch,
max(s.s_clustering_coef) s_clustering_coef_mean_max_stitch,
min(s.s_clustering_coef) s_clustering_coef_mean_min_stitch,
std(s.s_clustering_coef) s_clustering_coef_mean_std_stitch,

avg(s.s_betweenness) s_betweenness_mean_stitch,
max(s.s_betweenness) s_betweenness_max_stitch,
min(s.s_betweenness) s_betweenness_min_stitch,
std(s.s_betweenness) s_betweenness_std_stitch,

avg(s.s_closeness) s_closeness_mean_stitch,
max(s.s_closeness) s_closeness_max_stitch,
min(s.s_closeness) s_closeness_min_stitch,
std(s.s_closeness) s_closeness_std_stitch,

avg(s.s_neighbor_sharing) s_neighbor_sharing_mean_stitch,
max(s.s_neighbor_sharing) s_neighbor_sharing_max_stitch,
min(s.s_neighbor_sharing) s_neighbor_sharing_min_stitch,
std(s.s_neighbor_sharing) s_neighbor_sharing_std_stitch,

avg(s.s_shortest_path_length) s_shortest_path_length_mean_stitch,
min(s.s_shortest_path_length) s_shortest_path_length_min_stitch,

max(s_eb.s_eb_max) s_eb_max_max_stitch,
min(s_eb.s_eb_min) s_eb_min_min_stitch,
avg(s_eb.s_eb_mean) s_eb_mean_mean_stitch,
avg(s_eb.s_eb_fraction) s_eb_fraction_mean_stitch,

avg(s_nsp.s_number_of_shortest_paths) s_number_of_shortest_paths_mean_stitch,
max(s_nsp.s_number_of_shortest_paths) s_number_of_shortest_paths_max_stitch,
min(s_nsp.s_number_of_shortest_paths) s_number_of_shortest_paths_min_stitch,
std(s_nsp.s_number_of_shortest_paths) s_number_of_shortest_paths_std_stitch

-- from az_dream_2015.drug_pair_to_ensp_idx_pair p
from az_dream_2015.drug_pair_to_ensp_idx_pair_stitch p
left join chemical_interactions_v2.biogrid_topo b using (ensp_1, ensp_2)
left join chemical_interactions_v2.biogrid_topo_eb b_eb using (ensp_1, ensp_2)
left join chemical_interactions_v2.biogrid_topo_nsp b_nsp using (ensp_1, ensp_2)
left join chemical_interactions_v2.gene_coexpression gc using (ensp_1, ensp_2)
left join chemical_interactions_v2.gene_essentiality ge1 on (ge1.ensp = p.ensp_1)
left join chemical_interactions_v2.gene_essentiality ge2 on (ge2.ensp = p.ensp_1)
left join chemical_interactions_v2.getint_topo gi using (ensp_1, ensp_2)
left join chemical_interactions_v2.getint_topo_eb gi_eb using (ensp_1, ensp_2)
left join chemical_interactions_v2.getint_topo_nsp gi_snp using (ensp_1, ensp_2)
left join chemical_interactions_v2.go_all go_all using (ensp_1, ensp_2)
left join chemical_interactions_v2.go_bp go_bp using (ensp_1, ensp_2)
left join chemical_interactions_v2.go_cc go_cc using (ensp_1, ensp_2)
left join chemical_interactions_v2.go_mf go_mf using (ensp_1, ensp_2)
left join chemical_interactions_v2.go_slim go_slim using (ensp_1, ensp_2)
left join chemical_interactions_v2.phylo phylo using (ensp_1, ensp_2)
left join chemical_interactions_v2.string_topo s using (ensp_1, ensp_2)
left join chemical_interactions_v2.string_topo_eb s_eb using (ensp_1, ensp_2)
left join chemical_interactions_v2.string_topo_nsp s_nsp using (ensp_1, ensp_2)
group by unique_id
);

show warnings;

-- create index unique_id_idx on az_dream_2015_features.clare_gbdd (unique_id);
-- create index a on az_dream_2015_features.clare_gbdd (d_1, d_2);
-- create index b on az_dream_2015_features.clare_gbdd (d_2, d_1);

create index unique_id_idx on az_dream_2015_features.clare_gbdd_stitch (unique_id);
create index a on az_dream_2015_features.clare_gbdd_stitch (d_1, d_2);
create index b on az_dream_2015_features.clare_gbdd_stitch (d_2, d_1);

show warnings;
