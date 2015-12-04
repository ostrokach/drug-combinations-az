create database az_dream_target_features;

create table az_dream_target_features.original_target_features as (
select 
p.*,
b.bg_degree,
b.bg_clustering_coef,
b.bg_betweenness,
b.bg_closeness,
b.bg_neighbor_sharing,
b.bg_shortest_path_length,
b_eb.bg_eb_max,
b_eb.bg_eb_min,
b_eb.bg_eb_mean,
b_eb.bg_eb_fraction,
b_nsp.bg_number_of_shortest_paths,
gc.coexpression,
gi.gi_degree,
gi.gi_clustering_coef,
gi.gi_betweenness,
gi.gi_closeness,
gi.gi_neighbor_sharing,
gi.gi_shortest_path_length,
gi_eb.gi_eb_max,
gi_eb.gi_eb_min,
gi_eb.gi_eb_mean,
gi_eb.gi_eb_fraction,
gi_snp.gi_number_of_shortest_paths,
go_all.go_all_sem_sim,
go_bp.go_bp_sem_sim,
go_cc.go_cc_sem_sim,
go_mf.go_mf_sem_sim,
go_slim.go_slim_sem_sim,
phylo.phylogenic_similarity,
s.s_degree,
s.s_clustering_coef,
s.s_betweenness,
s.s_closeness,
s.s_neighbor_sharing,
s.s_shortest_path_length,
s_eb.eb_max,
s_eb.eb_min,
s_eb.eb_mean,
s_eb.eb_fraction,
s_nsp.s_number_of_shortest_paths
-- chemical_interactions_v2.gene_pair_essentiality(ge1.gene_essentiality, ge2.gene_essentiality) gene_pair_essentiality,
from az_dream.drug_pair_to_ensp_idx_pair p
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
