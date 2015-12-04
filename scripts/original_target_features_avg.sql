create table az_dream_target_features.original_target_features as (
select 
p.unique_id,
avg(b.bg_degree) bg_degree,
avg(b.bg_clustering_coef) bg_clustering_coef, 
avg(b.bg_betweenness) bg_betweenness,
avg(b.bg_closeness) bg_closeness,
avg(b.bg_neighbor_sharing) bg_neighbor_sharing,
avg(b.bg_shortest_path_length) bg_shortest_path_length,
avg(b_eb.bg_eb_max) bg_eb_max,
avg(b_eb.bg_eb_min) bg_eb_min,
avg(b_eb.bg_eb_mean) bg_eb_mean,
avg(b_eb.bg_eb_fraction) bg_eb_fraction,
avg(b_nsp.bg_number_of_shortest_paths) bg_number_of_shortest_paths,
avg(gc.coexpression) coexpression,
avg(avg(ge1.gene_essentiality), avg(ge2.gene_essentiality)) gene_essentiality,
avg(gi.gi_degree) gi_degree,
avg(gi.gi_clustering_coef) gi_clustering_coef,
avg(gi.gi_betweenness) gi_betweenness,
avg(gi.gi_closeness) gi_closeness,
avg(gi.gi_neighbor_sharing) gi_neighbor_sharing,
avg(gi.gi_shortest_path_length) gi_shortest_path_length,
avg(gi_eb.gi_eb_max) gi_eb_max,
avg(gi_eb.gi_eb_min) gi_eb_min,
avg(gi_eb.gi_eb_mean) gi_eb_mean,
avg(gi_eb.gi_eb_fraction) gi_eb_fraction,
avg(gi_snp.gi_number_of_shortest_paths) gi_number_of_shortest_paths,
avg(go_all.go_all_sem_sim) go_all_sem_sim,
avg(go_bp.go_bp_sem_sim) go_bp_sem_sim,
avg(go_cc.go_cc_sem_sim) go_cc_sem_sim,
avg(go_mf.go_mf_sem_sim) go_mf_sem_sim,
avg(go_slim.go_slim_sem_sim) go_slim_sem_sim,
avg(phylo.phylogenic_similarity) phylogenic_similarity,
avg(s.s_degree) s_degree,
avg(s.s_clustering_coef) s_clustering_coef,
avg(s.s_betweenness) s_betweenness,
avg(s.s_closeness) s_closeness,
avg(s.s_neighbor_sharing) s_neighbor_sharing,
avg(s.s_shortest_path_length) s_shortest_path_length,
avg(s_eb.s_eb_max) s_eb_max,
avg(s_eb.s_eb_min) s_eb_min,
avg(s_eb.s_eb_mean) s_eb_mean,
avg(s_eb.s_eb_fraction) s_eb_fraction,
avg(s_nsp.s_number_of_shortest_paths) s_number_of_shortest_paths
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
