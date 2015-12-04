-- 
drop table drug_pair_to_ensp_pair;

create table drug_pair_to_ensp_pair (
unique_id varchar(255), 
ensp_pair varchar(255),
ensp_1 varchar(255),
ensp_2 varchar(255) 
);

insert into drug_pair_to_ensp_pair 
SELECT 
unique_id,
ep.ensp_pair ensp_pair,
ep.ensp_1 ensp_1,
ep.ensp_2 ensp_2
FROM az_dream.all_drug_pairs dp
join az_dream.drug_to_target t1 on (t1.drug = dp.compound_a)
join az_dream.drug_to_target t2 on (t2.drug = dp.compound_b)
join az_dream.ensp_to_ensp_pair ep1 on (ep1.ensp = t1.ensp)
join az_dream.ensp_to_ensp_pair ep2 on (ep2.ensp = t2.ensp and ep1.ensp_pair = ep2.ensp_pair)
join az_dream.ensp_pairs ep on (ep.ensp_pair = ep2.ensp_pair);

create index a on drug_pair_to_ensp_pair (unique_id, ensp_pair);
create index b on drug_pair_to_ensp_pair (unique_id, ensp_1, ensp_2);
create index c on drug_pair_to_ensp_pair (ensp_pair, unique_id);
create index d on drug_pair_to_ensp_pair (ensp_1, ensp_2, unique_id);


-- 
drop table drug_pair_to_ensp_idx_pair;

create table drug_pair_to_ensp_idx_pair (
unique_id varchar(255), 
ensp_idx_pair varchar(255),
ensp_1 varchar(255),
ensp_2 varchar(255) 
);

insert into drug_pair_to_ensp_idx_pair 
SELECT 
unique_id,
ep.ensp_idx_pair ensp_idx_pair,
ep.ensp_1 ensp_1,
ep.ensp_2 ensp_2
FROM az_dream.all_drug_pairs dp
join az_dream.drug_to_target t1 on (t1.drug = dp.compound_a)
join az_dream.drug_to_target t2 on (t2.drug = dp.compound_b)
join az_dream.ensp_to_ensp_pair ep1 on (ep1.ensp = t1.ensp)
join az_dream.ensp_to_ensp_pair ep2 on (ep2.ensp = t2.ensp and ep1.ensp_idx_pair = ep2.ensp_idx_pair)
join az_dream.ensp_idx_pairs ep on (ep.ensp_idx_pair = ep2.ensp_idx_pair);

create index a on drug_pair_to_ensp_idx_pair (unique_id, ensp_idx_pair);
create index b on drug_pair_to_ensp_idx_pair (unique_id, ensp_1, ensp_2);
create index c on drug_pair_to_ensp_idx_pair (ensp_idx_pair, unique_id);
create index d on drug_pair_to_ensp_idx_pair (ensp_1, ensp_2, unique_id);


-- 
drop table drug_pair_to_gene_pair;

create table drug_pair_to_gene_pair (
unique_id varchar(255),
gene_pair varchar(255),
gene_name_1 varchar(255),
gene_name_2 varchar(255)
);

insert into drug_pair_to_gene_pair 
SELECT 
unique_id,
ep.gene_pair gene_pair,
ep.gene_name_1 gene_name_1,
ep.gene_name_2 gene_name_2
FROM az_dream.all_drug_pairs dp
join az_dream.drug_to_target t1 on (t1.drug = dp.compound_a)
join az_dream.drug_to_target t2 on (t2.drug = dp.compound_b)
join az_dream.hgnc_to_hgnc_pair ep1 on (ep1.gene_name = t1.hgnc_name)
join az_dream.hgnc_to_hgnc_pair ep2 on (ep2.gene_name = t2.hgnc_name and ep2.gene_pair = ep1.gene_pair)
join az_dream.hgnc_pairs ep on (ep.gene_pair = ep2.gene_pair);

create index a on drug_pair_to_gene_pair (unique_id, gene_pair);
create index b on drug_pair_to_gene_pair (unique_id, gene_name_1, gene_name_2);
create index c on drug_pair_to_gene_pair (gene_pair, unique_id);
create index d on drug_pair_to_gene_pair (gene_name_1, gene_name_2, unique_id);
