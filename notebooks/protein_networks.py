import pandas as pd
import sqlalchemy as sa
import common.graph_analysis

pd.options.mode.chained_assignment = None


def get_good_genes():
    sql_query = """\
    select distinct g
    from (
        select hgnc_name g
        from drug_to_hgnc_target
        union all
        select hgnc_name g
        from drug_to_hgnc_target_stitch
    ) t \
    """
    print(sql_query)
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/az_dream_2015')
    good_genes = set(pd.read_sql_query(sql_query, engine)['g'])
    print("Good genes: {}".format(len(good_genes)))
    return good_genes


def get_gex_bad_genes(cell_line):
    sql_query = """\
    select g
    from gex_gbgc
    where c = '{}' and gex_gbgc < 3.1 \
    """.format(cell_line)
    print(sql_query)
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/az_dream_2015_features')
    gex_bad_genes = set(pd.read_sql_query(sql_query, engine)['g'])
    print("Bad GEX genes: {}".format(len(gex_bad_genes)))
    return gex_bad_genes


def get_cnv_bad_genes(cell_line):
    sql_query = """\
    select g
    from cnv_gbgc
    where c = '{}' and (
        min_cn_gbgc < 2 or
        disruption_status_gbgc_max > 0
    ) \
    """.format(cell_line)
    print(sql_query)
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/az_dream_2015_features')
    cnv_bad_genes = set(pd.read_sql_query(sql_query, engine)['g'])
    print("Bad CNV genes: {}".format(len(cnv_bad_genes)))
    return cnv_bad_genes


def get_mutations_bad_genes(cell_line):
    sql_query = """\
    select g
    from mutations_gbgc
    where c = '{}' and (
        f_mutation_very_bad_density_gbgc > 0
        or f_mutation_maybe_bad_density_gbgc > 0.0005
        or f_mutation_density_gbgc > 0.005
    ) \
    """.format(cell_line)
    print(sql_query)
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/az_dream_2015_features')
    mutations_bad_genes = set(pd.read_sql_query(sql_query, engine)['g'])
    print("Mutated genes: {}".format(len(mutations_bad_genes)))
    return mutations_bad_genes


def get_ppi_data():
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/protein_networks')
    ppi_df = pd.read_sql_table('mentha', engine)
    ppi_df = (
        ppi_df[['gene_1', 'gene_2', 'weight']]
        .drop_duplicates(subset=['gene_1', 'gene_2'])  # note the duplicates here
        .rename(columns={'gene_1': 'id_1', 'gene_2': 'id_2'})
    )
    print("PPI nrows: {}".format(len(ppi_df)))
    return ppi_df


def get_string_data():
    engine = sa.create_engine(
        'mysql://biodata:kimlab-biodata@192.168.6.19:3306/string_10_0')
    string_df = pd.read_sql_query(
        "select * from string_10_0.hgnc_network where combined_score > 250",
        engine)
    string_df['weight'] = string_df['combined_score'] / 1000
    string_df = (
        string_df[['gene_1', 'gene_2', 'weight']]
        .rename(columns={'gene_1': 'id_1', 'gene_2': 'id_2'})
    )
    print("String nrows: {}".format(len(string_df)))
    return string_df


def gene_network_properties(df, bad_genes=None, good_genes=None, cell_line=None):
    """
    Parameters
    ----------
    good_genes : set
        Gene IDs that will be included in the final result.
    bad_genes : set
        Gene IDs that will be excluded from the network
        when calculating network properties.
    cell_line : str
        The cell line to include under the `c` column
        in the final DataFrame.
    """
    print("Calculating graph properties...")

    if bad_genes is not None:
        print("Rows before removing bad genes: {}".format(len(df)))
        df = (
            df[
                (~df['id_1'].isin(bad_genes)) &
                (~df['id_2'].isin(bad_genes))
            ]
        )
        print("Rows after removing bad genes: {}".format(len(df)))

    vertex_df, edge_df, all_edge_df = common.graph_analysis.main(df)

    vertex_df = (
        vertex_df
        .reset_index()
        .rename(columns={'id': 'g'})
        .drop('graph_idx', axis=1)
    )

    edge_df = (
        edge_df
        .rename(columns={'id_1': 'g_1', 'id_2': 'g_2'})
        .drop(pd.Index(['graph_idx_1', 'graph_idx_2']), axis=1)
    )

    all_edge_df = (
        all_edge_df
        .rename(columns={'id_1': 'g_1', 'id_2': 'g_2'})
    )

    if good_genes is not None:
        print("Removing all but the good genes...")
        vertex_df = vertex_df[
            vertex_df['g'].isin(good_genes)
        ]
        edge_df = edge_df[
            edge_df['g_1'].isin(good_genes) &
            edge_df['g_2'].isin(good_genes)
        ]
        all_edge_df = all_edge_df[
            all_edge_df['g_1'].isin(good_genes) &
            all_edge_df['g_2'].isin(good_genes)
        ]

    if cell_line is not None:
        vertex_df['c'] = cell_line
        edge_df['c'] = cell_line
        all_edge_df['c'] = cell_line

    return vertex_df, edge_df, all_edge_df


def save_data(vertex_df, edge_df, all_edge_df, prefix, cell_line=None):
    import csv2sql
    db = csv2sql.DataFrameToMySQL(
        'mysql://strokach:@192.168.6.19:3306/az_dream_2015_features',
        'protein_networks/{}'.format(cell_line if cell_line else ''),
        '192.168.6.8',
        echo=False
    )
    db.import_table(vertex_df, prefix + '_vertex_gbg', if_exists='append')
    db.import_table(edge_df, prefix + '_edge_gbgg', if_exists='append')
    db.import_table(all_edge_df, prefix + '_all_edge_gbgg', if_exists='append')


if __name__ == '__main__':
    import common
    import argparse

    common.configure_logging()

    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_line', type=str, default=None)
    args = parser.parse_args()

    print("Cell line: {}".format(args.cell_line))

    good_genes = get_good_genes()

    if args.cell_line is None:
        gene_sets = [('all', None), ]
    else:
        gex_bad_genes = get_gex_bad_genes(args.cell_line)
        cnv_bad_genes = get_cnv_bad_genes(args.cell_line)
        mutations_bad_genes = get_mutations_bad_genes(args.cell_line)
        cnv_mutation_bad_genes = cnv_bad_genes | mutations_bad_genes
        print("Bad CNV or mutated genes: {}".format(len(cnv_mutation_bad_genes)))
        gene_sets = [('gex', gex_bad_genes), ('cnv_mut', cnv_mutation_bad_genes)]

    ppi_df = get_ppi_data()
    for prefix, bad_genes in gene_sets:
        vertex_df, edge_df, all_edge_df = gene_network_properties(
            ppi_df, bad_genes, good_genes, args.cell_line)
        save_data(vertex_df, edge_df, all_edge_df, 'ppi_' + prefix, args.cell_line)
    del ppi_df, vertex_df, edge_df, all_edge_df

    # string_df = get_string_data()
    # for prefix, bad_genes in gene_sets:
    #    vertex_df, edge_df, all_edge_df = gene_network_properties(
    #        string_df, bad_genes, good_genes, args.cell_line)
    #    save_data(vertex_df, edge_df, all_edge_df, 'string_' + prefix, args.cell_line)
    # del string_df, vertex_df, edge_df, all_edge_df
    #
    print("Done!")
