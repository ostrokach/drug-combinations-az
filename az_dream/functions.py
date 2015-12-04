# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 18:42:55 2015

@author: strokach
"""
from IPython.display import display
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
import requests
from requests.utils import quote
import re


# %%
def get_weighted_stats(
        df, group_by_column, feature_column, weight_column,
        mean_range=(-100, 100), std_range=(0, 100)):
    """
    """
    def get_average(x):
        try:
            return np.average(x[feature_column].values, weights=x[weight_column])
        except ZeroDivisionError as e:
            print('The following error occured:\n{}'.format(e))
            return np.average(x[x[weight_column] > 0][feature_column].values)

    def get_std(x):
        try:
            return np.sqrt(np.average(
                (x[feature_column].values - get_average(x))**2, weights=x[weight_column]
            ))
        except ZeroDivisionError as e:
            print('The following error occured:\n{}'.format(e))
            return np.sqrt(np.average(
                (x[x[weight_column] > 0][feature_column].values - get_average(x))**2
            ))

    df_gp = df.groupby(group_by_column)

    df_2 = pd.DataFrame(
        data={
            'MEAN': df_gp.apply(get_average),
            'STD': df_gp.apply(get_std),
            'LEN': df_gp[weight_column].apply(np.sum),
            'MEAN_UNWEIGHTED': df_gp[feature_column].apply(np.mean),
            'STD_UNWEIGHTED': df_gp[feature_column].apply(np.std),
            'LEN_UNWEIGHTED': df_gp[feature_column].apply(len),
        },
    )
    df_2 = df_2.reset_index()

    # Normalize the mean
    df_2.loc[(df_2['MEAN'] <= mean_range[0]), 'MEAN'] = mean_range[0]
    df_2.loc[(df_2['MEAN'] >= mean_range[1]), 'MEAN'] = mean_range[1]

    # Normalize the standard deviation
    df_2['STE'] = df_2['STD']
    df_2.loc[(df_2['STD'] <= std_range[0]), 'STD'] = std_range[0]
    df_2.loc[(df_2['STD'] >= std_range[1]), 'STD'] = std_range[1]
    # Standard error
    df_2.loc[(df_2['STE'] <= std_range[0]), 'STE'] = np.nan
    df_2.loc[(df_2['STE'] >= std_range[1]), 'STE'] = np.nan
    temp_len = df_2['LEN'].apply(lambda x: max(x, 1))
    df_2['STE'] = df_2['STE'] / np.sqrt(temp_len)
    df_2.loc[pd.isnull(df_2['STE']), 'STE'] = np.nanmax(df_2['STE'])
    df_2['STE_UNWEIGHTED'] = df_2['STD_UNWEIGHTED'] / np.sqrt(df_2['LEN_UNWEIGHTED'])
    # Confidence score
    df_2['log_ste'] = np.log(df_2['STE'].values)
    df_2['norm_log_ste'] = (
        (df_2['log_ste'] - df_2['log_ste'].min()) /
        (df_2['log_ste'].max() - df_2['log_ste'].min())
    )
    # The higher the normalized standard error, the lower the confidence
    df_2['CONFIDENCE'] = 1 - df_2['norm_log_ste']

    df_2.drop(['log_ste', 'norm_log_ste'], axis=1, inplace=True)

    return df_2


# %%
def cross_validate(clf, df, x_columns, y_column):
    combination_ids = set(df['COMBINATION_ID'])
    for combination_id in combination_ids:
        print('Cross validating {}...'.format(combination_id))
        x_train = df[df['COMBINATION_ID'] != combination_id][x_columns].values
        y_train = df[df['COMBINATION_ID'] != combination_id][y_column].values
        clf.fit(x_train, y_train)
        x_test = df[df['COMBINATION_ID'] == combination_id][x_columns].values
        df.loc[df['COMBINATION_ID'] == combination_id, 'prediction'] = clf.predict(x_test)


# %%
class CrossValidate:

    def __init__(self, clf_options, df, x_columns, y_column):
        self.clf_options = clf_options
        self.df = df
        self.x_columns = x_columns
        self.y_column = y_column

    def run(self):
        self.clf_options['random_state'] = np.random.RandomState(42)
        clf = GradientBoostingRegressor(**self.clf_options)
        cross_validate(clf, self.df, self.x_columns, self.y_column)

    def save(self, filename):
        self.df['prediction'].to_pickle(filename)


# %%
pug = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
properties = [
    'MolecularFormula',
    'MolecularWeight',
    'CanonicalSMILES',
    'IsomericSMILES',
    'InChI',
    'InChIKey',
    'IUPACName',
    'XLogP',
    'ExactMass',
    'MonoisotopicMass',
    'TPSA',
    'Complexity',
    'Charge',
    'HBondDonorCount',
    'HBondAcceptorCount',
    'RotatableBondCount',
    'HeavyAtomCount',
    'IsotopeAtomCount',
    'AtomStereoCount',
    'DefinedAtomStereoCount',
    'UndefinedAtomStereoCount',
    'BondStereoCount',
    'DefinedBondStereoCount',
    'UndefinedBondStereoCount',
    'CovalentUnitCount',
    'Volume3D',
    'XStericQuadrupole3D',
    'YStericQuadrupole3D',
    'ZStericQuadrupole3D',
    'FeatureCount3D',
    'FeatureAcceptorCount3D',
    'FeatureDonorCount3D',
    'FeatureAnionCount3D',
    'FeatureCationCount3D',
    'FeatureRingCount3D',
    'FeatureHydrophobeCount3D',
    'ConformerModelRMSD3D',
    'EffectiveRotorCount3D',
    'ConformerCount3D',
    'Fingerprint2D'
]
properties_string = ','.join(properties)


class NotFoundError(Exception):
    pass


def _get_properties(qtype, key):
    """
    Parameters
    ----------
    qtype : str
        {'cid', 'smiles'}
    """
    print(key)
    url = (
        pug + 'compound/{}/{}/property/{}/JSON'
        .format(qtype, quote(key, safe=''), properties_string)
    )
    r = requests.get(url)
    # print(url)
    try:
        properties = r.json()
    except ValueError as e:
        print(e)
        print(r.reason)
        if r.reason not in ['Not Found']:
            print(url)
        raise NotFoundError
    try:
        properties = properties['PropertyTable']['Properties']
    except KeyError as e:
        print(e)
        print(properties)
        print(url)
        raise NotFoundError
    assert len(properties) == 1
    return properties[0]


def get_drug_features(df):
    """Query the PubChem REST API for drug properties.

    https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    result_list = []
    for idx, row in df.iterrows():
        key = row['SMILES or PubChem ID']
        print()
        print(key)
        print('-' * 80)
        if pd.isnull(key):
            result_list.append(row)
            continue
        elif key.isnumeric():
            keys = [key]
            qtype = 'cid'
        else:
            qtype = 'smiles'
            keys = key.split(';')

        for key in keys:
            try:
                properties = _get_properties(qtype, key)
            except NotFoundError:
                continue
            row = row.append(pd.Series(properties))
            break
        result_list.append(row)
    return pd.concat(result_list, axis=1).transpose()


# %%
def get_drug_features_2(cids):
    cids = {int(c) for c in cids if pd.notnull(c)}
    url_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON'
    gene_names = {cid: [] for cid in cids}
    atc_codes = {cid: [] for cid in cids}
    for cid in cids:
        url = url_template.format(cid)
        r = requests.get(url)
        result = r.json()
        try:
            records = result['Record']['Section']
        except KeyError:
            print(result.keys())
            print(r.reason)
            print(url)
            raise
        return records
        for record in records:
            print(record)
            if 'TOCHeading' not in record:
                print(record)
                break

            if record['TOCHeading'] == 'Gene Name':
                print(record)
                for info in record['Information']:
                    if info['Name'] == 'Gene Name':
                        gene_names[cid].append(info['StringValue'])

            if record['TOCHeading'] == 'ATC Code':
                print(record)
                for info in record['Information']:
                    if info['Name'] == 'ATC Code':
                        atc_code = info['StringValue'].split()[0]
                        atc_code_2 = (
                            info['URL'][
                                info['URL'].index('?code=') + 6:
                                info['URL'].index('&')
                            ]
                        )
                        assert atc_code == atc_code_2
                        atc_codes[cid].append(atc_code)
    return gene_names, atc_codes


# %%
def p3print(recs, n_spaces):
    for rec in recs:
        if type(rec) == dict:
            if 'TOCHeading' in rec:
                print(' ' * n_spaces + 'TOCHeading: {}'.format(rec['TOCHeading']))
            if 'Description' in rec:
                print(' ' * n_spaces + 'Description: {}'.format(rec['Description']))
            if 'Information' in rec:
                print(' ' * (n_spaces + 4) + 'Information')
                for info in rec['Information']:
                    p3print(info, n_spaces+4)
            if 'Section' in rec:
                print(' ' * (n_spaces + 4) + 'Section')
                for section in rec['Section']:
                    p3print(section, n_spaces+4)
        else:
            print(rec)


# %%
def get_cids(cid):
    """
    https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    similar_cids_url = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/cids/JSON'
        '?cids_type=same_connectivity'
    )
    url = similar_cids_url.format(cid)
    r = requests.get(url)
    all_cids = r.json()['IdentifierList']['CID']
    return all_cids


def get_cid_website(cid):
    """
    https://pubchem.ncbi.nlm.nih.gov/
    """
    cid_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/'
    url = cid_url.format(cid)
    r = requests.get(url)
    return r.json()['Record']['Section']


def get_sids(cids):
    similar_sids_url = (
        'http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/sids/JSON'
    )
    chunk_size = 20
    sids = []
    for i in list(range(0, len(cids), chunk_size))[:-1]:
        cids_chunk = cids[i:i + chunk_size]
        url = similar_sids_url.format(','.join(str(cid) for cid in cids_chunk))
        r = requests.get(url)
        for info in r.json()['InformationList']['Information']:
            try:
                cid = info['CID']
                for sid in info['SID']:
                    sids.append((cid, sid))
            except KeyError:
                if 'CID' not in info or len(info) > 1:
                    raise
    return sids


def get_sid_data(sids):
    dbs_to_keep = [
        'ChEMBL', 'DrugBank', 'BindingDB', 'Broad Institute', 'KEGG', 'IUPHAR-DB',
    ]
    chunk_size = 20
    result = []
    for i in list(range(0, len(sids), chunk_size))[:-1]:
        sids_chunk = sids[i:i + chunk_size]
        sid_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{}/JSON'
        url = sid_url.format(','.join(str(sid) for sid in sids_chunk))
        r = requests.get(url)
        for info in r.json()['PC_Substances']:
            sid = info['sid']['id']
            db_name = info['source']['db']['name']
            db_id = info['source']['db']['source_id']['str']
            if db_name in dbs_to_keep:
                result.append((sid, db_name, db_id))
    return result


def get_atcs(cids):
    names = []
    mesh_synonyms = []
    depositor_supplied_synonyms = []
    atcs = []
    targets = []
    general_functions = []
    #
    remove_links = re.compile('[<>\(\)]')
    for cid in cids:
        recs = get_cid_website(cid)
        for rec in recs:
            if rec.get('TOCHeading') == 'Names and Identifiers':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'Record Title':
                        for info in subrec['Information']:
                            if info.get('Name') == 'Record Title':
                                names.append((cid, info['StringValue']))
                    elif subrec.get('TOCHeading') == 'Synonyms':
                        for subsubrec in subrec['Section']:
                            for info in subsubrec['Information']:
                                if info.get('Name') == 'MeSH Synonyms':
                                    for mesh_synonym in info['StringValueList']:
                                        mesh_synonyms.append((cid, mesh_synonym))
                                if info.get('Name') == 'Depositor-Supplied Synonyms':
                                    for syn in info['StringValueList']:
                                        depositor_supplied_synonyms.append((cid, syn))
            if rec.get('TOCHeading') == 'Pharmacology and Biochemistry':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'ATC Code':
                        for info in subrec['Information']:
                            if info.get('Name') == 'ATC Code':
                                atc = info['StringValue'].split()[0]
                                atcs.append((cid, atc))
            if rec.get('TOCHeading') == 'Biomolecular Interactions and Pathways':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'DrugBank Interactions':
                        for subsubrec in subrec['Section']:
                            if subsubrec.get('TOCHeading') in [
                                    'Target', 'Enzyme', 'Transporter', 'General Function']:
                                for info in subsubrec['Information']:
                                    if info.get('Name') in ['Target', 'Enzyme', 'Transporter']:
                                        name = ''.join(
                                            remove_links.split(info['StringValue'])[::2])
                                        bioentity = info['URL'].split('/')[-1]
                                        targets.append((cid, info['Name'], name, bioentity))
                                    elif info.get('Name') == 'General Function':
                                        g_function = ''.join(
                                            remove_links.split(info['StringValue'])[::2])
                                        general_functions.append((cid, g_function))

    return names, mesh_synonyms, depositor_supplied_synonyms, atcs, targets, general_functions


# %%
if __name__ == '__main__':
    data = [
        ('CNC(=O)c1ccc(cc1F)N2C(=S)N(C(=O)C2(C)C)c3ccc(C#N)c(c3)C(F)(F)F;'
         'CC1(C(=O)N(C(=S)N1c2ccc(c(c2)F)C(=O)NC)c3ccc(c(c3)C(F)(F)F)C#N)C'),
        ('CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F'),

    ]
    columns = ['SMILES or PubChem ID']
    df = pd.DataFrame(data, columns=columns)
    get_drug_features(df)
