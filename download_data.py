# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 14:25:46 2015

@author: strokach
"""
import os
import os.path as op
from IPython.display import display

from synapseclient import Activity
from synapseclient import Entity, Project, Folder, File
from synapseclient import Evaluation, Submission, SubmissionStatus
from synapseclient import Wiki

import synapseclient


# %%
syn = synapseclient.Synapse()
syn.login('ostrokach', ''.join(reversed(['er', 'at', 'Sk'])))

try:
    BASE_DIR = op.abspath(op.dirname(__file__))
except NameError:
    BASE_DIR = os.getcwd()
    assert op.isfile(op.join(BASE_DIR, 'download_data.py'))


# %%
def download_data(entity_id, path):
    print('entity_id: {}'.format(entity_id))
    print('path: {}'.format(path))
    #
    os.makedirs(path, exist_ok=True)
    try:
        syn.get(entity_id, downloadLocation=path)
    except KeyError as e:
        print('The following error occurred: {}'.format(e))
    #
    results = syn.query('SELECT * FROM entity WHERE parentId=="{}"'.format(entity_id))['results']
    # display(pd.DataFrame(results))
    for result in results:
        new_path = op.join(path, result['entity.name'].lower().replace(' ', '_'))
        download_data(result['entity.id'], new_path)


# %%
entity_id = 'syn4231880'
path = BASE_DIR
download_data(entity_id, path)
