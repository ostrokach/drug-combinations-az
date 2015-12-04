# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 14:25:46 2015

@author: strokach


This file can only be ran using Python 2.7!!!


"""
import os
import os.path as op

#from synapseclient import Activity
#from synapseclient import Entity, Project, Folder, File
#from synapseclient import Evaluation, Submission, SubmissionStatus
#from synapseclient import Wiki

import synapseclient
syn = synapseclient.Synapse()

syn.login()


#%%
def _get_children(entity_id, path):
    print

    try:
        os.makedirs(path)
    except OSError:
        pass

    item = syn.get(entity_id, downloadLocation=path)
    print('results:', item)
    print(type(item))

    if isinstance(item, synapseclient.Folder):
        path = op.join(path, item.properties['name'].lower().replace(' ', '_'))
    print('path:', path)

    results = syn.query('SELECT id, name FROM entity WHERE parentId=="{}"'.format(entity_id))
    print('results:', results)

    if 'results' in results:
        for entity in results['results']:
            _get_children(entity['entity.id'], path)


#%%
_get_children('syn4923176', '.')
