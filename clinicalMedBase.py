#!/usr/bin/env python3.6
import clinicalMedBaseFunctions as cmb
import requests
import pandas as pd
import importlib as imp
import time
import numpy as np
import matplotlib.pyplot as plt
import textwrap as txtwrp
import os
imp.reload(cmb)

# os.chdir("/Users/Oliver/AnacondaProjects/clinicalMed/")
# os.getcwd()

# get the OAUTH2 token
cfg = cmb.loadJSON('config.json')
token_endpoint = 'https://icdaccessmanagement.who.int/connect/token'
client_id = cfg['client_id']
client_secret = cfg['client_secret']
scope = 'icdapi_access'
grant_type = 'client_credentials'

# set data to post
payload = {'client_id': client_id,
           'client_secret': client_secret,
           'scope': scope,
           'grant_type': grant_type}

# make request
sess = requests.Session()
sess.verify = '/anaconda3/envs/clinicalTree/ssl/cacert.pem'

token_req = sess.post(token_endpoint, data=payload).json()
token = token_req['access_token']

# uri = 'https://id.who.int/icd/entity/search?q={Neuro}'
ICD_uri = 'https://id.who.int/icd/entity/'
neuro_id = '1296093776'
neuro_uri = ICD_uri + neuro_id  # Neurological Disease
# HTTP header fields to set
headers = {'Authorization':  'Bearer '+token,
           'Accept': 'application/json',
           'Accept-Language': 'en'}

# make request for NeuroCatData (ncd)
neuro_dat = cmb.ICD11_result(json_dat=sess.get(neuro_uri, headers=headers).json())
neuro_dat.storage_dict.update({'generation': 0})

df = pd.DataFrame(columns=['id', 'title', 'parents_id', 'children_id', 'generation'])

# Initialise dataframe with matriarch
df.loc[neuro_dat.id] = neuro_dat.storage_dict
parent_stack = [[neuro_id, 0]]
generation_limit = 2
i = 0
parent_stack
while parent_stack and (i < 1000):
    parent = parent_stack.pop()
    parent_id = parent[0]
    parent_gen = parent[1]
    if parent_gen < generation_limit:  # If parent is not at generational limit
        parent_dat = cmb.ICD11_result(json_dat=sess.get(ICD_uri+parent_id, headers=headers).json())
        for child_id in parent_dat.children_id:
            if child_id not in df['id']:
                # If not already stored request from API
                child_dat = cmb.ICD11_result(json_dat=sess.get(ICD_uri+child_id,
                                                               headers=headers).json())
                # Add generation to storage_dict
                child_dat.storage_dict.update({'generation': parent_gen+1})
                # Store storage_dict in dataframe
                df.loc[child_dat.id] = child_dat.storage_dict
                # Add child to parent_stack
                parent_stack.append([child_id, parent_gen+1])
                i += 1
            # time.sleep(1+np.random.uniform())  # Add some time between API calls

parent_stack
adj_path = 'raw_data/ICD_adjacency_list.tsv'
df.to_csv(adj_path, sep='\t')
test_df = pd.DataFrame.from_csv(adj_path, sep='\t')  # From database rather than API.

x = []
y = []
s = []
for gen in df['generation'].unique():
    node_dat = df[df['generation'] == gen]
    n_node = len(node_dat)
    y = y + list([gen*-2]*n_node)
    x = x + list((np.linspace(0, n_node-1, n_node)-((n_node-1)/2))*1)
    s = s + ([txtwrp.fill(raw_str, width=20) for raw_str in node_dat['title']])
x
fig = plt.figure(num=1, figsize=(12, 2), dpi=100, frameon=False)
ax = fig.add_subplot(1, 1, 1)
for i in range(len(x)):
    ax.text(x[i], y[i], s[i], ha='center', fontsize=3, va='center')

ax.set_ylim([min(y)-1, max(y)+1])
ax.set_xlim([min(x)-1, max(x)+1])
fig
