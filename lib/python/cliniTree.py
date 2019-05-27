#!/usr/bin/env python3.6
import os
import cliniTreeFunctions as ctf
import requests
import pandas as pd
import importlib as imp
import time
import numpy as np
import matplotlib.pyplot as plt
import textwrap as txtwrp
import json
from datetime import datetime
import ast
import csv
import sys
from matplotlib.lines import Line2D
from scipy.interpolate import BSpline, UnivariateSpline, CubicSpline
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
import plotly
import inspect

imp.reload(ctf)

os.chdir("/Users/Oliver/AnacondaProjects/clinicalMed/lib/python")
# os.getcwd()


class ICD_data_sess:
    # Function assessing the health and leafiness of currently downloaded ICD Tree.
    # Input:
    # 'class_spec' = 1st Generation disease class specifier to reduce API burden
    def __init__(self):
        # get the OAUTH2 token
        self.cfg = ctf.loadJSON('../config.json')
        self.client_id = self.cfg['client_id']
        self.client_secret = self.cfg['client_secret']
        self.scope = self.cfg['scope']
        self.grant_type = self.cfg['grant_type']
        self.sess_start_time = datetime.now().timestamp()
        # set data to post
        payload = {'client_id': self.client_id,
                   'client_secret': self.client_secret,
                   'scope': self.scope,
                   'grant_type': self.grant_type}

        # Open requests HTTP session
        self.sess = requests.Session()
        # Add local SSL certificate to session.
        self.sess.verify = '/anaconda3/envs/clinicalTree/ssl/cacert.pem'
        # Request access token with creds.
        token_endpoint = 'https://icdaccessmanagement.who.int/connect/token'
        try:
            token_req = self.sess.post(token_endpoint, data=payload).json()
            token = token_req['access_token']

            # Set HTTP header fields for future gets
            self.headers = {'Authorization':  'Bearer '+token,
                            'Accept': 'application/json',
                            'Accept-Language': 'en'}
        except:
            print('Connection unsuccessful.')
            print('Offline session initilised.')
        pd.options.display.max_columns = None
        pd.options.display.max_rows = 100

    def load_local_dat(self):
        try:
            self.adj_list = pd.read_csv(self.cfg['adj_path'], sep='\t', index_col=0)
            self.adj_list.id = self.adj_list.id.astype(int)
            self.parent_queue = ctf.read_parent_tsv(self.cfg['queue_path'])
            print('Found local ICD adjacency, loading prior list...')
        except(FileNotFoundError):
            print('ICD adjacency list and parent queue not found')
            self.adj_list, self.parent_queue = self.init_adj_list()

    def req(self, uri):
        return self.sess.get(uri, headers=self.headers).json()

    def init_adj_list(self):
        # State relevant URIs for requsting different sections of the ICD or All.
        # uri = 'https://id.who.int/icd/entity/search?q={Neuro}'

        # make request for RootData (ncd)
        root_dat = ctf.ICD11_result(json_dat=self.req(self.cfg['ICD_root_uri']))
        root_dat.storage_dict.update({'generation': 0})  # Set as generation 0
        # Initialise adjacency dataframe with root.
        self.adj_list = pd.DataFrame(columns=['id', 'title', 'parents', 'children',
                                     'generation', 'date_created'])
        self.adj_list.loc[root_dat.id] = root_dat.storage_dict
        # Remove id due to use in URIs later, already used as index for tree building
        self.adj_list.at['entity', 'id'] = 0
        # Initialiase parent queue with root.
        self.parent_queue = [[root_dat.id, 0]]

    def plot_summary(self):
        time_created = datetime.now().strftime('%d/%m/%Y %H:%M')
        timestamps = self.adj_list['date_created']
        cumsum_nodes = []
        for stamp in timestamps:
            cumsum_nodes.append(sum(self.adj_list['date_created'] < stamp))
        figParams = {'num': 1,
                     'figsize': (10, 5),
                     'dpi': 100,
                     'frameon': False}
        fig = plt.figure(**figParams)
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(timestamps[50:len(timestamps)], cumsum_nodes[50:len(cumsum_nodes)],
                   s=5, marker='.')
        fig.suptitle('ICD Node Acquisition: '+time_created)
        fig.savefig('logs/ICD_acquisition.pdf', dpi=300,
                    format='pdf', pad_inches=0.1, bbox_inches='tight')
        plt.close()

    def save_local_dat(self):
        # Save as /tsv.s
        self.adj_list.to_csv(self.cfg['adj_path'], sep='\t')
        ctf.save_tsv(self.parent_queue, self.cfg['queue_path'])
        self.plot_summary()

    def acquire_from_parent_queue(self, attempt_limit=1000, gen_limit=10):
        # Iterate over parents until node count, generation limit or all parents have been added.
        self.load_local_dat()
        node_limit = 10000
        self.start_parent_no = len(self.parent_queue)
        self.start_node_no = len(self.adj_list)
        ICD_root_uri = self.cfg['ICD_root_uri']
        n_i = 0
        a_i = 0
        while self.parent_queue and (a_i < attempt_limit) and (n_i < node_limit):
            # Opo parent from list and allocate
            parent = self.parent_queue.pop()
            parent_id = parent[0]
            if parent_id == 'entity':  # Special case for ICD root
                parent_id = ''
            parent_gen = parent[1]
            # If parent not exceeding generation limit...
            if parent_gen < gen_limit:
                # Request
                parent_uri = ICD_root_uri+str(parent_id)
                parent_dat = ctf.ICD11_result(json_dat=self.req(parent_uri))
                for child_id in parent_dat.children:
                    # If not already stored, request from API..
                    if child_id not in self.adj_list['id']:
                        child_uri = ICD_root_uri+str(child_id)
                        child_dat = ctf.ICD11_result(json_dat=self.req(child_uri))
                        # Add generation to storage_dict and store in adjacency list...
                        child_dat.storage_dict.update({'generation': parent_gen+1})
                        self.adj_list.loc[child_dat.id] = child_dat.storage_dict
                        # Add child to parent_queue if any children present.
                        if child_dat.children:
                            self.parent_queue.append([child_id, parent_gen+1])
                        n_i += 1
                        # time.sleep(0.5+np.random.uniform())  # Add some time between API calls
            else:
                self.parent_queue.insert(0, parent)
            a_i += 1
        self.duplicate_no = self.adj_list.duplicated(subset='id').sum()
        pre_dup = len(self.adj_list)
        self.adj_list.drop_duplicates(subset='id', inplace=True)
        post_dup = len(self.adj_list)
        print(str(pre_dup-post_dup)+' duplicates removed.')
        self.end_parent_no = len(self.parent_queue)
        self.end_node_no = len(self.adj_list)
        print('Saving newly acquired data:')
        self.save_local_dat()
        print(str(self.end_node_no-self.start_node_no)+' nodes added to adjacency.')
        print(str(self.end_parent_no-self.start_parent_no)+' parents removed from parent queue.')
        print('Total no. of nodes in adjacency = '+str(len(self.adj_list)))
        print('Total no. of parents in queue = '+str(len(self.parent_queue)))

    def condition_search(self, search_str):
        search_idx = []
        for i, r in enumerate(self.adj_list.title):
            if search_str in r:
                search_idx.append(i)
        return search_idx

    # def count_descendants(self, root_id, baton):
    #     print(root_id)
    #     print('Root_id = ' + str(root_id))
    #     root = self.adj_list.loc[str(root_id)]
    #     # print(root)
    #     baton['N_children'] = baton['N_children']+len(root.children)
    #     for c in root.children:
    #         count_descendants(self, c, baton)
    #     return baton




node_dict = {'neuro': 1296093776,
             'cardio': 426429380,
             'haem': 1766440644,
             'gastro': 1256772020,
             'resp': 197934298,
             'derm': 1639304259,
             'parkinsonism': 2024168133,
             'movement_disorders': 384289569,
             'cardiac_arrhythmias': 1457291912,
             'onc': 1630407678,
             'immune': 1954798891,
             'endo': 21500692,
             'mental': 224423054,
             'sleep': 274880002,
             'opth': 868865928,
             'audio': 1218729044,
             'msk': 1473673350,
             'gum': 30659757,
             'sex': 577470983,
             'obs': 714000734,
             'perinate': 1306203631,
             'develop': 223744320,
             'injury_poison': 435227771,
             'external': 850137482,
             'traditional': 718687701}

dat_sess = ICD_data_sess()
dat_sess.load_local_dat()
baton = {}
baton['N_children'] = 0

def export_tree(root, adj_list, style='h'):
    root = ctf.node(root)
    node_tree = ctf.DrawTree(root, adj_list)
    dt = ctf.firstwalk(node_tree)
    min = ctf.secondwalk(dt)
    fn = root.title
    ctf.plot_tree(node_tree, '../output/'+fn+'_tree.pdf', style=style)

export_tree(root=dat_sess.adj_list.loc[str(node_dict['neuro'])],
            adj_list = dat_sess.adj_list,
            style='h')



fp = '../raw_data/taskPerf.tsv'
df = pd.read_csv(fp, sep='\t', index_col=0)
df = df[df.acc.notna()]
df.drop(columns=list(df.columns[9:len(df.columns)]), inplace=True)
df['task_diff_sum'] = [np.sum(ast.literal_eval(r)) for r in df.task_domain_ic_wm_cf]
df['task_diff_div'] = [3-ast.literal_eval(r).count(0) for r in df.task_domain_ic_wm_cf]
df['mean_cohort_age'] = [np.mean(ast.literal_eval(r)) for r in df.age_range]
df.acc+(100-(df.acc))*((df.task_diff_sum*df.task_diff_div)/27)

df['adj_acc'] = df.acc+(100-(df.acc))*((df.task_diff_sum*df.task_diff_div)/27)

def exp_fit(x, y, N=200):
    m, c = np.polyfit(np.log(x), y, 1)
    x_fit = np.linspace(min(x), max(x), N)
    y_fit = m*np.log(x_fit) + c
    return x_fit, y_fit

def lin_fit(x, y, N=200):
    tt = inspect.getfullargspec(lin_fit)
    m, c = np.polyfit(np.log(x), y, 1)
    x_fit = np.linspace(min(x), max(x), N)
    y_fit = m*x_fit + c
    return x_fit, y_fit, tt

studies = list(np.unique(df.study_ref))
[r in studies for r in df.study_ref]
plot_y_vars = ['adj_comp', 'acc']
plot_x_var = ['mean_cohort_age']

plot_dat = dash_scatter_dat(df, plot_x_var, plot_y_vars)

figParams = {'num': 1,
             'figsize': (10, 10),
             'dpi': 100,
             'frameon': False}
fig = plt.figure(**figParams)
ax = fig.add_subplot(1, 1, 1)

x, y, tt = lin_fit([1, 2, 3], [1, 2, 3], 100)

print(tt)
    ax.scatter(dat.age_range, dat.acc, marker='.', s=dat.cohort_N*5, alpha=0.4)


ax.set_xlabel('Age (years)')
ax.set_ylabel('Normalised Task Competence')
fig
