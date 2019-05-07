#!/usr/bin/env python3.6
import os
import buchheim as buc
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

imp.reload(buc)
imp.reload(ctf)

os.chdir("/Users/Oliver/AnacondaProjects/clinicalMed/")
# os.getcwd()


class ICD_data_sess:
    # Function assessing the health and leafiness of currently downloaded ICD Tree.
    # Input:
    # 'class_spec' = 1st Generation disease class specifier to reduce API burden
    def __init__(self):
        # get the OAUTH2 token
        self.cfg = ctf.loadJSON('config.json')
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


def refresh_adj_list():
    dat_sess = ICD_data_sess()
    dat_sess.load_local_dat()
    dat_sess.acquire_from_parent_queue()
    return dat_sess


class node:
    def __init__(self, dat):
        try:
            self.parents = ast.literal_eval(dat.parents)
        except(KeyError):
            self.parents = []
        self.children = ast.literal_eval(dat.children)
        self.id = dat.id
        self.title = dat.title

    # def __str__(self): return "%s: x=%s mod=%s" % (self.tree, self.x, self.mod)
    def __str__(self): return 'Node id:'+str(self.id)

    def __repr__(self): return self.__str__()


class ICD_adjacency:
    def __init__(self, adj_path):
        self.data_frame = pd.read_csv(adj_path, sep='\t', index_col=0)


class DrawTree(object):
    def __init__(self, tree, adj_df, parents=None, depth=0, number=1):
        self.id = tree.id
        self.title = tree.title
        self.x = -1.
        self.y = depth
        self.tree = tree
        child_array = []
        i = 0
        for c in tree.children:
            try:
                child_array.append(DrawTree(node(adj_df.loc[str(c)]), adj_df,
                                            self, depth+1, i+1))
                i += 1
            except(KeyError):
                child_array = []
        self.children = child_array
        self.parents = parents
        # if type(self.parents) is not DrawTree:
        #     print('damn')
        #     print(type(self.parents))
        #     print(self.tree.id)
        self.thread = None
        self.mod = 0
        self.ancestor = self
        self.change = self.shift = 0
        self.lmost_sibling = None
        # This is the number of the node in its group of siblings 1..n
        self.number = number

    def left(self):
        return self.thread or len(self.children) and self.children[0]

    def right(self):
        return self.thread or len(self.children) and self.children[-1]

    def lbrother(self):
        n = None
        if self.parents:
            for node in self.parents.children:
                if node == self:
                    return n
                else:
                    n = node
        return n

    def get_lmost_sibling(self):
        if not self.lmost_sibling and self.parents and self != \
                self.parents.children[0]:
                    self.lmost_sibling = self.parents.children[0]
        return self.lmost_sibling
    # lmost_sibling = property(get_lmost_sibling)

    # def __str__(self): return "%s: x=%s mod=%s" % (self.tree, self.x, self.mod)
    def __str__(self): return str(self.tree)+': x='+str(self.x)+' mod='+str(self.mod)

    def __repr__(self): return self.__str__()


def firstwalk(v, distance=1.):
    if len(v.children) == 0:
        if v.get_lmost_sibling():
            v.x = v.lbrother().x + distance
        else:
            v.x = 0.
    else:
        default_ancestor = v.children[0]
        for w in v.children:
            firstwalk(w)
            # print('v = '+str(v))
            # print('w = '+str(w))
            default_ancestor = apportion(w, default_ancestor, distance)
        # print("finished v =", v.tree, "children")
        execute_shifts(v)

        midpoint = (v.children[0].x + v.children[-1].x) / 2

        ell = v.children[0]
        arr = v.children[-1]
        w = v.lbrother()
        if w:
            v.x = w.x + distance
            v.mod = v.x - midpoint
        else:
            v.x = midpoint
    return v


def apportion(v, default_ancestor, distance):
    w = v.lbrother()
    if w is not None:
        # print('w = '+str(w))
        # in buchheim notation:
        # i == inner; o == outer; r == right; l == left; r = +; l = -
        vir = vor = v
        vil = w
        vol = v.get_lmost_sibling()
        sir = sor = v.mod
        sil = vil.mod
        sol = vol.mod
        while vil.right() and vir.left():
            vil = vil.right()
            vir = vir.left()
            vol = vol.left()
            vor = vor.right()
            vor.ancestor = v
            shift = (vil.x + sil) - (vir.x + sir) + distance
            if shift > 0:
                move_subtree(ancestor(vil, v, default_ancestor), v, shift)
                sir = sir + shift
                sor = sor + shift
            sil += vil.mod
            sir += vir.mod
            sol += vol.mod
            sor += vor.mod
        if vil.right() and not vor.right():
            vor.thread = vil.right()
            vor.mod += sil - sor
        else:
            if vir.left() and not vol.left():
                vol.thread = vir.left()
                vol.mod += sir - sol
            default_ancestor = v
    return default_ancestor


def move_subtree(wl, wr, shift):
    subtrees = wr.number - wl.number
    # print(str(wl.tree)+' is conflicted with '+str(wr.tree) +
    #       ' moving '+str(subtrees)+' shift '+str(shift))
    # print wl, wr, wr.number, wl.number, shift, subtrees, shift/subtrees
    # print(str(subtrees))
    wr.change -= shift / subtrees
    wr.shift += shift
    wl.change += shift / subtrees
    wr.x += shift
    wr.mod += shift


def execute_shifts(v):
    shift = change = 0
    for w in v.children[::-1]:
        # print("shift:", w, shift, w.change)
        w.x += shift
        w.mod += shift
        change += w.change
        shift += w.shift + change


def ancestor(vil, v, default_ancestor):
    # the relevant text is at the bottom of page 7 of
    # "Improving Walker's Algorithm to Run in Linear Time" by Buchheim et al, (2002)
    # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.16.8757&rep=rep1&type=pdf
    if vil.ancestor in v.parents.children:
        return vil.ancestor
    else:
        return default_ancestor


def secondwalk(v, m=0, depth=0, min=None):
    v.x += m
    v.y = depth

    if min is None or v.x < min:
        min = v.x

    for w in v.children:
        min = secondwalk(w, m + v.mod, depth+1, min)

    return min


def draw_bezier(node1, node2, curve_ydim=0, curviness=0.5):
    if node1[curve_ydim] < node2[curve_ydim]:
        start = np.array(node1)
        end = np.array(node2)
    else:
        start = np.array(node2)
        end = np.array(node1)

    h = end-start
    x_buf = curviness*(h[0])
    y_buf = curviness*(h[1])

    if not curve_ydim:
        start_buf = start + np.array([x_buf, 0])
        end_buf = end + np.array([-x_buf, 0])
    else:
        start_buf = start + np.array([0, y_buf])
        end_buf = end + np.array([0, -y_buf])

    verts = [start,
             start_buf,
             end_buf,
             end]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4]

    path = Path(verts, codes)
    return path


def plot_tree(node_tree, fp, style='v'):
    def node_draw(node, ax, style=style):
        text_params = {'fontsize': 2,
                       'color': 'xkcd:black'}
        if style == 'v':
            coords = [node.x, -node.y]
            text_params.update({'horizontalalignment': 'center',
                                'verticalalignment': 'bottom',
                                'rotation': 'vertical'})
        elif style == 'h':
            coords = [node.y, node.x]
            text_params.update({'horizontalalignment': 'right',
                                'verticalalignment': 'center'})
        elif style == 'r':
            print('Radial function not yet implemented')

        ax.scatter(coords[0], coords[1], s=0.5, marker='o', alpha=1, c='xkcd:black')
        s = txtwrp.fill(node.title, 30)
        t = ax.text(coords[0], coords[1], s, **text_params)
        t.set_path_effects([patheffects.Stroke(linewidth=0.3, foreground='xkcd:white'),
                            patheffects.Normal()])
        for child in node.children:
            node_draw(child, ax)

    def conn_draw(node, ax, lin_style='bezier', style=style):
        # print(str(node.x)+'  '+str(node.y))
        for child in node.children:
            conn_draw(child, ax)
            if lin_style == 'linear':
                ax.add_line(Line2D([child.x, node.x], [-child.y, -node.y], linewidth=0.4,
                            alpha=0.5, color='xkcd:black'))
            elif lin_style == 'bezier':
                if style == 'v':
                    path = draw_bezier([child.x, -child.y], [node.x, -node.y], curve_ydim=1)
                elif style == 'h':
                    path = draw_bezier([child.y, child.x], [node.y, node.x], curve_ydim=0)
                patch = patches.PathPatch(path, facecolor='none', lw=0.4, alpha=0.5)
                ax.add_patch(patch)

    # fig_size = 25/node_tree.depth

    if style == 'v':
        fig_params = {'num': 1,
                      'figsize': (10, 3),
                      'dpi': 200,
                      'frameon': False}
    elif style == 'h':
        fig_params = {'num': 1,
                      'figsize': (3, 10),
                      'dpi': 200,
                      'frameon': False}

    fig = plt.figure(**fig_params)
    ax = fig.add_subplot(1, 1, 1)
    node_draw(node_tree, ax)
    conn_draw(node_tree, ax)
    ax.set_ylabel('Generations')
    # ax.set_title('ICD Tree Test')
    ax.set_axis_off()
    fig.savefig(fp, dpi=100,
                format='pdf', pad_inches=0, bbox_inches='tight')
    plt.close()


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
# dat_sess.adj_list.loc[str(384289569)]
root = node(dat_sess.adj_list.loc[str(node_dict['cardiac_arrhythmias'])])
node_tree = DrawTree(root, dat_sess.adj_list)
dt = firstwalk(node_tree)
min = secondwalk(dt)

plot_tree(node_tree, 'output/verttest.pdf', style='v')
