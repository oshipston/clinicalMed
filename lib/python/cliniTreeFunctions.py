#!/usr/bin/env python3.6
import csv
import ast
import sys
from json import load as jsonLoad
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.patheffects as patheffects
import textwrap as txtwrp

def loadJSON(fname):
    # Load configuration
    f = open(fname)  # Open config file...
    cfg = jsonLoad(f)  # Load data...
    f.close()  # Close config file...
    return cfg


class ICD11_result:
    """A wrapper for a json return from the ICD 11 entity"""
    def __init__(self, json_dat):
        self.raw_json = json_dat
        self.title = json_dat['title']['@value']
        self.uri = json_dat['@id']
        self.datetime_created = round(datetime.now().timestamp(), -1)
        try:
            self.id = int(json_dat['@id'][json_dat['@id'].rfind('/')+1:len(json_dat['@id'])])
        except(ValueError):
            self.id = json_dat['@id'][json_dat['@id'].rfind('/')+1:len(json_dat['@id'])]
            # ?? make self.id = 0
        try:
            self.parents_uri = json_dat['parent']
            try:
                self.parents = [int(p[p.rfind('/')+1:len(p)]) for p in json_dat['parent']]
            except(ValueError):
                self.parents = [p[p.rfind('/')+1:len(p)] for p in json_dat['parent']]
                # ?? make self.parents = None/False
        except(KeyError):
            self.parents_uri = []
            self.parents = []

        try:
            self.children_uri = json_dat['child']
            self.children = [int(c[c.rfind('/')+1:len(c)]) for c in json_dat['child']]
        except(KeyError):
            self.children_uri = []
            self.children = []
        self.storage_dict = {'title': self.title,
                             'id': self.id,
                             'parents': self.parents,
                             'children': self.children,
                             'date_created': self.datetime_created}


def save_tsv(dat, fp):
    with open(fp, 'w', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(dat)


def read_parent_tsv(fp):
    with open(fp) as f_dat:
        read_dat = csv.reader(f_dat, delimiter="\t")
        output = [r for r in read_dat]
        output = output[0]
        output = [ast.literal_eval(r) for r in output]
    return output

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


def firstwalk(v, distance=2.):
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
