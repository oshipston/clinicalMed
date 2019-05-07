#!/usr/bin/env python3.6
from json import load as jsonLoad
from datetime import datetime
import pandas as pd
import csv
import ast


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
