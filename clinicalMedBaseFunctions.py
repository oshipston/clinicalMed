#!/usr/bin/env python3.6
from json import load as jsonLoad


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
        self.id = json_dat['@id'][json_dat['@id'].rfind('/')+1:len(json_dat['@id'])]
        self.parents_uri = json_dat['parent']
        self.parents_id = [p[p.rfind('/')+1:len(p)] for p in json_dat['parent']]
        try:
            self.children_uri = json_dat['child']
            self.children_id = [c[c.rfind('/')+1:len(c)] for c in json_dat['child']]
        except(KeyError):
            self.children_uri = []
            self.children_id = []
        self.storage_dict = {'title': self.title,
                             'id': self.id,
                             'parents_id': self.parents_id,
                             'children_id': self.children_id}
