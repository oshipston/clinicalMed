import os
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
# os.chdir("/Users/Oliver/AnacondaProjects/clinicalMed/lib/python")

def logit_sigmoid(x, arg):
  return np.array(1/(1 + np.exp(-arg[0]*(x-arg[1]))))

def linear(x, arg):
   return np.array(arg[0]*x + arg[1])

def cubic(x, arg):
    return np.array(arg[0]*(x**3) + arg[1]*(x**2) + arg[2]*(x**3) + arg[3])

def piecewise(x, funcDict, split=[0.5]):
    assert len(funcDict) == (len(split)+1), print('Too few functions in func dict for no. of splits')
    trace = np.array([])
    start = 0
    sIdx = [0, int(np.array(split)*len(x)), len(x)]
    for s, k in enumerate(funcDict):
        func = funcDict[k]['func']
        pwTrace = func(x[sIdx[s]:sIdx[s+1]], funcDict[k]['args'])
        print(len(pwTrace))
        trace = np.append(trace, np.array(pwTrace))
    print(len(trace))
    return np.array(trace)

age_range = [-40, 1040]
x = np.linspace(age_range[0], age_range[1], 5000)

traces = []
traces.append(logit_sigmoid(x, [0.02, 300]))
traces.append(linear(x, [0.0005, .2]))
traces.append(cubic((x-600)/400, [1, 1, 1, 0.5]))
pwDict = {'linear': {'func': linear, 'args': [0.5, 2]},
          'cubic': {'func': cubic, 'args': [0.5, 2, 1, 0]}}
traces.append(piecewise(x, pwDict))
tt = piecewise(x, pwDict)
figParams = {'num': 1,
             'figsize': (10, 5),
             'dpi': 100,
             'frameon': False}

fig = plt.figure(**figParams)
ax = fig.add_subplot(1, 1, 1)
for trace in traces:
    ax.plot(x, trace, linewidth=1)

ax.plot(x, (cubic((x-600)/400, [1, 1, 1, 0.5])+linear(x, [0.005, 2])), linewidth=1)

ax.set_xlim(age_range)
ax.set_ylim([0, 1])

fig
fig.suptitle('Potential Developmental Trajectories': '+time_created)
fig.savefig('../devTrajectories.pdf', dpi=300,
            format='pdf', pad_inches=0.1, bbox_inches='tight')
plt.close()
