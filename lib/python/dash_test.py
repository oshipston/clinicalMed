import dash
import flask
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import pandas as pd
import plotly.graph_objs as go
from inspect import getfullargspec

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import ast as ast
# import os
# os.chdir("/Users/Oliver/AnacondaProjects/clinicalMed/lib/python")
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, assets_url_path='/assets/')
# app.server
# dir(app)

fp = '../raw_data/taskPerf.tsv'
df = pd.read_csv(fp, sep='\t', index_col=0)
df = df[df.acc.notna()]
df.drop(columns=list(df.columns[9:len(df.columns)]), inplace=True)
df['task_diff_sum'] = [np.sum(ast.literal_eval(r)) for r in df.task_domain_ic_wm_cf]
df['task_diff_div'] = [3-ast.literal_eval(r).count(0) for r in df.task_domain_ic_wm_cf]
df['adj_acc'] = df.acc+(100-(df.acc))*((df.task_diff_sum*df.task_diff_div)/27)
df['mean_cohort_age'] = [np.mean(ast.literal_eval(r)) for r in df.age_range]

var_dict = {'acc': {'label': 'Task Accuracy',
                          'color': 'rgb(77, 124, 198)'},
            'adj_acc': {'label': 'Difficulty Adjusted Task Accuracy',
                         'color': 'rgb(155, 49, 152)'}}

class article:
    def __init__(self, text):
        self.id = []
        self.title = []
        self.title = []
        self.text = []
        self.plots = []
        self.dataframe = []

def autoscale(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))


def sig_fit(x, y, N=200):
    def sigmoid(p,x):
        x0,y0,c,k = p
        y = c / (1 + np.exp(-k*(x-x0))) + y0
        return y

    def residuals(p,x,y):
        return y - sigmoid(p,x)
    # x = autoscale(x)
    # y = autoscale(y)
    p_guess=(np.median(x), np.median(y), 1.0, 1.0)
    p_fit, _ = leastsq(
        residuals,p_guess,args=(x,y),full_output=0)
    m, c = np.polyfit(np.log(x), y, 1)
    x_fit = np.linspace(min(x), max(x), N)
    y_fit = sigmoid(p_fit, x_fit)
    return x_fit, y_fit


def exp_fit(x, y, N=200):
    m, c = np.polyfit(np.log(x), y, 1)
    x_fit = np.linspace(min(x), max(x), N)
    y_fit = m*np.log(x_fit) + c
    return x_fit, y_fit

def lin_fit(x, y, N=200):
    m, c = np.polyfit(np.log(x), y, 1)
    x_fit = np.linspace(min(x), max(x), N)
    y_fit = m*x_fit + c
    return x_fit, y_fit

def residuals(p,x,y):
    return y - sigmoid(p,x)

y_vars = ['adj_acc', 'acc']
y_var = 'adj_acc'
x_var = 'mean_cohort_age'

x = df[x_var]
y = df[y_var]

x_fit, y_fit = sig_fit(x, y)
fig = plt.figure(1)
ax = fig.add_subplot(111)


ax.scatter(x, y)
ax.plot(x_fit, y_fit)

fig


lin_traces = []
scatter_traces = []

y_var_dd = dcc.Dropdown(id='y_var_dd',
                        options=[{'label': 'Task Accuracy', 'value': 'acc'},
                                 {'label': 'Difficulty Adjusted Task Accuracy', 'value': 'adj_acc'}],
                        value=['adj_acc'],
                        multi=True)
study_var_cl = dcc.Checklist(id='study_var_cl',
                             className='checklist',
                             options=[dict(label=study, value=study) for study in np.unique(df.study_ref)],
                             values=list(np.unique(df.study_ref)))

def generate_table(dataframe, max_rows=10):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))]
    )

app.layout = html.Div(className='mainbody',
    children=[
        html.Div(className='container', children=[
            html.H1(children='Developing Executive Functions'),
            html.P('''Contrary to popular belief, Lorem Ipsum is not simply random text.
            It has roots in a piece of classical Latin literature from 45 BC, making it
            over 2000 years old. Richard McClintock, a Latin professor at Hampden-Sydney
            College in Virginia, looked up one of the more obscure Latin words, consectetur,
            from a Lorem Ipsum passage, and going through the cites of the word in classical
            literature, discovered the undoubtable source. Lorem Ipsum comes from sections 1.10.32
            and 1.10.33 of "de Finibus Bonorum et Malorum" (The Extremes of Good and Evil) by Cicero,
            written in 45 BC. This book is a treatise on the theory of ethics, very popular during
            the Renaissance. The first line of Lorem Ipsum, "Lorem ipsum dolor sit amet..", comes
            from a line in section 1.10.32.'''),
            # generate_table(df, max_rows=10),
            html.Div(className='plotContainer', children=[
                html.Div(className='ten columns', children=[
                    dcc.Graph(id='taskCompetence'),
                    y_var_dd
                ]),
                html.Div(className='two columns', children=[
                    study_var_cl
                ])
            ])
        ])
    ])

@app.callback(
    Output(component_id='taskCompetence', component_property='figure'),
    [Input(component_id='y_var_dd', component_property='value'),
     Input(component_id='study_var_cl', component_property='values')]
)
def update_figure(y_vars, studies):
    # print([arg for arg in getfullargspec(update_figure))
    traces = [go.Scatter(
        x=[-1,30],
        y=[50, 50],
        mode='lines',
        opacity=1,
        line={'color': 'rgb(200, 200, 200)',
              'width': 1,
              'dash': '5px'},
        name='Chance')]
    for y_var in y_vars:
        study_idx = [r in studies for r in df.study_ref]
        if sum(study_idx) !=0:
            x = df[x_var].values[study_idx]
            y = df[y_var].values[study_idx]
            ms = (df.cohort_N[study_idx]+2)/2
            try:
                x_exp_fit, y_exp_fit = exp_fit(x, y)
                traces.append(go.Scatter(
                    x=x_exp_fit,
                    y=y_exp_fit,
                    mode='lines',
                    opacity=1,
                    line = {'color': var_dict[y_var]['color']},
                    name='Exponential Fit: '+var_dict[y_var]['label']))
            except:
                print('Exponential fit not available')
                # raise Warning('Exponential fit not available')
            try:
                x_sig_fit, y_sig_fit = sig_fit(x, y)
                traces.append(go.Scatter(
                    x=x_sig_fit,
                    y=y_sig_fit,
                    mode='lines',
                    opacity=1,
                    line = {'color': var_dict[y_var]['color']},
                    name='Sigmoidal Fit: '+var_dict[y_var]['label']))
            except:
                print('Sigmoidal fit not available')
                # raise Warning('Sigmoidal fit not available')
            traces.append(go.Scatter(
                x=x,
                y=y,
                text=df.study_ref[study_idx]+'<br>'+df.task[study_idx],
                mode='markers',
                opacity=0.3,
                name=var_dict[y_var]['label'],
                marker={
                    'size': ms,
                    'color': var_dict[y_var]['color'],
                    'line': {'width': 0.5, 'color': 'white'}
                }))
    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'title': 'Age (Years)', 'range': [min(df.mean_cohort_age)-2, max(df.mean_cohort_age)+2]},
            yaxis={'title': 'Task Competence', 'range': [40, 110]},
            margin={'l': 60, 'b': 60, 't': 10, 'r': 10},
            legend={'x': 0.1, 'y': 0.05},
            hovermode='closest'
        )
    }

if __name__ == '__main__':
    app.run_server(debug=True, port=8081)
