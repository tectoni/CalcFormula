
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday, 17 Jan 2025, 14:48

@author: pappel
"""
import matplotlib.pyplot as plt
import mpltern
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools

def gt_tern(df, savepic):
    fig = go.Figure()   # create a plain plotly figure object

    # define a list with some plotly colors (some are CSS colors)
    colors = ['rgba(50, 200, 20, 0.6)', 'rgba(250, 50, 0, 0.6)', 'rgba(10, 200, 200, 0.6)', 'yellow', 'white', 'tomato', 'rgba(50,10,180, 0.5)', 'olive', 'olive', 'cyan']
    # define a list with some plotly markers (symbols)
    markers = ['diamond', 'circle', 'square', 'cross', 'x', 'circle-x-open', 'star-diamond', 'circle', 'square', 'triangle-se', 'star-triangle-down'   ]
    # this is to extract the sample name from the df, it relates to the first 5 characters in the Comment column
    list_of_samples = df['sample'].str[:5].unique()
    # cycle through colors and markers to prevent from running out of markers
    color_cycle = itertools.cycle(colors)
    marker_cycle = itertools.cycle(markers)

  # now create tuples of corresponding parameters and iterate them
    for s,c,m in zip(list_of_samples, color_cycle, marker_cycle):
        sub_df = df.loc[df['sample'].str[:5] == s]
        # Add scatter trace with medium sized markers
        # this plots the data with the markers
        fig.add_trace(
            go.Scatterternary(
                a=sub_df['Ca'],
                b=sub_df['Fe2'],
                c=sub_df['Mg'],
                mode='markers',
                marker=dict(
                    color=c,
                    size=6,
                    line=dict(
                        color='black',
                        width=1
                    ),
                    symbol=m
                ),
                showlegend=True,
                # the name is for the legend of the data series
                name=s
            )
        )
    # Add title and layout adjustments
    fig.update_layout(
        title='Garnet',
        ternary=dict(
            aaxis_title='Ca',
            baxis_title='Fe2',
            caxis_title='Mg'
        )
    )

    if savepic:
        fig.write_image('./plots/gtternary.pdf')
    fig.show()
    return

def fsp_tern(df, savepic):
    fig = px.scatter_ternary(df, a = df['K'], b=df['Na'], c=['Ca'])
    if savepic:
        fig.write_image('./plots/fspternary.pdf')
    fig.show()
    return

def px_tern(df):  # here using matplotlib
    ax = plt.subplot(projection='ternary')
    ax.scatter(df['Ca'], df['Fe'], df['Mg'])
    if savepic:
        plt.savefig('./plots/cpxtern.pdf')
    plt.show()
    return

def px_quad():
    return

def bt_xfe():
    return

def sodicAmphs(df, savepic):
    fig = go.Figure()   # create a plain plotly figure object

    # define a list with some plotly colors (some are CSS colors)
    colors = ['rgba(50, 200, 20, 0.6)', 'rgba(250, 50, 0, 0.6)', 'rgba(10, 200, 200, 0.6)', 'yellow', 'white', 'tomato', 'rgba(50,10,180, 0.5)', 'olive', 'olive', 'cyan']
    # define a list with some plotly markers (symbols)
    markers = ['diamond', 'circle', 'square', 'cross', 'x', 'circle-x-open', 'star-diamond', 'circle', 'square', 'triangle-se', 'star-triangle-down'   ]
    # this is to extract the sample name from the df, it relates to the first 5 characters in the Comment column
    list_of_samples = df['sample'].str[:5].unique()
    # cycle through colors and markers to prevent from running out of markers
    color_cycle = itertools.cycle(colors)
    marker_cycle = itertools.cycle(markers)

  # now create tuples of corresponding parameters and iterate them
    for s,c,m in zip(list_of_samples, color_cycle, marker_cycle):
        sub_df = df.loc[df['sample'].str[:5] == s]
        # Add scatter trace with medium sized markers
        # this plots the data with the markers
        fig.add_trace(
            go.Scatter(
                x=sub_df['Si'],
                y=(sub_df['Fe2'] + sub_df['Mg']) / sub_df['Mg'],
                mode='markers',
                marker=dict(
                    color=c,
                    size=6,
                    line=dict(
                        color='black',
                        width=1
                    ),
                    symbol=m
                ),
                showlegend=True,
                # the name is for the legend of the data series
                name=s
            )
        )
    # Add title and layout adjustments
    fig.update_layout(
        xaxis = dict(autorange="reversed"),
        title='Sodic Amphiboles',
        xaxis_title='Si',
        yaxis_title='(Fe2+Mg)/Mg'
    )
    

    if savepic:
        fig.write_image('./plots/sodicAmphs.pdf')
    fig.show()
    return

    
