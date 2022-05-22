#!/usr/bin/env python3

import click
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt

@click.command()
@click.option('-f', '--file', required=True, help='csv file containing the results to plot')
@click.option('--l2_norm',show_default=True, default=False, help='plot l2 Error Norm')
@click.option('--sup_norm',show_default=True, default=False, help='plot Supremum Error Norm')
@click.option('--zero_form',show_default=True, default=False, help='plot results for zero form')
@click.option('--one_form',show_default=True, default=False, help='plot results for one form')
@click.option('--two_form',show_default=True, default=False, help='plot results for two form')
@click.option('--scale',show_default=True,type=click.Choice(["linear", "semilog", "loglog"], case_sensitive=False),  default="linear", help='choose scales of the axis')
def plot(file, l2_norm, sup_norm, zero_form, one_form, two_form, scale):

    df = pd.read_csv(file)
    norm_names = ["L2ErrorZero", "L2ErrorOne", "L2ErrorTwo", "SupErrorZero", "SupErrorOne", "SupErrorTwo"]
    df_y = df.loc[:,norm_names]
    df_x = df[["hMax"]]
    x = df_x.to_numpy().flatten()

    df_x_label = "h_Max"

    kind_index = []
    if l2_norm :
        kind_index.append(0)
    if sup_norm :
        kind_index.append(1)

    form_index = []
    if zero_form :
        form_index.append(0)
    if one_form :
        form_index.append(1)
    if two_form :
        form_index.append(2)

    m = len(kind_index)
    n = len(form_index)

    # if non is selected select all
    if m <= 0:
        kind_index = range(2)
        m = 2

    if n <= 0:
        form_index = range(3)
        n = 3

    fig, axs = plt.subplots(m,n)
    for i in range(m):
        for j in range(n):
            if m == 1 or n == 1:
                curr_axs = axs[i + j]
            else:
                curr_axs = axs[i, j]
            df_col_idx = kind_index[i] * n + form_index[j]
            y = df_y.loc[:, norm_names[df_col_idx]].to_numpy()
            if scale.lower() == 'semilog':
                curr_axs.semilogy(x,y)
            if scale.lower() == 'loglog':
                curr_axs.loglog(x,y)
            if scale.lower() == 'linear':
                curr_axs.plot(x,y)
            curr_axs.set(ylabel=norm_names[df_col_idx], xlabel=df_x_label, xlim=(max(x), min(x)))

    plt.show()

if __name__ == '__main__':
    plot()
