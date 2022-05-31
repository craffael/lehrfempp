#!/usr/bin/env python3

import click
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

@click.command()
@click.option('-f', '--file', required=True, help='csv file containing the results to plot')
@click.option('--l2_norm',show_default=True, default=False, help='plot l2 Error Norm')
@click.option('--sup_norm',show_default=True, default=False, help='plot Supremum Error Norm')
@click.option('--zero_form',show_default=True, default=False, help='plot results for zero form')
@click.option('--one_form',show_default=True, default=False, help='plot results for one form')
@click.option('--two_form',show_default=True, default=False, help='plot results for two form')
@click.option('--scale',show_default=True,type=click.Choice(["linear", "semilog", "loglog"], case_sensitive=False),  default="loglog", help='choose scales of the axis')

def plot(file, l2_norm, sup_norm, zero_form, one_form, two_form, scale):

    df = pd.read_csv(file)
    df_y = df.iloc[1:,4:]
    # print table
    print(df.iloc[1:, 12::12].head())

    num_k = df_y.shape[1] // 12
    ks = df.iloc[0,4::12].to_numpy()

    df_x = df[["hMax"]]
    df_x = df_x.iloc[1:]
    x = df_x.to_numpy().flatten()
    df_x_label = "h_Max"

    names = ["Sup Error Zero Form", "Sup Error One Form", "Sup Error Two Form", "L2 Error Zero Form", "L2 Error One Form", "L2 Error Two Form"]

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

    # if non is selected select only l2
    if m <= 0:
        kind_index = [1]
        m = 1

    if n <= 0:
        form_index = range(3)
        n = 3

    fig, axs = plt.subplots(m,n)

    # loop over the chosen norms 
    for i in range(m):

        # loop over the chosen forms
        for j in range(n):
            if m == 1 or n == 1:
                curr_axs = axs[i + j]
            else:
                curr_axs = axs[i, j]

            # rounds the axis on two siginficant positions
            curr_axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            curr_axs.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            # loop over the k values in the formula -\Delta u + k^2 u = f
            for k in range(num_k):
                df_col_idx = 2 * (kind_index[i] * n + form_index[j]) + 2 * k * n * m
                y = df_y.iloc[:, df_col_idx].to_numpy()
                if scale.lower() == 'semilog':
                    curr_axs.semilogy(x,y)
                if scale.lower() == 'loglog':
                    curr_axs.loglog(x,y)
                if scale.lower() == 'linear':
                    curr_axs.plot(x,y)
            curr_axs.set(ylabel=names[kind_index[i] * n + form_index[j]], xlabel=df_x_label, xlim=(max(x), min(x)))
            curr_axs.legend(ks)

    plt.show()

if __name__ == '__main__':
    plot()
