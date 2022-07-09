#!/usr/bin/env python3

import click
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

@click.command()
@click.option('-f', '--file', required=True, help='csv file containing the results to plot')
@click.option('-s', '--stat', type=click.Choice(['mean', 'max', 'min', 'sd'], case_sensitive=False), help='How to determine the convergence rate')

def plot(file, stat):

    df = pd.read_csv(file)
    df_y = df.iloc[1:,4:]

    num_k = df_y.shape[1] // 12
    ks = df.iloc[0,4::12].to_numpy()

    # remove the sup norm
    df_y = df_y.iloc[:,1::2]
    pat = [False,False,False,True,True,True];
    pat = np.repeat([pat], num_k, axis = 0).flatten()
    df_y = df_y.iloc[1:,pat]
    
    # get only orders 
    if(stat == 'max'):
        mean_order = np.apply_along_axis(np.max, 0, df_y)
    elif(stat == 'min'):
        mean_order = np.apply_along_axis(np.min, 0, df_y)
    elif(stat == 'sd'):
        mean_order = np.apply_along_axis(np.std, 0, df_y)
    else:
        mean_order = np.apply_along_axis(np.mean, 0, df_y)

    fig, axs = plt.subplots(1,1)

    axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    axs.plot(ks, mean_order[0::3])
    axs.plot(ks, mean_order[1::3])
    axs.plot(ks, mean_order[2::3])

    # set labels and swap the x direction to get decreasing error
    axs.set(ylabel="algebraic order", xlabel="k")
    axs.legend(["zero component", "one component", "two component"])

    plt.savefig("eigenvalsdirac.png", dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    plot()
