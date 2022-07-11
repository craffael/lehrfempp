#!/usr/bin/env python3

import click
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

@click.command()
@click.option('-f', '--file', required=True, help='csv file containing the results to plot')
@click.option('-l', '--log', required=False, is_flag=True, help='Indicates if log log plot should be drawn')

def plot(file, log):

    plt.rcParams.update({
        "text.usetex": True,
        "font.family" : 'Palatino', 
        "font.size" : 11,
    })

    df = pd.read_csv(file)

    # exclude supremum norm and meta data
    df_y = df.iloc[:,4:]

    # get h values
    df_h = df[["hMax"]].iloc[1:]
    x = df_h.to_numpy().flatten()
    df_x_label = "h"

    # print table that is, print all l2 errors
    print(df_y)

    fig, axs = plt.subplots(1,1)

    # rounds the axis on two siginficant positions
    axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # add the three lines
    yzero = df_y.iloc[1:,1].to_numpy()
    yone = df_y.iloc[1:,3].to_numpy()
    ytwo = df_y.iloc[1:,5].to_numpy()
    inter = yzero[0]
    yref = x * x
    ylin = pow(x,1.5)

    # compute rates
    ratezero = np.log(yzero[1:6] / yzero[0:5]) / np.log(x[1:6] / x[0:5])
    rateone = np.log(yone[1:6] / yone[0:5]) / np.log(x[1:6] / x[0:5])
    ratetwo = np.log(ytwo[1:6] / ytwo[0:5]) / np.log(x[1:6] / x[0:5])
    
    print("Rate for zero component", ratezero)
    print("Rate for one component", rateone)
    print("Rate for two component", ratetwo)

    if(log):
        axs.loglog(x,yref)
        axs.loglog(x,ylin)
        axs.loglog(x,yzero, marker='+')
        axs.loglog(x,yone, marker='+')
        axs.loglog(x,ytwo, marker='+')
    else:
        axs.plot(x,yref)
        axs.loglog(x,ylin)
        axs.plot(x,yzero)
        axs.plot(x,yone)
        axs.plot(x,ytwo)


    # set labels and swap the x direction to get decreasing error
    axs.set(ylabel="$L^2$ error", xlabel=df_x_label, xlim=(max(x), min(x)))
    axs.legend(["$h^{2}$", "$h^{1.5}$","zero component", "one component", "two component"])

    #plt.title("$L^2$-Error of Hodge Laplacian Source Problems")
    plt.savefig("diracConvergenceDebug.png", dpi=300, bbox_inches='tight')
    plt.show()



if __name__ == '__main__':
    plot()
