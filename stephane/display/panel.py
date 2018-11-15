


import matplotlib.pyplot as plt
import stephane.display.graphes as graphes

def make(subplot,fignum=1,axis='off'):
    """
    From a dictionnary of subplot location, generate a panel
    """
    axes = {}
    fig = graphes.set_fig(fignum)
    for num in subplot:
        axes[num] = plot_axes(fig,num,axis=axis)
    return fig,axes

def draw(plots,axes):
    for key in axes.keys():
        sca(axes[key])
        plots[key]
        
def plot_axes(fig,num,axis='on'):
    
    ax = fig.add_subplot(num)
   # ax.set_aspect('equal', adjustable='box')
    ax.axis(axis)
    
    return ax    

def sca(ax):
    """
    Set current axis.
    INPUT
    -----
    ax : axes object
    """
    plt.sca(ax)
