# A few useful functions for plotting.
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import aplpy

def plot_casa(xlabel='DETX', xunit='None', ylabel='DETY', yunit='None', figsize=[8, 6], nrows=1, ncols=1, coords=False,
              wcs=None, dpi=100):
    params = {'legend.fontsize': 'x-large',
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'medium',
              'ytick.labelsize': 'medium'}
    plt.rcParams.update(params)
    if coords == False:
        fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)
    if coords == True:
        fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi, subplot_kw={'projection': wcs})
    if xunit == 'None':
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(xlabel + ' [' + xunit + ']')
    if yunit == 'None':
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel(ylabel + ' [' + yunit + ']')
    if coords == True:
        ax.coords.grid(color='white')
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
        ax.tick_params(axis='x', which='major', width=1.00, length=4, direction='out', labelsize='x-large', top=False,
                       bottom=True)
        ax.tick_params(axis='y', which='major', width=1.00, length=4, direction='out', labelsize='x-large', right=False,
                       left=True)
        ax.tick_params(axis='both', which='minor', length=4)

    if coords == False:
        ax.tick_params(axis='both', which='major', width=1.00, length=4, direction='out', labelsize=18)
        ax.tick_params(axis='both', which='minor', width=0.75, length=2, direction='out', labelsize=12)
    ax.grid(False)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.xaxis.label.set_size(20)
    # ax.yaxis.label.set_size(20)
    ax.set_xlim(100, 380)
    ax.set_ylim(140, 370)
    ax.yaxis.labelpad = 50

    return fig, ax


def set_colorbar(fig, im, title, labelpad=25, fraction=0.046, pad=0.04):
    cbar = fig.colorbar(im, fraction=fraction, pad=pad)
    cbar.set_label(title, rotation=270, labelpad=labelpad)
    cbar.ax.tick_params(labelsize='x-large')
    return cbar


def plot_ncasa(xlabel='DETX', xunit='None', ylabel='DETY', yunit='None', figsize=[12, 12], nrows=1, ncols=1,
               coords=False, wcs=None):
    if coords == False:
        fig, axarr = plt.subplots(nrows, ncols, figsize=figsize)
    if coords == True:
        fig, axarr = plt.subplots(nrows, ncols, figsize=figsize, subplot_kw={'projection': wcs})
    axarr = axarr.flatten()
    for i in range(nrows * ncols):

        if xunit == 'None':
            axarr[i].set_xlabel(xlabel)
        else:
            axarr[i].set_xlabel(xlabel + ' [' + xunit + ']')
        if yunit == 'None':
            axarr[i].set_ylabel(ylabel)
        else:
            axarr[i].set_ylabel(ylabel + ' [' + yunit + ']')
        if coords == True:
            axarr[i].coords.grid(color='white')
            axarr[i].coords['ra'].set_axislabel('RA (J2000)')
            axarr[i].coords['dec'].set_axislabel('DEC (J2000)')
            axarr[i].tick_params(axis='both', which='major', width=1.00, length=8, direction='in', labelsize=18)
            axarr[i].tick_params(axis='both', which='minor', length=4)
        if coords == False:
            axarr[i].tick_params(axis='both', which='major', width=1.00, length=8, direction='in', labelsize=18)
            axarr[i].tick_params(axis='both', which='minor', width=0.75, length=4, direction='in', labelsize=12)
        # axarr[i].set_title(title)
        axarr[i].grid(False)
        axarr[i].xaxis.set_ticks_position('both')
        axarr[i].yaxis.set_ticks_position('both')
        axarr[i].xaxis.label.set_size(20)
        axarr[i].yaxis.label.set_size(20)


# ------------------------------------------------------------------------------------------------
# save plot into png or pdf files
def save_figures(fig, axarr, save_name, dpi=128, ftype='png'):
    print("Saving Plot %s ...") % (save_name)
    plt.savefig("%s.%s" % (save_name, ftype), dpi=dpi, bbox_inches='tight')


def aplpy_plot_casa(file):
    fig = aplpy.FITSFigure(file)
    fig.recenter(350.86737, 58.816033, radius=0.05)
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Flux density (Jy)')
    # radio.set_title('helo world km.s$^{-1}$')
    fig.colorbar.set_font(size='xx-large', weight='medium', \
                          stretch='normal', family='sans-serif', \
                          style='normal', variant='normal')
    fig.colorbar.set_axis_label_font(size=16, weight='bold')
    fig.axis_labels.set_font(size='xx-large', weight='medium', \
                             stretch='normal', family='sans-serif', \
                             style='normal', variant='normal')
    fig.tick_labels.set_yformat('dd:mm')
    fig.tick_labels.set_font(size='xx-large', weight='medium', \
                             stretch='normal', family='sans-serif', \
                             style='normal', variant='normal')
    fig.set_nan_color('white')
    # fig.set_theme('pretty')
    fig.colorbar.set_width(0.2)
    fig.colorbar.set_pad(0.05)
