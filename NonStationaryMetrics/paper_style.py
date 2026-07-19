# -*- coding: utf-8 -*-
"""
Shared single-column style for the non-stationary brachistochrone paper.

COL   = width of one column in a two-column layout (inches).
DCOL  = full text width (two columns) for the few wide figures.
set_style() sets rcParams sized for figures reduced to column width.
savefig() writes both PDF (vector, for the paper) and PNG (preview).

All figure text is in English.
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COL = 3.4          # single-column width (inch); revtex/AAS ~ 3.35, MNRAS ~ 3.32
DCOL = 7.0         # full two-column width

def set_style():
    plt.rcParams.update({
        'figure.constrained_layout.use': True,
        'font.size': 8,
        'axes.titlesize': 8,
        'axes.labelsize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 6.5,
        'legend.handlelength': 1.6,
        'lines.linewidth': 1.2,
        'axes.linewidth': 0.7,
        'xtick.major.width': 0.7,
        'ytick.major.width': 0.7,
        'mathtext.fontset': 'dejavusans',
        'savefig.dpi': 300,
        'figure.dpi': 150,
    })

def savefig(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(out_dir, f'{name}.{ext}'), bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {name} (column-width)')
