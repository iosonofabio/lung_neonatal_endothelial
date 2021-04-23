# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/05/19
content:    Fig 9 of endo paper: hyperoxia and cell-cell interactions
'''
import os
import sys
import glob
import gzip
import pickle
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions


fig_fdn = '../../figures/endomese_share/endo_paper_figure_9/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(cellType == "endothelial") & (doublet == 0)',
        inplace=True)
    print('Total endothelial cells analyzed: {:}'.format(ds0.n_samples))

    vs = ds0.samplesheet[['embed_endo_1', 'embed_endo_2']].copy()
    vs.columns = ['dimension 1', 'dimension 2']

    ds7 = ds0.query_samples_by_metadata('Timepoint == "P7"')
    ds7p = ds7.split(['cellSubtype', 'Treatment'])

    if True:
        print('Cumulative plots of ligands and receptors in lymphatics')
        genes = [
            'Wnt2',
            'Wnt5a',
            'Fzd3',
            'Fzd4',
            'Ryk',
            'Lrp5',
            'Lrp6',
            'Notch1',
            'Notch2',
            'Jag1',
            'Jag2',
            'Dll4',
        ]
        ncells = {key: dsi.n_samples for key, dsi in ds7p.items()}
        cmap = {'normal': 'slateblue', 'hyperoxia': 'tomato'}
        lss = {'normal': '-', 'hyperoxia': '--'}
        fig, axs = plt.subplots(4, 3, figsize=(4.5, 4.5), sharey=True)
        axs = axs.ravel()
        for i, (gene, ax) in enumerate(zip(genes, axs)):
            for cond in ['normal', 'hyperoxia']:
                x = ds7p[('Lymphatic EC', cond)].counts.loc[gene].values.copy()
                x.sort()
                y = 1.0 - np.linspace(0, 1, len(x))
                x = np.insert(x, 0, 0)
                y = np.insert(y, 0, 1.0)
                x += 0.1
                nc = ncells[('Lymphatic EC', cond)]
                co = cond[0].upper()
                label = f'{co}, {nc}'
                ax.plot(x, y, lw=2, color=cmap[cond], ls=lss[cond], label=label)
            if i == 2:
                ax.legend(
                    title='Condition:',
                    loc='upper left',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )
            ax.grid(True)
            ax.set_xlim(0.05, 10**4 + 1)
            ax.set_yticks([0, 0.5, 1])
            ax.set_xscale('log')
            ax.set_xticks([1, 100])
            if i >= len(axs) - 3:
                ax.set_xticklabels(['$1$', '$100$'])
            else:
                ax.set_xticklabels(['', ''])
            ax.set_title(gene)
        fig.text(0.42, 0.02, 'Gene expression [cpm]', ha='center')
        fig.text(0.02, 0.52, 'Fraction of cells > x', va='center', rotation=90)
        fig.tight_layout(h_pad=0.1, rect=(0.03, 0.03, 1, 1))
        if True:
            fxf = fig_fdn+'cumulaives_Wnt_lymphatics'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    plt.ion()
    plt.show()
