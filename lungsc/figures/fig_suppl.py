# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Supplementary figures
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions


fig_fdn = '../../figures/endomese_share/endo_paper_supplementary/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(doublet == 0)',
        inplace=True)
    print('Total cells analyzed: {:}'.format(ds0.n_samples))
    ds0.samplesheet['TimepointHO'] = ds0.samplesheet['Timepoint'].copy()
    ds0.samplesheet.loc[ds0.samplesheet['Treatment'] == 'hyperoxia', 'TimepointHO'] = 'P7_HO'

    ds = ds0.query_samples_by_metadata('Treatment == "normal"')

    #vs = ds.samplesheet[['embed_endo_1', 'embed_endo_2']].copy()

    genes_imm = ['Vegfa', 'Vegfb', 'Fgf3', 'Fgf7', 'Dll3', 'Mcpt8', 'Mcpt4']
    immune_subtypes = [
            'Mac I',
            'Mac II',
            'Mac III',
            'Mac IV',
            'Mac V',
            'DC I',
            'DC II',
            'DC III',
            'mast cell',
            'basophil',
            'neutrophil',
            'B cell',
            'NK cell',
            'T cell',
            'IL cell',
        ]
    fig, axs = plt.subplots(1, 2, figsize=(4, 7), sharey=True)
    for ax, tp in zip(axs, ['E18.5', 'P1']):
        dstp = ds.query_samples_by_metadata('TimepointHO == @tp', local_dict=locals())
        dstp.plot.dot_plot(
            axis='samples',
            group_by='cellSubtype',
            group_order=immune_subtypes,
            plot_list=genes_imm,
            ax=ax,
            layout='horizontal',
            cmap='viridis',
            vmax='max_single',
            vmin='min_single',
            min_size=0,
        )
        ax.set_title(tp)
        ax.set_xticklabels(genes_imm, rotation=90)
    fig.tight_layout()

    genes = ['Col1a1', 'Col5a1', 'Col6a1', 'Col13a1', 'Col16a1']
    cell_subtypes = [
            #'Mac I',
            #'Mac II',
            #'Mac III',
            #'Mac IV',
            #'Mac V',
            #'DC I',
            #'DC II',
            #'DC III',
            #'mast cell',
            #'basophil',
            #'neutrophil',
            #'B cell',
            #'NK cell',
            #'T cell',
            #'IL cell',
            'Early Car4- capillaries',
            'Car4+ capillaries',
            'Proliferative EC',
            'Venous EC',
            'Arterial EC I',
            'Arterial EC II',
            'Lymphatic EC',
        ]
    fig, axs = plt.subplots(1, 2, figsize=(4, 7), sharey=True)
    for ax, tp in zip(axs, ['P7', 'P7_HO']):
        dstp = ds0.query_samples_by_metadata('TimepointHO == @tp', local_dict=locals())
        dstp.plot.dot_plot(
            axis='samples',
            group_by='cellSubtype',
            group_order=cell_subtypes,
            plot_list=genes,
            ax=ax,
            layout='horizontal',
            cmap='viridis',
            vmax='max_single',
            vmin='min_single',
            min_size=0,
        )
        ax.set_title(tp)
        ax.set_xticklabels(genes, rotation=90)
    fig.tight_layout()

