# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/20
content:    Fig2 for the endo paper: harmonization with Tabula Muris.
'''
import os
import sys
import glob
import gzip
import argparse
import subprocess as sp
import numpy as np
import pandas as pd
import loompy
import scanpy

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, concatenate


fdn_data = '../../data/sequencing/datasets/endo_final/'
fns_loom = {
        'ours': fdn_data+'endo_normoxia.loom',
        'tabulamuris': '../../data/tabulamuris/FACS_alltissues/endos.loom',
        }

fig_fdn = '../../figures/endomese_share/endo_paper_figure_2/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    print('Load our endo data')
    ds = Dataset(
            dataset={
                'path': fns_loom['ours'],
                'index_samples': '_index',
                'index_features': 'GeneName',
                },
        )
    ds.samplesheet['Tissue'] = 'Lung'
    ds.samplesheet.index.name = 'CellID'
    ds.counts.columns.name = 'CellID'
    ds.samplesheet['Dataset'] = 'ours'

    print('Load Tabula Muris FACS endothelial')
    ds_tm = Dataset(
            dataset={
                'path': fns_loom['tabulamuris'],
                'index_samples': 'CellID',
                'index_features': 'GeneName',
                },
        )
    ds_tm.samplesheet['Tissue'] = ds_tm.samplesheet['tissue']
    ds_tm.samplesheet['Mousename'] = 'TM' + ds_tm.samplesheet['mouse.id']
    ds_tm.samplesheet['Dataset'] = 'TabulaMuris'
    ds_tm.samplesheet['Timepoint'] = '3MO'

    # They are all endos, without further classification (ouch!)
    ds_tm.samplesheet['Cell Type'] = ds_tm.samplesheet['cell_ontology_class']
    ds_tm.samplesheet['cellSubtype'] = 'unknown (TM)'

    print('Restrict TM to lung')
    ds_tm0 = ds_tm
    ds_tm = ds_tm0.query_samples_by_metadata('Tissue == "Lung"')

    print('Concatenate datasets')
    dsme = concatenate([ds, ds_tm], missing='pad')

    print('Normalize')
    dsme.counts.normalize('counts_per_million', inplace=True)

    if False:
        print('Feature selection')
        features = dsme.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
        dsf = dsme.query_features_by_name(features)

        print('PCA')
        dsc = dsf.dimensionality.pca(n_dims=50, robust=False, return_dataset='samples')

        print('BBKNN')
        adata = dsme.to_AnnData()
        adata.obsm['X_pca'] = dsc.counts.values.T
        scanpy.external.pp.bbknn(
                adata,
                batch_key='Mousename',
                neighbors_within_batch=3,
                set_op_mix_ratio=0.8,
                )
        adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']  # BUGFIX

        print('UMAP')
        scanpy.tl.umap(adata)
        vs = pd.DataFrame(
                adata.obsm['X_umap'],
                index=dsme.samplenames,
                columns=['dimension 1', 'dimension 2'],
                )

        ##vs = dsc.dimensionality.umap()

    else:
        print('Load UMAP from file')
        vs = pd.read_csv(
            '../../data/sequencing/datasets/endo_final/umap_with_TM_lung.tsv',
            sep='\t', index_col=0)

    if False:
        print('Plot embedding of metadata and many genes')
        genes = [
            'Dataset', 'cellSubtype', 'Tissue', 'Timepoint',
            'Gja5', 'Bmx', 'Fn1', 'Ctsh', 'Kcne3', 'Cdh13', 'Car8',
            'Mmp16', 'Slc6a2', 'Thy1', 'Mmrn1', 'Ccl21a', 'Reln',
            'Neil3', 'Mki67', 'Aurkb', 'Ankrd37', 'Peg3', 'Mest',
            'Hpgd', 'Cd36', 'Car4', 'Sirpa', 'Fibin',
            ]
        cmaps = {
            'Data source': {'Tabula Muris': 'darkred', 'ours': 'steelblue'},
            'cellSubtype': {
                'Lymphatic EC': '#ffda4f',
                'Venous EC': '#4380c9',
                'Arterial EC I': '#bd4f6f',
                'Arterial EC II': '#ee4242',
                'Proliferative EC': '#664fbd',
                'Nonproliferative embryonic EC': '#7d7d7d',
                'Early Car4- capillaries': '#a6bd4f',
                'Late Car4- capillaries': '#4fbd9d',
                'Car4+ capillaries': '#3ab646',
                'Proliferative venous EC': '#1e1e1e',
                'unknown (TM)': 'deeppink',
            },
            'Timepoint': {
                'E18.5': 'navy',
                'P1': 'gold',
                'P7': 'tomato',
                'P21': 'firebrick',
                '3MO': '#1e1e1e',
                }
        }

        fig = plt.figure(figsize=(12, 10))
        gs = fig.add_gridspec(11, 12)
        axs = []
        axs.append(fig.add_subplot(gs[0:3, 0:3]))
        axs.append(fig.add_subplot(gs[0:3, 3:6]))
        axs.append(fig.add_subplot(gs[0:3, 6:9]))
        axs.append(fig.add_subplot(gs[0:3, 9:12]))
        for i in range(4):
            for j in range(6):
                axs.append(fig.add_subplot(
                    gs[3+i*2:3+(i+1)*2, j*2: (j+1)*2]),
                    )
        for i, (gene, ax) in enumerate(zip(genes, axs)):
            cmap = cmaps.get(gene, 'viridis')
            dsme.plot.scatter_reduced(
                    vs,
                    color_by=gene,
                    color_log=(i >= 4),
                    cmap=cmap,
                    ax=ax,
                    alpha=0.2 - (0.1 * (gene == 'tissue')),
                    s=15,
                    )
            ax.set_axis_off()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(gene)
            if i >= 4:
                continue
            handles, labels = [], []
            for key, color in ax._singlet_cmap.items():
                handles.append(ax.scatter([], [], color=color))
                labels.append(key)
            ax.legend(
                    handles, labels, loc='lower center',
                    fontsize=8, ncol=2,
                    bbox_to_anchor=(0.5, 1.2), bbox_transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()

    if True:
        print('Plot embedding of a few things only')
        genes = ['Dataset', 'cellSubtype', 'Timepoint']

        cmaps = {
            'Data source': {'Tabula Muris': 'darkred', 'ours': 'steelblue'},
            'cellSubtype': {
                'Lymphatic EC': '#ffda4f',
                'Venous EC': '#4380c9',
                'Arterial EC I': '#bd4f6f',
                'Arterial EC II': '#ee4242',
                'Proliferative EC': '#664fbd',
                'Nonproliferative embryonic EC': '#7d7d7d',
                'Early Car4- capillaries': '#a6bd4f',
                'Late Car4- capillaries': '#4fbd9d',
                'Car4+ capillaries': '#3ab646',
                'Proliferative venous EC': '#1e1e1e',
                'unknown (TM)': 'deeppink',
            },
            'Timepoint': {
                'E18.5': 'navy',
                'P1': 'gold',
                'P7': 'tomato',
                'P21': 'firebrick',
                '3MO': '#1e1e1e',
                }
        }

        fig, axs = plt.subplots(1, 3, figsize=(9, 4))
        for i, (gene, ax) in enumerate(zip(genes, axs)):
            cmap = cmaps.get(gene, 'viridis')
            dsme.plot.scatter_reduced(
                    vs,
                    color_by=gene,
                    color_log=False,
                    cmap=cmap,
                    ax=ax,
                    alpha=0.2,
                    s=15,
                    )
            ax.set_axis_off()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(gene)
            if i >= 4:
                continue
            handles, labels = [], []
            for key, color in ax._singlet_cmap.items():
                handles.append(ax.scatter([], [], color=color))
                labels.append(key)
            ax.legend(
                    handles, labels, loc='lower center',
                    fontsize=8, ncol=2,
                    bbox_to_anchor=(0.5, 1.2), bbox_transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'comparison_with_tabula_muris'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    plt.ion()
    plt.show()
