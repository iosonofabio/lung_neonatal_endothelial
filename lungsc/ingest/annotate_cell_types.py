# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/06/19
content:    Annotate cell types for the full dataset
'''
import os
import sys
import glob
import gzip
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# Ensure leidenalg is correct
sys.path.insert(0, os.path.abspath('../../packages'))
from lungsc.ingest.load_dataset import versions, DatasetLung


if __name__ == '__main__':

    version = versions[-1]
    ds = DatasetLung.load(version=version)

    print('Feature selection')
    features = ds.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
    dsf = ds.query_features_by_name(features)

    print('PCA')
    dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

    #print('UMAP')
    #vsu = dsc.dimensionality.umap()

    print('tSNE')
    vst = dsc.dimensionality.tsne(perplexity=20)

    print('knn graph')
    neighbors, _, _ = dsc.graph.knn(axis='samples', n_neighbors=10, return_sparse=False)

    print('Leiden clustering')
    edges = set()
    for nei in neighbors:
        for n in nei:
            edges.add(frozenset(n))
    edges = [tuple(e) for e in edges]

    ds.samplesheet['cluster'] = ds.cluster.leiden(
            axis='samples',
            edges=edges,
            resolution_parameter=0.0005,
            )

    print('Guess cell type for each cluster')
    di = {'CD45': 'immune', 'CD31': 'endothelial', 'mesench': 'mesenchymal'}
    ds.samplesheet['cellType'] = ''
    for cl in ds.samplesheet['cluster'].unique():
        ind = ds.samplesheet['cluster'] == cl
        st = ds.samplesheet.loc[ind, 'SortType']
        ds.samplesheet.loc[ind, 'cellType'] = di[st.value_counts().idxmax()]

    ds.samplesheet['cellRoughSubtype'] = ''

    def assign_ct(cluster_numbers, cell_type):
        if np.isscalar(cluster_numbers):
            cluster_numbers = [cluster_numbers]
        ind = ds.samplesheet['cluster'].isin(cluster_numbers)
        ds.samplesheet.loc[ind, 'cellType'] = cell_type

    print('Plot dimensionality reduction of dataset')
    vs = vst
    fig, axs = plt.subplots(3, 11, figsize=(17.3, 5.5), sharex=True, sharey=True)
    axs = axs.ravel()
    marker_genes = [
            ('Pecam1', 'Cd31'),
            ('Ptprc', 'Cd45'),
            'Col6a2',
            'Cd3e',
            'Gzma',
            ('Ms4a1', 'Cd20'),
            'Cpa3',
            'Mcpt4',
            'Plac8',
            'Cd68',
            'Batf3',
            'Cd209c',
            'Stfa2',
            'Gja5',
            'Maf',
            'Car8',
            'Car4',
            'Tgfbi',
            'Mcam',
            'Pdgfra',
            'Hhip',
            'Wnt2',
            'Trpc6',
            'Mki67',
            ]
    markers = [
            'SortType',
            ] + marker_genes[:-1]
    if 'Treatment' in ds.samplesheet.columns:
        markers.append('Treatment')
    else:
        markers.append(marker_genes[-1])
    markers += [
            ('number_of_genes_1plusreads', 'n genes'),
            'Timepoint',
            'Gender',
            'Mousename',
            'cluster',
            'cellType',
            'cellRoughSubtype',
            ]
    mgs = [x if isinstance(x, str) else x[0] for x in marker_genes]
    for ipl, (gene, ax) in enumerate(zip(markers, axs)):
        print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
        if isinstance(gene, str):
            gene, title = gene, gene
        else:
            gene, title = gene
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=10,
                alpha=0.04 + 0.1 * (gene not in ['annotated', 'cluster', 'cellType', 'Mousename', 'Treatment', 'Timepoint']),
                color_by=gene,
                color_log=(gene in mgs + ['number_of_genes_1plusreads']),
                )
        ax.grid(False)
        ax.set_title(title)

        if gene == 'cluster':
            for com in ds.samplesheet['cluster'].unique():
                vsc = vs.loc[ds.samplesheet[gene] == com]
                xc, yc = vsc.values.mean(axis=0)
                ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                ax.text(xc, yc, str(com), fontsize=8, ha='center', va='bottom')

        if gene in ('Treatment', 'Timepoint'):
            import matplotlib.lines as mlines
            d = ax._singlet_cmap
            handles = []
            labels = []
            for key, color in d.items():
                h = mlines.Line2D(
                    [], [], color=color, marker='o', lw=0,
                    markersize=5,
                    )
                handles.append(h)
                if gene == 'Treatment':
                    key = key[0]
                    loc = 'upper left'
                else:
                    loc = 'lower left'
                labels.append(key.upper())
            if gene == 'Timepoint':
                labels_old = list(labels)
                labels = ['E18.5', 'P1', 'P7', 'P21']
                handles = [handles[labels_old.index(li)] for li in labels]
                ncol = 2
            else:
                ncol = 1
            ax.legend(handles, labels, loc=loc, fontsize=6, ncol=ncol)

    fig.tight_layout()

    plt.ion()
    plt.show()

    if False:
        print('Export cell type annotations')
        import loompy
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        with loompy.connect(fn_good) as dsl:
            dsl.ca['cellType'] = ds.samplesheet['cellType'].values

    if False:
        print('Export cell subtype annotations')
        import loompy
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        with loompy.connect(fn_good) as dsl:
            dsl.ca['cellRoughSubtype'] = ds.samplesheet['cellRoughSubtype'].values
