# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 1 for the endo paper.
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_1/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(cellType == "endothelial") & (doublet == 0)',
        inplace=True)
    print('Total endothelial cells analyzed: {:}'.format(ds0.n_samples))
    ds = ds0.query_samples_by_metadata('Treatment == "normal"')

    #print('Load tSNE from file')
    #vs = pd.read_csv(
    #    '../../data/sequencing/datasets/all_{:}/tsne_with_hyperoxia_endo.tsv'.format(version),
    #    sep='\t',
    #    index_col=0,
    #    )
    #vs = vs.loc[ds.samplenames]
    vs = ds.samplesheet[['embed_endo_1', 'embed_endo_2']]
    vs.columns = ['dimension 1', 'dimension 2']
    # FIXME: location of lymphatic was tuned by hand (x + 9.5, y - 6)

    if False:
        print('Plot embedding with clusters, with legend')
        gene, title = ('cellSubtype', 'subtype')
        ln = ds.samplesheet[gene].str.len().max()
        cmap = {
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
        }
        fig, ax = plt.subplots(figsize=(4.8 + 0.07 * ln, 4.2))
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=15,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )

        d = ax._singlet_cmap
        handles = []
        labels = []
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Proliferative EC',
            'Nonproliferative embryonic EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            'Car4+ capillaries',
            'Lymphatic EC',
            #'Proliferative venous EC',
            ]
        for key in group_order:
            color = d[key]
            h = mlines.Line2D(
                [], [], color=color, marker='o', lw=0,
                markersize=5,
                )
            handles.append(h)
            labels.append(key.upper())
        ncol = 1
        fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )

        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}.{:}'.format(
                        gene, ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}.{:}'.format(
                        gene, ext))

    if False:
        print('Plot embedding with clusters, no legend')
        gene, title = ('cellSubtype', 'subtype')
        cmap = {
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
        }
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=15,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}_nolegend.{:}'.format(
                        gene, ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}_nolegend.{:}'.format(
                        gene, ext))

    if False:
        print('Plot timepoints')
        gene, title = 'Timepoint', 'timepoint'
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(4.8 + 0.07 * ln, 4.2))
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=15,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )

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
            labels.append(key.upper())
        labels_old = list(labels)
        labels = ['E18.5', 'P1', 'P7', 'P21']
        handles = [handles[labels_old.index(li)] for li in labels]
        ncol = 1
        fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )

        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}.{:}'.format(
                        gene, ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_tsne_metadata_{:}.{:}'.format(
                        gene, ext))

    if False:
        print('Plot embedding of genes')
        genes = [
            'Ankrd37',
            'Depp1',
            'Pde3a',
            'Hpgd',
            ]

        for gene in genes:
            title = gene
            fig, ax = plt.subplots(figsize=(4.8, 4.2))
            cmap = 'viridis'
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=15,
                    alpha=0.40,
                    cmap=cmap,
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            ax.set_title(gene)
            ax.set_axis_off()
            fig.tight_layout()

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_tsne_gene_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_tsne_gene_{:}.{:}'.format(
                            gene, ext))

    if True:
        print('Dotplot without hyperoxia')
        genes_dotplot = [
            'Gja5', 'Bmx',
            'Fn1',
            'Ctsh', 'Kcne3', 'Cdh13',
            'Car8', 'Mmp16', 'Slc6a2',
            'Thy1', 'Mmrn1', 'Ccl21a', 'Reln',
            'Neil3', 'Mki67', 'Aurkb',
            'Depp1', 'Ankrd37',
            'Peg3', 'Mest', 'Hpgd',
            'Cd36',
            'Car4', 'Sirpa', 'Fibin',
            ]
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            #'Proliferative venous EC',
            'Proliferative EC',
            'Nonproliferative embryonic EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            'Car4+ capillaries',
            ]

        cmap = {
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
        }
        fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 17]}, figsize=(9.3, 4))
        ax = axs[1]
        ds.plot.dot_plot(
            axis='samples',
            group_by='cellSubtype',
            group_order=group_order,
            plot_list=genes_dotplot,
            ax=ax,
            layout='horizontal',
            cmap='viridis',
            tight_layout=False,
            vmax='max_single',
            vmin='min_single',
            min_size=0,
            )
        #ax.set_xticklabels(group_order, rotation=90)
        ax.set_yticks([])
        for ict, ct in enumerate(group_order):
            # Semi-transparent background
            r = plt.Rectangle(
                (-0.5, ict - 0.5, -0.5), len(genes_dotplot), 1,
                facecolor=cmap[ct],
                edgecolor='none',
                alpha=0.1,
                zorder=0.5,
                )
            ax.add_patch(r)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)

        # Number of cells
        ncells = ds.samplesheet['cellSubtype'].value_counts()
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(ax.get_yticks())
        #ax2.set_xticklabels([str(ncells[x]) for x in group_order])

        # Legends
        sfun = ax._singlet_dotmap['fraction_size_map']
        handles1 = [
            ax.scatter([], [], color='grey', s=sfun(0)),
            ax.scatter([], [], color='grey', s=sfun(0.1)),
            ax.scatter([], [], color='grey', s=sfun(0.25)),
            ax.scatter([], [], color='grey', s=sfun(0.5)),
            ax.scatter([], [], color='grey', s=sfun(0.75)),
            ax.scatter([], [], color='grey', s=sfun(1)),
            ]
        labels1 = ['0%', '10%', '25%', '50%', '75%', '100%']
        leg1 = ax.legend(
                handles1, labels1,
                title='Fraction of\npopulation:    ',
                bbox_to_anchor=(1.01, 1.13),
                loc='upper left',
                )
        cfun = ax._singlet_dotmap['level_color_map']

        handles2 = [
            ax.scatter([], [], color=cfun(0.1), s=sfun(1)),
            ax.scatter([], [], color=cfun(0.6), s=sfun(1)),
            ax.scatter([], [], color=cfun(1.0), s=sfun(1)),
            ]
        labels2 = ['Low', 'Medium', 'High']
        leg2 = ax.legend(
                handles2, labels2,
                title='Expression\ncompared to\nhighest\npopulation:',
                bbox_to_anchor=(1.01, -0.02),
                loc='lower left',
                )
        ax.add_artist(leg1)

        ax = axs[0]
        ax.set_xticks([])
        ax.set_xlim(0, 1)
        ax.set_ylim(*axs[1].get_ylim())
        ax.set_yticks(np.arange(len(group_order)))
        for x in range(len(group_order)):
            ax.text(0.5, x, str(x+1), ha='center', va='center')
        for x in range(len(group_order) - 1):
            ax.axhline(x+0.5, lw=1, color='k')
        ax.set_yticklabels(group_order, rotation=0)


        if True:
            # Macrovascular
            r = plt.Rectangle(
                (-0.5, -0.5), len(genes_dotplot), 4,
                facecolor='none',
                edgecolor='black',
                alpha=0.5,
                zorder=10,
                lw=2,
                clip_on=False,
                )
            axs[1].add_patch(r)
            # Microvascular
            r = plt.Rectangle(
                (-0.5, 3.5), len(genes_dotplot), 5,
                facecolor='none',
                edgecolor='grey',
                alpha=0.5,
                zorder=10,
                lw=2,
                clip_on=False,
                )
            axs[1].add_patch(r)

            # Macrovascular
            r = plt.Rectangle(
                (0, -0.5), 1, 4,
                facecolor='none',
                edgecolor='black',
                alpha=0.5,
                zorder=10,
                lw=2,
                clip_on=False,
                )
            axs[0].add_patch(r)
            # Microvascular
            r = plt.Rectangle(
                (0, 3.5), 1, 5,
                facecolor='none',
                edgecolor='grey',
                alpha=0.5,
                zorder=10,
                lw=2,
                clip_on=False,
                )
            axs[0].add_patch(r)

        fig.tight_layout(rect=(0, 0, 1, 0.96), w_pad=-0.3)

        if True:
            for ext in ['svg', 'pdf', ['png', 300]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_dotplot.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_dotplot.{:}'.format(
                        ext))

        plt.ion()
        plt.show()

