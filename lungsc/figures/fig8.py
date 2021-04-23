# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/05/19
content:    Look for distinguishing features of the hyperoxia samples
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_8/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-2]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(cellType == "endothelial") & (doublet == 0)',
        inplace=True)
    print('Total endothelial cells analyzed: {:}'.format(ds0.n_samples))

    ds0.samplesheet['TimepointHO'] = ds0.samplesheet['Timepoint'].copy()
    ds0.samplesheet.loc[ds0.samplesheet['Treatment'] == 'hyperoxia', 'TimepointHO'] = 'P7_HO'
    ds0.samplesheet['cellSubtypeCoarse'] = ds0.samplesheet['cellSubtype'].replace({
            'Nonproliferative embryonic EC': 'Car4- capillaries',
            'Early Car4- capillaries': 'Car4- capillaries',
            'Late Car4- capillaries': 'Car4- capillaries',
            })
    ds = ds0.query_samples_by_metadata('cellType == "endothelial"')

    print('Read pregenerated DEG in normoxia vs hyperoxia')
    ct = 'endothelial'
    fn_comps = '../../data/gene_lists/comps_hyperoxia_normoxia_{:}.pkl'.format(ct)
    with open(fn_comps, 'rb') as f:
        comps = pickle.load(f)

    print('Figure out number of cells')
    ncells = {key: comp.iloc[0][['# cells normoxia P7', '# cells hyperoxia P7']] for key, comp in comps.items()}
    ncells = pd.DataFrame(ncells).T.astype(int)

    if False:
        print('Plot cell number changes')
        ln = ncells.index.str.len().max()
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            'Proliferative venous EC',
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
        fig, ax = plt.subplots(figsize=(5, 3))
        x = [0, 1]
        for cst in group_order:
            if cst not in ncells.index:
                continue
            row = ncells.loc[cst]
            y = row[['# cells normoxia P7', '# cells hyperoxia P7']].values
            ax.plot(x, y, 'o-', lw=2, label=cst, color=cmap[cst])
        ax.grid(True)
        ax.set_xlim(-0.2, 1.2)
        ax.set_xticks(x)
        ax.set_yscale('log')
        ax.set_ylabel('# cells')
        ax.legend(
                loc='upper left', title='Cell subtype:',
                bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                )
        ax.set_xticklabels(['Normoxia', 'Hyperoxia'])
        fig.tight_layout()
        if True:
            fxf = 'endo_hypoeroxia_ncells'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext))

    if False:
        print('Plot cell fraction changes')
        ln = ncells.index.str.len().max()
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            'Proliferative venous EC',
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
        frac_cells = 1.0 * ncells / ncells.sum(axis=0)
        fig, ax = plt.subplots(figsize=(5, 3))
        x = [0, 1]
        for cst in group_order:
            if cst not in ncells.index:
                continue
            row = frac_cells.loc[cst]
            y = row[['# cells normoxia P7', '# cells hyperoxia P7']].values
            ax.plot(x, y, 'o-', lw=2, label=cst, color=cmap[cst])
        ax.grid(True)
        ax.set_xlim(-0.2, 1.2)
        ax.set_xticks(x)
        ax.set_yscale('log')
        ax.set_ylabel('Fraction of endothelial cells')
        ax.legend(
                loc='upper left', title='Cell subtype:',
                bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                )
        ax.set_xticklabels(['Normoxia', 'Hyperoxia'])
        fig.tight_layout()
        if True:
            fxf = 'endo_hypoeroxia_frac_cells'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext))

    if False:
        print('Plot cell fraction changes, compact')
        ln = ncells.index.str.len().max()
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            'Proliferative venous EC',
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
        data = (1.0 * ncells / ncells.sum(axis=0)).loc[group_order]
        data.loc['_'] = 0
        data = data.loc[['_'] + list(data.index[:-1])]
        ys = 100 * (1.0 - data.cumsum(axis=0))

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(4.5, 2.9))
        x1 = 0.1
        x = [0, x1, 1 - x1, 1]
        xint = np.linspace(0, 1, 100)
        yints = []
        for i in range(len(ys)):
            y = [ys.iloc[i, 0]] * 2 + [ys.iloc[i, 1]] * 2
            yint = interpolate.pchip_interpolate(x, y, xint)
            yints.append(yint)
            ax.plot(xint, yint, lw=2, color='k', alpha=0.3)
            if i != 0:
                if data.iloc[i].max() > 0.01:
                    label = ys.index[i]
                else:
                    label = ''
                ax.fill_between(
                        xint, yint, yints[-2],
                        color=cmap[ys.index[i]],
                        label=label,
                        )
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Normoxia', 'Hyperoxia'])
        ax.xaxis.set_ticks_position('top')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylim(0, 100)
        ax.set_ylabel('Percent of\nendothelial cells')
        ax.legend(
                title='Cell subtype:',
                loc='upper left',
                bbox_to_anchor=(1.11, 1.01),
                bbox_transform=ax.transAxes,
                fontsize=10,
                )
        fig.tight_layout()
        if True:
            fxf = 'endo_hypoeroxia_ncells_compact'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext))

    if True:
        print('Find DEGs in P1 -> P7 and see what happens with HO')
        dsp = ds.split(['cellSubtype', 'TimepointHO'])
        csts = [
            'Arterial EC I',
            ]
        for cst in csts:
            print(cst)
            comp = dsp[(cst, 'P7')].compare(
                    dsp[(cst, 'P1')],
                    method='kolmogorov-smirnov-rich',
                    )

    if True:
        print('Wedge plot of normal, HO, and P1')
        ds.obs['cellSubtype_TimepointHO'] = ds.obs['cellSubtype'] + '_' + ds.obs['TimepointHO']
        csts = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
        ]
        xlabels2 = ['Art I', 'Art II', 'Veins']
        tps = ['P1', 'P7', 'P7_HO']
        xlabels3 = ['P1', 'P7', 'HO']
        group_order = []
        xlabels = []
        for cst in csts:
            group_order.extend([cst+'_'+tp for tp in tps])
            xlabels.extend(xlabels3)
        ds4 = ds.query_samples_by_metadata('cellSubtype in @csts', local_dict=locals())

        genes_dotplot = [
                'Adrb2',
                'Sox7',
                'Ntrk2',
                'Npr3',
                'Tslp',
                #'Npr1',
                #'Hpgd',
                'Ets1',
                'Gja4',
                #'Runx1',
                #'Cdh5',
                'Sept2',
                'Cd34',
        ]
        fnpl = ['P1_P7', 'P7_HO']
        figsize = (2, 1 + 0.4 * len(genes_dotplot))
        fig, ax = plt.subplots(figsize=figsize)
        ds4a = ds4.copy()
        ds4a.counts.log(inplace=True)
        ds4a = ds4a.average('samples', by='cellSubtype_TimepointHO')
        for ig, gene in enumerate(genes_dotplot):
            norm = ds4a.counts.loc[gene].max()
            for j, g in enumerate(csts):
                go = group_order[j * 3: (j + 1) * 3]
                x = j * 3 + np.arange(3)
                y = ds4a.counts.loc[gene, go].values
                y = (y + 1) / (norm + 1)

                # Check type
                if np.abs(y[1] - y[0]) > 0.1:
                    if (y[2] > y[1]) ^ (y[1] > y[0]):
                        color = 'royalblue'
                    else:
                        color = 'tomato'
                elif np.abs(y[2] - y[0]) > 0.1:
                    color = 'purple'
                else:
                    color = 'grey'

                y -= 0.5
                y *= 0.8
                y += len(genes_dotplot) - 1 - ig
                ax.plot(x, y, lw=2, color='k')
                ax.add_artist(
                        plt.Rectangle(
                            (j * 3 - 0.5, len(genes_dotplot) - ig - 1.5),
                            3, 1, edgecolor='none', facecolor=color, alpha=0.5,
                            zorder=0))
            if ig != 0:
                ax.axhline(ig - 0.5, lw=1, color='grey')
        ax.set_yticks(np.arange(len(genes_dotplot)))
        ax.set_yticklabels(genes_dotplot[::-1])
        ax.set_ylim(-0.5, len(genes_dotplot) - 0.5)
        ax.set_xlim(-0.5, len(group_order) - 0.5)
        ax.set_xticks(np.arange(len(group_order)))
        ax.set_xticklabels(xlabels, rotation=90)
        for i in range(len(csts) - 1):
            ax.axvline(i * 3 + 2.5, lw=1, color='grey')
        ax2 = ax.twiny()
        ax2.set_xticks(3 * np.arange(len(csts)) + 1)
        ax2.set_xticklabels(xlabels2, rotation=90)
        ax2.set_xlim(ax.get_xlim())
        fig.tight_layout()

        if True:
            fnf = fig_fdn+'hyperoxia_wedge_plots_endo_ArtI'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fnf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fnf, ext))

    if False:
        print('Select genes that are high in both fold change and fraction of cells')
        from collections import defaultdict, Counter
        th = defaultdict(lambda: 0.5)
        th['Early Car4- capillaries'] = 0.3
        th['Proliferative EC'] = 0.3
        th['Car4+ capillaries'] = 0.3
        # These have at least 5 cells per condition
        group_order = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            'Proliferative EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            'Car4+ capillaries',
            ]
        nshorts = {'up': Counter(), 'down': Counter()}
        fc_total = pd.Series(np.zeros(len(comps[group_order[0]])), index=comps[group_order[0]].index)
        for cst in group_order:
            comp = comps[cst]
            comp = comp.loc[comp['statistic'] > th[cst]]
            comp = comp.loc[np.abs(comp['log2_fold_change']) > 1]
            comp_up = comp.loc[comp['log2_fold_change'] > 0]
            genes_up_all = comp_up.index
            comp_down = comp.loc[comp['log2_fold_change'] < 0]
            genes_dw_all = comp_down.index

            # L1 metric (filtered and capped)
            fc_total.loc[genes_up_all] += np.minimum(3, comp.loc[genes_up_all, 'log2_fold_change'])
            fc_total.loc[genes_dw_all] += np.maximum(-3, comp.loc[genes_dw_all, 'log2_fold_change'])

    if False:
        print('Plot arrows')
        genes_common_up = list(fc_total.nlargest(8).index)
        genes_common_down = list(fc_total.nsmallest(8).index)
        cmap2 = {
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
        fig, axs = plt.subplots(2, 7, figsize=(13, 6.5), sharex=True, sharey=True)
        genes_common_both = [genes_common_up, genes_common_down]
        for k in range(2):
            for i in range(axs.shape[1]):
                gene = genes_common_both[k][i]
                ax = axs[k, i]
                ax.set_title(gene)
                x = np.arange(len(group_order))
                y1 = np.array([comps[cst].at[gene, 'avg_normoxia'] for cst in group_order])
                y2 = np.array([comps[cst].at[gene, 'avg_hyperoxia'] for cst in group_order])
                y1 = np.log10(0.1 + y1)
                y2 = np.log10(0.1 + y2)
                for j in range(len(x)):
                    ax.arrow(
                            x[j], y1[j], 0, y2[j] - y1[j],
                            length_includes_head=True,
                            width=0.25, head_length=min(0.2, np.abs(y2[j] - y1[j])),
                            facecolor=cmap2[group_order[j]],
                            edgecolor='none',
                            )
                ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)
                ax.set_xticks(x)
                ax.set_xticklabels(group_order, rotation=90)
                ax.set_yticks([-1, 0, 1, 2, 3, 4])
                ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$'])
                ax.set_ylim(-1, 4)
        axs[0, 0].set_ylabel('Up\nin HO\n[cpm]', rotation=0, labelpad=25)
        axs[1, 0].set_ylabel('Down\nin HO\n[cpm]', rotation=0, labelpad=25)
        fig.tight_layout()

        if True:
            fxf = 'endo_hypoeroxia_arrows'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'{:}.{:}'.format(fxf, ext))

    if False:
        print('Get longer lists for pathway analysis')
        genes_long_up = list(fc_total.nlargest(100).index)
        genes_long_down = list(fc_total.nsmallest(100).index)

        print('UP:')
        print(','.join(genes_long_up))

        print('DOWN:')
        print(','.join(genes_long_down))

        # METASCAPE HERE

    if False:
        print('Restrict to TFs in 2 interesting populations')
        subs = [
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            'Proliferative EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            'Car4+ capillaries',
            ]

        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()

        th = defaultdict(lambda: 0.5)
        th['Early Car4- capillaries'] = 0.3
        th['Proliferative EC'] = 0.3
        th['Car4+ capillaries'] = 0.3
        comp_tfs = {}
        for cst in subs:
            comp = comps[cst]
            comp['is_TF'] = comp.index.isin(tfs)
            comp = comp.loc[comp['statistic'] > th[cst]]
            comp = comp.loc[np.abs(comp['log2_fold_change']) >= 1]
            comp_tf = comp.loc[comp['is_TF']].sort_values(
                    'log2_fold_change', ascending=False)
            comp_tfs[cst] = comp_tf

        print('Output result to Excel')
        with pd.ExcelWriter(fig_fdn+'DEG_TFs.xlsx') as ex:
            for key, df in comp_tfs.items():
                df.to_excel(ex, sheet_name=key)


    plt.ion()
    plt.show()

