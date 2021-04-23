# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 7 for the endo paper: cell-cell interactions.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions


fig_fdn = '../../figures/endomese_share/endo_paper_figure_7/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-2]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.samplesheet['TimepointHO'] = ds0.samplesheet['Timepoint'].copy()
    ds0.samplesheet.loc[ds0.samplesheet['Treatment'] == 'hyperoxia', 'TimepointHO'] = 'P7_HO'

    ds_imm = ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0)',
        inplace=False)

    ds0.query_samples_by_metadata(
        '(cellType == "endothelial") & (doublet == 0)',
        inplace=True)
    print('Total endothelial cells analyzed: {:}'.format(ds0.n_samples))

    ds0.samplesheet['cellSubtypeCoarse'] = ds0.samplesheet['cellSubtype'].replace({
            'Nonproliferative embryonic EC': 'Car4- capillaries',
            'Early Car4- capillaries': 'Car4- capillaries',
            'Late Car4- capillaries': 'Car4- capillaries',
            })
    ds = ds0.query_samples_by_metadata('Treatment == "normal"')

    if False:
        print('Plot the cumulative distributions')
        focus = 'Car4+ capillaries'
        dsp = ds.split('cellSubtype')
        csts = [
            'Car4+ capillaries',
            'Nonproliferative embryonic EC',
            'Proliferative EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            'Arterial EC I',
            'Arterial EC II',
            'Venous EC',
            'Lymphatic EC',
            ]
        titles = [
            'Car4+ caps',
            'Embryonic',
            'Proliferative',
            'Early caps',
            'Late caps',
            'Art I',
            'Art II',
            'Venous',
            'Lymphatic',
            ]
        pops = list(zip(csts, titles))
        pairs = {
            'Car4+ capillaries': [('Apln', 'Aplnr')],
            }
        for ligand, receptor in pairs['focus']:
            fig, axs = plt.subplots(3, 3, figsize=(5, 5), sharex=True, sharey=True)
            axs = axs.ravel()
            for i in range(9):
                ax = axs[i]
                cst, title = pops[i]
                dsi = dsp[cst]
                x1 = np.sort(np.log10(dsi.counts.loc[ligand].values + 0.1))
                x2 = np.sort(np.log10(dsi.counts.loc[receptor].values + 0.1))
                ax.plot(x1, 1.0 - np.linspace(0, 1, len(x1)), color='tomato', label=ligand)
                ax.plot(x2, 1.0 - np.linspace(0, 1, len(x2)), color='steelblue', label=receptor)
                ax.set_yticks([0, 0.5, 1.0])
                ax.set_yticklabels(['0%', '50%', '100%'])
                ax.set_xticks([-1, 1, 3, 5])
                ax.set_xticklabels(['$0$', '$10$', '$10^3$', '$10^5$'])
                ax.grid(True)
                ax.set_xlim(-1, 6)
                ax.set_title(titles[i], fontsize=10)
                if i == 0:
                    ax.legend(fontsize=8, frameon=False)
            fig.text(0.56, 0.02, 'Gene expression [cpm]', fontsize=10, ha='center')
            fig.text(0.02, 0.52, 'Fraction of cells > x', fontsize=10, rotation=90, va='center')
            fig.tight_layout(rect=(0.04, 0.04, 1, 1))

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_cumulative_{:}_{:}.{:}'.format(
                            ligand, receptor, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_cumulative_{:}_{:}.{:}'.format(
                            ligand, receptor, ext))

    if False:
        print('Make arrow plot outwards of a certain population')
        focus = 'Car4+ capillaries'
        cmap = {
                'Lymphatic EC': 'greenyellow',
                'Venous EC': 'turquoise',
                'Arterial EC I': 'tomato',
                'Arterial EC II': 'orange',
                'Proliferative EC': 'gold',
                'Nonproliferative embryonic EC': 'silver',
                'Early Car4- capillaries': 'coral',
                'Late Car4- capillaries': 'chocolate',
                'Car4+ capillaries': 'peru',
                'Proliferative venous EC': 'steelblue',
            }
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

        interactions = {'Car4+ capillaries': [
            ('Apln-Aplnr', '+',
             {'high': [
                'Nonproliferative embryonic EC',
                'Proliferative EC',
                'Venous EC',
                'Early Car4- capillaries',
                ],
              'mid': [
                'Late Car4- capillaries',
                ],
              'low': [
                'Arterial EC I',
                'Arterial EC II',
                ],
              }),
            ('Fgf1-Fgfr3', '+',
             {'high': [
                'Arterial EC II',
                'Venous EC',
                ],
              'mid': [
                'Arterial EC I',
                'Car4+ capillaries',
                ],
              'low': [
                'Early Car4- capillaries',
                ],
              }),
            ('Itgb5-Fn1', '-',
             {'high': [
                'Arterial EC II',
                'Lymphatic EC',
                ],
              'mid': [
                ],
              'low': [
                'Arterial EC I',
                'Venous EC',
                ],
              }),
            ],
        }

        fig = plt.figure(figsize=(9, 6))
        gs = fig.add_gridspec(2, 3)
        axs = []
        axs.append(fig.add_subplot(gs[:, :2], projection='polar'))
        axs.append(fig.add_subplot(gs[0, 2], projection='polar'))
        axs.append(fig.add_subplot(gs[1, 2], projection='polar'))
        for iint, (label, direction, dic) in enumerate(interactions[focus]):
            ax = axs[iint]
            ax.scatter(
                0, 0, s=200,
                marker='s',
                zorder=10,
                facecolor=cmap['Car4+ capillaries'], alpha=0.95,
                )

            n = len(group_order)
            x = 2 * np.pi * np.arange(n) / n
            y = np.ones(n)
            colors = np.array([cmap[x] for x in group_order])

            ind = np.array([x != 'Car4+ capillaries' for x in group_order])
            ax.scatter(
                    x[ind], y[ind], c=colors[ind],
                    s=200, alpha=0.95, zorder=10, marker='o',
                    clip_on=False,
                    )

            ind = ~ind
            ax.scatter(
                    x[ind], y[ind], c=colors[ind],
                    s=200, alpha=0.95, zorder=10, marker='s',
                    clip_on=False,
                    )
            ax.grid(True, alpha=0.2)

            # Plot the arrows
            widths = {'low': 0.1, 'mid': 0.2, 'high': 0.4}
            for key, targets in dic.items():
                for target in targets:
                    xi = 2 * np.pi * group_order.index(target) / n
                    if direction == '+':
                        y0, dy = 0.1, 0.8
                    else:
                        y0, dy = 0.9, -0.8
                    ax.arrow(
                        xi, y0, 0, dy,
                        length_includes_head=True,
                        lw=3,
                        head_width=widths[key] * 0.6 * (1 + int(iint > 0)),
                        head_length=widths[key] * 0.6,
                        width=widths[key] * 0.1 * (1 + int(iint > 0)),
                        overhang=0.1,
                        edgecolor='none',
                        facecolor='k',
                        zorder=10,
                        )

            ax.set_rmin(0)
            ax.set_rmax(1)
            ax.set_rticks([])  # Less radial ticks
            #ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax.set_xticks(x)
            if iint == 0:
                ticks = [
                    'Arterial I',
                    '\nArterial II',
                    '\n\n\nVenous',
                    '\n\nProliferative',
                    'Nonproliferative   \nembryonic',
                    'Early Car4-\ncapillaries',
                    'Late Car4-\ncapillarie\n\n',
                    'Car4+ capillaries\n\n\n',
                    'Lymphatic\n',
                    #'Proliferative venous EC',
                    ]
                ax.set_xticklabels(ticks)
                ax.xaxis.set_tick_params(pad=35)
                ax.set_title(label+'\n\n')
            else:
                ax.set_xticklabels([])
                ax.set_title(label)
        fig.tight_layout(h_pad=4)

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_signaling_arrows_{:}.{:}'.format(
                        focus, ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_signaling_arrows_{:}.{:}'.format(
                        focus, ext))

    if False:
        print('Heatmap with number of interactions')
        pvals = {}
        for tp in ['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']:
            print('Timepoint: {:}'.format(tp))
            fn_res = '../../data/cellphonedb/endothelial/results/{:}/pvalues.txt'.format(tp)
            pval = pd.read_csv(
                fn_res,
                sep='\t',
                index_col=1,
                )
            pvals[tp] = pval

        print('Split by time')
        dsp = ds.split('Timepoint')

        csts = [
                'Arterial EC I',
                'Arterial EC II',
                'Venous EC',
                'Lymphatic EC',
                'Proliferative EC',
                'Nonproliferative embryonic EC',
                'Early Car4- capillaries',
                'Late Car4- capillaries',
                'Car4+ capillaries',
                #'Proliferative venous EC',
                ]
        alpha = 0.01
        n = len(csts)
        nints = np.ma.masked_all((6, n, n), np.int64)
        for i, ct1 in enumerate(csts):
            for j, ct2 in enumerate(csts):
                for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                    col = ct1+'|'+ct2
                    if col not in pvals[tp]:
                        col = ct2+'|'+ct1
                    if col not in pvals[tp]:
                        continue
                    nints[k, i, j] = (pvals[tp][col] < alpha).sum()
        nints_sym = nints + nints.swapaxes(1, 2)

        if False:
            print('Plot by time in separate figures, each normalized')
            for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                nints_plot = pd.DataFrame(nints_sym[k], index=csts, columns=csts)
                fig, ax = plt.subplots(figsize=(6, 6))
                sns.heatmap(nints_plot, ax=ax, cmap='plasma', vmin=0)
                ax.set_title(tp)
                fig.tight_layout()

                if True:
                    for ext in ['svg', 'pdf', ['png', 600]]:
                        if isinstance(ext, list):
                            ext, dpi = ext
                            fig.savefig(fig_fdn+'endo_signaling_heatmap_{:}.{:}'.format(
                                tp, ext),
                                dpi=dpi)
                        else:
                            fig.savefig(fig_fdn+'endo_signaling_heatmap_{:}.{:}'.format(
                                tp, ext))

        if False:
            print('Plot by time in single figure, each normalized')
            fig, axs = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
            axs = axs.ravel()
            for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                ax = axs[k]
                nints_plot = pd.DataFrame(nints_sym[k], index=csts, columns=csts)
                sns.heatmap(
                        nints_plot, ax=ax, cmap='plasma',
                        vmin=0,
                        )
                ax.set_title(tp)
            fig.suptitle('     Number of specific interactions (P < 0.01)')
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_signaling_heatmap.{:}'.format(
                            ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_signaling_heatmap.{:}'.format(
                            ext))

        if True:
            print('Plot by time in single figure, common normalization')
            fig = plt.figure(
                figsize=(9.2, 8),
                )
            gs = fig.add_gridspec(2, 3, width_ratios=[15] * 2 + [1])
            axs = []
            axs.append(fig.add_subplot(gs[0, 0]))
            axs.append(fig.add_subplot(gs[0, 1]))
            axs.append(fig.add_subplot(gs[1, 0]))
            axs.append(fig.add_subplot(gs[1, 1]))
            cax = fig.add_subplot(gs[:, -1])
            vmax = nints_sym[:4].max()
            ns = (ds.samplesheet
                    .loc[:, ['Timepoint', 'cellSubtype', 'cellType']]
                    .groupby(['Timepoint', 'cellSubtype'])
                    .count()
                    .loc[:, 'cellType']
                    .unstack().T.fillna(0)
                    .astype(int))
            for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21']):
                ax = axs[k]
                nints_plot = pd.DataFrame(nints_sym[k], index=csts, columns=csts)
                nk = ns[tp]
                for cst, n in nk.items():
                    if n < 3:
                        nints_plot.loc[cst] = np.nan
                        nints_plot.loc[:, cst] = np.nan
                sns.heatmap(
                        nints_plot, ax=ax, cmap='plasma',
                        vmin=0,
                        vmax=vmax,
                        cbar=False,
                        )
                ax.set_title(tp)
                xlim = ax.get_xlim()
                ax.set_xlim(xlim[0], xlim[1] + 1)
                yts = ax.get_yticks()
                ytls = axs[0].get_yticklabels()
                for yt, ytl in zip(yts, ytls):
                    n = ns.loc[ytl.get_text(), tp]
                    ax.text(xlim[1] + 0.5, yt, str(n), ha='center', va='center',
                            fontsize=9)
                if k < 2:
                    ax.set_xticklabels([])
                if k in (1, 3):
                    ax.set_yticklabels([])

                if k >= 2:
                    xticks = list(ax.get_xticks())
                    ax.set_xticks(xticks + [xticks[-1] + 1])
                    ax.set_xticklabels(list(csts) + ['number of cells'])

            norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='plasma'), cax=cax)
            cb.set_label('number of interactions', labelpad=15)
            fig.suptitle('     Number of specific interactions (P < 0.01)')
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            if False:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_abs.{:}'.format(
                            ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_abs.{:}'.format(
                            ext))

        if False:
            print('Focus on the effect of hyperoxia')
            tps = ['P7', 'P7_HO']
            idx = [2, 5]
            ddi = (nints_sym[5] - nints_sym[2])
            ddif = (nints_sym[5] - nints_sym[2]) / (0.5 * (nints_sym[5] + nints_sym[2]))
            fig, ax = plt.subplots(1, 1, figsize=(6.2, 5.2))
            ax.set_xlim(-0.5, len(csts) - 0.5)
            ax.set_ylim(len(csts) - 0.5, -0.5)
            xmax = nints_sym[2].max()

            def color_fun(x, cmap='plasma'):
                cmap = plt.cm.get_cmap(cmap)
                v = 1.0 * x / 175.0
                return cmap(v)

            def size_fun(x):
                return 500. * np.abs(x) / xmax

            for i, csti in enumerate(csts):
                for j, cstj in enumerate(csts):
                    x = nints_sym[2, i, j]
                    dx = ddi[i, j]
                    color = color_fun(x)
                    s = size_fun(dx)
                    marker = '^' if dx > 0 else 'v'
                    ax.scatter([i], [j], s=s, color=color, marker=marker)
            ax.set_xticks(np.arange(len(csts)))
            ax.set_yticks(np.arange(len(csts)))
            ax.set_xticklabels(csts, rotation=90)
            ax.set_yticklabels(csts)
            ax.set_title('Effect of hyperoxia on\ninteraction numbers')

            vals = [0, 50, 100, 175]
            labels = [str(x) for x in vals]
            handles = [ax.scatter([], [], marker='s', color=color_fun(x)) for x in vals]
            ax.legend(
                    handles, labels, loc='upper left', title='# interactions\nin normoxia',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )

            vals = [40, 20, 10, 5, 0, -5, -10, -20, -40, -80]
            labels = [str(x) for x in vals]
            handles = [
                ax.scatter([], [], s=size_fun(x),
                           marker='^' if x > 0 else 'v',
                           c='grey') for x in vals]
            ax2 = ax.twinx()
            ax2.set_axis_off()
            ax2.legend(
                    handles, labels, loc='upper left', title='Difference in\n# interactions\nin hyperoxia',
                    bbox_to_anchor=(1.01, 0.41), bbox_transform=ax2.transAxes,
                    )

            fig.tight_layout()

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_hyperoxia.{:}'.format(
                            ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_hyperoxia.{:}'.format(
                            ext))

    if False:
        print('Heatmap with number of interactions, coarse')
        pvals = {}
        for tp in ['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']:
            print('Timepoint: {:}'.format(tp))
            fn_res = '../../data/cellphonedb/endothelial_coarse/results/{:}/pvalues.txt'.format(tp)
            pval = pd.read_csv(
                fn_res,
                sep='\t',
                index_col=1,
                )
            pvals[tp] = pval

        print('Split by time')
        dsp = ds.split('Timepoint')

        csts = [
                'Lymphatic EC',
                'Arterial EC II',
                'Venous EC',
                'Arterial EC I',
                'Proliferative EC',
                'Car4- capillaries',
                'Car4+ capillaries',
                ]
        alpha = 0.01
        n = len(csts)
        nints = np.ma.masked_all((6, n, n), np.int64)
        for i, ct1 in enumerate(csts):
            for j, ct2 in enumerate(csts):
                for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                    col = ct1+'|'+ct2
                    if col not in pvals[tp]:
                        col = ct2+'|'+ct1
                    if col not in pvals[tp]:
                        continue
                    nints[k, i, j] = (pvals[tp][col] < alpha).sum()
        nints_sym = nints + nints.swapaxes(1, 2)

        if False:
            print('Plot by time in single figure, common normalization')
            fig = plt.figure(
                figsize=(9.2, 7),
                )
            gs = fig.add_gridspec(2, 4, width_ratios=[15] * 3 + [1])
            axs = []
            axs.append(fig.add_subplot(gs[0, 0]))
            axs.append(fig.add_subplot(gs[0, 1]))
            axs.append(fig.add_subplot(gs[0, 2]))
            axs.append(fig.add_subplot(gs[1, 0]))
            axs.append(fig.add_subplot(gs[1, 1]))
            axs.append(fig.add_subplot(gs[1, 2]))
            cax = fig.add_subplot(gs[:, -1])
            vmax = nints_sym.max()
            ns = (ds0.samplesheet
                     .loc[:, ['TimepointHO', 'cellSubtypeCoarse', 'cellType']]
                     .groupby(['TimepointHO', 'cellSubtypeCoarse'])
                     .count()
                     .loc[:, 'cellType']
                     .unstack().T.fillna(0)
                     .astype(int))
            ns['all'] = ns[['E18.5', 'P1', 'P7', 'P21']].sum(axis=1)
            ns = ns[['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']]
            ns = ns.loc[csts]
            for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                ax = axs[k]
                nints_plot = pd.DataFrame(nints_sym[k], index=csts, columns=csts)
                nk = ns[tp]
                for cst, n in nk.items():
                    if n < 1:
                        nints_plot.loc[cst] = 0
                        nints_plot.loc[:, cst] = 0
                sns.heatmap(
                        nints_plot, ax=ax, cmap='plasma',
                        vmin=0,
                        vmax=vmax,
                        cbar=False,
                        )
                ax.set_title(tp)
                xlim = ax.get_xlim()
                ax.set_xlim(xlim[0], xlim[1] + 1)
                yts = ax.get_yticks()
                ytls = axs[0].get_yticklabels()
                for yt, ytl in zip(yts, ytls):
                    n = ns.loc[ytl.get_text(), tp]
                    ax.text(xlim[1] + 0.5, yt, str(n), ha='center', va='center',
                            fontsize=9)
                if k < 3:
                    ax.set_xticklabels([])
                if k in (1, 2, 4, 5):
                    ax.set_yticklabels([])

                if k >= 3:
                    xticks = list(ax.get_xticks())
                    ax.set_xticks(xticks + [xticks[-1] + 1])
                    ax.set_xticklabels(list(csts) + ['number of cells'])

            norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='plasma'), cax=cax)
            cb.set_label('number of interactions', labelpad=15)
            fig.suptitle('     Number of specific interactions (P < 0.01)')
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_abs_coarse.{:}'.format(
                            ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_abs_coarse.{:}'.format(
                            ext))

        if True:
            print('Focus on the effect of hyperoxia, coarse')
            tps = ['P7', 'P7_HO']
            idx = [2, 5]
            ddi = (nints_sym[5] - nints_sym[2])
            ddif = (nints_sym[5] - nints_sym[2]) / (0.5 * (nints_sym[5] + nints_sym[2]))
            fig, ax = plt.subplots(1, 1, figsize=(5.6, 4.6))
            ax.set_xlim(-0.5, len(csts) - 0.5)
            ax.set_ylim(len(csts) - 0.5, -0.5)
            xmax = nints_sym[2].max()

            def color_fun(x, cmap='plasma'):
                cmap = plt.cm.get_cmap(cmap)
                v = 1.0 * x / 175.0
                return cmap(v)

            def size_fun(x):
                return 500. * np.abs(x) / xmax

            for i, csti in enumerate(csts):
                for j, cstj in enumerate(csts):
                    x = nints_sym[2, i, j]
                    dx = ddi[i, j]
                    color = color_fun(x)
                    s = size_fun(dx)
                    marker = '^' if dx > 0 else 'v'
                    ax.scatter([i], [j], s=s, color=color, marker=marker)
            ax.set_xticks(np.arange(len(csts)))
            ax.set_yticks(np.arange(len(csts)))
            ax.set_xticklabels(csts, rotation=90)
            ax.set_yticklabels(csts)
            ax.set_title('Effect of hyperoxia on\ninteraction numbers')

            vals = [0, 50, 100, 175]
            labels = [str(x) for x in vals]
            handles = [ax.scatter([], [], marker='s', color=color_fun(x)) for x in vals]
            ax.legend(
                    handles, labels, loc='upper left', title='# interactions\nin normoxia',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )

            vals = [40, 20, 10, 5, 0, -5, -10, -20, -40, -80]
            labels = [str(x) for x in vals]
            handles = [
                ax.scatter([], [], s=size_fun(x),
                           marker='^' if x > 0 else 'v',
                           c='grey') for x in vals]
            ax2 = ax.twinx()
            ax2.set_axis_off()
            ax2.legend(
                    handles, labels, loc='upper left', title='Difference in\n# interactions\nin hyperoxia',
                    bbox_to_anchor=(1.01, 0.51), bbox_transform=ax2.transAxes,
                    )

            fig.tight_layout()

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_hyperoxia_coarse.{:}'.format(
                            ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_signaling_heatmap_hyperoxia_coarse.{:}'.format(
                            ext))


    if False:
        print('Plot marginals')
        from scipy import interpolate
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
        fig, ax = plt.subplots(figsize=(6.5, 3))
        # Last one is hyperoxia
        nints_marg = nints_sym[[0, 1, 2, 3]].sum(axis=1).T
        for i, cst in enumerate(csts):
            x = np.array([0, 1, 2, 3])
            y = nints_marg[i]
            x = x[~(y.mask)]
            y = y[~(y.mask)]
            outx = np.linspace(x[0], x[-1], 100)
            outy = interpolate.pchip_interpolate(x, y, outx)
            ax.scatter(x, y, color=cmap[cst], label=cst)
            ax.plot(outx, outy, color=cmap[cst], lw=2)

        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(['E18.5', 'P1', 'P7', 'P21'])
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('Number of interactions')
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
            title='Cell subtype:',
            )
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal.{:}'.format(
                        ext))

    if False:
        print('Plot marginals, coarse')
        pvals = {}
        for tp in ['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']:
            print('Timepoint: {:}'.format(tp))
            fn_res = '../../data/cellphonedb/endothelial_coarse/results/{:}/pvalues.txt'.format(tp)
            pval = pd.read_csv(
                fn_res,
                sep='\t',
                index_col=1,
                )
            pvals[tp] = pval

        print('Split by time')
        dsp = ds.split('Timepoint')

        csts = [
                'Lymphatic EC',
                'Arterial EC II',
                'Venous EC',
                'Arterial EC I',
                'Proliferative EC',
                'Car4- capillaries',
                'Car4+ capillaries',
                ]
        alpha = 0.01
        n = len(csts)
        nints = np.ma.masked_all((6, n, n), np.int64)
        for i, ct1 in enumerate(csts):
            for j, ct2 in enumerate(csts):
                for k, tp in enumerate(['E18.5', 'P1', 'P7', 'P21', 'all', 'P7_HO']):
                    col = ct1+'|'+ct2
                    if col not in pvals[tp]:
                        col = ct2+'|'+ct1
                    if col not in pvals[tp]:
                        continue
                    nints[k, i, j] = (pvals[tp][col] < alpha).sum()
        nints_sym = nints + nints.swapaxes(1, 2)

        from scipy import interpolate
        cmap = {
            'Lymphatic EC': '#ffda4f',
            'Venous EC': '#4380c9',
            'Arterial EC I': '#bd4f6f',
            'Arterial EC II': '#ee4242',
            'Proliferative EC': '#664fbd',
            'Car4- capillaries': '#4fbd9d',
            'Car4+ capillaries': '#3ab646',
        }
        csts = [
            'Lymphatic EC',
            'Arterial EC II',
            'Venous EC',
            'Arterial EC I',
            'Proliferative EC',
            'Car4- capillaries',
            'Car4+ capillaries',
        ]
        fig, axs = plt.subplots(figsize=(6.5, 3))
        # Last one is hyperoxia
        nints_marg = nints_sym[[0, 1, 2, 3, 4]].sum(axis=1).T.data
        for i, cst in enumerate(csts):
            x = np.array([0, 1, 2, 3])
            y = nints_marg[i][x]
            outx = np.linspace(x[0], x[-1], 100)
            outy = interpolate.pchip_interpolate(x, y, outx)
            ax.scatter(x, y, color=cmap[cst], label=cst)
            ax.plot(outx, outy, color=cmap[cst], lw=2)

        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(['E18.5', 'P1', 'P7', 'P21'])
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('Number of interactions')
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
            title='Cell subtype:',
            )
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal_coarse.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal_coarse.{:}'.format(
                        ext))

    if False:
        print('Plot marginals with hyperoxia')
        from scipy import interpolate
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
        # Last one is hyperoxia
        nints_marg = nints_sym[[2, 5]].sum(axis=1).T
        x = np.array([0, 1])
        for i, cst in enumerate(csts):
            y = nints_marg[i]
            x = x[~(y.mask)]
            y = y[~(y.mask)]
            ax.scatter(x, y, color=cmap[cst], label=cst)
            ax.plot(x, y, color=cmap[cst], lw=2)

        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(['normoxia', 'hyperoxia'])
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('Number of interactions')
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
            title='Cell subtype:',
            )
        ax.set_xlim(-0.2, 1.2)
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal_hyperoxia.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_number_interactions_by_timepoint_marginal_hyperoxia.{:}'.format(
                        ext))

    if False:
        print('Distributions of Cxcl12 in Arterial cells, Cxcr4 in Basophils, Ackr3 in Venous')
        gene_cells = [
            ('Cxcl12', 'Arterial EC'),
            ('Cxcl12', 'Car4- capillaries'),
            #('Cxcr4', 'basophil'),
            #('Ackr3', 'Venous EC'),
        ]
        dsf = ds.query_features_by_name(['Cxcl12', 'Cxcr4', 'Ackr3'])
        dsf.obs['cellSubtypeCoarse'] = dsf.obs['cellSubtypeCoarse'].replace({
            'Arterial EC I': 'Arterial EC',
            'Arterial EC II': 'Arterial EC',
        })
        dsf_imm = ds_imm.query_features_by_name(['Cxcl12', 'Cxcr4', 'Ackr3'])
        tps = ['E18.5', 'P1', 'P7', 'P21']
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
        }
        fig, axs = plt.subplots(1, 2, figsize=(5, 2.3), sharex=True, sharey=True)
        for tp in tps:
            for (gene, ct), ax in zip(gene_cells, axs):
                if ct == 'basophil':
                    dsi = dsf_imm
                else:
                    dsi = dsf
                dsict = dsi.query_samples_by_metadata(
                        '(TimepointHO == @tp) & (cellSubtypeCoarse == @ct)', local_dict=locals(),
                        )
                if False:
                    x = 0.1 + np.sort(dsict.counts.loc[gene].values)
                    y = 1.0 - np.linspace(0, 1, len(x))
                else:
                    from scipy.stats import gaussian_kde
                    xpoints = dsict.counts.loc[gene].values + 0.1
                    if len(xpoints) <= 1:
                        continue
                    x = np.logspace(-1, 4)
                    y = gaussian_kde(np.log10(xpoints))(np.log10(x))
                ax.plot(x, y, lw=2, color=cmap[tp], label=tp)
                ax.grid(True)
                ax.set_xscale('log')
                if tp == tps[0]:
                    genet = gene.capitalize()
                    ax.set_title(f'{ct}')
                    if (gene, ct) == gene_cells[0]:
                        ax.set_ylabel('Density of cells')
        axs[0].legend(fontsize=9)
        #axs[1].arrow(0.04, 0.6, 0, -0.4, color='k', head_width=0.05, transform=ax.transAxes)
        #axs[1].arrow(0.7, 0.1, 0.05, 0.3, color='k', head_width=0.05, transform=ax.transAxes)
        fig.text(0.52, 0.03, 'Cxcl12 expression [cpm]', ha='center')
        fig.tight_layout(rect=(0, 0.05, 1, 1))

        if True:
            fxf = f'{fig_fdn}/distributions_Cxcl12'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(f'{fxf}.{ext}', dpi=dpi)
                else:
                    fig.savefig(f'{fxf}.{ext}')

    if False:
        print('Add panel showing basophils and endos')
        ct1 = 'basophil'
        #tp = 'E18.5'
        gene_cells = [
            ('Rspo3', 'Lgr5', 'Proliferative EC'),
            ('Ccl4', 'Ackr2', 'Proliferative EC'),
            ('Kit', 'Kitl', 'Proliferative EC'),
        ]
        #dsi1 = ds_imm.query_samples_by_metadata('(cellSubtype == @ct1) & (TimepointHO == @tp)', local_dict=locals())
        dsi1 = ds_imm.query_samples_by_metadata('(cellSubtype == @ct1)', local_dict=locals())
        fig, axs = plt.subplots(1, len(gene_cells), figsize=(3 * len(gene_cells), 3), sharey=True, sharex=True)
        for (gene1, gene2, ct2), ax in zip(gene_cells, axs):
            x1 = 0.1 + np.sort(dsi1.counts.loc[gene1].values)
            y1 = 1.0 - np.linspace(0, 1, len(x1))
            #dsi2 = ds0.query_samples_by_metadata('(cellSubtype == @ct2) & (TimepointHO == @tp)', local_dict=locals())
            dsi2 = ds0.query_samples_by_metadata('(cellSubtype == @ct2)', local_dict=locals())
            x2 = 0.1 + np.sort(dsi2.counts.loc[gene2].values)
            y2 = 1.0 - np.linspace(0, 1, len(x2))
            ax.plot(x1, y1, color='grey', lw=2, label=f'{gene1} in {ct1}')
            ax.plot(x2, y2, color='tomato', lw=2, label=f'{gene2} in {ct2}')
            ax.legend()
            ax.grid()
            ax.set_xscale('log')
            ax.set_yscale('log')
        fig.tight_layout()

    if False:
        print('Lgr5 in Proliferative ECs')
        gene_cells = [
            ('Lgr5', 'Proliferative EC'),
        ]
        dsf = ds0.query_features_by_name([g for g, _ in gene_cells])
        tps = ['E18.5', 'P1', 'P7', 'P7_HO']
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            'P7_HO': 'purple'
        }
        fig, axs = plt.subplots(1, len(gene_cells), figsize=(3 * len(gene_cells), 3), sharex=True, sharey=True)
        if len(gene_cells) == 1:
            axs = [axs]
        for tp in tps:
            for (gene, ct), ax in zip(gene_cells, axs):
                dsict = dsf.query_samples_by_metadata(
                        '(TimepointHO == @tp) & (cellSubtype == @ct)', local_dict=locals(),
                        )
                x = 0.1 + np.sort(dsict.counts.loc[gene].values)
                y = 1.0 - np.linspace(0, 1, len(x))
                ax.plot(x, y, lw=2, color=cmap[tp], label=tp)
                ax.grid(True)
                ax.set_xscale('log')
                ax.set_yscale('log')
                if tp == tps[0]:
                    genet = gene.capitalize()
                    ax.set_title(f'{genet} in {ct}')
                    if (gene, ct) == gene_cells[0]:
                        ax.set_ylabel('Fraction of cell\nexpressing > x')
        axs[0].legend(fontsize=9)
        fig.text(0.52, 0.03, 'Gene expression [cpm]', ha='center')
        fig.tight_layout(rect=(0, 0.05, 1, 1))

    if True:
        print('Horizontal bar plots')
        tps = ['E18.5', 'P1', 'P7', 'P21']
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            'P7_HO': 'purple'
        }

        # Cxcr4, immune
        gene = 'Cxcr4'
        dsf_imm = ds_imm.query_features_by_name([gene])
        dsf_imm.counts.log(inplace=True)
        dsp = dsf_imm.split(['TimepointHO', 'cellSubtype'])

        fig, axs = plt.subplots(4, 1, figsize=(2.8, 6), sharex=True)
        for tp, ax in zip(tps, axs):
            ax.set_title(tp)
            dspi = {k[1]: v for k, v in dsp.items() if k[0] == tp}
            dspi = {k: v for k, v in dspi.items() if v.n_samples >= 3}
            avgs = {k: v.counts.loc[gene].mean() for k, v in dspi.items()}
            stds = {k: v.counts.loc[gene].std() for k, v in dspi.items()}
            df = pd.DataFrame([avgs, stds], index=['avg', 'std']).T
            df['ncells'] = pd.Series({k: v.n_samples for k, v in dspi.items()})
            df['serr'] = df['std'] / np.sqrt(df['ncells'] - 1)
            top = df.nlargest(5, 'avg')
            y = len(top) - np.arange(len(top))
            x = top['avg'].values
            ax.barh(y, x + 1, height=0.9, color=cmap[tp])
            dx = top['serr'].values
            for i in range(len(top)):
                ax.plot([x[i] - dx[i] + 1, x[i] + dx[i] + 1], [y[i]] * 2, lw=2, color='k')
            ax.set_yticks(y)
            ax.set_yticklabels(top.index)
            ax.grid(True, axis='x')
            ax.set_xlim(0, 2.8)
            ax.set_xticks([0, 1, 2])
        axs[-1].set_xticklabels(['$0$', '$1$', '$10$'])
        axs[-1].set_xlabel(f'{gene} expr [cpm]')
        fig.tight_layout()

        if True:
            fxf = f'{fig_fdn}/bar_horizontal_{gene}'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(f'{fxf}.{ext}', dpi=dpi)
                else:
                    fig.savefig(f'{fxf}.{ext}')

        # Ackr3, endothelial
        gene = 'Ackr3'
        dsf = ds.query_features_by_name([gene])
        dsf.counts.log(inplace=True)
        dsp = dsf.split(['TimepointHO', 'cellSubtype'])

        fig, axs = plt.subplots(4, 1, figsize=(2.8, 3.6), sharex=True)
        for tp, ax in zip(tps, axs):
            ax.set_title(tp)
            dspi = {k[1]: v for k, v in dsp.items() if k[0] == tp}
            dspi = {k: v for k, v in dspi.items() if v.n_samples >= 3}
            avgs = {k: v.counts.loc[gene].mean() for k, v in dspi.items()}
            stds = {k: v.counts.loc[gene].std() for k, v in dspi.items()}
            df = pd.DataFrame([avgs, stds], index=['avg', 'std']).T
            df['ncells'] = pd.Series({k: v.n_samples for k, v in dspi.items()})
            df['serr'] = df['std'] / np.sqrt(df['ncells'] - 1)
            top = df.nlargest(2, 'avg')
            y = len(top) - np.arange(len(top))
            x = top['avg'].values
            ax.barh(y, x + 1, height=0.9, color=cmap[tp])
            dx = top['serr'].values
            for i in range(len(top)):
                ax.plot([x[i] - dx[i] + 1, x[i] + dx[i] + 1], [y[i]] * 2, lw=2, color='k')
            ax.set_yticks(y)
            ax.set_yticklabels(top.index)
            ax.grid(True, axis='x')
            ax.set_xlim(0, 4.5)
            ax.set_xticks([0, 1, 2, 3, 4])
        axs[-1].set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$'])
        axs[-1].set_xlabel(f'{gene} expr [cpm]')
        fig.tight_layout()

        if True:
            fxf = f'{fig_fdn}/bar_horizontal_{gene}'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(f'{fxf}.{ext}', dpi=dpi)
                else:
                    fig.savefig(f'{fxf}.{ext}')

    plt.ion()
    plt.show()
