# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 3 for the endo paper: macrovasculature
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_3/'


if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(cellType == "endothelial") & (doublet == 0)',
        inplace=True)
    print('Total endothelial cells analyzed: {:}'.format(ds0.n_samples))
    ds = ds0.query_samples_by_metadata('Treatment == "normal"')

    vs = ds.samplesheet[['embed_endo_1', 'embed_endo_2']].copy()

    if False:
        print('Plot embedding of genes')
        genes = [
            'Thy1',
            'Prox1',
            'Ccl21a',
            'Vwf',
            'Car8',
            'Mmp2',
            'Gja5',
            'Bmx',
            'Mgp',
            'Dkk2',
            'Kcne3',
            ]

        for gene in genes:
            title = gene
            fig, ax = plt.subplots(figsize=(2.5, 2.4))
            cmap = 'viridis'
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.40,
                    cmap=cmap,
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            ax.text(0.95, 0.95, gene, ha='right', va='top',
                    fontsize=12,
                    transform=ax.transAxes)
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

    if False:
        print('Find markers for Arterial I vs II')
        ds1 = ds.query_samples_by_metadata('cellSubtype == "Arterial EC I"')
        ds2 = ds.query_samples_by_metadata('cellSubtype == "Arterial EC II"')
        comp = ds2.compare(ds1)
        comp.rename(
                columns={'avg_self': 'avg_II', 'avg_other': 'avg_I'},
                inplace=True,
                )
        deg = comp.loc[comp['statistic'] >= 0.4].sort_values('log2_fold_change')
        deg.to_csv(fig_fdn+'DEGs_arterial_II_I.tsv', sep='\t', index=True)

    if False:
        print('Pathway analysis on DEGs between Art I and II')
        deg = pd.read_csv(fig_fdn+'DEGs_arterial_II_I.tsv', sep='\t', index_col=0)
        print('Up in II')
        print(','.join(deg.loc[deg['log2_fold_change'] > 0].index))
        print('Down in II')
        print(','.join(deg.loc[deg['log2_fold_change'] < 0].index))

    if False:
        print('Plot time sensitive and cell subtype specific genes')
        dsi = ds.query_samples_by_metadata('cellSubtype in ("Arterial EC I", "Arterial EC II")')
        dsip = dsi.split(['cellSubtype', 'Timepoint'])

        #deg = pd.read_csv(
        #        fig_fdn+'DEGs_arterial_II_I.tsv',
        #        sep='\t',
        #        index_col=0,
        #        )
        ## Top Art II markers
        #deg1 = deg.index[:24]
        ## Top Art I markers
        #deg2 = deg.index[-24:][::-1]

        deg1 = ['Gpihbp1', 'Tbx3', 'Mest', 'Kit']
        deg2 = ['Cytl1', 'Efna5', 'Polr3h', 'Ctsh', 'Serpinf1', 'Dkk2',
                'Lmo1', 'Clu']
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }

        for ig, genes in enumerate([deg1, deg2]):
            cst = 'Arterial EC '+('I' * (ig + 1))
            nrows = 1 + 1 * (len(genes) > 4)
            fig, axs = plt.subplots(
                    nrows, 4,
                    figsize=(7, 0.2 + 1.7 * nrows), sharex=True, sharey=True)
            axs = axs.ravel()
            for i in range(len(axs)):
                gene = genes[i]
                ax = axs[i]
                for tp in ['E18.5', 'P1', 'P7', 'P21']:
                    x = sorted(dsip[(cst, tp)].counts.loc[gene] + 0.1)
                    y = list(1.0 - np.linspace(0, 1, len(x)))
                    x.insert(0, 0.1)
                    y.insert(0, 1)
                    x.append(x[-1])
                    y.append(0)
                    if i == 0:
                        label = tp
                    else:
                        label = ''
                    ax.plot(x, y, lw=2, label=label, color=cmap[tp])
                ax.set_xlim(0.09, 1e5)
                ax.set_xscale('log')
                ax.set_ylim(-0.02, 1.02)
                ax.set_title(gene)
                ax.grid(True)
            fig.legend(
                    loc='upper left', fontsize=10, title='Timepoint:',
                    bbox_to_anchor=(0.83, 0.92),
                    bbox_transform=fig.transFigure)
            #fig.suptitle('Upregulated in {:}'.format(cst))
            fig.text(0.45, 0.02, 'Gene expression [cpm]', ha='center', va='bottom')
            fig.tight_layout(rect=(0, 0.06, 0.83, 1))

            if True:
                inf = 'I' * (ig + 1)
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_Art{:}_cum_time.{:}'.format(
                            inf, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_Art{:}_cum_time.{:}'.format(
                                inf, ext))


    if False:
        print('Figure out what gene of the Art II/I are only at E18.5')
        deg = pd.read_csv(
                fig_fdn+'DEGs_arterial_II_I.tsv',
                sep='\t',
                index_col=0,
                )

        ds1e = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC I") & (Timepoint in ("E18.5", "P1"))')
        ds1l = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC I") & (Timepoint in ("P7", "P21"))')
        ds2e = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC II") & (Timepoint in ("E18.5", "P1"))')
        ds2l = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC II") & (Timepoint in ("P7", "P21"))')

        # Top Art II markers
        deg1 = deg.index[:24]
        # Top Art I markers
        deg2 = deg.index[-24:][::-1]

        for ig, genes in enumerate([deg1, deg2]):
            fig, axs = plt.subplots(3, 8, figsize=(14, 6), sharex=True, sharey=True)
            axs = axs.ravel()
            for i in range(len(axs)):
                gene = genes[i]
                ax = axs[i]
                x = np.sort(ds1e.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art I Early')
                x = np.sort(ds1l.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art I Late')
                x = np.sort(ds2e.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art II Early')
                x = np.sort(ds2l.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art II Late')
                ax.set_xlim(0.09, 1e5)
                ax.set_xscale('log')
                ax.set_ylim(-0.02, 1.02)
                ax.set_title(gene)
                ax.grid(True)
                if (ig == 1) & (i == 0):
                    ax.legend(fontsize=6)
            fig.suptitle('Upregulated in Arterial {:}'.format('I'*(ig+1)))
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_ArtI_II_early_late_{:}.{:}'.format(
                            ig+1, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_ArtI_II_early_late_{:}.{:}'.format(
                            ig+1, ext))

    if False:
        print('Figure out what genes of the Art II/I are P1/P7')
        deg = pd.read_csv(
                fig_fdn+'DEGs_arterial_II_I.tsv',
                sep='\t',
                index_col=0,
                )

        ds1e = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC I") & (Timepoint == "P1")')
        ds1l = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC I") & (Timepoint == "P7")')
        ds2e = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC II") & (Timepoint == "P1")')
        ds2l = ds.query_samples_by_metadata(
            '(cellSubtype == "Arterial EC II") & (Timepoint == "P7")')

        # Top Art II markers
        deg1 = deg.index[:24]
        # Top Art I markers
        deg2 = deg.index[-24:][::-1]

        for ig, genes in enumerate([deg1, deg2]):
            fig, axs = plt.subplots(3, 8, figsize=(14, 6), sharex=True, sharey=True)
            axs = axs.ravel()
            for i in range(len(axs)):
                gene = genes[i]
                ax = axs[i]
                x = np.sort(ds1e.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art I P1')
                x = np.sort(ds1l.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art I P7')
                x = np.sort(ds2e.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art II P1')
                x = np.sort(ds2l.counts.loc[gene]) + 0.1
                ax.plot(x, 1.0 - np.linspace(0, 1, len(x)), lw=2, label='Art II P7')
                ax.set_xlim(0.09, 1e5)
                ax.set_xscale('log')
                ax.set_ylim(-0.02, 1.02)
                ax.set_title(gene)
                ax.grid(True)
                if (ig == 1) & (i == 0):
                    ax.legend(fontsize=6)
            fig.suptitle('Upregulated in Arterial {:}'.format('I'*(ig+1)))
            fig.tight_layout(rect=(0, 0, 1, 0.95))

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'endo_ArtI_II_P1_P7_{:}.{:}'.format(
                            ig+1, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_ArtI_II_P1_P7_{:}.{:}'.format(
                            ig+1, ext))

    if False:
        print('Heatmap Art I vs II')
        gene_table = pd.read_excel(
            #'../../data/gene_lists/endo_Arterial EC II_vs_Arterial EC I heat Map.xlsx',
            '../../data/gene_lists/endo_Arterial EC II_vs_Arterial EC I heat map revised.xlsx',
            #'endo_Arterial EC II_vs_Arterial',
            'Sheet1',
            index_col='GeneName',
            squeeze=True,
            )
        # NOTE: Cristina reclassified Fn1, ask whether it has dual function
        #gene_table['Fn1'] = 'Extracellular Matrix'
        # Reorder: this is simplest
        order = [
            'Plvap', 'Ace', 'Aqp1', 'Bmpr2', 'Foxf1', 'Ets1', 'Cd36', 'Robo4', 'Kit', 'Emp2',
            'Calcrl', 'Ly6e', 'Adgrf5', 'Tmem204', 'Tbx3', 'Car2', 'Psen1',
            'Jag1', 'Igf2', 'Fn1', 'Serpinf1', 'Cdh13',
            'Col4a1', 'Col4a2', 'Lama3',
            'Col5a2', 'Col18a1', 'Bgn', 'Spint2', 'Kazald1',
            'Fbln5', 'Fbn1', 'Eln', 'Fbln2', 'Lox',
            'Cldn5', 'Ltbp4', 'Ltbp3', 'Htra1', 'Prdm16',# 'Fos',
            'Cd63', 'Cd59a', 'Itgb4', 'Efna5', 'Lamb2',
            ]
        gene_table = gene_table.loc[order]
        dsi = ds.query_samples_by_metadata('cellSubtype in ("Arterial EC I", "Arterial EC II")')
        dsi.counts.log(inplace=True)
        dsia = dsi.average('samples', 'cellSubtype')
        dsi.counts.unlog(inplace=True)

        genes = gene_table.index.tolist()

        data = dsia.counts.loc[genes].T
        vmin = -1
        vmax = data.values.max()

        fig, ax = plt.subplots(figsize=(10, 3))
        sns.heatmap(
            data,
            ax=ax,
            cmap='plasma',
            vmin=vmin,
            vmax=vmax,
            fmt='.1f',
            xticklabels=True,
            yticklabels=True,
            cbar=False,
            )
        for tk in ax.get_yticklabels():
            tk.set_rotation(0)
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)
        ax.set_xlim(0, len(genes))
        ax.set_ylabel('')
        y0, y1 = ax.get_ylim()

        pathway = None
        cuts = []
        mids = []
        for i, gene in enumerate(genes):
            pw = gene_table[gene]
            if pathway is None:
                pathway = pw
                continue
            if pw != pathway:
                if len(cuts) == 0:
                    mids.append(0.5 * i)
                else:
                    mids.append(0.5 * (i + cuts[-1]))
                cuts.append(i)
                pathway = pw
        mids.append(0.5 * (i + 1 + cuts[-1]))
        for cut in cuts:
            ax.plot([cut] * 2, [y0, y1 + (y1 - y0) * 1.5], lw=2, color='grey', clip_on=False)

        fig.tight_layout(rect=(0, 0, 0.9, 0.55))

        vals = [-1, 0, 1, 2, 3]
        labels = ['$0$', '$1$', '$10$', '$10^2$', '$10^3$']
        handles = []
        cmap = plt.cm.get_cmap('plasma')
        for val in vals:
            color = cmap(1.0 * (val - vmin) / (vmax - vmin))
            r = plt.Rectangle((0, 0), 0, 0, facecolor=color)
            handles.append(r)
        ax.legend(
                handles, labels,
                loc='upper left', title='Gene\nexpr\n[cpm]',
                bbox_to_anchor=(1.01, 2.5), bbox_transform=ax.transAxes,
                )

        for mid in mids:
            txt = gene_table.values[int(mid)]
            fields = txt.split()
            txt = ''
            for field in fields:
                if len(txt.split('\n')[-1]+' '+field) > 15:
                    txt += '\n'+field
                else:
                    txt += ' '+field
            ax.text(mid, y1 + 0.1 * (y1 - y0), txt, rotation=90, clip_on=False, ha='center', va='bottom')

        if True:
            fxf = 'arterial_heatmap'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'{:}.{:}'.format(
                        fxf, ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'{:}.{:}'.format(
                        fxf, ext))

    if False:
        print('Plot abundances in time of macrovascular')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)
        frac = frac.loc[['Venous EC', 'Lymphatic EC', 'Arterial EC I', 'Arterial EC II']]
        sts = list(frac.index)

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(4, 2))
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
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        labeld = {
            'Arterial EC I': 'Art I',
            'Arterial EC II': 'Art II',
            'Venous EC': 'Venous',
            'Lymphatic EC': 'Lymphatic',
            }
        for ist, st in enumerate(sts):
            y = frac.loc[st, xorder]

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
            out = np.vstack([outx, outy])
            if ist >= 10:
                ls = '--'
                m = 's'
            else:
                ls = '-'
                m = 'o'

            ax.scatter(
                    x, y + 0.1,
                    marker=m,
                    lw=2, alpha=0.8,
                    edgecolor=cmap[st],
                    facecolor='none',
                    zorder=10 - ist,
                    )
            ax.plot(out[0], out[1] + 0.1, lw=2,
                    alpha=0.4,
                    color=cmap[st], label=labeld[st],
                    zorder=10 - ist,
                    ls=ls,
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.legend(
                title='Cell subtype:',
                bbox_to_anchor=(1.04, 1.02),
                bbox_transform=ax.transAxes,
                loc='upper left',
                fontsize=9,
                ncol=1,
                )
        ax.set_ylim(0.1, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of\nendothelial cells')
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(['0.1%', '1%', '10%', '100%'])

        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'macrovascular_abundances.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'macrovascular_abundances.{:}'.format(
                        ext))

    if True:
        print('Plot abundances in time of macrovascular (among themselves)')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac = frac.loc[['Venous EC', 'Lymphatic EC', 'Arterial EC I', 'Arterial EC II']]
        frac *= 100 / frac.sum(axis=0)
        sts = list(frac.index)

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(4, 2))
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
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        labeld = {
            'Arterial EC I': 'Art I',
            'Arterial EC II': 'Art II',
            'Venous EC': 'Venous',
            'Lymphatic EC': 'Lymphatic',
            }
        for ist, st in enumerate(sts):
            y = frac.loc[st, xorder]

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
            out = np.vstack([outx, outy])
            if ist >= 10:
                ls = '--'
                m = 's'
            else:
                ls = '-'
                m = 'o'

            ax.scatter(
                    x, y + 0.1,
                    marker=m,
                    # Lymphatic are gold, need a little more contrast
                    lw=2, alpha=0.8 + 0.2 * ('Lymphatic' in st),
                    edgecolor=cmap[st],
                    facecolor='none',
                    zorder=10 - ist - 3 * ('Lymphatic' in st),
                    )
            ax.plot(out[0], out[1] + 0.1, lw=2,
                    # Lymphatic are gold, need a little more contrast
                    alpha=0.6 + 0.4 * ('Lymphatic' in st),
                    color=cmap[st], label=labeld[st],
                    zorder=10 - ist - 3 * ('Lymphatic' in st),
                    ls=ls,
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.legend(
                title='Cell subtype:',
                bbox_to_anchor=(1.04, 1.02),
                bbox_transform=ax.transAxes,
                loc='upper left',
                fontsize=9,
                ncol=1,
                )
        ax.set_ylim(2, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of\nmacrovascular cells')
        ax.set_yticks([3, 10, 100])
        ax.set_yticklabels(['3%', '10%', '100%'])

        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'macrovascular_abundances_among.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'macrovascular_abundances_among.{:}'.format(
                        ext))

    if False:
        print('Plot abundances in time of Art I and II, with bars')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)
        frac = frac.loc[['Arterial EC I', 'Arterial EC II']]
        frac.loc['Other'] = 100 - frac.sum(axis=0)
        sts = list(frac.index)

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(3, 2))
        colors = sns.color_palette('muted', n_colors=len(sts))
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        labeld = {'Arterial EC I': 'Art I', 'Arterial EC II': 'Art II'}
        for ist, st in enumerate(sts):
            y = frac.loc[st, xorder]

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
            out = np.vstack([outx, outy])
            if ist >= 10:
                ls = '--'
                m = 's'
            else:
                ls = '-'
                m = 'o'

            ax.scatter(
                    x, y + 0.1,
                    marker=m,
                    lw=2, alpha=0.8,
                    edgecolor=colors[ist],
                    facecolor='none',
                    zorder=10 - ist,
                    )
            ax.plot(out[0], out[1] + 0.1, lw=2,
                    alpha=0.4,
                    color=colors[ist], label=labeld[st],
                    zorder=10 - ist,
                    ls=ls,
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.legend(
                #title='Cell subtype:',
                #bbox_to_anchor=(1.04, 1.02),
                #loc="upper left",
                loc='upper right',
                fontsize=9,
                ncol=2,
                )
        ax.set_ylim(0.1, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of\nendothelial cells')
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(['0.1%', '1%', '10%', '100%'])

        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'arterial_abundances.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'arterial_abundances.{:}'.format(
                        ext))

    if False:
        print('Plot RNA velocity')
        csts = [
            'Arterial EC I',
            'Arterial EC II',
        ]
        ds12 = ds.query_samples_by_metadata(
                'cellSubtype in @csts',
                local_dict=locals(),
                inplace=False)

        print('Load combined velocity file')
        fn_combined = '../../data/sequencing/datasets/all_{:}/velocity_endo.loom'.format(version)
        import scvelo as scv
        adata_endo = scv.read(fn_combined, cache=True)
        adata_endo.var_names_make_unique()

        print('Restrict to subtypes')
        adata = adata_endo[ds12.samplenames]

        print('Follow tutorial')
        # show proportions of spliced/unspliced abundances
        #scv.utils.show_proportions(adata)

        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=25, n_neighbors=10)

        scv.tl.velocity(adata)

        scv.tl.velocity_graph(adata)

        # Use embedding from the paper
        adata.obsm['X_tsne'] = ds12.samplesheet.loc[adata.obs_names, ['embed_endo_1', 'embed_endo_2']].values

        cmap = {
            'Arterial EC I': '#bd4f6f',
            'Arterial EC II': '#ee4242',
        }
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        ds12.plot.scatter_reduced(
            ('embed_endo_1', 'embed_endo_2'),
            color_by='cellSubtype',
            cmap=cmap,
            color_log=False,
            ax=ax,
            s=30,
            alpha=0.7,
            )
        scv.pl.velocity_embedding_stream(
            adata,
            basis='tsne',
            size=0,
            ax=ax,
            alpha=0.5,
            )
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(top=2.5)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'arterial_velocity'
            for ext in ['svg', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    if False:
        print('Plot embedding of time')
        csts = [
            'Arterial EC I',
            'Arterial EC II',
        ]
        ds12 = ds.query_samples_by_metadata(
                'cellSubtype in @csts',
                local_dict=locals(),
                inplace=False)
        gene, title = 'Timepoint', 'timepoint'
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        fig, ax = plt.subplots(figsize=(2.5, 2.4))
        ds12.plot.scatter_reduced_samples(
                ('embed_endo_1', 'embed_endo_2'),
                ax=ax,
                s=12,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.grid(False)
        ax.text(0.95, 0.95, gene, ha='right', va='top',
                fontsize=12,
                transform=ax.transAxes)
        ax.set_axis_off()
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'endo_embedding_timepoint'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    plt.ion()
    plt.show()
