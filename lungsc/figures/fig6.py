# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 6 for the endo paper: proliferative ECs.
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_5/'

if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds = DatasetLung.load(preprocess=True, include_hyperoxia=False, version=version)
    ds.query_samples_by_metadata(
            '(cellType == "endothelial") & (Treatment == "normal") & (doublet == 0)',
            inplace=True)

    ds.samplesheet['is_prolif'] = ds.samplesheet['cellSubtype'] == 'Proliferative EC'

    #print('Load tSNE from file')
    #vs = pd.read_csv(
    #    '../../data/sequencing/datasets/all_{:}/tsne_with_hyperoxia_endo.tsv'.format(version),
    #    sep='\t',
    #    index_col=0,
    #    )
    #vs = vs.loc[ds.samplenames]
    vs = ds.samplesheet[['embed_endo_1', 'embed_endo_2']].copy()
    vs.columns = ['dimension 1', 'dimension 2']

    if True:
        print('Plot embedding of interesting genes')
        genes = [
            'Mki67',
            'Hmmr',
            'Neil3',
            'Mybl2',
            'Depp1',
            'Sparcl1',
            'Rgcc',
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
                        fig.savefig(fig_fdn+'endo_embedding_gene_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_embedding_gene_{:}.{:}'.format(
                            gene, ext))

    if False:
        print('Count proliferative fraction at each time')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        print('Plot fraction trajectories with emphasis on proliferative cells')
        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(3.5, 2.7))
        st = 'Proliferative EC'
        color = 'r'
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        y = frac.loc[st, xorder]

        outx = np.linspace(0, 3, 100)
        outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
        out = np.vstack([outx, outy])

        ax.scatter(
                x, y + 0.1, lw=2, alpha=0.8,
                edgecolor=color,
                facecolor='none',
                zorder=10,
                )
        ax.plot(out[0], out[1] + 0.1, lw=2,
                alpha=0.4 + 0.3 * (st == 'proliferative'),
                color=color, label=st,
                zorder=10,
                )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of EC')
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(['0.1%', '1%', '10%', '100%'])
        ax.set_ylim(0.8, 101)
        xmin, xmax = ax.get_xlim()
        ax.plot([-1, 5], [10**0.5] * 2, lw=1, color='k', alpha=0.2)
        ax.plot([-1, 5], [10**1.5] * 2, lw=1, color='k', alpha=0.2)
        ax.set_xlim(xmin, xmax)

        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_proliferative_abundances.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_proliferative_abundances.{:}'.format(
                        ext))

    if False:
        print('Count all cell type fractions at each time')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        print('Plot fraction trajectories with emphasis on proliferative cells')
        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(7, 2.7))
        sts = frac.index.tolist()
        if '' in sts:
            sts.remove('')
        sts.remove('Proliferative EC')
        sts.insert(0, 'Proliferative EC')
        colors = ['r'] + sns.color_palette('muted', n_colors=len(sts) - 1)
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        for ist, st in enumerate(sts):
            y = frac.loc[st, xorder]

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
            out = np.vstack([outx, outy])

            ax.scatter(
                    x, y + 0.1, lw=2, alpha=0.8,
                    edgecolor=colors[ist],
                    facecolor='none',
                    zorder=10 - ist,
                    )
            ax.plot(out[0], out[1] + 0.1, lw=2,
                    alpha=0.4 + 0.3 * (st == 'proliferative'),
                    color=colors[ist], label=st,
                    zorder=10 - ist,
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.legend(title='Cell subtype:', bbox_to_anchor=(1.04, 1), loc="upper left")
        ax.set_ylim(0.1, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of EC')
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(['0.1%', '1%', '10%', '100%'])

        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_subtype_abundances.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_subtype_abundances.{:}'.format(
                        ext))

    if False:
        print('Find TFs that are expressed in nonproliferative cells')
        dsp = ds.split('is_prolif')
        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()
        fns = ds.featurenames
        tfs = [x for x in tfs if x in fns]
        dstf = {key: val.query_features_by_name(tfs) for key, val in dsp.items()}

        if dstf[False].n_samples > 300:
            dstf[False].subsample(300, inplace=True)

        exp = dstf[False].counts.mean(axis=1).to_frame()
        exp.rename(columns={0: 'expression'}, inplace=True)
        exp['Family'] = ''
        for tf in exp.index:
            row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
            exp.at[tf, 'Family'] = row['Family']
        exp.sort_values('expression', ascending=False, inplace=True)
        exp['rank'] = np.arange(len(tfs))
        exp = exp.loc[tfs]
        megs = exp.nlargest(50, 'expression')

        print('Find TFs that are differentially expressed in proliferative cells')
        comp = dstf[True].compare(dstf[False])
        comp['Family'] = ''
        for tf in comp.index:
            row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
            comp.at[tf, 'Family'] = row['Family']

        comp.sort_values(['statistic', 'log2_fold_change'], ascending=[False, False], inplace=True)
        comp['rank'] = np.arange(len(tfs))
        comp = comp.loc[tfs]

        ndegs = 50
        degs = comp.nlargest(ndegs, 'statistic')
        comp_sorted = comp.sort_values('statistic', ascending=False)

        # Count families of TFs
        fam_counts = comp['Family'].value_counts()
        fam_counts['Others'] = fam_counts.pop('Others')
        fam_counts_null = megs['Family'].value_counts()
        fam_counts_deg = degs['Family'].value_counts()
        for fam in fam_counts.index:
            if fam not in fam_counts_null:
                fam_counts_null[fam] = 0
            if fam not in fam_counts_deg:
                fam_counts_deg[fam] = 0
        fam_counts_null = fam_counts_null[fam_counts.index]
        fam_counts_deg = fam_counts_deg[fam_counts.index]
        fam_frac = 1.0 * fam_counts / fam_counts.sum()
        fam_frac_deg = 1.0 * fam_counts_deg / fam_counts_deg.sum()
        fam_frac_null = 1.0 * fam_counts_null / fam_counts_null.sum()

        if False:
            print('Pie charts of TF families')
            fig, ax = plt.subplots()
            size = 0.3
            cmap = plt.get_cmap("tab20c")
            colors = cmap(np.linspace(0, 1, len(fam_counts)))
            ax.pie(
                fam_frac.values, radius=1, colors=colors,
                wedgeprops=dict(width=size, edgecolor='w'))
            wedges, _ = ax.pie(
                fam_frac_deg.values, radius=1-size, colors=colors,
                wedgeprops=dict(width=size, edgecolor='w'))
            handles = []
            texts = []
            for iw, wedge in enumerate(wedges):
                tf = inc.index.values[iw]
                print(tf)
                if fam_counts_deg[tf] > 0:
                    handles.append(wedge)
                    texts.append('{:} ({:.1f})'.format(tf, inc[tf]))
            ax.legend(
                    handles, texts, loc='upper left', title='TF family\n(fold enrichment):',
                    bbox_to_anchor=(1.01, 1.01),
                    bbox_transform=ax.transAxes,
                    )
            fig.tight_layout(rect=(0, 0, 0.93, 1))

        print('Sankey plots of TF families')
        from scipy.interpolate import pchip_interpolate
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.5))
        axs = axs.ravel()
        data = [{
            'fam_frac1': fam_frac, 'fam_frac2': fam_frac_null,
            'fam_counts1': fam_counts, 'fam_counts2': fam_counts_null},
            {
            'fam_frac1': fam_frac, 'fam_frac2': fam_frac_deg,
            'fam_counts1': fam_counts, 'fam_counts2': fam_counts_deg},
            ]
        titles = ['Expressed in nonproliferative ECs', 'Upregulated in proliferative ECs']
        for i in range(2):
            ax = axs[i]
            datum = data[i]

            fam_frac1 = datum['fam_frac1']
            fam_counts1 = datum['fam_counts1']
            fam_frac2 = datum['fam_frac2']
            fam_counts2 = datum['fam_counts2']

            inc = fam_frac2 / fam_frac1
            fam_highlight = fam_counts2.index[(fam_counts2 > 1) & (inc > 2)]
            colors = sns.color_palette('husl', n_colors=len(fam_highlight))
            cdict = dict(zip(fam_highlight, colors))

            xbase = -0.1
            ngreys = 0
            texts = []
            handles = []
            yleft = 1.0
            yright = 1.0
            for tf in fam_frac1.index:
                fr = fam_frac1[tf]
                frdeg = fam_frac2[tf]

                yb_left = yleft - fr
                yb_right = yright - frdeg
                if tf in cdict:
                    color = cdict[tf]
                else:
                    greys = ['grey', 'lightgrey', 'darkgrey']
                    color = greys[ngreys % len(greys)]
                    ngreys += 1
                rect = plt.Rectangle(
                        (xbase, yb_left), 0.1, fr, facecolor=color, edgecolor='none',
                        )
                ax.add_artist(rect)
                rect = plt.Rectangle(
                        (1, yb_right), 0.1, frdeg, facecolor=color, edgecolor='none',
                        )
                ax.add_artist(rect)

                if tf in fam_highlight:
                    texts.append('{:} ({:.1f})'.format(tf, inc[tf]))
                    handles.append(rect)

                    # Line connecting top and bottom
                    x1, x2 = xbase + 0.2, 0.9
                    x = np.linspace(0, 1, 100)
                    y1, y2 = pchip_interpolate(
                            [0, x1, x2, 1],
                            [[yleft, yleft, yright, yright],
                             [yb_left, yb_left, yb_right, yb_right]],
                            x,
                            axis=1,
                            )
                    ax.plot(x, y1, lw=1.5, color=color)
                    ax.plot(x, y2, lw=1.5, color=color)
                    ax.fill_between(x, y1, y2, color=color, alpha=0.3)

                yleft -= fr
                yright -= frdeg

            ax.set_ylabel('Fraction of TFs')
            ax.set_xlim(xbase, 1.1)
            ax.set_xticks([])
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 0.5, 1])
            ax.legend(
                    handles, texts, loc='upper left',
                    title='TF family\n(fold enrichment):',
                    bbox_to_anchor=(1.01, 1.01),
                    bbox_transform=ax.transAxes,
                    )
            ax.set_title(titles[i])
        fig.tight_layout()

        if True:
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig(fig_fdn+'endo_TFF_proliferative.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'endo_TFF_proliferative.{:}'.format(
                        ext))

    if False:
        print('Save differentially expressed TFs between P7 and E18.5')
        dsp = ds.split(['Timepoint', 'is_prolif'])
        comp = dsp[('P7', True)].compare(dsp[('E18.5', True)])
        comp.rename(columns={'avg_self': 'avg_P7', 'avg_other': 'avg_E18.5'}, inplace=True)
        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()
        fns = ds.featurenames
        tfs = [x for x in tfs if x in fns]
        comptf = comp.loc[tfs]
        comptf['Family'] = ''
        for tf in comptf.index:
            row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
            comptf.at[tf, 'Family'] = row['Family']

        deglist1 = comptf.sort_values('statistic', ascending=False)
        deglist2 = comptf.loc[comptf['Family'].isin(['MYB', 'E2F'])].sort_values('statistic', ascending=False)
        deglist3 = deglist1.loc[deglist1['statistic'] >= 0.2].sort_values('log2_fold_change', ascending=False)
        if True:
            deglist1.to_csv(
                fig_fdn+'proliferative_DEG_P7_vs_E18.5_TFs.tsv',
                sep='\t',
                index=True,
                )
            deglist2.to_csv(
                fig_fdn+'proliferative_DEG_P7_vs_E18.5_TFs_MYBandE2F.tsv',
                sep='\t',
                index=True,
                )
            deglist3.to_csv(
                fig_fdn+'proliferative_DEG_P7_vs_E18.5_TFs_statistic_0.2+.tsv',
                sep='\t',
                index=True,
                )



    if False:
        print('Same TF family analysis, but separately for E18.5 and P7')
        dsp = ds.split(['Timepoint', 'is_prolif'])
        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()
        fns = ds.featurenames
        tfs = [x for x in tfs if x in fns]
        dstf = {key: val.query_features_by_name(tfs) for key, val in dsp.items()}

        fam_counts = {}
        fam_count_degs = {}
        fam_fracs = {}
        fam_frac_degs = {}
        incs = {}
        tps = ['E18.5', 'P7']
        for tp in tps:
            if dstf[(tp, False)].n_samples > 300:
                dstf[(tp, False)].subsample(300, inplace=True)
            comp = dstf[(tp, True)].compare(dstf[(tp, False)])
            comp['Family'] = ''
            for tf in comp.index:
                row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
                comp.at[tf, 'Family'] = row['Family']

            comp.sort_values(['statistic', 'log2_fold_change'], ascending=[False, False], inplace=True)
            comp['rank'] = np.arange(len(tfs))
            comp = comp.loc[tfs]

            ndegs = 50
            degs = comp.nlargest(ndegs, 'statistic')
            comp_sorted = comp.sort_values('statistic', ascending=False)

            # Count families of TFs
            fam_count = comp['Family'].value_counts()
            fam_count['Others'] = fam_count.pop('Others')
            fam_count_deg = degs['Family'].value_counts()
            for fam in fam_count.index:
                if fam not in fam_count_deg:
                    fam_count_deg[fam] = 0
            fam_count_deg = fam_count_deg[fam_count.index]
            fam_frac = 1.0 * fam_count / fam_count.sum()
            fam_frac_deg = 1.0 * fam_count_deg / fam_count_deg.sum()
            inc = fam_frac_deg / fam_frac

            fam_counts[tp] = fam_count
            fam_count_degs[tp] = fam_count_deg
            fam_fracs[tp] = fam_frac
            fam_frac_degs[tp] = fam_frac_deg
            incs[tp] = inc

        fam_highlight = []
        for tp in tps:
            tmp = list(fam_count_degs[tp].index[(fam_count_degs[tp] > 1) & (incs[tp] > 1)])
            fam_highlight.extend(tmp)
        fam_highlight = list(set(fam_highlight))

        from scipy.interpolate import pchip_interpolate
        fig, axs = plt.subplots(1, 2, figsize=(8, 3), sharey=True, sharex=True)
        axs = axs.ravel()
        colors = sns.color_palette('husl', n_colors=len(fam_highlight))
        cdict = dict(zip(fam_highlight, colors))

        texts = []
        handles = []
        tps = ['E18.5', 'P7']
        for i in range(2):
            ax = axs[i]
            tp = tps[i]
            fam_frac = fam_fracs[tp]
            fam_frac_deg = fam_frac_degs[tp]

            xbase = -0.1
            ngreys = 0
            yleft = 1.0
            yright = 1.0
            for tf in fam_frac.index:
                fr = fam_frac[tf]
                frdeg = fam_frac_deg[tf]

                yb_left = yleft - fr
                yb_right = yright - frdeg
                if tf in cdict:
                    color = cdict[tf]
                else:
                    greys = ['grey', 'lightgrey', 'darkgrey']
                    color = greys[ngreys % len(greys)]
                    ngreys += 1
                rect = plt.Rectangle(
                        (xbase, yb_left), 0.1, fr, facecolor=color, edgecolor='none',
                        )
                ax.add_artist(rect)
                rect = plt.Rectangle(
                        (1, yb_right), 0.1, frdeg, facecolor=color, edgecolor='none',
                        )
                ax.add_artist(rect)

                if tf in fam_highlight:
                    if i == 0:
                        texts.append('{:} ({:.1f})'.format(tf, inc[tf]))
                        handles.append(rect)

                    # Line connecting top and bottom
                    x1, x2 = xbase + 0.2, 0.9
                    x = np.linspace(0, 1, 100)
                    y1, y2 = pchip_interpolate(
                            [0, x1, x2, 1],
                            [[yleft, yleft, yright, yright],
                             [yb_left, yb_left, yb_right, yb_right]],
                            x,
                            axis=1,
                            )
                    ax.plot(x, y1, lw=1.5, color=color)
                    ax.plot(x, y2, lw=1.5, color=color)
                    ax.fill_between(x, y1, y2, color=color, alpha=0.3)

                yleft -= fr
                yright -= frdeg

            ax.set_xlim(xbase, 1.1)
            ax.set_xticks([])
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 0.5, 1])
            ax.set_title(tp)
        axs[0].set_ylabel('Fraction of TFs')
        axs[1].legend(
                handles, texts, loc='upper left', title='TF family\n(fold enrichment):',
                bbox_to_anchor=(1.01, 1.01),
                bbox_transform=ax.transAxes,
                )
        fig.tight_layout()

    if False:
        print('Differential expression between E18.5 and P7 proliferative cells')
        dsi = dsp[True]
        dsie = dsi.query_samples_by_metadata('Timepoint == "E18.5"')
        dsi7 = dsi.query_samples_by_metadata('Timepoint == "P7"')
        comp = dsi7.compare(dsie)
        comp['log2_fc'] = np.log2(dsi7.counts.mean(axis=1) + 0.1) - np.log2(dsie.counts.mean(axis=1) + 0.1)

    if False:
        print('DEGs for proliferative cells, only TFs')
        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()

        fns = ds.featurenames
        tfs = [x for x in tfs if x in fns]

        dstf = {key: val.query_features_by_name(tfs) for key, val in dsp.items()}
        dstf[False].subsample(500, inplace=True)
        comp = dstf[True].compare(dstf[False])
        comp['Family'] = ''
        for tf in comp.index:
            row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
            comp.at[tf, 'Family'] = row['Family']
        degs = comp.nsmallest(50, 'P-value')

        fam_counts = comp['Family'].value_counts()
        fam_counts_deg = degs['Family'].value_counts()
        for fam in fam_counts.index:
            if fam not in fam_counts_deg:
                fam_counts_deg[fam] = 0
        fam_counts_deg = fam_counts_deg[fam_counts.index]

        fam_frac = 1.0 * fam_counts / fam_counts.sum()
        fam_frac_deg = 1.0 * fam_counts_deg / fam_counts_deg.sum()
        inc = fam_frac_deg / fam_frac
        inc_top = inc.nlargest(10)

        print('Pie charts of TF families')
        fig, ax = plt.subplots()
        size = 0.3
        cmap = plt.get_cmap("tab20c")
        colors = cmap(np.linspace(0, 1, len(fam_counts)))
        ax.pie(
            fam_frac.values, radius=1, colors=colors,
            wedgeprops=dict(width=size, edgecolor='w'))
        wedges, _ = ax.pie(
            fam_frac_deg.values, radius=1-size, colors=colors,
            wedgeprops=dict(width=size, edgecolor='w'))
        handles = []
        texts = []
        for iw, wedge in enumerate(wedges):
            tf = inc.index.values[iw]
            print(tf)
            if fam_counts_deg[tf] > 0:
                handles.append(wedge)
                texts.append('{:} ({:.1f})'.format(tf, inc[tf]))
        ax.legend(
                handles, texts, loc='upper left', title='TF family\n(fold enrichment):',
                bbox_to_anchor=(1.01, 1.01),
                bbox_transform=ax.transAxes,
                )
        fig.tight_layout(rect=(0, 0, 0.93, 1))

        degs.to_csv('../../data/gene_lists/endo_DEGS_proliferative_TFonly.tsv', sep='\t', index=True)
        fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families.png')
        fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families.pdf')

        print('Try and get a sense of stats via looking at random TFs')
        n = len(degs)
        nrep = 100
        enrichs = pd.DataFrame([], columns=inc_top.index)
        for ir in range(nrep):
            print('rep #{:}'.format(ir+1))
            ind = np.arange(len(comp))
            np.random.shuffle(ind)
            ind = ind[:n]
            fake = comp.iloc[ind]['Family']
            fam_counts_fake = fake.value_counts()
            for fam in fam_counts.index:
                if fam not in fam_counts_fake:
                    fam_counts_fake[fam] = 0
            fam_frac_fake = 1.0 * fam_counts_fake / fam_counts_fake.sum()
            inc_fake = fam_frac_fake / fam_frac
            enrichs.loc[ir] = [inc_fake[tf] for tf in enrichs.columns]

        fig, axs = plt.subplots(2, 5, figsize=(10, 4))
        axs = axs.ravel()
        for i, ax in enumerate(axs):
            tf = enrichs.columns[i]
            data = enrichs[tf].values
            x = np.sort(data)
            y = 100 * (1.0 - np.linspace(0, 1, len(x)))
            ax.plot(x, y, lw=2, color='grey')
            ax.grid(True)
            ax.set_title(tf)
            xe = inc_top[tf]
            ax.plot([xe] * 2, [0, 100], color='red', lw=1)
        fig.text(0.52, 0.02, 'TF family enrichment', ha='center')
        fig.text(0.02, 0.52, 'Percent of simulations with enrichment > x', ha='center', va='center', rotation=90)
        fig.tight_layout(rect=(0.03, 0.05, 1, 1))
        fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families_enrichments.png')
        fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families_enrichments.pdf')

    if False:
        print('Analyze separately the E18.5 and the P7')
        ds.samplesheet['is_prolif'] = ds.samplesheet['cellSubtype'] == 'Proliferative EC'
        dsp2 = ds.split(['Timepoint', 'is_prolif'])

        print('DEGs for proliferative cells, only TFs')
        tf_table = pd.read_csv('../../data/gene_lists/Mus_musculus_TF.txt', sep='\t', index_col=2)
        tfs = tf_table['Symbol'].unique()
        fns = ds.featurenames
        tfs = [x for x in tfs if x in fns]
        dstf = {key: val.query_features_by_name(tfs) for key, val in dsp2.items()}
        degs_tp = {}
        for tp in ['E18.5', 'P7']:
            print(tp)
            if dstf[(tp, False)].n_samples > 300:
                dstf[(tp, False)].subsample(300, inplace=True)
            comp = dstf[(tp, True)].compare(dstf[(tp, False)])
            comp['Family'] = ''
            for tf in comp.index:
                row = tf_table.loc[tf_table['Symbol'] == tf].iloc[0]
                comp.at[tf, 'Family'] = row['Family']

            comp.sort_values(['statistic', 'log2_fold_change'], ascending=[False, False], inplace=True)
            comp['rank'] = np.arange(len(tfs))
            comp = comp.loc[tfs]

            degs = comp.nlargest(50, 'statistic')
            degs.to_csv('../../data/gene_lists/endo_DEGS_proliferative_TFonly_{:}.tsv'.format(tp.replace('.', '')), sep='\t', index=True)

            comp.sort_values('statistic', ascending=False).to_csv('../../data/gene_lists/endo_DEGS_proliferative_TFonly_{:}_all.tsv'.format(tp.replace('.', '')), sep='\t', index=True)

            degs_tp[tp] = comp

    if False:
        print('Pie charts of TF families')
        for tp in ['E18.5', 'P7']:
            print(tp)
            comp = degs_tp[tp]
            degs = comp.nlargest(50, 'statistic')
            fam_counts = comp['Family'].value_counts()
            fam_counts_deg = degs['Family'].value_counts()
            for fam in fam_counts.index:
                if fam not in fam_counts_deg:
                    fam_counts_deg[fam] = 0
            fam_counts_deg = fam_counts_deg[fam_counts.index]

            fam_frac = 1.0 * fam_counts / fam_counts.sum()
            fam_frac_deg = 1.0 * fam_counts_deg / fam_counts_deg.sum()
            inc = fam_frac_deg / fam_frac
            inc_top = inc.nlargest(10)

            fig, ax = plt.subplots()
            size = 0.3
            cmap = plt.get_cmap("tab20c")
            colors = cmap(np.linspace(0, 1, len(fam_counts)))
            ax.pie(
                fam_frac.values, radius=1, colors=colors,
                wedgeprops=dict(width=size, edgecolor='w'))
            wedges, _ = ax.pie(
                fam_frac_deg.values, radius=1-size, colors=colors,
                wedgeprops=dict(width=size, edgecolor='w'))
            handles = []
            texts = []
            for iw, wedge in enumerate(wedges):
                tf = inc.index.values[iw]
                print(tf)
                if fam_counts_deg[tf] > 0:
                    handles.append(wedge)
                    texts.append('{:} ({:.1f})'.format(tf, inc[tf]))
            ax.legend(
                    handles, texts, loc='upper left', title='TF family\n(fold enrichment):',
                    bbox_to_anchor=(1.01, 1.01),
                    bbox_transform=ax.transAxes,
                    )
            ax.set_title(tp)
            fig.tight_layout(rect=(0, 0, 0.93, 1))
            fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families_{:}.png'.format(tp.replace('.', '')))
            fig.savefig('../../figures/endomese_share/pies/endo_proliferative_TF_families_{:}.pdf'.format(tp.replace('.', '')))


    if False:
        fig, ax = plt.subplots(1, 1, figsize=(6, 8))
        n = 50
        yticklabels = list(degs_tp['E18.5'].nlargest(n, 'statistic').index.values)
        yticklabels2 = list(degs_tp['P7'].nlargest(n, 'statistic').index.values)
        ax.set_xlim(-0.05, 1.05)
        ax.set_yticks(np.arange(n))
        ax.set_yticklabels(yticklabels)
        ax2 = ax.twinx()
        ax2.set_yticks(np.arange(n))
        ax2.set_yticklabels(yticklabels2)
        for i, gene in enumerate(yticklabels):
            if gene not in yticklabels2:
                ax.text(0.03, i, '?', ha='center', va='center')
                continue
            i2 = yticklabels2.index(gene)
            ax.arrow(0, i, 1, i2 - i, length_includes_head=True)
        for i, gene in enumerate(yticklabels2):
            if gene not in yticklabels:
                ax.text(0.97, i, '?', ha='center', va='center')

        ax.set_ylim(n - 0.5, -0.5)
        ax2.set_ylim(n - 0.5, -0.5)
        ax.set_title('DEGs for proliferation')
        ax.set_xticks([-0.03, 1.03])
        ax.set_xticklabels(['E18.5', 'P7'])
        ax.xaxis.tick_top()

        import matplotlib.lines as mlines
        from collections import Counter
        tf_plot = np.union1d(yticklabels, yticklabels2)
        families = Counter()
        for tf in tf_plot:
            families[tf_table.loc[tf_table['Symbol'] == tf, 'Family'].iloc[0]] += 1

        labels = [x[0] for x in families.most_common() if x[0] != 'Others'] + ['Others']
        colors = sns.color_palette('husl', len(labels))
        cmap = dict(zip(labels, colors))
        handles = []
        for key, color in cmap.items():
            h = mlines.Line2D(
                [], [], color=color, marker='o', lw=0,
                markersize=5,
                )
            handles.append(h)
            labels.append(key.upper())
        ax.legend(
                handles,
                labels,
                title='TF family:',
                loc='upper left',
                fontsize=10,
                ncol=1,
                bbox_to_anchor=(1.31, 1.01),
                bbox_transform=ax.transAxes,
                )

        for i, tf in enumerate(yticklabels):
            fam = tf_table.loc[tf_table['Symbol'] == tf, 'Family'].iloc[0]
            ax.plot([-0.03] * 2, [i-0.5, i+0.5], lw=5, color=cmap[fam])

        for i, tf in enumerate(yticklabels2):
            fam = tf_table.loc[tf_table['Symbol'] == tf, 'Family'].iloc[0]
            ax.plot([1.03] * 2, [i-0.5, i+0.5], lw=5, color=cmap[fam])

        fig.tight_layout()

    plt.ion()
    plt.show()

