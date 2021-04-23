# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 4 for the endo paper: Car4- capillaries
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_4/'


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
    vs.columns = ['dimension 1', 'dimension 2']

    if False:
        print('Plot embedding of genes')
        genes = [
            'Mest',
            'Peg3',
            'Sparcl1',
            'Depp1',
            'Rgcc',
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
        print('Plot abundances in time of microvascular')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)
        frac = frac.loc[[
            'Nonproliferative embryonic EC',
            'Early Car4- capillaries',
            'Late Car4- capillaries',
            #'Car4+ capillaries',
            ]]
        sts = list(frac.index)

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(5, 2))
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        labeld = {
                'Early Car4- capillaries': 'Early Capillaries',
                'Late Car4- capillaries': 'Late Capillaries',
                'Nonproliferative embryonic EC': 'Embryonic Capillaries',
                'Car4+ capillaries': 'Car4+ Capillaries',
            }
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
        for ist, st in enumerate(sts):
            y = frac.loc[st, xorder]

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
            out = np.vstack([outx, outy])
            if ist >= 2:
                ls = '-'
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
                    alpha=0.7,
                    color=cmap[st],
                    zorder=10 - ist,
                    ls=ls,
                    )
            ax.plot([], [], m+ls,
                    lw=2,
                    alpha=0.7,
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
                    fig.savefig(fig_fdn+'microvascular_abundances.{:}'.format(
                        ext),
                        dpi=dpi)
                else:
                    fig.savefig(fig_fdn+'microvascular_abundances.{:}'.format(
                        ext))

    if False:
        print('Analyze pseudotime of Car4-')
        cst = None # TODO TODO
        #cst = 'Car4+ capillaries'
        dsi = ds.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())
        vsi = vs.loc[dsi.samplenames]

        print('Find initial cell, e.g. highest cell cycle')
        stem = vsi['dimension 1'].idxmax()
        stem_idx = dsi.samplenames.tolist().index(stem)

        fn_pt = '../../data/sequencing/datasets/all_{:}/{:}_pseudotime_Car4+_nobranch.tsv'.format(version, cst)
        if not os.path.isfile(fn_pt):
            import anndata
            import scanpy

            print('Feature selection')
            features = dsi.feature_selection.overdispersed_within_groups(
                    'Mousename',
                    n_features=500,
                    inplace=False,
                    )
            dsf = dsi.query_features_by_name(features)

            print('PCA')
            dsc = dsf.dimensionality.pca(
                n_dims=25,
                robust=False,
                return_dataset='samples',
                )

            print('Convert to AnnData')
            adata = anndata.AnnData(
                dsc.counts.T,
                obs=dsc.samplesheet,
                var=dsc.featuresheet,
                )

            print('Perform a bunch of operations before pseudotime')
            adata.uns['iroot'] = stem_idx
            scanpy.pp.neighbors(
                    adata,
                    n_neighbors=15,
                    n_pcs=0,
                    use_rep='X',
                    knn=True,
                    random_state=0,
                    method='umap',
                    metric='correlation',
                    metric_kwds={},
                    copy=False,
                    )
            scanpy.tl.diffmap(
                    adata,
                    n_comps=15,
                    copy=False,
                    )

            print('Compute pseudotime')
            scanpy.tl.dpt(
                    adata,
                    n_dcs=10,
                    n_branchings=0,
                    min_group_size=0.01,
                    allow_kendall_tau_shift=True,
                    copy=False,
                    )
            dsi.samplesheet['pseudotime'] = adata.obs['dpt_pseudotime']
            dsi.samplesheet['pseudotime'].to_csv(fn_pt, sep='\t', index=True)
        else:
            pt = pd.read_csv(fn_pt, sep='\t', index_col=0, squeeze=True)
            pt = pt.loc[dsi.samplenames]
            dsi.samplesheet['pseudotime'] = pt

        print('Plot embedding with pseudotime')
        fig, axs = plt.subplots(1, 2, figsize=(6.2, 3))
        axs = axs.ravel()
        ax, gene = axs[0], 'Timepoint'
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        dsi.plot.scatter_reduced_samples(
                vsi,
                ax=ax,
                s=12,
                alpha=0.30,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.set_title(gene)
        ax.grid(False)
        ax.set_axis_off()

        ax = axs[1]
        dsi.plot.scatter_reduced_samples(
                vsi,
                ax=ax,
                s=12,
                alpha=0.30,
                cmap='plasma',
                color_by='pseudotime',
                color_log=False,
                )
        xstem, ystem = vsi.loc[stem].values
        ax.scatter(
            [xstem], [ystem],
            s=200,
            marker='*',
            edgecolor='k',
            facecolor='none',
            lw=2,
            )
        ax.grid(False)
        ax.set_axis_off()
        ax.set_title('Pseudotime')

        pt = dsi.samplesheet['pseudotime']
        times = np.linspace(0, 1, 10)
        traj = [{'pseudotime': 0, 'x': xstem, 'y': ystem}]
        for it in range(len(times) - 1):
            ind = (pt >= times[it])
            if it != len(times) - 1:
                ind &= (pt < times[it + 1])
            if ind.sum() == 0:
                continue
            xm, ym = vsi.loc[ind].mean(axis=0)
            tm = 0.5 * (times[it] + times[it + 1])
            traj.append({
                'pseudotime': tm,
                'x': xm,
                'y': ym,
                })
        traj = pd.DataFrame(traj)
        traj.sort_values('pseudotime', inplace=True)

        from scipy.interpolate import splprep, splev
        points = traj[['x', 'y']].values.T
        tck, u = splprep(points, s=0)
        new_points = splev(np.linspace(0, traj['pseudotime'].iloc[-2], 100), tck)
        traji = pd.DataFrame(data=new_points, index=['x', 'y']).T
        traj_colors = sns.color_palette('plasma', n_colors=len(traji) - 1)
        ax = axs[1]
        for i in range(len(traji) - 1):
            ax.plot(
                traji['x'].values[i:i+2], traji['y'].values[i:i+2],
                lw=2, color=traj_colors[i],
                )
        xi = traj['x'].values[-2]
        yi = traj['y'].values[-2]
        dx = traj['x'].values[-1] - xi
        dy = traj['y'].values[-1] - yi
        ax.arrow(
            xi, yi, dx, dy,
            color=traj_colors[-1],
            lw=2,
            length_includes_head=True,
            head_width=0.15,
            head_length=0.25,
            overhang=0.1,
            )

        fig.suptitle('{:} pseudotiming'.format(cst))
        fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=0.9, w_pad=0.01)

        if True:
            fxf = fig_fdn+'pseudotime_Car4+'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    if False:
        dsi = ds.query_samples_by_metadata(
            'cellSubtype == "Nonproliferative embryonic EC"',
            )

        genes = {
            'pro': [
                'Gas5',
                'Npm1',
            ],
            'anti': [
                'Sparcl1',
                'Rgcc',
                'Mest',
            ]}

        # correlates
        genes['pro'].extend([
            'Atp5c1',
            #'Atp5c1-ps',
            'Cd81',
            'Ifitm2',
            'S100a16',
            'Ppp1ca',
            'Gpihbp1',
            'Nrep',
            'Cavin1',
            'Gnb1',
            'Ywhaz',
            ])
        genes['anti'].extend([
            'Pecam1',
            'Col4a1',
            'Col4a2',
            'Dpysl2',
            'Tjp1',
            'Itgb1',
            ])

        genes_all = genes['pro'] + genes['anti']

        corr = dsi.correlation.correlate_features_features(
            genes_all, 'all',
            )

        print('Hierarchical clustering')
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage, leaves_list

        pdis = squareform(np.round(1.0 - corr[corr.index].values, 5))
        z = linkage(pdis, optimal_ordering=True)
        ll = leaves_list(z)[::-1]

        print('Plot heatmap')
        data = corr[corr.index].iloc[ll].T.iloc[ll].T
        for i in range(data.shape[0]):
            data.iloc[i, i] = np.nan
        vmax = np.nanmax((data.abs() - np.eye(data.shape[0])).values)
        fig, ax = plt.subplots(figsize=(5.1, 4.1))
        sns.heatmap(
                data,
                ax=ax,
                cmap='coolwarm',
                vmin=-vmax,
                vmax=vmax,
                cbar_kws={'label': 'Spearman $\\rho$ in\nembryonic miVEC'},
                )
        for tk in ax.get_yticklabels():
            tk.set_rotation(0)
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'pro_anti_angiogenic_heatmap'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    if False:
        print('Analyze RNA velocity')
        csts = [
            'Early Car4- capillaries',
            'Late Car4- capillaries',
        ]
        dsi = ds.query_samples_by_metadata(
                'cellSubtype in @csts',
                local_dict=locals(),
                inplace=False)

        print('Load combined velocity file')
        fn_combined = '../../data/sequencing/datasets/all_{:}/velocity_endo.loom'.format(version)
        import scvelo as scv
        adata_endo = scv.read(fn_combined, cache=True)
        adata_endo.var_names_make_unique()

        print('Restrict to subtypes')
        adata = adata_endo[dsi.samplenames]

        print('Follow tutorial')
        # show proportions of spliced/unspliced abundances
        #scv.utils.show_proportions(adata)

        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=25, n_neighbors=10)

        scv.tl.velocity(adata)

        scv.tl.velocity_graph(adata)

        # Use embedding from the paper
        adata.obsm['X_tsne'] = ds.samplesheet.loc[adata.obs_names, ['embed_endo_1', 'embed_endo_2']].values

        cmap = {
            'Early Car4- capillaries': '#a6bd4f',
            'Late Car4- capillaries': '#4fbd9d',
        }
        fig, ax = plt.subplots(figsize=(4.2, 4))
        ds.plot.scatter_reduced(
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
            fxf = fig_fdn+'Car4-_velocity'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    if True:
        print('Heatmap of Car4- changes in time, with pathways')
        ctms = ['Early Car4- capillaries', 'Late Car4- capillaries']

        #fn_heat = '../../data/gene_lists/endo_early_late_caps_heatmap.xlsx'
        fn_heat = '../../data/gene_lists/Heat Map Early versus Late Car4 v3.0.xlsx'
        pathways_df = pd.read_excel(
                fn_heat,
                engine='openpyxl',
                ).set_index('GeneName')['Pathway'].iloc[:-1]
        genes = pathways_df.index


        ## Try for fun a few more genes
        #add_genes = [
        #     'Hpgd',
        #     #'Rsrp1',
        #     #'Ifitm3',
        #     #'Ptprb',
        #     #'Ehd4',
        #     'Fmo1',
        #     'H2-K1',
        #     'Dennd2d',
        #     #'Itga1',
        #     'S1pr1',
        #     'Ak7',
        #     'Sertad2',
        #     #'Adgrf5',
        #     'Tmem209',
        #     #'Tspan7',
        #     #'Cd36',
        #     #'Itm2b',
        #     #'Malat1',
        #     #'B2m',
        #     ]

        #for gene in add_genes:
        #    pathways_df.loc[gene] = 'NA'
        #genes = list(genes) + add_genes
        genes = list(genes)

        # Limit to genes in the counts table
        fns = ds.featurenames
        for gene in genes:
            if gene not in fns:
                print(gene)

        from collections import defaultdict
        pathways = defaultdict(list)
        for gene, pw in pathways_df.items():
            pathways[pw].append(gene)
        pathways = list(pathways.items())

        data_geo = pd.DataFrame([], index=genes)
        data_art = data_geo.copy()
        for cst in ctms:
            dsi = ds.query_samples_by_metadata(
                'cellSubtype == @cst',
                local_dict=locals(),
                )
            data_geo[cst] = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data_art[cst] = np.log10(0.1 + dsi.counts.loc[genes].mean(axis=1))

        data = data_geo

        ## Normalize by max expression of that gene
        #data += 1
        #data = (data.T / data.max(axis=1)).T

        fig, axs = plt.subplots(
                2, 1, figsize=(11, 2.6), sharex=True,
                gridspec_kw={'height_ratios': [1, 7]})
        sns.heatmap(
                data.T,
                ax=axs[1],
                cmap='plasma',
                vmin=-1,
                vmax=data.values.max(),
                fmt='.1f',
                xticklabels=True,
                yticklabels=True,
                cbar=False,
            )
        for i, (gene, exp) in enumerate(data.iterrows()):
            xm = i + 0.5
            ym = 1
            m = (exp.iloc[-1] - exp.iloc[0])
            theta = np.arctan(m)
            xa = xm - 0.3 * np.cos(theta)
            ya = ym + 0.3 * np.sin(theta)
            dx = 0.6 * np.cos(theta)
            dy = - 0.6 * np.sin(theta)
            #axs[1].arrow(
            #        xa, ya, dx, dy, color='k',
            #        head_width=0.2, head_length=0.3, overhang=0.07,
            #        length_includes_head=True,
            #        alpha=0.5,
            #        )

        for tk in axs[1].get_yticklabels():
            tk.set_rotation(0)
        for tk in axs[1].get_xticklabels():
            tk.set_rotation(90)
            tk.set_fontsize(8)
        axs[1].set_ylim(data.shape[1], 0)
        axs[1].set_xlim(0, len(genes))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            if i != 0:
                axs[1].plot([i] * 2, [0, len(genes)], lw=2, color='lightgrey', alpha=0.9)
            i += len(gns)

        # Legend
        #labels = ['none', 'low', 'mid', 'high']
        labels = ['$0$', '$1$', '$10$', '$100$', '$1000$', '$2000$']
        sfun = lambda x: plt.cm.plasma((np.log10(x+0.1) + 1) / (data.values.max() + 1))
        handles = []
        for x in labels:
            h = axs[1].scatter([], [], marker='s', s=50, color=sfun(float(x[1:-1])))
            handles.append(h)
        leg = axs[1].legend(
                handles, labels,
                title='Expression [cpm]:',
                bbox_to_anchor=(1.01, 0.3),
                loc='center left',
                )

        axs[0].set_ylim(0, 1)
        axs[0].set_xlim(0, len(genes))
        color_d = dict(zip(
            (x[0] for x in pathways),
            sns.color_palette('muted', n_colors=len(pathways)),
            ))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            w = len(gns)
            rect = plt.Rectangle(
                    (i, 0), w, 1,
                    facecolor=color_d[pw],
                    edgecolor='none',
                    lw=0,
                    )
            axs[0].add_artist(rect)

            wt = i + 0.5 * w
            ht = 2 + 1.5 * (ipw % 2)
            axs[0].text(
                    wt, ht, pw, ha='center', va='bottom',
                    fontsize=10,
                    clip_on=False,
                    )

            if ipw % 2:
                axs[0].plot(
                        [wt] * 2, [ht - 0.2, 1.2], lw=1, color='k',
                        clip_on=False,
                        )

            i += w
        axs[0].set_axis_off()
        fig.tight_layout(h_pad=0.01)

        if True:
            fxf = fig_fdn+'early_late_heatmap_pathways'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    if False:
        print('Dot plot of Car4- changes in time, with pathways')
        ctms = ['Early Car4- capillaries', 'Late Car4- capillaries']

        #fn_heat = '../../data/gene_lists/endo_early_late_caps_heatmap.xlsx'
        fn_heat = '../../data/gene_lists/Heat Map Early versus Late Car4 v3.0.xlsx'
        pathways_df = pd.read_excel(
                fn_heat,
                engine='openpyxl',
                ).set_index('GeneName')['Pathway'].iloc[:-1]
        genes = pathways_df.index

        # Limit to genes in the counts table
        fns = ds.featurenames
        for gene in genes:
            if gene not in fns:
                print(gene)

        from collections import defaultdict
        pathways = defaultdict(list)
        for gene, pw in pathways_df.items():
            pathways[pw].append(gene)
        pathways = list(pathways.items())

        dsi = ds.query_samples_by_metadata(
                'cellSubtype in @ctms', local_dict=locals())
        dsi.query_features_by_name(np.unique(genes), inplace=True)

        fig, axs = plt.subplots(
                2, 1, figsize=(11, 2.4), sharex=True,
                gridspec_kw={'height_ratios': [1, 5]})
        dsi.plot.dot_plot(
            group_by='cellSubtype',
            group_order=ctms,
            plot_list=genes,
            ax=axs[1],
            xoffset=0.5,
            yoffset=0.5,
            vmin=0,
            vmax=1000,
            )
        axs[1].grid(False)
        #sns.heatmap(
        #        data.T,
        #        ax=axs[1],
        #        cmap='plasma',
        #        vmin=-1,
        #        vmax=data.values.max(),
        #        fmt='.1f',
        #        xticklabels=True,
        #        yticklabels=True,
        #        cbar=False,
        #    )
        #for i, (gene, exp) in enumerate(data.iterrows()):
        #    xm = i + 0.5
        #    ym = 1
        #    m = (exp.iloc[-1] - exp.iloc[0])
        #    theta = np.arctan(m)
        #    xa = xm - 0.3 * np.cos(theta)
        #    ya = ym + 0.3 * np.sin(theta)
        #    dx = 0.6 * np.cos(theta)
        #    dy = - 0.6 * np.sin(theta)
        #    axs[1].arrow(
        #            xa, ya, dx, dy, color='k',
        #            head_width=0.2, head_length=0.3, overhang=0.07,
        #            length_includes_head=True,
        #            alpha=0.5,
        #            )

        for tk in axs[1].get_yticklabels():
            tk.set_rotation(0)
        for tk in axs[1].get_xticklabels():
            tk.set_rotation(90)
            tk.set_fontsize(8)
        axs[1].set_ylim(len(ctms), 0)
        axs[1].set_xlim(0, len(genes))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            if i != 0:
                axs[1].plot([i] * 2, [0, len(genes)], lw=2, color='lightgrey', alpha=0.9)
            i += len(gns)

        # Legend
        labels = ['none', 'low', 'mid', 'high']
        labels = ['$0$', '$1$', '$10$', '$100$', '$1000$', '$2000$']
        sfun = lambda x: axs[1]._singlet_dotmap['level_color_map'](x / 1000.0)
        handles = []
        for x in labels:
            h = axs[1].scatter([], [], marker='s', s=50, color=sfun(float(x[1:-1])))
            handles.append(h)
        leg = axs[1].legend(
                handles, labels,
                title='Expression [cpm]:',
                bbox_to_anchor=(1.01, 0.3),
                loc='center left',
                )

        axs[0].set_ylim(0, 1)
        axs[0].set_xlim(0, len(genes))
        color_d = dict(zip(
            (x[0] for x in pathways),
            sns.color_palette('muted', n_colors=len(pathways)),
            ))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            w = len(gns)
            rect = plt.Rectangle(
                    (i, 0), w, 1,
                    facecolor=color_d[pw],
                    edgecolor='none',
                    lw=0,
                    )
            axs[0].add_artist(rect)

            wt = i + 0.5 * w
            ht = 2 + 1.5 * (ipw % 2)
            axs[0].text(
                    wt, ht, pw, ha='center', va='bottom',
                    fontsize=10,
                    clip_on=False,
                    )

            if ipw % 2:
                axs[0].plot(
                        [wt] * 2, [ht - 0.2, 1.2], lw=1, color='k',
                        clip_on=False,
                        )

            i += w
        axs[0].set_axis_off()
        fig.tight_layout(h_pad=0.01)

        if False:
            fxf = fig_fdn+'early_late_dotplot_pathways'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    if False:
        print('Fraction of genes increased/decreased between early and late')
        ctms = ['Early Car4- capillaries', 'Late Car4- capillaries']
        dsp = {
            key.split()[0].lower(): ds.query_samples_by_metadata('cellSubtype == @key', local_dict=locals())
            for key in ctms
            }
        reps = 1000
        comps = np.zeros((ds.n_features, reps))
        for irep in range(reps):
            i0 = np.random.randint(dsp['early'].n_samples)
            i1 = np.random.randint(dsp['late'].n_samples)
            comps[:, irep] = np.log2(dsp['late'].counts.iloc[:, i1] + 0.1) - \
                np.log2(dsp['early'].counts.iloc[:, i0] + 0.1)
        comps = pd.DataFrame(
                comps,
                index=ds.featurenames,
                )

        thresholds = np.linspace(0.2, 0.7, 30)
        ngenes = {'up': [], 'down': []}
        for th in thresholds:
            nplus = ((comps > 0).mean(axis=1) >= th).sum()
            nminus = ((comps < 0).mean(axis=1) >= th).sum()
            ngenes['up'].append(nplus)
            ngenes['down'].append(nminus)
        fig, ax = plt.subplots(figsize=(2.5, 2))
        ax.plot(thresholds, ngenes['down'], lw=2, label='down', color='dodgerblue')
        ax.plot(thresholds, ngenes['up'], lw=2, label='up', color='tomato')
        ax.set_xlabel('Minimal fraction of\n2-cell comparisons')
        ax.set_ylabel('N genes')
        ax.grid(True)
        ax.set_yscale('log')
        ax.legend(loc='upper left', bbox_to_anchor=(0.85, 1.05), bbox_transform=ax.transAxes)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'early_late_n_genes_each'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


        ngenes = {
            'up': (comps > 0).sum(axis=0),
            'down': (comps < 0).sum(axis=0),
            }
        cmap = {'up': 'tomato', 'down': 'dodgerblue'}
        fig, ax = plt.subplots(figsize=(2.3, 2))
        for key in ['down', 'up']:
            ngene = ngenes[key]
            ax.plot(np.sort(ngene), 1.0 - np.linspace(0, 1, reps), lw=2, label=key, color=cmap[key])
        ax.set_xlabel('N genes')
        ax.set_ylabel('Fraction of comps\nwith n genes > x')
        ax.grid(True)
        ax.legend(loc='upper left', bbox_to_anchor=(0.4, 1.05), bbox_transform=ax.transAxes)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'early_late_n_genes_cumulative'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

        print('Just count the expressed genes')
        n_early = (dsp['early'].counts > 0).sum(axis=0)
        n_late = (dsp['late'].counts > 0).sum(axis=0)
        cmap = {
            'early': '#a6bd4f',
            'late': '#4fbd9d',
            }
        fig, ax = plt.subplots(figsize=(2.5, 2))
        ax.plot(np.sort(n_early), 1.0 - np.linspace(0, 1, dsp['early'].n_samples), lw=2, color=cmap['early'], label='Early')
        ax.plot(np.sort(n_late), 1.0 - np.linspace(0, 1, dsp['late'].n_samples), lw=2, color=cmap['late'], label='Late')
        ax.set_xlabel('N genes expressed')
        ax.set_ylabel('Fraction of cells\nexpressing > x genes')
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_xlim(left=100)
        ax.legend(loc='lower left')
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'early_late_n_genes_expressed'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    plt.ion()
    plt.show()
