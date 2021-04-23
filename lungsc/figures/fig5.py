# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/08/19
content:    Fig 5 for the endo paper: Car4+ cells
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
            'Cd36',
            'Hpgd',
            'Sox7',
            'Sox17',
            'Sox18',
            'Car4',
            'Tbx2',
            'Fibin',
            'Sirpa',
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
                        fig.savefig(fig_fdn+'endo_embedding_gene_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'endo_embedding_gene_{:}.{:}'.format(
                            gene, ext))

    if True:
        print('Analyze pseudotime of Car4+')
        cst = 'Car4+ capillaries'
        ds4 = ds.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())
        vsi = vs.loc[ds4.samplenames]

        print('Find initial cell, e.g. highest cell cycle')
        stem = vsi['dimension 1'].idxmax()
        stem_idx = ds4.samplenames.tolist().index(stem)

        fn_pt = '../../data/sequencing/datasets/all_{:}/{:}_pseudotime_Car4+_nobranch.tsv'.format(version, cst)
        if not os.path.isfile(fn_pt):
            import anndata
            import scanpy

            print('Feature selection')
            features = ds4.feature_selection.overdispersed_within_groups(
                    'Mousename',
                    n_features=500,
                    inplace=False,
                    )
            dsf = ds4.query_features_by_name(features)

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
            ds4.samplesheet['pseudotime'] = adata.obs['dpt_pseudotime']
            ds4.samplesheet['pseudotime'].to_csv(fn_pt, sep='\t', index=True)
        else:
            pt = pd.read_csv(fn_pt, sep='\t', index_col=0, squeeze=True)
            pt = pt.loc[ds4.samplenames]
            ds4.samplesheet['pseudotime'] = pt

    if False:
        print('Plot embedding with pseudotime')
        fig, axs = plt.subplots(1, 2, figsize=(4, 2))
        axs = axs.ravel()
        ax, gene = axs[0], 'Timepoint'
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        ds4.plot.scatter_reduced_samples(
                vsi,
                ax=ax,
                s=40,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.set_title(gene)
        ax.grid(False)
        ax.set_axis_off()

        ax = axs[1]
        ds4.plot.scatter_reduced_samples(
                vsi,
                ax=ax,
                s=40,
                alpha=0.40,
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

        pt = ds4.samplesheet['pseudotime']
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
        fig.tight_layout(h_pad=0.9, w_pad=0.01)

        if True:
            fxf = fig_fdn+'pseudotime_Car4+'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    if False:
        print('Analyze correlates of pseudotime in Car4+ capillaries')
        cst = 'Car4+ capillaries'
        ds4 = ds.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())
        pt = pd.read_csv(fn_pt, sep='\t', index_col=0, squeeze=True)
        pt = pt.loc[ds4.samplenames]
        ds4.samplesheet['pseudotime'] = pt

        corr = ds4.correlation.correlate_features_phenotypes(
                'pseudotime',
                method='pearson',
                fillna=0,
                )

        fig, ax = plt.subplots(figsize=(3.3, 2.7))
        sns.kdeplot(corr, ax=ax, color='grey', zorder=10, shade=True, alpha=0.8, legend=False)
        nout = 8
        ax.scatter(corr.nlargest(nout).values, [0] * nout, color='darkred', s=20, lw=2, zorder=12)
        ax.scatter(corr.nsmallest(nout).values, [0] * nout, color='darkblue', s=20, lw=2, zorder=12)
        ax.text(0.97, 0.95, '\n'.join(corr.nlargest(nout).index), ha='right', va='top', fontsize=10, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='darkred', alpha=0.9))
        ax.text(0.03, 0.95, '\n'.join(corr.nsmallest(nout).index), ha='left', va='top', fontsize=10, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='darkblue', alpha=0.9))
        ax.plot([0] * 2, list(ax.get_ylim()), lw=2, ls='--', color='k', alpha=0.5, zorder=14)
        ax.grid(True)
        ax.set_xlabel('Spearman $\\rho$ with pseudotime\nin Car4+ capillaries')
        ax.set_ylabel('Density')
        ax.set_xlim(-1, 1)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'correlation_with_pseudotime_Car4+'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

        nsmooth = 5
        ds_sor = dsi.sort_by_metadata('pseudotime', axis='samples')
        ds_sor.counts.smoothen(nsmooth, inplace=True)
        pt_sor = ds_sor.samplesheet['pseudotime']


        print('Differential expression between final and initial pseudotime')
        nqua = 6
        ds_sor.samplesheet['pseudotime_quantile'] = pd.qcut(pt_sor, nqua, labels=False)

        #ds_sora = ds_sor.average('samples', 'pseudotime_quantile')
        #comp = (ds_sora.counts.iloc[:, -1] + 0.1) / (ds_sora.counts.iloc[:, 0] + 0.1)
        ds_sorp = ds_sor.split('pseudotime_quantile')
        comp = ds_sorp[nqua-1].compare(ds_sorp[0])
        comp['-log2_fold_change'] = -comp['log2_fold_change']


        corr = ds_sor.correlation.correlate_features_phenotypes(
                'pseudotime',
                method='spearman',
                fillna=0,
                )

        nout = 10
        genes_both = [
            corr.nlargest(nout).index,
            corr.nsmallest(nout).index,
        ]
        #genes_both = [
        #    comp.loc[corr > 0.3].nlargest(nout, 'log2_fold_change').index,
        #    comp.loc[corr < -0.55].nsmallest(nout, 'log2_fold_change').index,
        #]
        #genes_both = [
        #    comp.nlargest(nout).index,
        #    comp.nsmallest(nout).index,
        #]
        #genes_both = [
        #    comp.loc[comp['statistic'] > 0.8].nlargest(nout, 'log2_fold_change').index,
        #    comp.loc[comp['statistic'] > 0.8].nsmallest(nout, 'log2_fold_change').index,
        #]
        #genes_both = [
        #    comp.loc[comp['log2_fold_change'] > 0].nlargest(nout, ['statistic', 'log2_fold_change']).index,
        #    comp.loc[comp['log2_fold_change'] < 0].nlargest(nout, ['statistic', '-log2_fold_change']).index,
        #]

        from scipy import interpolate
        fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
        x = pt_sor
        yss = np.log10(ds_sor.counts + 0.1)

        colors = sns.color_palette('husl', n_colors=nout)
        titles = ['Increasing with pseudotime', 'Decreasing with pseudotime']
        for iax, (ax, genes) in enumerate(zip(axs, genes_both)):
            for i in range(nout):
                y = yss.loc[genes[i]]

                bins = np.linspace(0, 1, 7)
                yb = np.zeros(len(y), int)
                for bi in bins[:-1]:
                    yb[x >= bi] += 1
                yb -= 1
                xm = 0.5 * (bins[1:] + bins[:-1])
                ym = np.array([np.mean(y[yb == i]) for i in range(len(bins) - 1)])

                xint = np.linspace(xm[0], xm[-1], 100)
                yint = interpolate.pchip_interpolate(xm, ym, xint)
                xm, ym = xint, yint

                norm = y.values.max()
                y = (y + 1) / (norm + 1)
                ym = (ym + 1) / (norm + 1)

                ax.plot(xm, ym, color=colors[i], alpha=0.8, lw=2)
                ax.scatter(x, y, color=colors[i], alpha=0.2, label=genes[i])

            ax.set_title(titles[iax])
            ax.grid(True)
            ax.set_xlabel('Pseudotime in Car4+ capillaries')
        axs[0].set_ylabel('$\\log_{10}$ (gene expression)\n[normalized to max cell]')
        axs[0].legend(loc='lower right')
        axs[1].legend(loc='lower left')
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'normalized_expression_along_pseudotime_Car4+'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    if True:
        print('Heatmap of genes changing with pseudotime')
        genes, pathways = pd.read_excel(
            '../../data/gene_lists/Heat Map Car4+ Developmental.xlsx').values.T

        pt_bins = [
            [0, 0.225],
            [0.175, 0.425],
            [0.375, 0.625],
            [0.575, 0.825],
            [0.775, 1.01],
            ]

        data = pd.DataFrame([], index=genes)
        for (bl, br) in pt_bins:
            dsi = ds4.query_samples_by_metadata(
                '(pseudotime >= @bl) & (pseudotime < @br)',
                local_dict=locals(),
                )
            mat = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data[(bl, br)] = mat

        ## Normalize by max expression of that gene
        #data += 1
        #data = (data.T / data.max(axis=1)).T

        # Normalize 0-1 for each gene
        data = (data.T - data.min(axis=1)).T
        data = (data.T / data.max(axis=1)).T

        fig, ((tax, ax_), (ax, cax)) = plt.subplots(
                2, 2, figsize=(8.4, 3.2),
                gridspec_kw={'height_ratios': [1, 8], 'width_ratios': [20, 1]},
                )
        ax_.set_axis_off()
        sns.heatmap(
            data.T,
            ax=ax,
            cmap='plasma',
            vmin=0,
            vmax=1,
            fmt='.1f',
            xticklabels=True,
            yticklabels=True,
            cbar=True,
            cbar_ax=cax,
            cbar_kws={'ticks': [0, 0.33, 0.67, 1], 'label': 'Relative expression'},
            )
        cax.set_yticklabels(['none', 'low', 'mid', 'high'])

        pwu = np.unique(pathways)
        pwd = {}
        for pw in pwu:
            tmp = (pathways == pw).nonzero()[0]
            pwd[pw] = {'min': tmp.min(), 'max': tmp.max(), 'mean': tmp.mean()}
        cmap = dict(zip(pwu, sns.color_palette('husl', len(pwu))))
        txposs = []
        txlabels = []
        for pw, d in pwd.items():
            x0 = d['min']
            dx = d['max'] + 1 - x0
            tax.add_artist(plt.Rectangle((x0, 0), dx, 1, color=cmap[pw]))
            txposs.append(d['mean'] + 0.5)
            txlabels.append(pw.replace(' ', '\n'))
            if x0 > 0:
                ax.axvline(x0, lw=2, color='grey')
                tax.axvline(x0, lw=2, color='grey')
        idx = np.argsort(txposs)
        txlabels = np.array(txlabels)[idx]
        txposs = np.array(txposs)[idx]
        tax.set_xticks(txposs)
        tax.set_xticklabels(txlabels)
        tax.xaxis.tick_top()
        tax.set_yticks([])
        tax.set_xlim(0, len(genes))

        #ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_ylabel('Pseudotime', rotation=270, labelpad=33)
        ax.arrow(
                -0.023, 0.84, 0, -0.6, color='k', lw=1.5,
                head_width=0.02,
                head_length=0.08,
                overhang=0.2,
                clip_on=False,
                transform=ax.transAxes,
                )
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)
        ax.set_xlim(0, len(genes))
        ax.set_ylim(len(pt_bins), 0)
        fig.tight_layout(h_pad=0.5)

        if True:
            fxf = fig_fdn+'heatmap_pseudotime_Car4+'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    if False:
        print('Differential expression between P1-P7 and P7-P21')
        cst = 'Car4+ capillaries'
        dsi = ds.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())
        dsip = dsi.split('Timepoint')

        comp1 = dsip['P7'].compare(dsip['P1'])
        comp2 = dsip['P21'].compare(dsip['P7'])

        deg1 = comp1.loc[comp1['statistic'] > 0.3].sort_values('log2_fold_change', ascending=False)
        deg1.rename(columns={'avg_self': 'avg_P7', 'avg_other': 'avg_P1'}, inplace=True)
        deg2 = comp2.loc[comp2['statistic'] > 0.3].sort_values('log2_fold_change', ascending=False)
        deg2.rename(columns={'avg_self': 'avg_P21', 'avg_other': 'avg_P7'}, inplace=True)

        deg1.to_csv(fig_fdn+'endo_DEG_Car4+_P1_vs_P7.tsv', sep='\t')
        print('Output result to Excel')
        with pd.ExcelWriter(fig_fdn+'endo_DEG_Car4+_P1_P7_P21.xlsx') as ex:
            deg1.to_excel(ex, sheet_name='P7 vs P1')
            deg2.to_excel(ex, sheet_name='P21 vs P7')

    if False:
        print('Pathway analysis between the various timepoints')
        comps = [
            ['P1', 'P7'],
            ['P7', 'P21'],
            ]
        npws = 4
        xmax = 13
        for tp1, tp2 in comps:
            fns = [f'{fig_fdn}DEG_{tp2}_vs_{tp1}_up_pathways.xlsx',
                   f'{fig_fdn}DEG_{tp2}_vs_{tp1}_down_pathways.xlsx']
            cmaps = [plt.cm.get_cmap('Reds'), plt.cm.get_cmap('Blues')]
            fig, ax = plt.subplots(1, 1, figsize=(5.5, 2.5))

            # UP/DOWN
            lab_ln = 0
            xmaxi = 0
            yticks = []
            yticklabels = []
            for ii, fn in enumerate(fns):
                data = pd.read_excel(fn, sheet_name='Enrichment')
                # Ignore subcategories
                data = data.loc[data['GroupID'].str.contains('Summary')]
                data.sort_values('LogP', inplace=True)
                data = data.iloc[:npws]
                for j, (_, row) in enumerate(data.iterrows()):
                    y = j + ii * (npws + 0.5)
                    x = -row['LogP']
                    color = cmaps[ii](1.0 * x / xmax)
                    ax.barh([y], [x], color=color)
                    yticks.append(y)
                for lab in data['Description'].values:
                    lab = lab[0].upper()+lab[1:]
                    yticklabels.append(lab)
                lab_ln = max(lab_ln, max(len(l) for l in yticklabels))
            ax.grid(axis='x')
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels, fontsize=8)
            ax.plot([0, ax.get_xlim()[1]], [npws - 0.25] * 2, lw=1, color='grey')
            ax.set_ylim(2 * len(data), -0.5)
            ax.yaxis.set_ticks_position('right')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylabel(
                    f'{tp2} > {tp1}\n\n\n\n\n\n{tp2} < {tp1}',
                    ha='center', va='center', rotation=0,
                    labelpad=30)
            if ax.get_xlim()[1] > 10:
                ax.set_xticks([0, 4, 8, 12])
            elif ax.get_xlim()[1] > 6:
                ax.set_xticks([0, 2, 4, 6])
            else:
                ax.set_xticks([0, 2, 4])
            fig.set_size_inches(2.3 + 0.06 * lab_ln, 2.5)

            fig.tight_layout(h_pad=0.5)

            if True:
                fxf = fig_fdn+f'pathway_barplots_{tp1}_{tp2}'
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                    else:
                        fig.savefig('{:}.{:}'.format(fxf, ext))


    if False:
        print('RNA velocity')
        cst = 'Car4+ capillaries'
        dsi = ds.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())

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
        adata.obsm['X_umap'] = ds.samplesheet.loc[adata.obs_names, ['embed_endo_1', 'embed_endo_2']].values

        cmap = {
            'Car4+ capillaries': '#3ab646',
        }
        fig, ax = plt.subplots(figsize=(2.6, 2.4))
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
            basis='umap',
            size=0,
            ax=ax,
            alpha=0.5,
            )
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(top=-6.4)
        fig.tight_layout()

        if True:
            fxf = fig_fdn+'Car4+_velocity'
            for ext in ['svg', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    plt.ion()
    plt.show()

