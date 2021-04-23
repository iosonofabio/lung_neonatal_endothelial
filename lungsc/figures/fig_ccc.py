# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/10/20
content:    Figure 6 for mese paper: cell communications.
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


fig_fdn = '../../figures/endomese_share/endo_paper_figure_ccc/'


ct_groups = {
        'Pericyte': ['Pericyte', 'Proliferating pericyte'],
        'Car4- caps': [
            #'Proliferative EC',
            "Nonproliferative embryonic EC",
            "Early Car4- capillaries",
            "Late Car4- capillaries",
        ],
        'Car4+ caps': [
            #'Proliferative EC',
            "Nonproliferative embryonic EC",
            "Car4+ capillaries",
        ],
        'Car4+ caps no emb': ["Car4+ capillaries"],
        'Arterial': ['Arterial EC I', 'Arterial EC II'],
        'Venous': ['Venous EC'],
        'Lymph': ['Lymphatic EC'],
        'Epi': ['AT1', 'AT2'],
        'DC': ['DC I', 'DC II', 'DC III'],
        'Mono': ['Mac V'],
        'AlvM': ['Mac I', 'Mac II', 'Mac III'],
        'IntM': ['Mac IV'],
        'Lympho': ['B cell', 'T cell', 'NK cell', 'IL cell'],
        'Neutro': ['neutrophil'],
        'Baso': ['basophil'],
        'Mast': ['mast cell'],
        'AT1': ['AT1'],
        'AT2': ['AT2'],
        'Ciliated': ['Ciliated'],
        'Clara': ['Clara'],
}


def get_interactions(organism, fns=None):
    interactions_fn = '../../data/cellphonedb/database/interaction_unpacked_mouse.tsv'
    pairs = pd.read_csv(interactions_fn, sep='\t', index_col=0)

    if organism == 'mouse':
        pairs = pairs[['gene_name_mouse_a', 'gene_name_mouse_b']]
    else:
        pairs = pairs[['gene_name_a', 'gene_name_b']]

    pairs = (pairs.iloc[:, 0] + '_' + pairs.iloc[:, 1]).value_counts().index
    pairs = pairs.sort_values()
    tmp = []
    for p in pairs:
        if '_' in p:
            tmp.append(p.split('_'))
        else:
            tmp.append([p, p])
    pairs = pd.DataFrame(tmp, index=pairs, columns=['Gene1', 'Gene2'])

    if fns is not None:
        idx = pairs.isin(fns).all(axis=1)
        pairs = pairs.loc[idx]

    return pairs


def load_epi_development(average='geometric'):
    from singlet import Dataset

    data_fn = '../../data/Cohen_et_al_2018/cohen2018_lung_epi_cpm.loom'
    dsepi = Dataset(dataset={
        'path': data_fn,
        'index_samples': 'obs_names',
        'index_features': 'var_names',
        'bit_precision': 32,
    })

    dsepi.counts.normalize(inplace=True)

    if not average:
        return dsepi

    if average == 'geometric':
        dsepi.counts.log(inplace=True)

    dsepi_av = dsepi.average('samples', by=['cellSubtype', 'Timepoint'])

    return dsepi_av


def load_epi_TMS(average='geometric'):
    from singlet import Dataset

    data_fn = '../../data/tabulamurissenis/lung_epi_droplet.h5ad'
    dsepi = Dataset(dataset={
        'path': data_fn,
        'index_samples': 'obs_names',
        'index_features': 'var_names',
        'bit_precision': 32,
    })
    dsepi.samplesheet['Timepoint'] = 'Adult'
    dsepi.samplesheet['cellSubtype'] = dsepi.samplesheet['free_annotation'].replace({
        'Alveolar Epithelial Type 2': 'AT2',
        'Alveolar Epithelial Type 1': 'AT1',

        })

    dsepi.counts.normalize(inplace=True)

    if not average:
        return dsepi

    if average == 'geometric':
        dsepi.counts.log(inplace=True)

    dsepi_av = dsepi.average('samples', by=['cellSubtype', 'Timepoint'])

    return dsepi_av


def load_epi(average='geometric'):
    import singlet

    ds_dev = load_epi_development(average=False)
    ds_adu = load_epi_TMS(average=False)
    dsepi = singlet.concatenate([ds_dev, ds_adu], missing='pad')

    if not average:
        return dsepi

    if average == 'geometric':
        dsepi.counts.log(inplace=True)

    dsepi_av = dsepi.average('samples', by=['cellSubtype', 'Timepoint'])

    return dsepi_av



if __name__ == '__main__':

    os.makedirs(fig_fdn, exist_ok=True)

    version = versions[-1]
    ds0 = DatasetLung.load(preprocess=True, version=version, include_hyperoxia=True)
    ds0.query_samples_by_metadata(
        '(doublet == 0) & (cellSubtype not in ("emt II", "", "low-quality", "Striated muscle"))',
        inplace=True)

    ds0.samplesheet['TimepointHO'] = ds0.samplesheet['Timepoint'].copy()
    ds0.samplesheet.loc[ds0.samplesheet['Treatment'] == 'hyperoxia', 'TimepointHO'] = 'P7_HO'
    

    ct_parents = {x['cellSubtype']: x['cellType'] for _, x in ds0.obs.iterrows()}
    ct_parents['AT1'] = 'epithelial'
    ct_parents['AT2'] = 'epithelial'
    ct_parents['Clara'] = 'epithelial'
    ct_parents['Ciliated'] = 'epithelial'

    ds = ds0.query_samples_by_metadata('cellType == "endothelial"')
    print('Total endothelial cells analyzed: {:}'.format(ds.n_samples))

    dsepi = load_epi(average=False)
    tps_epi = {
        'E18.5': 'E18.5',
        'P1': 'P1.25',
        'P7': 'P7',
        'P21': 'Adult',
    }
    dsp = ds0.split('TimepointHO')
    dsepip = dsepi.split('Timepoint')

    if True:
        print('Cxcl12/Cxcr4/Cxcr7')
        import singlet
        dsa = []
        for g1, csts1 in ct_groups.items():
            if g1 == 'Epi':
                continue
            dsia = ds0.query_samples_by_metadata(
                    'cellSubtype in @csts1',
                    local_dict=locals()).average('samples', by='TimepointHO')
            dsia.obs['new_name'] = [(g1, tp) for tp in dsia.samplenames]
            dsia.rename('samples', 'new_name', inplace=True)
            dsa.append(dsia)
        dsa = singlet.concatenate(dsa)

        dsepia = dsepi.average('samples', by=['cellSubtype', 'Timepoint'])

        genes = ['Cxcl12', 'Cxcr4', 'Ackr3']

        df = dsa.counts.loc[genes]
        df.columns = pd.MultiIndex.from_arrays(
                [[x[0] for x in df.columns],
                 [x[1] for x in df.columns]],
                )
        for (cst, tp)in dsepia.samplenames:
            if tp not in tps_epi:
                continue
            tmp = pd.Series(np.zeros(len(genes)), index=genes)
            for gene in genes:
                if gene in dsepia.featurenames:
                    tmp[gene] = dsepia.counts.loc[gene, (cst, tp)]
            df[(cst, tps_epi[tp])] = tmp
        df = df.T
        df.loc[('Car4+ caps', 'E18.5')] = 0

        # Select only some groups
        groups_plot = [
            #'Pericyte',
            'Car4- caps', 'Car4+ caps',
            #'VSM',
            'Arterial',
            'Venous', 'Lymph',
            #'AdvF', 'AlvF', 'MyoF/ASM',
            'DC', 'Mono',
            'AlvM', 'IntM', 'Lympho', 'Neutro', 'Baso', 'Mast',
            #'AT1', 'AT2', 'Ciliated', 'Clara',
            ]
        df = df.loc[groups_plot]

        name_conv = {'Car4+ caps no emb': 'Car4+ caps'}
        cmap = {'mesenchymal': 'forestgreen', 'endothelial': 'tomato',
                'immune': 'black', 'epithelial': 'royalblue'}
        for tp in ['E18.5', 'P1', 'P7', 'P21']:  # FIXME
            datum = df.loc[df.index.get_loc_level(tp, 1)[0]]

            fig, ax = plt.subplots(figsize=(3, 2.5))
            dat1 = datum[genes[0]].sort_values(ascending=False)
            th = 0.1 * dat1.iloc[0]
            dat1 = dat1[dat1 > th]
            norm = dat1.max()
            yticks1 = []
            yticklabels1 = []
            for ict, ((cst, _), val) in enumerate(dat1.items()):
                ct = ct_parents[ct_groups[cst][0]]
                y = 0.5 * len(dat1) - ict
                x0 = 0
                x = val / norm
                ax.barh(y, x, left=x0, height=0.8, color=cmap[ct])
                yticks1.append(y)
                yticklabels1.append(name_conv.get(cst, cst))
            ax.set_yticks(yticks1)
            ax.set_yticklabels(yticklabels1)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.plot([0, 0], [yticks1[0] + 0.5, yticks1[-1] - 0.5], lw=1, color='k')

            ax2 = ax.twinx()
            yticks2 = []
            yticklabels2 = []

            dat2 = datum[genes[1]].sort_values(ascending=False)
            th = 0.1 * dat2.iloc[0]
            dat2 = dat2[dat2 > th]
            dat2 = dat2.iloc[:5]
            norm = dat2.max()
            for ict, ((cst, _), val) in enumerate(dat2.items()):
                ct = ct_parents[ct_groups[cst][0]]
                y = 6.5 - ict
                x0 = 2
                x = x0 - val / norm
                ax2.barh(y, x0, left=x, height=0.8, color=cmap[ct])
                yticks2.append(y)
                yticklabels2.append(name_conv.get(cst, cst))
            ax2.plot([2] * 2, [yticks2[0] + 0.5, yticks2[-1] - 0.5], lw=1, color='k')

            dat2 = datum[genes[2]].sort_values(ascending=False)
            th = 0.1 * dat2.iloc[0]
            dat2 = dat2[dat2 > th]
            dat2 = dat2.iloc[:5]
            norm = dat2.max()
            for ict, ((cst, _), val) in enumerate(dat2.items()):
                ct = ct_parents[ct_groups[cst][0]]
                y = -2 - ict
                x0 = 2
                x = x0 - val / norm
                ax2.barh(y, x0, left=x, height=0.8, color=cmap[ct])
                yticks2.append(y)
                yticklabels2.append(name_conv.get(cst, cst))
            ax2.plot([2] * 2, [-1.5, yticks2[-1] - 0.5], lw=1, color='k')

            ax2.set_yticks(yticks2)
            ax2.set_yticklabels(yticklabels2)
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            ax2.set_xticks([])
            ax.set_xlim(0, 2)
            ax.set_ylim(-6, 8)
            ax2.set_ylim(-6, 8)

            ax.text(0, 0.5 * len(dat1) + 1.5, genes[0], fontsize=14, ha='center', clip_on=False)
            ax.text(2, 8, genes[1], fontsize=14, ha='center', clip_on=False)
            ax.text(2, -0.5, genes[2], fontsize=14, ha='center', clip_on=False)
            fig.text(0.025, 0.95, tp, ha='left', va='top',
                     transform=fig.transFigure, fontsize=14)
            fig.tight_layout()

            if True:
                gene1, gene2, gene3 = genes
                fnf = fig_fdn+f'bars_{gene1}_{gene2}_{gene3}_{tp}'
                for ext in ['svg', 'pdf', ['png', 300]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig('{:}.{:}'.format(fnf, ext), dpi=dpi)
                    else:
                        fig.savefig('{:}.{:}'.format(fnf, ext))

        print('Average expression increasing over time, Cxcl12')
        pairs = [
            ('Cxcl12', 'Arterial'),
            ('Ackr3', 'Venous'),
            ('Cxcr4', 'Baso'),
            ('Cxcr4', 'Neutro'),
            ('Cxcr4', 'Mono'),
            ('Cxcr4', 'AlvM'),
            ]
        genes, cts = list(zip(*pairs))
        cmap = {
            'Arterial': 'tomato',
            'Venous': 'royalblue',
            'Baso': 'grey',
            'Neutro': 'black',
            'Mono': 'turquoise',
            'AlvM': 'greenyellow',
            }
        mmap = {'Cxcl12': 'o', 'Cxcr4': 's', 'Ackr3': '^'}

        tps = ['E18.5', 'P1', 'P7', 'P21']
        from scipy.interpolate import pchip_interpolate
        fig, axs = plt.subplots(
            1, 2, figsize=(4.5, 2.5), sharey=True,
            gridspec_kw={'width_ratios': [4, 1]},
            )
        ax = axs[0]
        x = np.arange(len(tps))
        xint = np.linspace(x[0], x[-1], 100)
        for ig, (gene, ct) in enumerate(zip(genes, cts)):
            idx = [(ct, tp) for tp in tps]
            y = dsa.counts.loc[gene][idx]
            yint = pchip_interpolate(x, y, xint)
            if (ig == 0) or (genes[ig - 1] != gene):
                label = f'{gene} in {ct}'
            else:
                label = f'    "     in {ct}'
            ax.scatter(
                x, y, color=cmap[ct], label=label, zorder=10,
                marker=mmap[gene],
                alpha=0.8,
                )
            ax.plot(xint, yint, lw=2, color=cmap[ct], zorder=10)
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(tps, rotation=90)
        ax.set_yscale('log')
        ax.set_ylabel('Expr [cpm]')
        ax = axs[1]
        x = [0, 1]
        for ig, (gene, ct) in enumerate(zip(genes, cts)):
            idx = [(ct, tp) for tp in ['P7', 'P7_HO']]
            y = dsa.counts.loc[gene][idx]
            if (ig == 0) or (genes[ig - 1] != gene):
                label = f'{gene} in {ct}'
            else:
                label = f'    "     in {ct}'
            ax.plot(
                x, y, color=cmap[ct], label=label, zorder=10,
                marker=mmap[gene],
                alpha=0.8,
                ls='--',
                )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(['normoxia', 'hyperoxia'], rotation=90)
        ax.set_xlim(-0.3, 1.3)

        ax.legend(
            loc='upper left', bbox_to_anchor=(1, 1.1),
            bbox_transform=ax.transAxes,
            title='Gene and cell type:',
            )
        ax.set_ylim(3, 1e4)
        fig.tight_layout(h_pad=0.1)

        if True:
            gene1, gene2, gene3 = genes[:3]
            fnf = fig_fdn+'avg_over_time_Cxcl12_pathway'
            for ext in ['svg', 'pdf', ['png', 300]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fnf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fnf, ext))


        print('Average expression increasing over time, Apln')
        pairs_all = [[
            ('Apln', 'Car4+ caps'),
            ('Aplnr', 'Car4- caps'),
            ], [
            ('Kitl', 'Car4+ caps'),
            ('Kit', 'Car4- caps'),
            ]]
        for pairs in pairs_all:
            genes, cts = list(zip(*pairs))
            cmap = {
                'Car4- caps': '#4fbd9d',
                'Car4+ caps': '#3ab646',
                }
            mmap = {'Apln': 'o', 'Aplnr': 's', 'Kitl': 'o', 'Kit': 's'}

            tps = ['E18.5', 'P1', 'P7', 'P21']
            from scipy.interpolate import pchip_interpolate
            fig, axs = plt.subplots(
                    1, 2, figsize=(4.5, 2.1), sharey=True,
                    gridspec_kw={'width_ratios': [4, 1]},
                )
            ax = axs[0]
            x = np.arange(len(tps))
            xint = np.linspace(x[0], x[-1], 100)
            for ig, (gene, ct) in enumerate(zip(genes, cts)):
                idx = [(ct, tp) for tp in tps]
                y = dsa.counts.loc[gene][idx]
                yint = pchip_interpolate(x, y, xint)
                if (ig == 0) or (genes[ig - 1] != gene):
                    label = f'{gene} in {ct}'
                else:
                    label = f'    "     in {ct}'
                ax.scatter(
                    x, y, color=cmap[ct], label=label, zorder=10,
                    marker=mmap[gene],
                    alpha=0.8,
                    )
                ax.plot(xint, yint, lw=2, color=cmap[ct], zorder=10)
            ax.grid(True)
            ax.set_xticks(x)
            ax.set_xticklabels(tps, rotation=90)
            ax.set_yscale('log')
            ax.set_ylabel('Expr [cpm]')
            ax = axs[1]
            x = [0, 1]
            for ig, (gene, ct) in enumerate(zip(genes, cts)):
                idx = [(ct, tp) for tp in ['P7', 'P7_HO']]
                y = dsa.counts.loc[gene][idx]
                if (ig == 0) or (genes[ig - 1] != gene):
                    label = f'{gene} in {ct}'
                else:
                    label = f'    "     in {ct}'
                ax.plot(
                    x, y, color=cmap[ct], label=label, zorder=10,
                    marker=mmap[gene],
                    alpha=0.8,
                    ls='--',
                    )
            ax.grid(True)
            ax.set_xticks(x)
            ax.set_xticklabels(['normoxia', 'hyperoxia'], rotation=90)
            ax.set_xlim(-0.3, 1.3)

            ax.legend(
                loc='upper left', bbox_to_anchor=(1, 1.1),
                bbox_transform=ax.transAxes,
                title='Gene and cell type:',
                )
            ax.set_ylim(80, 1.5e3)
            fig.tight_layout(h_pad=0.1)

            if True:
                gene1, gene2 = genes[:2]
                fnf = fig_fdn+f'avg_over_time_{gene1}_pathway'
                for ext in ['svg', 'pdf', ['png', 300]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig('{:}.{:}'.format(fnf, ext), dpi=dpi)
                    else:
                        fig.savefig('{:}.{:}'.format(fnf, ext))



    if False:
        print('Plot interaction marginals, coarse')
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
        fig, axs = plt.subplots(
                1, 2, figsize=(5, 3), gridspec_kw={'width_ratios': [4, 1]},
                sharey=True,
            )
        ax = axs[0]
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

        ax = axs[1]
        for i, cst in enumerate(csts):
            x = np.array([0, 1])
            y = nints_marg[i][[2, 4]]
            ax.plot(x, y, color=cmap[cst], lw=2, label=cst, marker='o',
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(['normoxia', 'hyperoxia'], rotation=90)
        ax.set_xlim(-0.3, 1.3)
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
            print('Focus on the effect of hyperoxia, coarse')
            tps = ['P7', 'P7_HO']
            idx = [2, 5]
            ddi = (nints_sym[5] - nints_sym[2])
            ddif = (nints_sym[5] - nints_sym[2]) / (0.5 * (nints_sym[5] + nints_sym[2]))
            fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.5))
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

            vals = [0, 50, 100, 175]
            labels = [str(x) for x in vals]
            handles = [ax.scatter([], [], marker='s', color=color_fun(x)) for x in vals]
            ax.legend(
                    handles, labels, loc='upper left', title='# interactions\nin normoxia',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )

            vals = [40, 20, 10, -10, -20, -40]
            labels = [str(x) for x in vals]
            handles = [
                ax.scatter([], [], s=size_fun(x),
                           marker='^' if x > 0 else 'v',
                           c='grey') for x in vals]
            ax2 = ax.twinx()
            ax2.set_axis_off()
            ax2.legend(
                    handles, labels, loc='upper left', title='Difference in\n# interactions\nin hyperoxia',
                    bbox_to_anchor=(1.01, 0.31), bbox_transform=ax2.transAxes,
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
        print('Cumulative plots of interactions in pre/postnatal lymphatics')
        ct = 'Arterial'
        genesd = {
            'Lymph': [
                'Cd36',
                'Osmr',
                'Lefty1',
                #'Pf4',
                'Pdgfd',
                #'Rarres2',
                #'Col4a6',
                #'Mst1',
            ],
            'Venous': [
                #'Ackr1',
                #'Itga9',
                #'Pthlh',
                #'Ager',
                #'Fgf1',
                'Thbs1',
                'Ntf3',
                'Col3a1',
                'Vgf',
                'Sirpa',
            ],
            'Arterial': [
                'Pthlh',
                'Cxcl16',
                'Inhbb',
                'Ddr1',
                'Col3a1',
            ],
        }
        genes = genesd[ct]
        dspp = ds.split(['cellSubtype', 'TimepointHO'])
        cmap = {'E18.5': 'slateblue', 'P1': 'tomato'}
        lss = {'E18.5': '-', 'P1': '--'}
        fig, axs = plt.subplots(1, len(genes), figsize=(6, 1.8), sharey=True)
        axs = axs.ravel()
        for i, (gene, ax) in enumerate(zip(genes, axs)):
            for cond in ['E18.5', 'P1']:
                cstl = ct_groups[ct]
                dspi = dsp[cond].query_samples_by_metadata(
                    'cellSubtype in @cstl', local_dict=locals(),
                )
                x = dspi.counts.loc[gene].values.copy()
                x.sort()
                y = 1.0 - np.linspace(0, 1, len(x))
                x = np.insert(x, 0, 0)
                y = np.insert(y, 0, 1.0)
                x += 0.1
                label = f'{cond}'
                ax.plot(x, y, lw=2, color=cmap[cond], ls=lss[cond], label=label)
            if i == len(axs) - 1:
                ax.legend(
                    title='Condition:',
                    loc='upper left',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )
            ax.grid(True)
            ax.set_xlim(0.05, 10**4 + 1)
            ax.set_yticks([0, 0.5, 1])
            ax.set_xscale('log')
            ax.set_xticks([1, 100])
            ax.set_xticklabels(['$1$', '$100$'])
            ax.set_title(gene)
        fig.text(0.42, 0.02, 'Gene expression in lymphatic EC [cpm]', ha='center')
        fig.text(0.02, 0.52, 'Fraction of cells > x', va='center', rotation=90)
        fig.tight_layout(h_pad=0.1, rect=(0.03, 0.03, 1, 1))

        if False:
            fnd = {'Lymph': 'lymphatics', 'Venous': 'veins', 'Arterial': 'arteries'}
            fndi = fnd[ct]
            fxf = fig_fdn+f'cumulaives_pre_postbirth_{fndi}'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))

    if False:
        print('Cumulative plots of ligands and receptors in lymphatics')
        ds7 = ds0.query_samples_by_metadata('Timepoint == "P7"')
        ds7p = ds7.split(['cellSubtype', 'Treatment'])
        genes = [
            'Wnt2',
            'Wnt5a',
            #'Fzd3',
            'Fzd4',
            'Ryk',
            #'Lrp5',
            #'Lrp6',
            #'Notch1',
            #'Notch2',
            #'Jag1',
            #'Jag2',
            #'Dll4',
        ]
        ncells = {key: dsi.n_samples for key, dsi in ds7p.items()}
        cmap = {'normal': 'slateblue', 'hyperoxia': 'tomato'}
        lss = {'normal': '-', 'hyperoxia': '--'}
        condl = {'normal': 'normoxia', 'hyperoxia': 'hyperoxia'}
        fig, axs = plt.subplots(1, len(genes), figsize=(6.5, 1.8), sharey=True)
        axs = axs.ravel()
        for i, (gene, ax) in enumerate(zip(genes, axs)):
            for cond in ['normal', 'hyperoxia']:
                x = ds7p[('Lymphatic EC', cond)].counts.loc[gene].values.copy()
                x.sort()
                y = 1.0 - np.linspace(0, 1, len(x))
                x = np.insert(x, 0, 0)
                y = np.insert(y, 0, 1.0)
                x += 0.1
                nc = ncells[('Lymphatic EC', cond)]
                co = condl[cond]
                label = f'{co}'
                ax.plot(x, y, lw=2, color=cmap[cond], ls=lss[cond], label=label)
            if i == len(axs) - 1:
                ax.legend(
                    title='Condition:',
                    loc='upper left',
                    bbox_to_anchor=(1.01, 1.01), bbox_transform=ax.transAxes,
                    )
            ax.grid(True)
            ax.set_xlim(0.05, 10**4 + 1)
            ax.set_yticks([0, 0.5, 1])
            ax.set_xscale('log')
            ax.set_xticks([1, 100])
            ax.set_xticklabels(['$1$', '$100$'])
            ax.set_title(gene)
        fig.text(0.42, 0.02, 'Gene expression in lymphatic EC [cpm]', ha='center')
        fig.text(0.02, 0.52, 'Fraction of cells > x', va='center', rotation=90)
        fig.tight_layout(h_pad=0.1, rect=(0.03, 0.03, 1, 1))

        if True:
            fxf = fig_fdn+'cumulaives_Wnt_lymphatics'
            for ext in ['svg', 'pdf', ['png', 600]]:
                if isinstance(ext, list):
                    ext, dpi = ext
                    fig.savefig('{:}.{:}'.format(fxf, ext), dpi=dpi)
                else:
                    fig.savefig('{:}.{:}'.format(fxf, ext))


    # TODO: run CellPhoneDB with EndoCoarse + immune to make new heatmaps
    # like current panel B (for supplmentary probably)
