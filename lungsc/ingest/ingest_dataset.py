# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/03/19
content:    Load dataset using singlet
'''
import os
import sys
import glob
import gzip
from collections import defaultdict
import numpy as np
import pandas as pd

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, CountsTable, FeatureSheet, CountsTableSparse, SampleSheet


versions = ('20190325', '20190513', '20190620', '20190828')


# Some genes we just want them no matter what
safe_genes = [
        # Basophils
        'Mcpt4',
        'Mcpt2',
        'Mcpt9',
        'Cma2',
        'Tpsb2',
        'Tpsab1',
        'Cma1',
        'Adamts9',
    ]


def make_samplesheet(cellnames):
    df = pd.DataFrame(data=[], index=cellnames)
    df['Timepoint'] = ''
    df['Mousename'] = ''
    df['Gender'] = ''
    df['SortType'] = ''
    df['Time [days]'] = 99
    df['Treatment'] = ''
    df['Well'] = ''
    for cn in cellnames:
        # SortType
        if ('Mesench' in cn) or ('mesench' in cn):
            df.at[cn, 'SortType'] = 'mesenchymal'
        elif 'CD45' in cn:
            df.at[cn, 'SortType'] = 'immune'
        elif 'CD31' in cn:
            df.at[cn, 'SortType'] = 'endothelial'
        else:
            raise ValueError('SortType not found: {:}'.format(cn))

        # Gender
        if ('Cornfield_F' in cn) or ('_F_' in cn):
            df.at[cn, 'Gender'] = 'F'
        elif ('Cornfield_M' in cn) or ('_M_' in cn):
            df.at[cn, 'Gender'] = 'M'
        else:
            raise ValueError('Gender not found: {:}'.format(cn))

        # Timepoint
        if ('_FE18' in cn) or ('_ME18' in cn) or ('E18_M_' in cn) or ('E18_F_' in cn):
            df.at[cn, 'Timepoint'] = 'E18.5'
            df.at[cn, 'Time [days]'] = -1
        elif ('_FP1' in cn) or ('_MP1' in cn) or ('P1_M_' in cn) or ('P1_F_' in cn):
            df.at[cn, 'Timepoint'] = 'P1'
            df.at[cn, 'Time [days]'] = 1
        elif ('_FP7' in cn) or ('_MP7' in cn) or ('hyperox' in cn) or ('P7_M_' in cn) or ('P7_F_' in cn):
            df.at[cn, 'Timepoint'] = 'P7'
            df.at[cn, 'Time [days]'] = 7
        elif ('_FP21' in cn) or ('_MP21' in cn) or ('P21_M_' in cn) or ('P21_F_' in cn):
            df.at[cn, 'Timepoint'] = 'P21'
            df.at[cn, 'Time [days]'] = 21
        else:
            raise ValueError('Timepoint not found: {:}'.format(cn))

        # Treatment
        if 'hyperox' in cn:
            df.at[cn, 'Treatment'] = 'hyperoxia'
        else:
            df.at[cn, 'Treatment'] = 'normal'

        # Mousename
        df.at[cn, 'Mousename'] = '{:}_{:}'.format(
                df.at[cn, 'Gender'],
                df.at[cn, 'Timepoint'])

        # Well
        df.at[cn, 'Well'] = cn.split('_')[-1]

    return df


if __name__ == '__main__':

    version = versions[-1]

    print('Dataset version: {:}'.format(version))

    # BLOCK 1: make raw loom file
    if False:
        print('Make and save samplesheet')
        fn_co = '../../data/sequencing/datasets/all_{:}/counts.tsv'.format(version)
        with open(fn_co, 'rt') as f:
            cellnames = f.readline().rstrip('\n').split('\t')[1:]
        samplesheet = make_samplesheet(cellnames)
        fn_ss = '../../data/sequencing/datasets/all_{:}/samplesheet.tsv'.format(version)
        samplesheet.to_csv(fn_ss, sep='\t', index=True)

        print('Load counts')
        fn_co = '../../data/sequencing/datasets/all_{:}/counts.tsv'.format(version)
        with open(fn_co, 'rt') as f:
            fields = f.readline().rstrip('\n').split('\t')
            dtypes = {key: np.float32 for key in fields[1:]}
            dtypes[fields[0]] = str
        counts_table = CountsTable(pd.read_csv(
                fn_co, sep='\t', index_col=0, dtype=dtypes
                ))

        print('Set spikeins and other features')
        counts_table._spikeins = tuple(counts_table.index[counts_table.index.str.startswith('ERCC-')].tolist())
        counts_table._otherfeatures = tuple(counts_table.index[counts_table.index.str.startswith('__')].tolist())

        print('Load featuresheet')
        fn_fs = '../../data/sequencing/datasets/all_{:}/featuresheet.tsv'.format(version)
        featuresheet = FeatureSheet(pd.read_csv(fn_fs, sep='\t', index_col=0))

        print('Calculate coverage')
        cov_uniquely_mapped = counts_table.iloc[len(counts_table._spikeins): -len(counts_table._otherfeatures)].sum(axis=0)
        samplesheet['coverage'] = cov_uniquely_mapped

        print('Separate out ERCCs')
        ercc = counts_table.loc[list(counts_table._spikeins)]

        print('Rename trashing the gene: prefix')
        gns = []
        for gn in counts_table.index:
            if gn.startswith('gene:'):
                gn = gn[5:]
            gns.append(gn)
        gns = pd.Index(gns, name=counts_table.index.name)
        counts_table.index = gns
        featuresheet.index = gns

        print('Calculate number of genes')
        n_genes = (counts_table.iloc[len(counts_table._spikeins): -len(counts_table._otherfeatures)] >= 1).sum(axis=0)
        samplesheet['number_of_genes_1plusreads'] = n_genes

        print('Exclude slashes for loom')
        featuresheet.rename(columns={'Chromosome/scaffold name': 'Chromosome'}, inplace=True)

        print('Add index name for loom')
        samplesheet.index.name = 'CellID'
        counts_table.columns.name = 'CellID'

        ds = Dataset(
                counts_table=counts_table,
                featuresheet=featuresheet,
                samplesheet=samplesheet,
                )

        print('Save raw dataset as loom file')
        fn_raw = '../../data/sequencing/datasets/all_{:}/raw.loom'.format(version)
        ds.to_dataset_file(fn_raw, fmt='loom')

    # BLOCK 2: filter genes and cells
    if True:
        print('Load raw loom file')
        import loompy
        fn_raw = '../../data/sequencing/datasets/all_{:}/raw.loom'.format(version)
        with loompy.connect(fn_raw) as dsl:

            print('Load samplesheet')
            samplesheet = pd.DataFrame([], index=dsl.ca['CellID'])
            for col in dsl.ca.keys():
                if col != 'CellID':
                    samplesheet[col] = dsl.ca[col]
            N = len(samplesheet)

            print('Set features')
            features = pd.Index(dsl.ra['GeneName'])
            L = len(features)

            print('Load featuresheet')
            fn_fs = '../../data/sequencing/datasets/all_{:}/featuresheet.tsv'.format(version)
            featuresheet = pd.read_csv(fn_fs, sep='\t', index_col=0)

            if len(featuresheet.index) != L:
                raise ValueError('Featuresheet has the wrong features')
            if not (featuresheet.index == features).all():
                raise ValueError('Featuresheet has the wrong features')

            print('Total cells before filtering: {:}'.format(N))
            print('Total features before filtering: {:}'.format(L))

            print('Set spikeins and other features')
            spikeins = tuple(features[features.str.startswith('ERCC-')].tolist())
            otherfeatures = tuple(features[features.str.startswith('__')].tolist())
            ind_mapped = np.ones(L, bool)
            ind_mapped[features.isin(spikeins)] = False
            ind_mapped[features.isin(otherfeatures)] = False

            #print('Set coverage and n genes')
            #coverage = np.zeros(N, np.float32)
            #ngenes1p = np.zeros(N, np.float32)
            #gsize = 400
            #ngroups = (N // gsize) + bool(N % gsize)
            #submat = np.zeros((ind_mapped.sum(), gsize), np.float32)
            #for ig in range(ngroups):
            #    ist, ien = ig * gsize, (ig + 1) * gsize
            #    if N < ien:
            #        ien = N
            #    submat = submat[:, :(ien - ist)]
            #    submat[:, :] = dsl[ind_mapped, ist: ien]
            #    print(
            #        'Group n: {:4d}, cells {:} to {:}'.format(ig, ist, ien),
            #        end='\r')
            #    coverage[ist: ien] = submat.sum(axis=0)
            #    ngenes1p[ist: ien] = (submat >= 1).sum(axis=0)
            #print()
            #samplesheet['coverage'] = coverage
            #samplesheet['number_of_genes_1plusreads'] = ngenes1p
            #dsl.ca['coverage'] = coverage
            #dsl.ca['number_of_genes_1plusreads'] = ngenes1p

            print('Select decently expressing cells')
            good_cells = samplesheet['coverage'] >= 50000
            good_cells &= samplesheet['number_of_genes_1plusreads'] >= 400
            good_cells = good_cells.index[good_cells]
            print('Number of decent cells: {:}'.format(len(good_cells)))
            with open('../../data/sequencing/datasets/all_{:}/good_cells.tsv'.format(version), 'wt') as f:
                f.write('\n'.join(good_cells))

            print('Select decently expressed genes')
            with open('../../data/sequencing/datasets/all_{:}/good_cells.tsv'.format(version), 'rt') as f:
                good_cells = f.read().split('\n')
            good_cells_ind = samplesheet.index.isin(good_cells).nonzero()[0]
            ncellsexp = np.zeros(L, np.int64)
            gsize = 400
            ngroups = (N // gsize) + bool(N % gsize)
            for ig in range(ngroups):
                ist, ien = ig * gsize, (ig + 1) * gsize
                if N < ien:
                    ien = N
                print(
                    'Group n: {:4d}, cells {:} to {:}'.format(ig, ist, ien),
                    end='\r')
                submat = dsl[:, ist: ien]
                indi = good_cells_ind[(good_cells_ind >= ist) & (good_cells_ind < ien)] - ist
                submat = submat[:, indi]
                # NOTE: this is expression before normalization (raw reads)
                ncellsexp += (submat >= 5).sum(axis=1)
            print()
            good_genes = pd.Series(ind_mapped & (ncellsexp >= 10), index=features)
            good_genes[featuresheet['GeneName'].isin(safe_genes)] = True
            good_genes = features[good_genes]
            print('Number of decent genes: {:}'.format(len(good_genes)))
            with open('../../data/sequencing/datasets/all_{:}/good_genes.tsv'.format(version), 'wt') as f:
                f.write('\n'.join(good_genes))

            print('Restrict and normalize data')
            with open('../../data/sequencing/datasets/all_{:}/good_genes.tsv'.format(version), 'rt') as f:
                good_genes = f.read().split('\n')
            with open('../../data/sequencing/datasets/all_{:}/good_cells.tsv'.format(version), 'rt') as f:
                good_cells = f.read().split('\n')
            L2 = len(good_genes)
            N2 = len(good_cells)
            good_genes_ind = featuresheet.index.isin(good_genes).nonzero()[0]
            good_cells_ind = samplesheet.index.isin(good_cells).nonzero()[0]
            coverage = samplesheet['coverage'].values
            gsize = 400
            ngroups = (N // gsize) + bool(N % gsize)
            matrix = np.zeros((L2, N2), np.float32)
            ii = 0
            for ig in range(ngroups):
                ist, ien = ig * gsize, (ig + 1) * gsize
                if N < ien:
                    ien = N
                print(
                    'Group n: {:4d}, cells {:} to {:}'.format(ig, ist, ien),
                    end='\r')
                submat = dsl[:, ist: ien]
                subcoverage = coverage[ist: ien]
                indi = good_cells_ind[(good_cells_ind >= ist) & (good_cells_ind < ien)] - ist
                ni = len(indi)
                if ni == 0:
                    continue
                submat = submat[good_genes_ind]
                matrix[:, ii: ii + ni] = 1e6 * submat[:, indi] / subcoverage[indi]
                ii += ni
            print()

        print('Filter genes')
        featuresheet_good = featuresheet.loc[good_genes]
        print('Exclude genes without a gene name')
        ind_discard = featuresheet_good['GeneName'].astype(str) == 'nan'
        print('For genes with more than one EnsemblID, sum over them')
        gn_counts = featuresheet_good['GeneName'].value_counts()
        genes_nonunique = gn_counts.index[gn_counts > 1]
        for gn in genes_nonunique:
            inds = (featuresheet_good['GeneName'] == gn).values.nonzero()[0]
            if gn == 'nan':
                ind_discard[inds] = True
            else:
                for i in inds[1:]:
                    matrix[inds[0]] += matrix[i]
                ind_discard[inds[1:]] = True
        matrix = matrix[~ind_discard]

        print('Finalize list of genes and cells')
        featuresheet_good = featuresheet_good.loc[~ind_discard]
        samplesheet_good = samplesheet.loc[good_cells]

        print('Save filtered dataset as loom file')
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        col_attrs = {'_index': samplesheet_good.index.values}
        for col in samplesheet_good.columns:
            col_attrs[col] = samplesheet_good[col].values
        row_attrs = {}
        for col in featuresheet_good.columns:
            if '/' in col:
                colnew = col[:col.find('/')]
            else:
                colnew = col
            row_attrs[colnew] = featuresheet_good[col].values
        loompy.create(
            fn_good,
            layers={'': matrix},
            row_attrs=row_attrs,
            col_attrs=col_attrs,
            )

    # BLOCK 3: check that it all worked
    if False:
        print('Load filtered loom file')
        import loompy
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        with loompy.connect(fn_good) as dsl:
            print('Load samplesheet')
            samplesheet = pd.DataFrame([], index=dsl.ca['CellID'])
            for col in dsl.ca.keys():
                if col != 'CellID':
                    samplesheet[col] = dsl.ca[col]

            print('Load featuresheet')
            featuresheet = pd.DataFrame([], index=dsl.ra['GeneName'])
            for col in dsl.ra.keys():
                featuresheet[col] = dsl.ra[col]

            print('Load counts')
            counts_table = pd.DataFrame(
                    dsl[:, :],
                    index=featuresheet.index,
                    columns=samplesheet.index,
                    )

        print('Convert to singlet data structures')
        samplesheet = SampleSheet(samplesheet)
        featuresheet = FeatureSheet(featuresheet)
        counts_table = CountsTable(counts_table)

        print('Set spikeins and other features')
        counts_table._spikeins = tuple(counts_table.index[counts_table.index.str.startswith('ERCC-')].tolist())
        counts_table._otherfeatures = tuple(counts_table.index[counts_table.index.str.startswith('__')].tolist())

        ds = Dataset(
                counts_table=counts_table,
                featuresheet=featuresheet,
                samplesheet=samplesheet,
                )

    # BLOCK 4: filter velocity
    if False:
        import loompy
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        with loompy.connect(fn_good) as dsl:
            print('Load samplesheet')
            samplesheet = pd.DataFrame([], index=dsl.ca['_index'])
            for col in dsl.ca.keys():
                if col != '_index':
                    samplesheet[col] = dsl.ca[col]

            print('Load featuresheet')
            featuresheet = pd.DataFrame([], index=dsl.ra['GeneName'])
            for col in dsl.ra.keys():
                featuresheet[col] = dsl.ra[col]

        fn_velo_raw = '../../data/sequencing/datasets/all_{:}/rnavelocity.loom'.format(version)
        with loompy.connect(fn_velo_raw) as dsl:
            snames = dsl.ca['CellID']

            # NOTE: velocity and htseq use different naming schemes, so we need
            # an adaptor here.
            ind_good_cells = []
            n_missing = 0
            for sn in samplesheet.index:
                for i, s in enumerate(snames):
                    if sn in s:
                        ind_good_cells.append(i)
                        break
                else:
                    n_missing += 1
                    print('Cell not found in velocity: {:}. Total # cells missing: {:}'.format(sn, n_missing))

            ind_good_cells = np.sort(ind_good_cells)
            n_good_cells = len(ind_good_cells)
            print('# good cells: {:}'.format(n_good_cells))

            gnames_count = featuresheet['GeneName'].values
            gnames = list(dsl.ra['Gene'])
            ind_good_genes = np.array([i for i, g in enumerate(gnames) if g in gnames_count])
            n_good_genes = len(ind_good_genes)
            print('# good genes: {:}'.format(n_good_genes))

            if n_good_cells < 100:
                raise ValueError('# cells too low')

            if n_good_genes < 100:
                raise ValueError('# genes too low')

            fn_velo_good = '../../data/sequencing/datasets/all_{:}/rnavelocity_good.loom'.format(version)

            print('Setting samplesheet')
            ca = {key: dsl.ca[key][ind_good_cells] for key in dsl.ca.keys()}
            print('Setting featuresheet')
            ra = {key: dsl.ra[key][ind_good_genes] for key in dsl.ra.keys()}

            matrices = {}
            for il, layer in enumerate(dsl.layers.keys()):
                print('Setting layer: {:}'.format(layer))
                x = dsl.layers[layer][:, :]
                matrices[layer] = x[ind_good_genes][:, ind_good_cells]

        print('Saving to file')
        loompy.create(fn_velo_good, matrices, ra, ca)
