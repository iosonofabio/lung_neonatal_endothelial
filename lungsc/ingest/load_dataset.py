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
import numpy as np
import pandas as pd

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, CountsTable, FeatureSheet, SampleSheet


versions = ('20190325', '20190513', '20190620', '20190828', '20191204', '20200129', '20200722', '20201006')


class DatasetLung(Dataset):
    @classmethod
    def load(
            cls,
            version=versions[-1],
            preprocess=False,
            include_hyperoxia=False,
            include_doublets=False,
            ):
        import loompy
        fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
        with loompy.connect(fn_good) as ds:

            print('Load samplesheet')
            if '_index' in ds.ca:
                ind_col = '_index'
            else:
                ind_col = 'CellID'
            samplesheet = pd.DataFrame([], index=ds.ca[ind_col])
            for col in ds.ca.keys():
                if col != ind_col:
                    samplesheet[col] = ds.ca[col]

            if 'Time [days]' not in samplesheet.columns:
                samplesheet['Time [days]'] = ''
                samplesheet.loc[samplesheet['Timepoint'] == 'E18.5', 'Time [days]'] = -1
                samplesheet.loc[samplesheet['Timepoint'] == 'P1', 'Time [days]'] = 1
                samplesheet.loc[samplesheet['Timepoint'] == 'P7', 'Time [days]'] = 7
                samplesheet.loc[samplesheet['Timepoint'] == 'P21', 'Time [days]'] = 21

            if 'Treatment' not in samplesheet.columns:
                tmp = np.array(['normal', 'hyperoxia'])
                samplesheet['Treatment'] = tmp[[int('hyperoxia' in x) for x in samplesheet.index]]

            if 'TimepointHO' not in samplesheet.columns:
                samplesheet['TimepointHO'] = samplesheet['Timepoint'].copy()
                samplesheet.loc[
                        samplesheet['Treatment'] == 'hyperoxia',
                        'TimepointHO'] = 'P7_HO'

            if 'Mousename' not in samplesheet.columns:
                tt = samplesheet['Treatment'].copy()
                tt.loc[tt == 'normal'] = ''
                samplesheet['Mousename'] = samplesheet['Gender'] + samplesheet['Timepoint']+tt

            print('Load featuresheet')
            featuresheet = pd.DataFrame([], index=ds.ra['GeneName'])
            for col in ds.ra.keys():
                featuresheet[col] = ds.ra[col]

            print('Load normalized counts')
            counts_table = pd.DataFrame(
                    ds[:, :],
                    index=featuresheet.index,
                    columns=samplesheet.index,
                    )

        print('Convert to singlet data structures')
        samplesheet = SampleSheet(samplesheet)
        featuresheet = FeatureSheet(featuresheet)
        counts_table = CountsTable(counts_table)
        counts_table._normalized = 'counts_per_million'

        print('Set spikeins and other features')
        counts_table._spikeins = tuple(counts_table.index[counts_table.index.str.startswith('ERCC-')].tolist())
        counts_table._otherfeatures = tuple(counts_table.index[counts_table.index.str.startswith('__')].tolist())

        self = cls(
                counts_table=counts_table,
                featuresheet=featuresheet,
                samplesheet=samplesheet,
                )

        if not include_hyperoxia:
            self.query_samples_by_metadata('Treatment == "normal"', inplace=True)

        if preprocess:
            feas = self.featurenames
            good_genes = pd.Series(
                data=np.ones(self.n_features, bool),
                index=self.featurenames)
            good_genes[feas.str.startswith('mt-')] = False
            good_genes[feas.str.startswith('Rpl')] = False
            good_genes[feas.str.startswith('Rps')] = False
            for i in range(10):
                good_genes[feas.str.startswith('Gm{:}'.format(i))] = False
            good_genes = good_genes.index[good_genes]
            self.query_features_by_name(good_genes, inplace=True)

        if not include_doublets:
            self.query_samples_by_metadata('doublet == 0', inplace=True)

        return self


if __name__ == '__main__':

    ds = DatasetLung.load()
