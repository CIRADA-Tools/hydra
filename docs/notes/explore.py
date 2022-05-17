import re
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.table import QTable
from astropy.table import vstack


class Clump(object):
    def __init__(self,cluster_cat,injected_cat,injected_cat_label,bins,depth='deep'):
        self.cluster = QTable.read(cluster_cat)
        self.injected = QTable.read(injected_cat)
        self.iname = injected_cat_label
        self.bins = bins
        self.depth= 'deep' if depth == 'shallow_deep' else depth

        # add s/n columns
        self.injected.add_column(
            name = 'sn_deep',
            col = self.injected['flux_total'].to(u.mJy)/self.injected['rms_noise_bane_deep'].to(u.mJy)
        )
        self.injected.add_column(
            name = 'sn_shallow',
            col = self.injected['flux_total'].to(u.mJy)/self.injected['rms_noise_bane_shallow'].to(u.mJy)
        )

    def get_dsc(self):
        return self.dsc

    def __sigfigs(self,x, sig=2):
        return np.round(x, sig-np.int(np.floor(np.log10(np.abs(x))))-1)

    def get_injected(self):
        return self.injected

    def get_cluster(self,bin_no=None):
        if not bin_no is None:
            qt = self.cluster[((self.cluster['image_type']==self.depth)&(self.cluster['image_type']==self.depth))]
            tol = 1.0e-05
            if bin_no == 1:
               qt=qt[(((self.bins[bin_no-1]-tol)<=qt[f"sn_bane_{self.depth}"])&(qt[f"sn_bane_{self.depth}"]<self.bins[bin_no]))] 
            elif bin_no == (len(self.bins)-1):
               qt=qt[((self.bins[bin_no-1]<=qt[f"sn_bane_{self.depth}"])&(qt[f"sn_bane_{self.depth}"]<=(self.bins[bin_no]+tol)))] 
            else:
               qt=qt[((self.bins[bin_no-1]<=qt[f"sn_bane_{self.depth}"])&(qt[f"sn_bane_{self.depth}"]<self.bins[bin_no]))] 
            return qt
        return self.cluster

    def get_bins(self):
        return self.bins

    def get_flux(self,bin_no):
        return {
            'flux_min': self.get_bins()[bin_no-1],
            'flux_max': self.get_bins()[bin_no],
            'flux_avg': (self.get_bins()[bin_no-1]+self.get_bins()[bin_no])/2.0,
            'flux_err': (self.get_bins()[bin_no]-self.get_bins()[bin_no-1])/2.0,
        }

    def get_sn(self,bin_no):
        return (self.bins[bin_no-1]+self.bins[bin_no])/2.0

    def get_bin(self,bin_no): 
        tol = 1.0e-05
        if bin_no == 1:
           r=self.injected[(((self.bins[bin_no-1]-tol)<=self.injected[f"sn_{self.depth}"])&(self.injected[f"sn_{self.depth}"]<self.bins[bin_no]))] 
        elif bin_no == (len(self.bins)-1):
           r=self.injected[((self.bins[bin_no-1]<=self.injected[f"sn_{self.depth}"])&(self.injected[f"sn_{self.depth}"]<=(self.bins[bin_no]+tol)))] 
        else:
           r=self.injected[((self.bins[bin_no-1]<=self.injected[f"sn_{self.depth}"])&(self.injected[f"sn_{self.depth}"]<self.bins[bin_no]))] 
        r.sort(f"sn_{self.depth}") 
        return r 

    def get_cat_ids(self,bin_no):
        return np.unique(self.get_bin(bin_no)['id'])

    def get_unmatched_ids(self,bin_no):
        # usage: used for extracting reliability plot info, that is not in the bin -- defects
        qt = self.cluster[(self.cluster['image_type']==self.depth)]
        for match_id in np.unique(qt[(qt['source_finder']==self.iname)]['match_id']):
            qt = qt[(qt['match_id']!=match_id)]
        flx = self.get_flux(bin_no)
        qt = qt[((flx['flux_min']<=qt['sn_bane'])&(qt['sn_bane']<flx['flux_max']))]
        qt=qt['ref_id','cat_id','clump_id','match_id','source_finder','image_type','sn_bane','flux_total']
        return qt

    def get_deep_ids(self,bin_no):
        qt = self.cluster[(self.cluster['source_finder']!=self.iname)]
        flx = self.get_flux(bin_no)
        qt_deep = qt[((flx['flux_min']<=qt['sn_bane_shallow'])&(qt['sn_bane_shallow']<flx['flux_max']))]
        deep_ids = list(np.unique(qt_deep[(qt_deep['image_type']=='deep')]['match_id']))
        shallow_ids = list()
        qt_shallow = qt[(qt['image_type']=='shallow')]
        for deep_id in deep_ids:
            shallow_ids.extend(list(np.unique(qt_shallow[(qt_shallow['match_id']==deep_id)]['match_id'])))
        complement_ids = list()
        for deep_id in deep_ids:
            if not deep_id in shallow_ids:
                complement_ids.append(deep_id)
        return {
            'inj_deep': [{'c': self.cluster[(self.cluster['match_id']==m)]['clump_id'][0], 'm': m} for m in deep_ids],
            'det_shallow': [{'c': self.cluster[(self.cluster['match_id']==m)]['clump_id'][0], 'm': m} for m in shallow_ids],
            'cmp_shallow': [{'c': self.cluster[(self.cluster['match_id']==m)]['clump_id'][0], 'm': m} for m in complement_ids],
            'bin_info': {
                'bin_no': bin_no,
                'sn_min': flx['flux_min'],
                'sn_max': flx['flux_max'],
                'sn_avg': flx['flux_avg'],
            }
        }

    def get_clump_ids(self,bin_no):
        # usage: used for analysis completed plots
        # preamble
        source_finders = ['aegean','caesar','profound','pybdsf','selavy']
        qt = self.cluster[(self.cluster['image_type']==self.depth)]
        def filter(mids,sf=None):
            v_mids = list()
            for mid in mids:
                match = qt[(qt['match_id']==mid)]
                for source_finder in np.unique(match['source_finder']):
                    if source_finder in (source_finders if sf is None else [sf]):
                        v_mids.append(mid)
                        break
            return list(np.unique(v_mids))
        def translate(mids,is_detailed=True):
            def get_info(m):
                match=qt[(qt['match_id']==m)]
                if is_detailed:
                    return {
                        'c': match['clump_id'][0],
                        'm': match['match_id'][0],
                        # TODO: Needs some work...
                        #'sn': self.__sigfigs(match['sn_bane'][0].value,4),
                        #'flux_mJy': self.__sigfigs(match['flux_total'][0].value,4),
                    }
                return match['clump_id'][0]
            if isinstance(mids,list):
                cids = list()
                for mid in mids:
                    cids.append(get_info(mid))
                #cids=list(np.unique(cids))
            elif isinstance(mids,dict):
                cids = dict()
                for key in mids:
                    cids[key]=list()
                    for mid in mids[key]:
                        cids[key].append(get_info(mid))
                    #cids[key]=list(np.unique(cids[key]))
            else:
                cids=get_info(mid)
            return cids

        # get catalogue match ids
        cat_ids = self.get_cat_ids(bin_no)
        match_ids = list()
        for cat_id in cat_ids:
            match_ids.extend(list(qt[((qt['cat_id']==cat_id)&(qt['source_finder']==self.iname))]['match_id']))
        match_ids = list(np.unique(match_ids))

        # break it down
        breakdowns = {'detected': filter(match_ids)}
        missed = list()
        for match_id in match_ids:
            if not match_id in breakdowns['detected']:
                missed.append(match_id)
        breakdowns['missed']=list(np.unique(missed))
        for source_finder in source_finders:
            breakdowns[source_finder]=filter(match_ids,source_finder)

        #return translate(breakdowns)
        return {
            'data': translate(breakdowns),
            'summary': translate(breakdowns,False),
            'stats': {
                'detected': len(breakdowns['detected']),
                'injected': len(breakdowns['detected'])+len(breakdowns['missed']),
            },
        }

    def get_bin_stats(self,bin_no_a,bin_no_b):
        def get_binfo(bin_no):
            datum = self.get_clump_ids(bin_no)
            ds = dict()
            for sf in ['aegean','caesar','profound','pybdsf','selavy']:
                ds[sf] = len(datum['data'][sf])
            return {
                'finders': ds,
                'summary': {
                    'detected': len(datum['data']['detected']),
                    'missed':   len(datum['data']['missed']),
                    'total':    len(datum['data']['detected'])+len(datum['data']['missed']),
                },
            }
        qt = QTable(names=['sn','sn_error','total','detected','fraction_detected','missed'],dtype=[np.float,np.float,np.int,np.int,np.float,np.int])
        for bino in range(bin_no_a,bin_no_b+1):
            df = get_binfo(bino)['summary']
            ds = self.get_flux(bino)
            qt.add_row([ds['flux_avg'],ds['flux_err'],df['total'],df['detected'],df['detected']/df['total'],df['missed']])
        return qt

    def get_clump_missed(self,clump_id):
        qt = QTable(names=['bin_no','sn','sn_error','clump_id','match_id','missed'],dtype=[np.int,np.float,np.float,np.int,np.int,np.int])
        for bino in range(1,len(self.get_bins())):
            df = self.get_clump_ids(bino)['data']
            for d in df['missed']:
                if d['c']==clump_id:
                    flx = self.get_flux(bino)
                    datum = [bino,flx['flux_avg'],flx['flux_err'],clump_id,d['m'],1]
                    print(datum)
                    qt.add_row(datum)
            for d in df['detected']:
                if d['c']==clump_id:
                    flx = self.get_flux(bino)
                    datum = [bino,flx['flux_avg'],flx['flux_err'],clump_id,d['m'],0]
                    print(datum)
                    qt.add_row(datum)
        return qt


class sdCompleteness(Clump):
    def __init__(self,completenenss_cat,cluster_cat,injected_cat,injected_cat_label,bins,depth='shallow_deep'):
        super(sdCompleteness,self).__init__(cluster_cat,injected_cat,injected_cat_label,bins)
        self.completeness = QTable.read(completenenss_cat)
        self.completeness['bin_no'] = self.completeness['bin_no']+1
        self.completeness.add_column(
            name = 'clump_id',
            col = self.completeness['clump_id_deep'],
            index = 2
        )
        self.completeness.add_column(
            name = 'match_id_deep',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['match_id'][0] if rid>-1 else -1 for rid in self.completeness['ref_id_deep']],
            index = 3
        )
        self.completeness.add_column(
            name = 'match_id_shallow',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['match_id'][0] if rid>-1 else -1 for rid in self.completeness['ref_id_shallow']],
            index = 4
        )

    def get_completeness(self,bin_no=None,source_finder=None):
        qt = self.completeness
        if not bin_no is None:
            flx = self.get_flux(bin_no)
            qt = qt[((flx['flux_min']<=qt['sn_avg'])&(qt['sn_avg']<flx['flux_max']))]
        if not source_finder is None:
            qt = qt[(qt['source_finder']==source_finder)]
        return qt[
            'bin_no',
            'ref_id_deep',
            'match_id_deep',
            'ref_id_shallow',
            'match_id_shallow',
            'clump_id',
            'sn_avg',
            'source_finder',
            'completeness',
            'sn_deep',
            'sn_shallow',
            'flux_total_deep',
            'flux_total_shallow',
            'clump_size_deep',
            'clump_size_shallow',
            #'rms_noise_bane_deep',
            #'rms_noise_bane_shallow',
        ]

    def get_matched(self,bin_no,source_finder=None):
        qt = self.get_completeness(bin_no,source_finder)
        return qt[(qt['ref_id_shallow']>-1)]

    def get_unmatched(self,bin_no,source_finder=None):
        qt = self.get_completeness(bin_no,source_finder)
        return qt[(qt['ref_id_shallow']<0)]

    def get_stats(self,bin_no):
        qt = self.get_completeness(bin_no)
        stats = dict()
        stats['bin'] = {'no': bin_no,'sn': self.get_flux(bin_no)['flux_avg']}
        for source_finder in np.unique(qt['source_finder']):
            stats[source_finder] = {
                'detecthed': len(qt[((qt['source_finder']==source_finder)&(qt['ref_id_shallow']>-1))]),
                'missed':    len(qt[((qt['source_finder']==source_finder)&(qt['ref_id_shallow']<0))]),
                'total':     len(qt[(qt['source_finder']==source_finder)])
            }
        stats['global'] = {
            'detecthed': len(qt[(qt['ref_id_shallow']>-1)]),
            'missed':    len(qt[(qt['ref_id_shallow']<0)]),
            'total': len(qt)
        }
        return stats


class sdReliability(Clump):
    def __init__(self,reliability_cat,cluster_cat,injected_cat,injected_cat_label,bins,depth='shallow_deep'):
        super(sdReliability,self).__init__(cluster_cat,injected_cat,injected_cat_label,bins)
        self.reliability = QTable.read(reliability_cat)
        self.reliability['bin_no'] = self.reliability['bin_no']+1
        self.reliability.add_column(
            name = 'clump_id',
            col = self.reliability['clump_id_shallow'],
            index = 2
        )
        self.reliability.add_column(
            name = 'match_id_shallow',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['match_id'][0] if rid>-1 else -1 for rid in self.reliability['ref_id_shallow']],
            index = 3
        )
        self.reliability.add_column(
            name = 'match_id_deep',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['match_id'][0] if rid>-1 else -1 for rid in self.reliability['ref_id_deep']],
            index = 4
        )

    def get_reliability(self,bin_no=None,source_finder=None):
        qt = self.reliability
        if not bin_no is None:
            flx = self.get_flux(bin_no)
            qt = qt[((flx['flux_min']<=qt['sn_avg'])&(qt['sn_avg']<flx['flux_max']))]
        if not source_finder is None:
            qt = qt[(qt['source_finder']==source_finder)]
        return qt[
            'bin_no',
            'ref_id_shallow',
            'match_id_shallow',
            'ref_id_deep',
            'match_id_deep',
            'clump_id',
            'sn_avg',
            'source_finder',
            'reliability',
            'sn_shallow',
            'sn_deep',
            'flux_total_shallow',
            'flux_total_deep',
            'clump_size_shallow',
            'clump_size_deep',
            #'rms_noise_bane_shallow',
            #'rms_noise_bane_deep',
        ]

    def get_matched(self,bin_no,source_finder=None):
        qt = self.get_reliability(bin_no,source_finder)
        return qt[(qt['ref_id_deep']>-1)]

    def get_unmatched(self,bin_no,source_finder=None):
        qt = self.get_reliability(bin_no,source_finder)
        return qt[(qt['ref_id_deep']<0)]

    def get_stats(self,bin_no):
        qt = self.get_reliability(bin_no)
        stats = dict()
        stats['bin'] = {'no': bin_no,'sn': self.get_flux(bin_no)['flux_avg']}
        for source_finder in np.unique(qt['source_finder']):
            stats[source_finder] = {
                'detecthed': len(qt[((qt['source_finder']==source_finder)&(qt['ref_id_deep']>-1))]),
                'missed':    len(qt[((qt['source_finder']==source_finder)&(qt['ref_id_deep']<0))]),
                'total':     len(qt[(qt['source_finder']==source_finder)])
            }
        stats['global'] = {
            'detecthed': len(qt[(qt['ref_id_deep']>-1)]),
            'missed':    len(qt[(qt['ref_id_deep']<0)]),
            'total': len(qt)
        }
        return stats


class fRatio(object):
    def __init__(self,injected_cat_label='EMUSim2x2'):
        # load base data
        self.cluster = QTable().read("sim_2x2/xtables/emu_simulated_2x2.hydra.cluster_catalogue.fits")
        self.rratios = {
            'deep': {
                'aegean':   QTable().read("sim_2x2/xratios/aegean_deep_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'caesar':   QTable().read("sim_2x2/xratios/caesar_deep_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'profound': QTable().read("sim_2x2/xratios/profound_deep_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'pybdsf':   QTable().read("sim_2x2/xratios/pybdsf_deep_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'selavy':   QTable().read("sim_2x2/xratios/selavy_deep_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
            },
            'shallow': {
                'aegean':   QTable().read("sim_2x2/xratios/aegean_shallow_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'caesar':   QTable().read("sim_2x2/xratios/caesar_shallow_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'profound': QTable().read("sim_2x2/xratios/profound_shallow_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'pybdsf':   QTable().read("sim_2x2/xratios/pybdsf_shallow_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
                'selavy':   QTable().read("sim_2x2/xratios/selavy_shallow_to_injected_flux_ratio_vs_sn_plot_data_table.fits"),
            },
            'shallow_deep': {
                'aegean':   QTable().read("sim_2x2/xratios/aegean_shallow_to_deep_flux_ratio_vs_sn_plot_data_table.fits"),
                'caesar':   QTable().read("sim_2x2/xratios/caesar_shallow_to_deep_flux_ratio_vs_sn_plot_data_table.fits"),
                'profound': QTable().read("sim_2x2/xratios/profound_shallow_to_deep_flux_ratio_vs_sn_plot_data_table.fits"),
                'pybdsf':   QTable().read("sim_2x2/xratios/pybdsf_shallow_to_deep_flux_ratio_vs_sn_plot_data_table.fits"),
                'selavy':   QTable().read("sim_2x2/xratios/selavy_shallow_to_deep_flux_ratio_vs_sn_plot_data_table.fits"),
            },
        }
        self.iname = injected_cat_label

        # compilate deep
        def get_flux_in(rid,depth):
            match_id = self.cluster[(self.cluster['ref_id']==rid)]['match_id'][0]
            flux_in = self.cluster[((self.cluster['match_id']==match_id)&(self.cluster['source_finder']==self.iname)&(self.cluster['image_type']==depth))]['flux_total'][0]
            return flux_in
        self.rdeep = QTable()
        for source_finder in self.rratios['deep']:
            self.rdeep = self.rratios['deep'][source_finder] if len(self.rdeep)<1 else vstack([self.rdeep,self.rratios['deep'][source_finder]])
        if 'fluxes' in self.rdeep.colnames:
            self.rdeep.rename_column('fluxes','sn_out')
        if 'in_flux' in self.rdeep.colnames:
            self.rdeep.rename_column('in_flux','sn_in')
        self.rdeep.remove_columns(['minus_1_sigma','frac_1_sigma','plus_1_sigma','minus_2_sigma','frac_2_sigma','plus_2_sigma'])
        self.rdeep.add_column(
            name = 'flux_total_in',
            col = [get_flux_in(rid,'deep') for rid in self.rdeep['ref_id']]
        )
        self.rdeep.add_column(
            name = 'flux_total_out',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['flux_total'][0] for rid in self.rdeep['ref_id']]
        )

        # compilate shallow
        self.rshallow = QTable()
        for source_finder in self.rratios['shallow']:
            self.rshallow = self.rratios['shallow'][source_finder] if len(self.rshallow)<1 else vstack([self.rshallow,self.rratios['shallow'][source_finder]])
        if 'fluxes' in self.rshallow.colnames:
            self.rshallow.rename_column('fluxes','sn_out')
        if 'in_flux' in self.rshallow.colnames:
            self.rshallow.rename_column('in_flux','sn_in')
        self.rshallow.remove_columns(['minus_1_sigma','frac_1_sigma','plus_1_sigma','minus_2_sigma','frac_2_sigma','plus_2_sigma'])
        self.rshallow.add_column(
            name = 'flux_total_in',
            col = [get_flux_in(rid,'shallow') for rid in self.rshallow['ref_id']]
        )
        self.rshallow.add_column(
            name = 'flux_total_out',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['flux_total'][0] for rid in self.rshallow['ref_id']]
        )

        # compilate shallow-deep
        self.rshallow_deep = QTable()
        for source_finder in self.rratios['shallow_deep']:
            self.rshallow_deep = self.rratios['shallow_deep'][source_finder] if len(self.rshallow_deep)<1 else vstack([self.rshallow_deep,self.rratios['shallow_deep'][source_finder]])
        if 'fluxes' in self.rshallow_deep.colnames:
            self.rshallow_deep.rename_column('fluxes','sn_out')
        if 'in_flux' in self.rshallow_deep.colnames:
            self.rshallow_deep.rename_column('in_flux','sn_in')
        self.rshallow_deep.remove_columns(['minus_1_sigma','frac_1_sigma','plus_1_sigma','minus_2_sigma','frac_2_sigma','plus_2_sigma'])
        self.rshallow_deep.add_column(
            name = 'flux_total_in',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['flux_total'][0] for rid in self.rshallow_deep['in_ref_id']]
        )
        self.rshallow_deep.add_column(
            name = 'flux_total_out',
            col = [self.cluster[(self.cluster['ref_id']==rid)]['flux_total'][0] for rid in self.rshallow_deep['ref_id']]
        )

        # do 3sigma cuts
        self.sdeep_cut = QTable()
        for row in self.rdeep:
            if row['flux_ratio']<=row['minus_3_sigma'] or row['plus_3_sigma']<=row['flux_ratio']:
                self.sdeep_cut = row if len(self.sdeep_cut)<1 else vstack([self.sdeep_cut,row])
        self.sshallow_cut = QTable()
        for row in self.rshallow:
            if row['flux_ratio']<=row['minus_3_sigma'] or row['plus_3_sigma']<=row['flux_ratio']:
                self.sshallow_cut = row if len(self.sshallow_cut)<1 else vstack([self.sshallow_cut,row])
        self.sshallow_deep_cut = QTable()
        for row in self.rshallow_deep:
            if row['flux_ratio']<=row['minus_3_sigma'] or row['plus_3_sigma']<=row['flux_ratio']:
                self.sshallow_deep_cut = row if len(self.sshallow_deep_cut)<1 else vstack([self.sshallow_deep_cut,row])

    def get_sdeep(self,sn_cut=False):
        return self.sdeep_cut if sn_cut else self.rdeep

    def get_sdeep_cut(self,s_in_min=None,s_in_max=None):
        sdeep = self.get_sdeep(True)
        if not s_in_min is None and not s_in_max is None:
            return sdeep[((s_in_min<=sdeep['sn_in'])&(sdeep['sn_in']<=s_in_max))]
        elif not s_in_min is None and s_in_max is None:
            return sdeep[(s_in_min<=sdeep['sn_in'])]
        elif s_in_min is None and not s_in_max is None:
            return sdeep[(sdeep['sn_in']<=s_in_max)]
        return sdeep

    def get_sshallow(self,sn_cut=False):
        return self.sshallow_cut if sn_cut else self.rshallow

    def get_sshallow_cut(self,s_in_min=None,s_in_max=None):
        sshallow = self.get_sshallow(True)
        if not s_in_min is None and not s_in_max is None:
            return sshallow[((s_in_min<=sshallow['sn_in'])&(sshallow['sn_in']<=s_in_max))]
        elif not s_in_min is None and s_in_max is None:
            return sshallow[(s_in_min<=sshallow['sn_in'])]
        elif s_in_min is None and not s_in_max is None:
            return sshallow[(sshallow['sn_in']<=s_in_max)]
        return sshallow

    def get_sshallow_deep(self,sn_cut=False):
        return self.sshallow_deep_cut if sn_cut else self.rshallow_deep

    def get_sshallow_deep_cut(self,s_in_min=None,s_in_max=None):
        sshallow_deep = self.get_sshallow_deep(True)
        if not s_in_min is None and not s_in_max is None:
            return sshallow_deep[((s_in_min<=sshallow_deep['sn_in'])&(sshallow_deep['sn_in']<=s_in_max))]
        elif not s_in_min is None and s_in_max is None:
            return sshallow_deep[(s_in_min<=sshallow_deep['sn_in'])]
        elif s_in_min is None and not s_in_max is None:
            return sshallow_deep[(sshallow_deep['sn_in']<=s_in_max)]
        return sshallow_deep



class Scatter:
    def __init__(self,fRatio_file):
        qt = QTable().read(fRatio_file)
        def r_3s(flx):
            qs = qt[(qt['in_flux']>=flx)]
            qb = qs[(qs['minus_3_sigma']<=qs['flux_ratio'])&(qs['flux_ratio']<=qs['plus_3_sigma'])]
            return len(qb)/len(qs)
        #qt.add_column(
        #    name = 'r_3sigma',
        #    col = [r_3s(in_flux) for in_flux in qt['in_flux']]
        #)
        qt.add_column(
            name = 'scatter_3s',
            #col = 1 - qt['r_3sigma']
            col = 1 - qt['frac_3_sigma']
        )
        self.qt = qt.copy()
        self.qt.sort('in_flux')

    def get_qt(self):
        return self.qt

    def get_1_to_10(self):
        qt = self.qt.copy()
        qt = qt[(1<=qt['in_flux'])&(qt['in_flux']<=10)]
        stats = {
            'avg': np.mean(qt['scatter_3s']),
            'min': abs(np.mean(qt['scatter_3s'])-min(qt['scatter_3s'])),
            'max': abs(max(qt['scatter_3s'])-np.mean(qt['scatter_3s'])),
            'std': np.std(qt['scatter_3s']),
        }
        print(f"Scatter: {stats['avg']*100:0.1f}+{stats['max']*100:0.1f}-{stats['min']*100:0.1f}%")
        return stats

    def get_n(self,flx):
        qt = self.qt.copy()
        qs = qt[(qt['in_flux']>=flx)]
        #qb = qs[(qs['minus_3_sigma']<=qs['flux_ratio'])&(qs['flux_ratio']<=qs['plus_3_sigma'])]
        #return 1 - len(qb)/len(qs)
        sm = np.abs(qs['in_flux'][0]-flx)
        scatter=qs['scatter_3s'][0]
        for rw in qt:
            if sm > np.abs(rw['in_flux']-flx):
                sm = np.abs(rw['in_flux']-flx)
                scatter= rw['scatter_3s']
        return scatter


def getno(sn,metric='c'): 
     bin_no = -1 
     if metric != 'r': # default to completess
         print("MODE: COMPLETENESS")
         bins = [0.150008,0.179516,0.214829,0.257088,0.307659,0.368178,0.440603,0.527273,0.630993,0.755115,0.903653,1.081410,1.294133,1.548701,1.853345,2.217915,2.654199,3.176305,3.801114,4.548828,5.443624,6.514436,7.795885,9.329408,11.164589,13.360767,15.988954,19.134129,22.897989,27.402235,32.792508,39.243098,46.962579,56.200553,67.255722,80.485545,96.317796,115.264398,137.937972,165.071648,197.542769,236.401260,282.903576,338.553327,405.149900,484.846635,580.220456,694.355189,830.941281,994.395121,1190.001845]
     else: # else reliability
         print("MODE: RELIABILITY")
         bins = [1.875198,2.133573,2.427548,2.762028,3.142594,3.575597,4.068262,4.628808,5.266589,5.992247,6.817889,7.757294,8.826134,10.042244,11.425916,13.000238,14.791478,16.829525,19.148384,21.786748,24.788638,28.204145,32.090259,36.511821,41.542608,47.266564,53.779196,61.189172,69.620133,79.212756,90.127100,102.545278,116.674498,132.750514,151.041567,171.852857,195.531633,222.472995,253.126476,288.003553,327.686174,372.836472,424.207812,482.657361,549.160391,624.826554,710.918393,808.872410,920.323040,1047.129915,1191.408897]
     for i in range(len(bins)-1): 
         if bins[i]<= sn and sn <=bins[i+1]: 
             bin_no = i+1 
             print(f"S/N: {(bins[i+1]+bins[i])/2.0}")
             break 
     return bin_no 


# debug function...
def translate(cat_id,finder,depth):
    src_root = f"/Users/susy/cirada/emu_pipeline/data/real/emu/trial4/emu_pilot_sample_2x2deg.hydra_dir"
    dst_root = f"/Users/susy/cirada/emu_pipeline/data/real/emu/trial4/debug/emu_pilot_sample_2x2deg.hydra_dir"
    src_cat = f"{src_root}/catalogues/{depth}/emu_pilot_sample_2x2deg.hydra.{finder}.{depth}.fits"
    dst_cat = f"{dst_root}/catalogues/{depth}/emu_pilot_sample_2x2deg.hydra.{finder}.{depth}.fits"
    cluster = f"{dst_root}/emu_pilot_sample_2x2deg.hydra.cluster_catalogue.fits"
    sf = Table().read(src_cat)
    df = Table().read(dst_cat)
    cf = QTable().read(cluster)
    sf_row = sf[(sf['id']==cat_id)]['id','ra','dec']
    df_row = df[(df['ra']==sf_row['ra'])&(df['dec']==sf_row['dec'])]['id','ra','dec']
    cf_row = cf[(cf['cat_id']==df_row['id'])&(cf['source_finder']==finder)&(cf['image_type']==depth)]['ref_id','clump_id','match_id','cat_id','source_finder','image_type','ra','dec','extent_semimajor','extent_semiminor','extent_angle']
    print(f"Translating: ['cat_id': {cat_id}, 'finder': {finder}, 'depth': {depth}]")
    print(f"Source> {src_root}")
    print(f"Destin> {dst_root}")
    print(f"catgen> emu_pilot_sample_2x2deg.hydra.{finder}.{depth}.fits:")
    print(f"clustr> emu_pilot_sample_2x2deg.hydra.cluster_catalogue.fits")
    print(f"")
    print(f"Source Row:")
    print(sf_row)
    print(f"")
    print(f"Destin Row:")
    print(df_row)
    print(f"")
    print(f"Cluster Table Row:")
    print(cf_row)
    print(f"")
    print(f"[DONE]")
    print(f"")
    print(f"CLUMP_ID: {cf_row['clump_id'][0]}")
    print(f"")
    return cf_row



def confc():
    qu = None
    qt = Table().read("emu_2x2/xtables/emu_pilot_sample_2x2deg.hydra.cluster_catalogue.fits")
    qt = qt[(qt['image_type']=='deep')]
    for mid in np.unique(qt["match_id"]):
        df = qt[(qt['match_id']==mid)].copy()
        if len(df)==1:
            qu = df if qu is None else vstack([qu,df])
    return qu











