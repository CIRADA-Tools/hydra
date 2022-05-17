import re
import sys
import time
import math
import yaml as yml
from astropy.table import QTable
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from progressbar import *

#import pyximport
#pyximport.install()
#from .clusterizer import clusterize
#
#def clusterize(pairings):
#    return clusterize(pairings)



def time_string(seconds):
    def f(n):
        return math.floor(n)
    t_str = ""
    if seconds > 0:
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        ms = int(1000.0*(seconds % 1))
        us = round(1000.0*((1000.0*seconds) % 1))
        t_str += ("%dh" % h) if f(h) > 0 else ""
        t_str += ("%dm" % m) if f(m) > 0 else ""
        t_str += ("%ds" % s) if f(s) > 0 else ""
        t_str += ("%dms" % ms) if f(ms) > 0 else ""
        t_str += ("%dus" % us) if f(us) > 0 else ""
    else:
        t_str += "0s"
    return t_str

def prt_time(timer):
    sys.stderr.write("[{0}]\n".format(time_string(time.clock() - timer)))


class Cluster(object):
    def __init__(self,catalogues,match_on_self=True):
        # tabulate a list of surveys from the catalgoues...
        self.surveys = [c["survey"] for c in catalogues]
        self.signature_dypes = {
            'ref_id':   {'dtype': np.int64, 'units': None},
            'clump_id': {'dtype': np.int64, 'units': None},
            'survey': {'dtype': np.dtype(f"S{np.max([len(c['survey']) for c in catalogues])}"), 'units': None},
            'cat_no': {'dtype': np.int64, 'units': None},
            'ra':  {'dtype': np.float64, 'units': catalogues[0]['coords'].ra.unit },
            'dec': {'dtype': np.float64, 'units': catalogues[0]['coords'].dec.unit},
            'extent_semimajor': {'dtype': np.float64, 'units': catalogues[0]['a_extents'].unit},
            'extent_semiminor': {'dtype': np.float64, 'units': catalogues[0]['b_extents'].unit},
            'extent_angle':     {'dtype': np.float64, 'units': catalogues[0]['t_extents'].unit}
        }
        self.signature = list(self.signature_dypes.keys())

        # time to clusterize!
        sys.stderr.write("Clusterizing...\n")
        timer = time.clock()
        self.clumps = self.__make_clumps(catalogues,match_on_self)
        prt_time(timer)


    def __metric(self,cat_1,idx_1,cat_2,idx_2):
        d_ra  = (cat_2['coords'][idx_2].ra  - cat_1['coords'][idx_1].ra)*np.cos(np.pi*cat_1['coords'][idx_1].dec.to(u.deg).value/180.0)
        d_dec = cat_2['coords'][idx_2].dec - cat_1['coords'][idx_1].dec
        eta = -np.arctan2(d_dec.value,d_ra.value)
        beta_1 = np.pi*cat_1['t_extents'][idx_1].value/180.0 - eta
        beta_2 = np.pi*cat_2['t_extents'][idx_2].value/180.0 - eta
        numerator = cat_1['a_extents'][idx_1]*cat_1['b_extents'][idx_1]
        denominator = np.sqrt((cat_1['a_extents'][idx_1]*np.cos(beta_1))**2+(cat_1['b_extents'][idx_1]*np.sin(beta_1))**2)
        r_1 = np.divide(numerator,denominator,out=np.zeros_like(numerator),where=numerator!=0)
        numerator = cat_2['a_extents'][idx_2]*cat_2['b_extents'][idx_2]
        denominator = np.sqrt((cat_2['a_extents'][idx_2]*np.cos(beta_2))**2+(cat_2['b_extents'][idx_2]*np.sin(beta_2))**2)
        r_2 = np.divide(numerator,denominator,out=np.zeros_like(numerator),where=numerator!=0)
        return (r_1 + r_2)


    def __get_pairings(self,catalogues,match_on_self=False):
        rows = sum([len(catalogue['ref_ids']) for catalogue in catalogues])
        timer = time.clock()
        sys.stderr.write("Analysing {:0,} rows...\n".format(rows))
        num_catalogues = len(catalogues)
        pairings = list()
        for i in range(num_catalogues):
            for j in range(i,num_catalogues):
                if (catalogues[i]['survey'] == catalogues[j]['survey']) and (not match_on_self):
                    continue
                sys.stderr.write("> Matching {0} and {1} components /w separation <= {2}...\n".format(catalogues[i]['survey'],catalogues[j]['survey'],catalogues[i]['max_extent']+catalogues[j]['max_extent']))
                idxc, idxcatalog, d2d, d3d = catalogues[j]['coords'].search_around_sky(catalogues[i]['coords'],catalogues[i]['max_extent']+catalogues[j]['max_extent'])
                filter = d2d <= self.__metric(catalogues[i],idxc,catalogues[j],idxcatalog)
                idxc = idxc[filter]
                idxcatalog = idxcatalog[filter]
                del d2d
                del d3d
                sys.stderr.write("> {0} ref_ids for [{1},{2}]...\n".format("Reducing" if i==j else "Collating",catalogues[i]['survey'],catalogues[j]['survey']))
                if i == j:
                    matches = np.column_stack((catalogues[i]['ref_ids'][idxc],catalogues[j]['ref_ids'][idxcatalog]))[catalogues[i]['ref_ids'][idxc] < catalogues[j]['ref_ids'][idxcatalog]]
                else:
                    matches = np.column_stack((catalogues[i]['ref_ids'][idxc],catalogues[j]['ref_ids'][idxcatalog]))
                pairings.append(matches)
                del idxc
                del idxcatalog
                sys.stderr.write("> Pairings: {:0,}\n".format(len(matches)))
        sys.stderr.write("> Tabulating...\n")
        lhs = list()
        rhs = list()
        for matches in pairings:
           lhs.extend(matches[:,0])
           rhs.extend(matches[:,1])
        pairings = np.column_stack((lhs,rhs))
        sys.stderr.write("> Sorting...\n")
        if len(pairings) > 0:
            icps = np.array([[idx,pairings[idx][0]] for idx in range(len(pairings))],dtype=np.int64)
            pairings = np.array([pairings[i] for i in icps[icps[:,1].argsort()][:,0]],dtype=np.int64)
        sys.stderr.write("> Done!\n")
        sys.stderr.write("Total Pairings: {:0,}\n".format(len(pairings)))
        prt_time(timer)

        return pairings


    def __clusterize(self,pairings):
        ##return self.__clusterize2(pairings)
        #import pyximport
        #pyximport.install()
        #import clusterizer
        #return clusterizer.clusterize(pairings)
        import pyximport
        #pyximport.install()
        pyximport.install(language_level=3)
        from .clusterizer import clusterize
        return clusterize(pairings)


    def __transpose(self,catalogues):
        sources = list()
        for catalogue in catalogues:
            ref_id = catalogue["ref_ids"]
            survey = catalogue["survey"]
            cat_no = catalogue["cat_nos"]
            ra  = catalogue["coords"].ra.value
            dec = catalogue["coords"].dec.value
            e_major = catalogue["a_extents"].value
            e_minor = catalogue["b_extents"].value
            e_angle = catalogue["t_extents"].value
            for i in range(len(ref_id)):
                sources.append([ref_id[i],None,survey,cat_no[i],ra[i],dec[i],e_major[i],e_minor[i],e_angle[i]])
        return sources


    def __merge(self,clusters,catalogues):
        sources = self.__transpose(catalogues)

        clumps = list()
        rows = len(clusters)
        row = 0
        sys.stderr.write("Max Counts: ")
        timer = time.clock()
        max_counts = max([s[0] for s in sources])+1
        counts = np.zeros(max_counts,dtype=np.int64)
        sys.stderr.write("{0} [{1}]\n".format(max_counts,time_string(time.clock() - timer)))
        sys.stderr.write("Re-integrating clusters with catalogue info...\n")
        timer = time.clock()
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
        for cluster in clusters:
            clump = list()
            for component in cluster:
                clump.append(sources[component-1])
                if sources[component-1][0] != component:
                    sys.stderr.write("Whoops! {0} != {1}\n".format(sources[component-1][0],component))
                    exit()
                counts[component] += 1
            clumps.append(clump)
            row += 1
            pbar.update(row if row < rows else rows)
        sys.stderr.write("\n")
        prt_time(timer)

        rows = len(counts)
        row = 1
        singles = list()
        sys.stderr.write("Grabbing non-clustered components from catalogues...\n")
        timer = time.clock()
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
        while row < rows:
            if counts[row] == 0:
                singles.append([sources[row-1]])
            row += 1
            pbar.update(row if row < rows else rows)
        sys.stderr.write("\n")
        prt_time(timer)

        timer = time.clock()
        sys.stderr.write("Clusters: %d\n" % len(clumps))
        sys.stderr.write("Singles: %d\n" % len(singles))
        sys.stderr.write("Merging and sorting...\n")
        clumps.extend(singles)
        icps = np.array([[idx,clumps[idx][0][0]] for idx in range(len(clumps))],dtype=np.int64)
        clumps = [clumps[i] for i in icps[icps[:,1].argsort()][:,0]]
        sys.stderr.write("Clumps: %d\n" % len(clumps))
        prt_time(timer)

        timer = time.clock()
        sys.stderr.write("Renomalizing ref_id's and adding clump indexes...\n")
        ref_idx = 1
        clump_idx = 1
        for clump in clumps:
            for component in clump:
                component[0] = ref_idx
                component[1] = clump_idx
                ref_idx += 1
            clump_idx += 1
        prt_time(timer)

        return clumps

    def __make_clumps(self,catalogues,match_on_self):
        return self.__merge(self.__clusterize(self.__get_pairings(catalogues,match_on_self)),catalogues)

    def get_clumps(self):
        return self.clumps

    def prt_clumps(self):
        print("clumps = [")
        for clump in self.clumps:
            print(" [")
            for component in clump:
                print("  {0},".format(component))
            print(" ],")
        print("]")

    def get_catalogue(self):
        catalogue = list()
        for clump in self.clumps:
            for component in clump:
                catalogue.append(component)
        return catalogue

    def prt_catalogue(self):
        catalogue = self.get_catalogue()
        print("catalogue = [")
        for item in catalogue:
            print(" {0},".format(item))
        print("]")

    def get_signature(self,is_mysql=False):
        if is_mysql:
            return ",".join(['`%s`' % s for s in self.signature])
        return ",".join(self.signature)

    def get_qtable(self):
        catalogue = self.get_catalogue()
        t = QTable(
            names=list(self.signature_dypes.keys()),
            dtype=[v['dtype'] for v in self.signature_dypes.values()]
        )
        for row in catalogue:
            t.add_row(row)
        for field in self.signature_dypes:
            if not self.signature_dypes[field]['units'] is None:
                t[field] = t[field]*self.signature_dypes[field]['units']
        return t

    def prt_csv(self):
        catalogue = self.get_catalogue()
        print(self.get_signature())
        for row in catalogue:
            # NB: Old code for printing strings in quotes...
            #str = ""
            #for col in row:
            #    # TODO: If changing to phthon3, will need to replace basestring with str...
            #    #    Notes: https://stackoverflow.com/questions/4843173/how-to-check-if-type-of-a-variable-is-string
            #    str+="," if col is None else "{0},".format("'%s'" % col if isinstance(col,basestring) else col)
            #print re.sub(r",$","",str)
            print(re.sub(r"(\[|\])","","{0}".format(row)))

