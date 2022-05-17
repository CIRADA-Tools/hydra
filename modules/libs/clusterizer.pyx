# Notes:
# [1] Tutorial
#    https://www.quora.com/What-is-the-difference-between-Cython-and-CPython
#    http://docs.cython.org/en/latest/src/tutorial/cython_tutorial.html
#    https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html#memoryviews
#    https://cython.readthedocs.io/en/latest/src/tutorial/memory_allocation.html
#    Arrays: https://cython.readthedocs.io/en/latest/src/tutorial/array.html
# [3] A la, numpy...
#    https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html#numpy-tutorial
# [2] Re: Python.h errors...
#    https://stackoverflow.com/questions/21530577/fatal-error-python-h-no-such-file-or-directory
# [3] Trees...
#     https://www.analyticsvidhya.com/blog/2016/04/complete-tutorial-tree-based-modeling-scratch-in-python/
#     https://stackoverflow.com/questions/34964878/python-generate-a-dictionarytree-from-a-list-of-tuples
#     https://en.wikipedia.org/wiki/K-d_tree
# [4] Clustering...
#     HDBSCAN: https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html
# [5] sklearn...
#     https://scikit-learn.org/stable/modules/clustering.html
#         https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.DistanceMetric.html
#             https://en.wikipedia.org/wiki/Great-circle_distance
#             http://docs.astropy.org/en/stable/coordinates/matchsep.html
#     https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering
#     https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.connected_components.html
# [6] GIL
#     https://realpython.com/python-gil/
# [7] Vectorizing
#     https://stackoverflow.com/questions/45133276/passing-c-vector-to-numpy-through-cython-without-copying-and-taking-care-of-me

import sys
import numpy as np
#cimport numpy as np
from cython cimport view
from progressbar import ProgressBar, Percentage, Bar

cdef list flatten(clumps):
    cdef list clump 
    cdef list flattened
    if isinstance(clumps, np.ndarray):
        flattened = clumps.flatten().tolist()
    else:
        flattened = list()
        for clump in clumps:
            flattened += clump
    return flattened

def deltas(pairings,clusters):
    pairings_set = np.unique(pairings.flatten())
    clusters_set = flatten(clusters)
    reduced_clusters_set = np.unique(clusters_set)
    delta_pairings_rcluster = len(pairings_set) - len(reduced_clusters_set)
    delta_clusters_rcluster = len(clusters_set) - len(reduced_clusters_set)
    return delta_pairings_rcluster,delta_clusters_rcluster

#cdef int is_intersecting(list clump_a,list clump_b):
#    cdef long a_i
#    cdef long b_i
#    if not (clump_b[len(clump_b)-1] < clump_a[0]):
#        for a_i in reversed(clump_a):
#            for b_i in clump_b:
#              if a_i < b_i:
#                  break
#              elif a_i == b_i:
#                  return 1
#    return 0
cdef int is_intersecting(list clump_a,list clump_b):
    # notes: https://stackoverflow.com/questions/3170055/test-if-lists-share-any-items-in-python
    return not set(clump_a).isdisjoint(clump_b)
#cdef int is_intersecting(list clump_a, list clump_b):
#    cdef long max_a = len(clump_a)-1
#    cdef long max_b = len(clump_b)-1
#    if max_a < 0 or max_b < 0 or clump_a[max_a] < clump_b[0] or clump_b[max_b] < clump_a[0]:
#        return 0
#    cdef long i = max_a
#    cdef long j = max_b
#    while i > -1:
#        while j > -1:
#            if clump_a[i] < clump_b[j]:
#                j -= 1
#            elif clump_a[i] == clump_b[j]:
#                return 1
#            else:
#                j = max_b
#                break
#        i -= 1
#    return 0

cdef class connections:
    cdef long max_index
    cdef long[:] counts
 
    def __init__(self, list clumps, no_counts = False):
        cdef list flattened
        if len(clumps) > 0:
            flattened = self.__flatten(clumps)
            self.max_index = np.max(np.array(flattened,dtype=np.int64))+1
            self.counts = np.ascontiguousarray(np.zeros(self.max_index,dtype=np.int64))
            if not no_counts:
                self.update_counts(flattened)
        else:
            self.max_index = -1
            self.counts = np.array(list(),dtype=np.int64)
 
    cdef list __flatten(self, list clumps):
        cdef list clump 
        cdef list flattened = list()
        for clump in clumps:
            flattened += clump
        return flattened

    cdef int update_counts(self, list clump):
        cdef int success = 0
        cdef long component
        if np.max(np.array(clump,dtype=np.int64)) <= self.max_index:
            for component in clump:
                self.counts[component] += 1
            success = 1
        return success
 
    cdef int is_connected(self, list clump):
        cdef long component
        for component in clump:
            if self.counts[component] > 1:
                return 1
        return 0

    cdef int is_in(self, list clump):
        cdef long component
        for component in clump:
            if self.counts[component] > 0:
                return 1
        return 0

cdef long clipping_off = 1
cdef object pfpbar = None
cdef long progress = 0
cdef list isolated_clumps = list()
cdef list prune(list clumps):
    global isolated_clumps
    cdef len_init = len(isolated_clumps)
    cdef Py_ssize_t rows = len(clumps)
    cdef Py_ssize_t row = 0

    if clipping_off:
        sys.stderr.write("Pruning {:0,} rows...\n".format(rows))
    #print "# Pruning {:0,} rows...".format(rows)

    cn = connections(clumps)
    cdef list clump
    cdef list connected_clumps = list()
    global progress
    global pfpbar
    if clipping_off:
        pfpbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
    while row < rows:
        clump = clumps[row]
        if cn.is_connected(clump):
            connected_clumps.append(clump)
        else:
            isolated_clumps.append(clump)
            if not clipping_off:
                progress += 1
                pfpbar.update(progress)
        row += 1
        if clipping_off:
            pfpbar.update(row)
    if clipping_off:
        sys.stderr.write("\n")

    #sys.stderr.write("Isolated Clumps: {:0,} rows...\n".format(len(isolated_clumps)))
    #print "# Isolated Pairings: {:0,} rows...".format(len(isolated_clumps))

    #sys.stderr.write("Connected Clumps: {:0,} rows...\n".format(len(connected_clumps)))
    #print "# Connected Pairings: {:0,} rows...".format(len(connected_clumps))


    #sys.stderr.write("Sanity Check: {:0,}\n".format(rows+len_init-len(isolated_clumps)-len(connected_clumps)))
    #print "# Sanity Check: {:0,}\n".format(rows-(len(isolated_clumps)-len_init)-len(connected_clumps))

    return connected_clumps

cdef long get_index(list clump,long index, bc_periodic = True):
    cdef long i_max = len(clump)
    if bc_periodic:
        return index % i_max
    else:
        i_max -= 1
        if index > i_max:
            return i_max
    return index

cdef list get_index_component_pairs(list clumps,long clump_index):
    cdef list index_component_pairs = list()
    cdef long rows = len(clumps)
    cdef long row = 0
    while row < rows:
        index_component_pairs.append([row,clumps[row][get_index(clumps[row],clump_index)]])
        row += 1
    return index_component_pairs

cdef list sort_clumps_on_index(list clumps,long index):
    if len(clumps) == 0:
        return list()
    icps = np.array(get_index_component_pairs(clumps,index),dtype=np.int64)
    return [clumps[i] for i in icps[icps[:,1].argsort()][:,0]]

#cdef long get_max_clump_size(list clumps):
#    cdef long max_clump_size = 0
#    cdef long i = 0
#    cdef long i_max = len(clumps)
#    while i < i_max:
#        if max_clump_size < len(clumps[i]):
#            max_clump_size = len(clumps[i])
#        i += 1
#    return max_clump_size

#cdef long get_number_of_clumps_at_size(list clumps,long clump_size):
#    cdef long n_clumps = 0
#    cdef long i = 0
#    cdef long i_max = len(clumps)
#    while i < i_max:
#        if clump_size == len(clumps[i]):
#            n_clumps += 1
#        i += 1
#    return n_clumps

cdef list fold(list clumps,long on_index = 0):
    cdef list sorted_clumps = sort_clumps_on_index(clumps,on_index)
    cdef Py_ssize_t rows = len(clumps)
    cdef Py_ssize_t row = 0

    if clipping_off:
        sys.stderr.write("Folding {:0,} rows...\n".format(rows))
    #print "# Folding: {:0,}".format(rows)

    cdef list clump
    cdef list folded_clumps = list()
    global progress
    global pfpbar
    if clipping_off:
        pfpbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
    while row < rows:
        clump = list(sorted_clumps[row])
        while row+1 < rows:
            if is_intersecting(clump,sorted_clumps[row+1]):
                clump = np.unique(clump+sorted_clumps[row+1]).tolist()
                if not clipping_off:
                    progress += 1
                    pfpbar.update(progress)
            else:
                break
            row += 1
        folded_clumps.append(clump)
        row += 1 if row < rows else 0
        if clipping_off:
            pfpbar.update(row)
    if clipping_off:
        sys.stderr.write("\n")

    return folded_clumps

cdef list clip(list hedges, long init = 1):
    global isolated_clumps
    if init:
        isolated_clumps = list()

    cdef long initial_isolated_clumps = len(isolated_clumps)
    def get_outbox_size():
        return len(isolated_clumps)-initial_isolated_clumps

    cdef long rows = len(hedges)
    sys.stderr.write("Clipping: {:0,} rows...\n".format(rows))
    #print "# Clipping: {0} rows...".format(rows)

    global pfpbar
    pfpbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
    global progress
    progress = 0
    global clipping_off
    clipping_off = 0

    cdef long idx = 0
    cdef long prune_level = 30001
    cdef long delta_prune = 1
    cdef long max_idx = 100
    cdef list clumps = prune(hedges)
    cdef long outbox_size = get_outbox_size()
    pfpbar.update(get_outbox_size())
    rows = len(clumps)
    while rows > prune_level and idx < max_idx:
        clumps = fold(clumps,idx)
        pfpbar.update(get_outbox_size()+rows-len(clumps))
        clumps = prune(clumps)
        new_outbox_size = get_outbox_size()
        pfpbar.update(new_outbox_size)
        if new_outbox_size - outbox_size < delta_prune:
            break
        else:
            outbox_size = new_outbox_size
        rows = len(clumps)
        idx += 1
    clumps = sort_clumps_on_index(prune(fold(prune(fold(clumps,0)),1)),0)
    pfpbar.update(get_outbox_size())
    sys.stderr.write("\n")
    ##print "# idx_max:",idx
    ##print "# rows:",rows,(">" if rows > prune_level else "<="),prune_level
    ##print "# delta_prune:",delta_prune

    clipping_off = 1
    return clumps

cdef list merge(list clumps):
    cdef Py_ssize_t rows = len(clumps)
    cdef Py_ssize_t row = 0

    sys.stderr.write("Merging {:0,} rows...\n".format(rows))
    #print "# Merging: {:0,}".format(rows)

    nc = connections(clumps,True)
    cdef long i
    cdef list merged_clumps = list()
    cdef list clump
    cdef int is_merged
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
    while row < rows:
        is_merged = 0
        clump = list(clumps[row])
        i = len(merged_clumps) - 1
        if nc.is_in(clump):
            while i > -1:
                if is_intersecting(clump,merged_clumps[i]):
                    merged_clumps[i] = np.unique(clump+merged_clumps[i]).tolist()
                    is_merged = 1
                    break
                i -= 1
        if is_merged == 0:
            merged_clumps.append(clump)
        nc.update_counts(clump)
        row += 1 if row < rows else 0
        pbar.update(row)
    pbar.update(rows)
    sys.stderr.write("\n")

    return merged_clumps 

def clusterize(pairings):
    global isolated_clumps
    isolated_clumps

    #cdef list clumps = sort_clumps_on_index(prune(merge(clip(pairings.tolist()))),0)
    cdef list clumps = clip(pairings.tolist())
    cdef long iters = 10
    cdef long iter = 0
    cdef long clump_size = len(clumps)
    cdef long new_clump_size
    while iter < iters and clump_size > 0:
        clumps = clip(merge(clumps),0)
        new_clump_size = len(clumps)
        if clump_size == new_clump_size:
            break
        clump_size = new_clump_size
        iter += 1
    clumps = sort_clumps_on_index(clumps,0)

    cdef Py_ssize_t rows = len(clumps)
    cdef Py_ssize_t row = 0

    sys.stderr.write("Clustering {:0,} rows...\n".format(rows))

    nc = connections(clumps,True)
    cdef long i, i_max
    cdef list clusters = list()
    cdef list cluster
    if rows > 0:
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=rows).start()
        while row < rows:
            cluster = clumps[row]
            if nc.is_in(cluster):
                while row+1 < rows: # pre-merge loop
                    if is_intersecting(clumps[row+1],cluster):
                        cluster = np.unique(cluster+clumps[row+1]).tolist()
                        row += 1
                    else:
                        break
                i = 0
                i_max = len(clusters)
                while i < i_max: # cluster loop
                    if is_intersecting(cluster,clusters[i]):
                        cluster = np.unique(cluster+clusters[i]).tolist()
                        clusters.pop(i)
                        i_max -= 1
                    else:
                        i += 1
            clusters.append(cluster)
            nc.update_counts(cluster)
            row += 1
            pbar.update(row)
            pbar.update(rows)
    sys.stderr.write("\n")

    rows = len(clusters)
    sys.stderr.write("Clusters > 2: {:0,}\n".format(rows))
    #print "# Clusters: {:0,}".format(rows)
    clusters.extend(isolated_clumps)

    rows = len(clusters)
    sys.stderr.write("Total Clusters: {:0,}\n".format(rows))
    #print "# Clusters: {:0,}".format(rows)

    delta_pc, delta_cc = deltas(pairings,clusters)
 
    sys.stderr.write("Delta: {:0,}\n".format(delta_pc))
    #print "# Delta: {:0,}".format(delta_pc)

    sys.stderr.write("Missing: {:0,}\n".format(delta_cc))
    #print "# Missing: {:0,}".format(delta_cc)

    return sort_clumps_on_index(clusters,0)
