import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.table import QTable
from astropy.units.quantity import Quantity as unit_type

# KEYS: ['ref_id', 'clump_id', 'subclump_id', 'match_id', 'cat_id', 'source_finder', 'source_finder_type', 'image_type', 'ra', 'dec', 'extent_semimajor', 'extent_semiminor', 'extent_angle', 'rms_noise_bane', 'sn_bane', 'rms_noise_bane_deep', 'sn_bane_deep', 'rms_noise_bane_shallow', 'sn_bane_shallow', 'residual_rms', 'residual_madfm', 'residual_sumsq', 'flux_peak', 'flux_peak_err', 'flux_total', 'flux_total_err']

rnmap = {
    'match_id': {
        'field': 'match_id',
        'title': ['Match', 'ID'],
        'align': 'c',
    },
    'source_finder': {
        'field': 'finder',
        'title': ['Source', 'Finder'],
        'align': 'l',
    },
    'image_type': {
        'field': 'depth',
        'title': ['Image', 'Depth'],
        'align': 'l',
    },
    'ra': {
        'field': 'ra',
        'title': ['RA', '($^\\circ$)'],
        #'title': ['RA ($^\\circ$)', ''],
        'align': 'r',
        'round_to': 3,
    },
    'dec': {
        'field': 'dec',
        'title': ['Dec', '($^\\circ$)'],
        #'title': ['Dec ($^\\circ$)', ''],
        'align': 'r',
        'round_to': 3,
    },
    'extent_semimajor': {
        'field': 'major',
        'mult_by': 2,
        'convert_to': u.arcsec,
        'title': ['$a$', '($^{\\prime\\prime}$)'],
        #'title': ['$a$ ($^{\\prime\\prime}$)', ''],
        'align': 'r',
        'round_to': 2,
    },
    'extent_semiminor': {
        'field': 'minor',
        'mult_by': 2,
        'convert_to': u.arcsec,
        'title': ['$b$', '($^{\\prime\\prime}$)'],
        #'title': ['$b$ ($^{\\prime\\prime}$)', ''],
        'align': 'r',
        'round_to': 2,
    },
    'extent_angle': {
        'field': 'pa',
        'title': ['$\\theta$', '($^\\circ$)'],
        #'title': ['$\\theta$ ($^\\circ$)', ''],
        'align': 'r',
        'round_to': 3,
    },
    'flux_total': {
        'field': 'flux',
        'title': ['Flux', '($m$Jy)'],
        #'title': ['Flux ($m$Jy)', ''],
        'align': 'r',
        'round_to': 4,
    },
    'sn_bane_shallow': {
        'field': 'sn',
        'title': ['S/N', ''],
        'align': 'r',
        'round_to': 4,
    },
    #'residual_rms': {
    #    'field': 'rms',
    #    'title': ['RMS', '$\\frac{\\displaystyle m\\mbox{Jy}}{\\displaystyle^{\\prime^2}\\mbox{beam}}$'],
    #    'align': 'c',
    #    'round_to': 7,
    #},
    'residual_madfm': {
        'field': 'madfm',
        #'title': ['MADFM', '$\\frac{\\displaystyle m\\mbox{Jy}}{\\displaystyle^{\\prime^2}\\mbox{beam}}$'],
        'title': ['MADFM', '$m\\mbox{Jy}/(\\prime^2\\mbox{beam})$'],
        'align': 'c',
        'round_to': 4,
        'is_exponent': None,
    },
    #'residual_sumsq': {
    #    'field': 'sumsq',
    #    'title': ['$\\Sigma I^2$', '$\\left(\\!\\frac{\\displaystyle m\\mbox{Jy}}{\\displaystyle^{\\prime}\\mbox{beam}}\\!\\right)^2$'],
    #    'align': 'c',
    #    'round_to': 3,
    #},
}

labels = {
    'aegean':   'Aegean',
    'caesar':   'Caesar',
    'profound': 'ProFound',
    'pybdsf':   'PyBDSF',
    'selavy':   'Selavy',
}

def bent():
    # init
    bent_tail_id = 2293
    compact_source_id = 2294
    root = "/Users/susy/cirada/emu_pipeline/docs/notes/emu_2x2"
    qt = QTable.read(f"{root}/xtables/emu_pilot_sample_2x2deg.hydra.cluster_catalogue.fits")
    qt.sort(['match_id','source_finder','image_type'])

    # extract the bent tail object
    qt = qt[(qt['clump_id']==bent_tail_id)|(qt['clump_id']==compact_source_id)]

    # filter on keys
    qt = qt[list(rnmap.keys())]

    # conversions
    for k in rnmap:
        if 'mult_by' in rnmap[k]:
            qt[k] *= 2
        if 'convert_to' in rnmap[k]:
            qt[k] = qt[k].to(rnmap[k]['convert_to'])
    qt['source_finder'] = [labels[s] for s in qt['source_finder']]
    qt['image_type'] = [depth.capitalize() for depth in qt['image_type']]

    # rename columns
    #qt.rename_columns(list(rnmap.keys()),[rnmap[k]['field'] for k in rnmap])
    fields =  [rnmap[k]['field'] for k in rnmap]
    qt.rename_columns(list(rnmap.keys()),fields)

    # print clump info
    print(qt)

    print("\n***\n")

    # create latex table
    def rnd(field,value):
        for key in rnmap:
            if field==rnmap[key]['field'] and 'round_to' in rnmap[key]:
                if not 'is_exponent' in rnmap[key]:
                    return f"{value.value if isinstance(value,unit_type) else value:0.{rnmap[key]['round_to']}f}"
                else:
                    return f"{value.value if isinstance(value,unit_type) else value:0.{rnmap[key]['round_to']}e}"
        return value if not isinstance(value,unit_type) else value.value
    tbs = list() 
    tbs.append("\\begin{table*}[hbt!]")
    tbs.append("\\caption{Bent tail example (\\texttt{clump\\_id} 2293).}")
    tbs.append("\\centering")
    #tbs.append("{\\footnotesize")
    tbs.append("\\begin{tabular}{@{\\;}%s@{\\;}}" % "@{\\;\\;\\;}".join([rnmap[k]['align'] for k in rnmap]))
    tbs.append("\\hline\\hline")
    for i in [0,1]:
        tbs.append("%s\\\\" % " & ".join(["\\multicolumn{1}{c}{"+f"{rnmap[k]['title'][i]}"+"}" for k in rnmap]))
    tbs.append("\\hline")
    for r in qt:
        tmp = list()
        for field in fields:
            tmp.append(f"{rnd(field,r[field])}")
        tbs.append("%s\\\\" % " & " .join(tmp))
    tbs.append("\\hline\\hline")
    tbs.append("\\end{tabular}")
    #tbs.append("}")
    tbs.append("\\label{tb:emu_2x2_clump_id_2293}")
    tbs.append("\\end{table*}")
    tbs.append("")

    print("\n".join(tbs))

    return qt









