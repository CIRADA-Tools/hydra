# TO-DO:
# [1] Handle empty catalogue case, re., metrics table
# [2] use multiprocessing.cpu_count() (and memory?
#     ... re., container sizes) to paralize, Hydra,
#     Typhon, and Cerberus.
# [3] See rms_noise_bane is expressed as mJy and not
#     mJy/beam in the deep PyBDSF cataloge (re., 2x2
#     simulated) -- fix...
import re
import gc
import os
import sys
import time
import glob
import aplpy
import shutil
import tarfile
sys.path.insert(0,"modules")
from libs import Typhon
from libs import homados
from libs import print_ok
from libs import print_warning
from libs import HydraIOError
from libs import Config
from libs import Cluster 
from libs import HydraRenderingError 
from libs import load_bane_rms_image
from libs import render
import click
import yaml as yml
yml.warnings({'YAMLLoadWarning': False})
from libs import FITS 
from astropy.io import fits
from astropy.table.table import QTable as QTableClass
from astropy import units as u
from astropy.table import Table
from astropy.table import QTable
from astropy.table import Column
from astropy.table import MaskedColumn
from astropy.table import hstack
from astropy.table import vstack
from astropy.coordinates import SkyCoord
from astropy.table.column import Column as ColumnType
from astropy.nddata.utils import Cutout2D
import numpy as np
import pandas as pd

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# notes: https://docs.python.org/3/library/multiprocessing.html
import multiprocessing

# define cache directory
def get_this_source_file_directory():
    def path(fullpath):
        return re.sub(r"(.*/).*$",r"\1",fullpath)
    return path(os.path.realpath(__file__))
escritoire = os.path.abspath(get_this_source_file_directory()+"/.escritoire")
cache = os.path.abspath(f"{escritoire}/data/hydra")

# define cache directory
def get_this_source_file_directory():
    def path(fullpath):
        return re.sub(r"(.*/).*$",r"\1",fullpath)
    return path(os.path.realpath(__file__))
escritoire = os.path.abspath(get_this_source_file_directory()+"/.escritoire")
cache = os.path.abspath(f"{escritoire}/data/hydra")

class HydraCroppingTool(FITS):
    def __init__(self,fits_image_file):
        super(HydraCroppingTool,self).__init__(fits_image_file)

    def __get_dim(self,header):
        dim = 0
        for i in range(1,5):
            if f"CTYPE{i}" in header:
                dim += 1
        #print_warning(f" ===> DIM: {dim}")
        return dim

    def get_residual_stats(self,ra,dec,size):
        position = SkyCoord(ra,dec)
        stamp = Cutout2D(self.get_data(),position=position,size=size,wcs=self.get_wcs(),mode='partial',copy=True)
        img = FITS(fits.PrimaryHDU(stamp.data, header=stamp.wcs.to_header()))
        return {'rms': img.get_rms(), 'madfm': img.get_madfm(), 'sumsq': img.get_sumsq()}

    def crop_image(self,ra,dec,size,output_file=None):
        stamp = Cutout2D(self.get_data(),position=SkyCoord(ra,dec),size=size,wcs=self.get_wcs(),mode='partial',copy=True)
        trimmed = fits.PrimaryHDU(stamp.data, header=stamp.wcs.to_header())
        if not output_file is None:
            trimmed.writeto(output_file,overwrite=True)
        return trimmed

    def create_residual_image(self,catalogue,residual_file,model_file):
        print(f"> Creating Residual Image...")

        # get some header info
        wcs = self.get_wcs()
        CDELT1 = np.abs(wcs.wcs.cdelt[0])
        CDELT2 = np.abs(wcs.wcs.cdelt[1])

        # build our model image
        scale = 8
        image = self.get_data()
        n_x = image.shape[1]-1
        n_y = image.shape[0]-1
        model = np.zeros(image.shape)
        def make_gaussian(ra_0,dec_0,smaj,smin,pa,peak):
            a = ((np.cos(pa)/smin)**2+(np.sin(pa)/smaj)**2)/2.0
            b = (np.sin(2.0*pa)/(smin)**2-np.sin(2.0*pa)/(smaj)**2)/2.0
            c = ((np.sin(pa)/smin)**2+(np.cos(pa)/smaj)**2)/2.0
            def gaussian(ra,dec):
                return peak*np.exp(-a*(ra-ra_0)**2-b*(ra-ra_0)*(dec-dec_0)-c*(dec-dec_0)**2)
            return gaussian
        for row in catalogue:
            # compute (ra,dec) in pixel units
            ra  = row['ra'].to(u.deg).value
            dec = row['dec'].to(u.deg).value
            r_pix = wcs.all_world2pix([[ra,dec]],0,ra_dec_order=True)[0]
            ra_0  = r_pix[0]
            dec_0 = r_pix[1]

            # compute semi-major/minor axes in pixel units
            pa = row['extent_angle'].to(u.rad).value
            a  = 2.0*row['extent_semimajor'].to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
            b  = 2.0*row['extent_semiminor'].to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
            smaj = a*np.sqrt((np.sin(pa)/CDELT1)**2+(np.cos(pa)/CDELT2)**2)
            smin = b*np.sqrt((np.cos(pa)/CDELT1)**2+(np.sin(pa)/CDELT2)**2)

            # get peak flux
            pf = row['flux_peak'].to(u.Jy/u.beam).value
            
            # create the gaussian function
            gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,pf)

            # make residuals
            ra_min  = max(0,int(np.floor(ra_0-scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
            ra_max  = min(n_x,int(np.ceil(ra_0+scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
            dec_min = max(0,int(np.floor(dec_0-scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
            dec_max = min(n_y,int(np.ceil(dec_0+scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
            i,j = np.mgrid[ra_min:ra_max+1,dec_min:dec_max+1]
            model[j,i] += gaussian(i,j)

        # create residual image
        image -= model

        ## output model image
        print(f">> Writing Model: {model_file}")
        fits.PrimaryHDU(model,header=self.get_header()).writeto(model_file,overwrite=True)

        ## output residual image
        print(f">> Writing Residaul: {residual_file}")
        fits.PrimaryHDU(image,header=self.get_header()).writeto(residual_file,overwrite=True)
        print(f"> [Done]")


# notes: syntax sync fromstart => crtl-l
class Hydra:
    def __init__(self,fits_image,use=None,catalogues_yml=None,use_archive=True):
        # define tar file
        self.tar_gz_file   = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".hydra.tar.gz",fits_image)

        # define deep/shallow info
        self.image_types  = ['deep','shallow']
        self.deep_image    = fits_image
        self.shallow_image = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".shallow.fits",self.deep_image)
        self.max_limit = 3.0*u.arcmin

        # we need to apply an agressive percent real deteictions cut, as we are dealing with 
        # low sample statics for out optimization process: i.e., C. Hale, et al., uses 98%, 
        # we're gonna do 90%!
        self.percent_real_detections_cut = 90.0
        #self.percent_real_detections_cut = 95.0

        # define clustering parameters
        #self.skirt_factor = 1.5
        self.skirt_factor = 1.0

        # define gc 
        self.default_bins = 50
        #self.default_bins = 100

        # paper plot limits
        #self.is_paper = True
        self.is_paper = False
        self.completeness_and_reliability_sn_plot_limits = (1,10**3)

        # ok, processing time
        if use_archive and os.path.isfile(self.tar_gz_file): # load archive
            def strip(fname):
                return re.sub(r"\.[Ff][Ii][Tt](|[Ss])$","",re.sub(r"^(.*?/)+","",fname))

            # define tar files
            self.deep_tar_gz_file    = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".typhon.tar.gz",self.deep_image)
            self.shallow_tar_gz_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".typhon.tar.gz",self.shallow_image) 
            archive_dir = re.sub("^(.*?/)+","",re.sub("\.tar\.gz$","_dir",self.tar_gz_file))

            # set location of noise files in archive, don't load
            self.rms_imgs = {
                'deep':    f"{archive_dir}/noise/deep/{strip(fits_image)}.hydra.bane.deep.noise.fits", 
                'shallow': f"{archive_dir}/noise/shallow/{strip(fits_image)}.hydra.bane.shallow.noise.fits" 
            }

            # load catalogue info
            self.catalogues = { # template
                'metrics': {
                    'rms_box_statistics': {
                        'deep': Table(),
                        'shallow': Table(),
                    },
                    'statistics': {
                        'deep': Table(),
                        'shallow': Table(),
                    },
                    'parameters': QTable(),
                },
                'modules': dict(),
                'cluster': QTable(),
                'clump': QTable(),
            }
            with tarfile.open(self.tar_gz_file,'r:gz') as fd:
                print_ok(f"LOADING ARCHIVE: {self.tar_gz_file}")
                def get_module_name(fname):
                    return re.sub(r"^.*?(\.hydra)\.(\w+?)\..*$",r"\2",fname)
                re_cluster                    = re.compile("^"+re.escape(f"{archive_dir}/{strip(fits_image)}.hydra.cluster_catalogue.fits"))
                re_clump                      = re.compile("^"+re.escape(f"{archive_dir}/{strip(fits_image)}.hydra.clump_catalogue.fits"))
                re_deep_rms_box_statistics    = re.compile("^"+re.escape(f"{archive_dir}/statistics/{strip(fits_image)}.hydra.deep_rms_box_statistics.fits"))
                re_shallow_rms_box_statistics = re.compile("^"+re.escape(f"{archive_dir}/statistics/{strip(fits_image)}.hydra.shallow_rms_box_statistics.fits"))
                re_deep_statistics            = re.compile("^"+re.escape(f"{archive_dir}/statistics/{strip(fits_image)}.hydra.deep_statistics.fits"))
                re_shallow_statistics         = re.compile("^"+re.escape(f"{archive_dir}/statistics/{strip(fits_image)}.hydra.shallow_statistics.fits"))
                re_prameters                  = re.compile("^"+re.escape(f"{archive_dir}/statistics/{strip(fits_image)}.hydra.global_metrics.fits"))
                re_extern                     = re.compile("^"+re.escape(f"{archive_dir}/injected/{strip(fits_image)}.hydra.")+"\w+?"+re.escape(".external_catalogue.fits"))
                re_cats_dir                   = re.compile("^%s/catalogues/" % re.escape(archive_dir))
                re_cats                       = {
                    'deep':    re.compile("^%s/catalogues/deep/" % re.escape(archive_dir)),
                    'shallow': re.compile("^%s/catalogues/shallow/" % re.escape(archive_dir))
                }
                re_gc                 = re.compile("^"+re.escape(f"{archive_dir}/plots/graphics_context.yml"))
                for member in fd.getmembers():
                    if re_cats_dir.search(member.name):
                        mname = get_module_name(member.name)
                        if not mname in self.catalogues['modules']:
                            self.catalogues['modules'][mname] = {
                                'file_names': {'deep': "", 'shallow': ""},
                                'tables': {'deep': QTable(), 'shallow': QTable()},
                            }
                        for depth in re_cats:
                            if re_cats[depth].search(member.name):
                                self.catalogues['modules'][mname]['file_names'][depth] = member.name
                                self.catalogues['modules'][mname]['tables'][depth] = Table.read(fd.extractfile(member),format='fits')
                    elif re_cluster.search(member.name):
                        self.catalogues['cluster'] = QTable(Table.read(fd.extractfile(member),format='fits'))
                    elif re_clump.search(member.name):
                        self.catalogues['clump'] = QTable(Table.read(fd.extractfile(member),format='fits'))
                    elif re_prameters.search(member.name):
                        self.catalogues['metrics']['parameters'] = QTable(Table.read(fd.extractfile(member),format='fits'))
                    elif re_deep_rms_box_statistics.search(member.name):
                        self.catalogues['metrics']['rms_box_statistics']['deep'] = QTable.read(fd.extractfile(member),format='fits')
                    elif re_shallow_rms_box_statistics.search(member.name):
                        self.catalogues['metrics']['rms_box_statistics']['shallow'] = QTable.read(fd.extractfile(member),format='fits')
                    elif re_deep_statistics.search(member.name):
                        self.catalogues['metrics']['statistics']['deep'] = Table.read(fd.extractfile(member),format='fits')
                    elif re_shallow_statistics.search(member.name):
                        self.catalogues['metrics']['statistics']['shallow'] = Table.read(fd.extractfile(member),format='fits')
                    elif re_extern.search(member.name):
                        cname = re.sub(r"^.*?\.hydra\.(\w+?)\.external_catalogue\.fits",r"\1",member.name)
                        if not 'extern' in self.catalogues:
                            self.catalogues['extern'] = dict()
                        self.catalogues['extern'][cname] = {'table': QTable(Table.read(fd.extractfile(member),format='fits'))}
                    elif re_gc.search(member.name):
                        self.graphics_context = yml.load(fd.extractfile(member).read())

            # let's do some house cleaning...
            if 'extern' in self.catalogues:
                # we need to get the proper case-sensitive names (i.e., dict-keys) for our external catalogues...
                pnames = [list(m.keys())[0] for m in self.graphics_context['plot']['labels']]
                for cname in self.catalogues['extern']:
                    for pname in pnames:
                        if cname.lower()==pname.lower():
                            self.catalogues['extern'][pname] = self.catalogues['extern'].pop(cname)
                            break

            # create source finder list
            self.source_finders = np.unique(self.catalogues['cluster']['source_finder'])
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    self.source_finders = self.source_finders[(self.source_finders!=catalogue)]

            # set plot info
            self.plot_colors = {list(m.keys())[0]:m[list(m.keys())[0]] for m in self.graphics_context['plot']['colors']}
            self.plot_labels = {list(m.keys())[0]:m[list(m.keys())[0]] for m in self.graphics_context['plot']['labels']}
            self.plot_collection = None
            #self.plot_collection = self.__collect_plot_info() # dev/debug (keep)
            self.residuals = None
            self.cutout_files = None
            
            ## TO-DO: temp
            #self.print_clump_table()
            #clump_sizes = list()
            #for source_finder in self.source_finders:
            #    clump_sizes.append(self.catalogues['clump'][f"n_{source_finder}"])
            #clump_sizes = list(np.unique(np.array(clump_sizes)[(np.array(clump_sizes)>0)]))
            #print_warning(f"clump_sizes: {clump_sizes}")
            #for clump_size in clump_sizes:
            #    clump_ids = list()
            #    for source_finder in self.source_finders:
            #        clump_ids.extend(list(self.catalogues['clump'][(self.catalogues['clump'][f"n_{source_finder}"]==clump_size)]['clump_id']))
            #    #print_warning(f"{clump_size} => {clump_ids}")
            #    clump_ids = list(np.unique(clump_ids))
            #    print_ok(f"{clump_size} => {clump_ids}")
            #exit()

        else: # create archive
            print_ok(f"CREATING ARCHIVE: {self.tar_gz_file}")
            # make sure we have the fits image file
            if not os.path.isfile(self.deep_image):
                print_warning(f"ERROR: File '{fits_image}' does not exist!")
                print_ok("Bye!")
                raise HydraIOError

            # make sure we have a shallow image; create, if doesn't exist
            if not os.path.isfile(self.shallow_image):
                print_ok(f"Creating Shallow image file: {self.shallow_image}")
                homados(["shallow",self.deep_image,"--no-noise","--verbose"],standalone_mode=False)
                print_ok("[Done]")

            # rms_box optimization flag
            is_optimize_rms_box = True

            # diagnostics flag
            is_diagnostics = False

            # load deep catalogues from deep_image.typhon.tar.gz archive file; create, if doesn't exist
            Deep_Archive = Typhon(
                fits_image_file     = self.deep_image,
                processing_cache    = cache,
                use                 = use,
                is_diagnostics      = is_diagnostics,
                is_fits_catalogue   = True,
                is_residual         = True,
                is_optimize_rms_box = is_optimize_rms_box,
                percent_real_detections_cut = self.percent_real_detections_cut
            )
            self.deep_tar_gz_file = Deep_Archive.get_archive_name()
            print_ok(f"> Loading: {self.deep_tar_gz_file}")
            for catalogue in Deep_Archive.get_catalogue_names():
                print_ok(f"> o {catalogue}")
            t_deep = Deep_Archive.load_catalogues()

            # load shallow catalogues from shallow_image.typhon.tar.gz archive file; create, if doesn't exist
            Shallow_Archive = Typhon(
                fits_image_file     = self.shallow_image,
                processing_cache    = cache,
                use                 = use,
                is_diagnostics      = is_diagnostics,
                is_fits_catalogue   = True,
                is_residual         = True,
                is_optimize_rms_box = is_optimize_rms_box,
                percent_real_detections_cut = self.percent_real_detections_cut
            )
            self.shallow_tar_gz_file = Shallow_Archive.get_archive_name() 
            print_ok(f"> Loading: {self.shallow_tar_gz_file}")
            for catalogue in Shallow_Archive.get_catalogue_names():
                print_ok(f"> o {catalogue}")
            t_shallow = Shallow_Archive.load_catalogues()

            ## debug
            #print_ok("***<TEST I>***")
            #print_warning("DEEP TEST")
            #qd = t_deep['aegean']['astropy_table']
            #print(qd[(321.14<=qd['ra'])&(qd['ra']<=321.15)&(-55.53<=qd['dec'])&(qd['dec']<=-55.52)]['ra','dec','a','b','pa'])
            #print_warning("SHALLOW TEST")
            #qs = t_shallow['aegean']['astropy_table']
            #print(qs[(321.14<=qs['ra'])&(qs['ra']<=321.15)&(-55.53<=qs['dec'])&(qs['dec']<=-55.52)]['ra','dec','a','b','pa'])
            ##print_warning("BYE!")
            #print_warning("[DONE]")
            #print_ok("***")
            ##exit()


            # get configuration
            config = Config()

            # initialize internal catalogues dict
            self.catalogues = dict()

            # create metrics entry
            self.catalogues['metrics'] = dict()

            # add typhon run rms box statistics
            self.catalogues['metrics']['rms_box_statistics'] = {
                'deep': Deep_Archive.load_rms_box_optimization_table(),
                'shallow': Shallow_Archive.load_rms_box_optimization_table()
            }

            # add typhon run statistics
            self.catalogues['metrics']['statistics'] = {
                'deep': Deep_Archive.load_optimization_table(),
                'shallow': Shallow_Archive.load_optimization_table()
            }

            # add typhon source-finder parameter optimizations
            self.catalogues['metrics']['optimizations'] = {
                'deep': Deep_Archive.get_optimizations(),
                'shallow': Shallow_Archive.get_optimizations()
            }

            # run bane to get deep and shallow rms images
            if is_optimize_rms_box:
                print_ok(f"Getting Typhon RMS box optimization pars...")
                # Get optimization pars, using try-except to save a lot of coding
                # nb: This only works because all source-finders uses a bane optimization,
                #     for their box_size and step_size input parameters (if supported).
                #     Also Aegean is assumed have these options, it works together with bane.
                try:
                    rms_box = {
                        'deep': {
                            'boxsize': self.catalogues['metrics']['optimizations']['deep']['aegean']['boxsize'],
                            'stepsize': self.catalogues['metrics']['optimizations']['deep']['aegean']['stepsize'],
                        },
                        'shallow': {
                            'boxsize': self.catalogues['metrics']['optimizations']['shallow']['aegean']['boxsize'],
                            'stepsize': self.catalogues['metrics']['optimizations']['shallow']['aegean']['stepsize'],
                        },
                    }
                except Exception as e:
                    rms_box = {'deep': None, 'shallow': None}
                    print_warning(f"> Whoops!")
                    print_warning(f">> Can't determine RMS box pars... this shouldn't happen!")
                    print_warning(f">> ERROR: {e}")
                    print_warning(f">> Setting: rms_box = {rms_box}")
                    print_warning(f">> We'll brute force it!")
                print_ok(f"> rms_box: {rms_box}")
                print_ok(f"[Done]")
            else:
                rms_box = {'deep': None, 'shallow': None}
            self.rms_imgs = {
                'deep':    load_bane_rms_image(self.deep_image,cache,rms_box['deep'],is_optimize_rms_box),
                'shallow': load_bane_rms_image(self.shallow_image,cache,rms_box['shallow'],is_optimize_rms_box)
            }

            # add typhon run parameters  metrics table
            def fetch_residual_image(module,depth):
                image = None
                work_dir = get_this_source_file_directory()
                def flush_cache(member):
                    os.chdir(cache)
                    subdir = re.sub(r"^(.*?/).*$",r"\1",member.name)
                    if os.path.isdir(subdir):
                        shutil.rmtree(subdir)
                    os.chdir(work_dir)
                archive = self.shallow_tar_gz_file if depth == 'shallow' else self.deep_tar_gz_file
                res = re.compile(r"(|\.shallow)\.%s\.residual\.fits$" % module)
                with tarfile.open(archive,"r") as fd:
                    for member in fd.getmembers():
                        if res.search(member.name):
                            fname = f"{cache}/{member.name}"
                            fd.extract(member,cache)
                            #image = FITS(fname)
                            image = HydraCroppingTool(fname)
                            flush_cache(member)
                return image
            def filter_pars(pars):
                def extract_params(ds_pars,module,params):
                    new_pars = dict()
                    for pr in ds_pars[module]:
                        if pr in params:
                            new_pars[pr] = ds_pars[module][pr]
                    ds_pars[module] = new_pars
                    return ds_pars
                variables = config.get_variables()
                for module in variables:
                    params = list()
                    for datum in variables[module]:
                        for param in datum:
                            if 'class' in datum[param] and datum[param]['class'] in ['rms_parameter','island_parameter']:
                                params.append(param)
                    pars = extract_params(pars,module,params)
                return pars
            qt = QTable()
            d_pars = filter_pars(Deep_Archive.optimize(use_existing=True))
            s_pars = filter_pars(Shallow_Archive.optimize(use_existing=True))
            rms_pars = [list(v['parameters'].keys())[0] if isinstance(list(v['parameters'].values())[0],str) else list(v['parameters'].keys())[1] for v in config.get_constraints().values()]
            def tr_name(pars,is_rms):
                keys = list(pars.keys())
                return keys[0 if is_rms else 1] if keys[0] in rms_pars else keys[1 if is_rms else 0]
            def tr_value(pars,is_rms):
                return pars[tr_name(pars,is_rms)]
            qt.add_column(
                name = 'source_finder',
                col  = list(d_pars.keys())+list(s_pars.keys())
            )
            qt.add_column(
                name = 'image_type',
                col = ['deep' for _ in d_pars]+['shallow' for _ in s_pars]
            )
            qt.add_column(
                name = 'rms_parameter_name',
                col  = [tr_name(d_pars[s],True) for s in d_pars]+[tr_name(s_pars[s],True) for s in s_pars]
            )
            qt.add_column(
                name = 'rms_parameter_value',
                col  = [tr_value(d_pars[s],True) for s in d_pars]+[tr_value(s_pars[s],True) for s in s_pars]
            )
            qt.add_column(
                name = 'isl_parameter_name',
                col  = [tr_name(d_pars[s],False) for s in d_pars]+[tr_name(s_pars[s],False) for s in s_pars]
            )
            qt.add_column(
                name = 'isl_parameter_value',
                col  = [tr_value(d_pars[s],False) for s in d_pars]+[tr_value(s_pars[s],False) for s in s_pars]
            )
            self.residuals = dict()
            for module in np.unique(qt['source_finder']):
                self.residuals[module] = {
                    'deep': fetch_residual_image(module,'deep'),
                    'shallow': fetch_residual_image(module,'shallow'),
                }
            ## TO-DO: --use flag bug. We probably need to vett the --use flag, as it's not been used for a while. That said, it seems to work for Typhon.
            ##
            ## The following print statements:
            ##
            #print_warning(f"self.residuals: {self.residuals}")
            #for t in qt:
            #    print_ok(f"t['source_finder']: {t['source_finder']}")
            #    print_ok(f"t['image_type']: {t['image_type']}")
            #    print_ok(f"self.residuals[{t['source_finder']}][{t['image_type']}]: {self.residuals[t['source_finder']][t['image_type']]}")
            #    print_ok(f"self.residuals[{t['source_finder']}][{t['image_type']}].get_rms(): {self.residuals[t['source_finder']][t['image_type']].get_rms()}")
            #    print_warning("***")
            #print_ok("OK")
            ##
            ## Run Produces:
            ##
            ##   (base) Continuum:emu_pipeline susy$ python hydra.py data/test/vlass/J085542+112459_s3arcmin_VLASS.fits --use data/test/vlass/J165209+212528_s3arcmin_VLASS.typhon.tar.gz --bypass-archived
            ##   ...
            ##   self.residuals: {'aegean': {'deep': None, 'shallow': <__main__.HydraCroppingTool object at 0x7fd0fdfcc990>}, 'caesar': {'deep': None, 'shallow': <__main__.HydraCroppingTool object at 0x7fd0fd4bcbd0>}, 'profound': {'deep': None, 'shallow': <__main__.HydraCroppingTool object at 0x7fd0fe0ad750>}, 'pybdsf': {'deep': None, 'shallow': <__main__.HydraCroppingTool object at 0x7fd0fdbf8990>}, 'selavy': {'deep': None, 'shallow': <__main__.HydraCroppingTool object at 0x7fd0fe0a5d50>}}
            ##   t['source_finder']: aegean
            ##   t['image_type']: deep
            ##   self.residuals[aegean][deep]: None
            ##   Traceback (most recent call last):
            ##     File "hydra.py", line 3895, in <module>
            ##       hydra()
            ##     File "/Users/susy/opt/anaconda3/lib/python3.7/site-packages/click/core.py", line 829, in __call__
            ##       return self.main(*args, **kwargs)
            ##     File "/Users/susy/opt/anaconda3/lib/python3.7/site-packages/click/core.py", line 782, in main
            ##       rv = self.invoke(ctx)
            ##     File "/Users/susy/opt/anaconda3/lib/python3.7/site-packages/click/core.py", line 1066, in invoke
            ##       return ctx.invoke(self.callback, **ctx.params)
            ##     File "/Users/susy/opt/anaconda3/lib/python3.7/site-packages/click/core.py", line 610, in invoke
            ##       return callback(*args, **kwargs)
            ##     File "hydra.py", line 3878, in hydra
            ##       hydra = Hydra(fits_file,use,catalogues_yml,not bypass_archived) \
            ##     File "hydra.py", line 420, in __init__
            ##       print_ok(f"self.residuals[{t['source_finder']}][{t['image_type']}].get_rms(): {self.residuals[t['source_finder']][t['image_type']].get_rms()}")
            ##   AttributeError: 'NoneType' object has no attribute 'get_rms'
            ##   (base) Continuum:emu_pipeline susy$ 
            ##
            ##  In short, fetch_residual_image(module,'deep') above is returning None -- INVESTIGATE!
            ##
            qt.add_column(
                name = 'residual_rms',
                col  = [self.residuals[t['source_finder']][t['image_type']].get_rms() for t in qt] 
            )
            qt.add_column(
                name = 'residual_madfm',
                col  = [self.residuals[t['source_finder']][t['image_type']].get_madfm() for t in qt] 
            )
            qt.add_column(
                name = 'residual_sumsq',
                col  = [self.residuals[t['source_finder']][t['image_type']].get_sumsq() for t in qt] 
            )
            print_ok(qt)
            self.catalogues['metrics']['parameters'] = qt

            # conglomerate deep and shallow catalogue information
            def translate(table,meta):
                cmap = dict()
                for datum in meta:
                    hname = list(datum.keys())[0]
                    fname = datum[hname]['header']
                    if hname != fname:
                        cmap[fname] = hname
                table.rename_columns(list(cmap.keys()),list(cmap.values()))
                table.add_column([r for r in range(1,len(table)+1)],name='id',index=0)
                return table
            meta = config.get_catalogue_meta()
            self.catalogues['modules'] = dict()
            for module in meta:
                self.catalogues['modules'][module] = {
                    'file_names': {
                        'deep':    t_deep[module]['catalogue'],
                        'shallow': t_shallow[module]['catalogue'],
                    },
                    'tables': {
                        'deep':    translate(t_deep[module]['astropy_table'],meta[module]),
                        'shallow': translate(t_shallow[module]['astropy_table'],meta[module])
                    }
                }
            # free up some memory
            del t_deep
            del t_shallow

            ## debug
            #print_warning("***<TEST II>***")
            #print_ok("DEEP TEST")
            #qd = self.catalogues['modules']['aegean']['tables']['deep']
            #print(qd.colnames)
            #print(qd[(321.14<=qd['ra'])&(qd['ra']<=321.15)&(-55.53<=qd['dec'])&(qd['dec']<=-55.52)]['id','ra','dec','semimajor','semiminor','pa'])
            #print_ok("SHALLOW TEST")
            #qs = self.catalogues['modules']['aegean']['tables']['shallow']
            #print(qs[(321.14<=qs['ra'])&(qs['ra']<=321.15)&(-55.53<=qs['dec'])&(qs['dec']<=-55.52)]['id','ra','dec','semimajor','semiminor','pa'])
            #print_ok("[DONE]")
            #print_warning("***")
            ##exit()

            # add catalogue rows to metrics
            self.catalogues['metrics']['parameters'].add_column(
                index = 6,
                name = 'components',
                col = [len(self.catalogues['modules'][t['source_finder']]['tables'][t['image_type']]) for t in self.catalogues['metrics']['parameters']]
            )

            # rms noise (a la bane) extraction tool
            def create_noise_column(qt,image_depth=None,is_suffix=False):
                wcs  = {
                    'deep':    self.rms_imgs['deep'].get_wcs(),
                    'shallow': self.rms_imgs['shallow'].get_wcs()
                }
                data = {
                    'deep':    self.rms_imgs['deep'].get_data(),
                    'shallow': self.rms_imgs['shallow'].get_data()
                }
                hdr =  {
                    'deep':    self.rms_imgs['deep'].get_header(),
                    'shallow': self.rms_imgs['shallow'].get_header()
                }
                noise = list()
                def get_noise(pos,depth):
                    r_pix = wcs[depth].all_world2pix([[pos.ra.to(u.deg).value,pos.dec.to(u.deg).value]],0,ra_dec_order=True)[0]
                    x_pix = min(int(np.round(r_pix[0])),int(hdr[depth]['NAXIS1'])-1)
                    y_pix = min(int(np.round(r_pix[1])),int(hdr[depth]['NAXIS2'])-1)
                    return data[depth][y_pix,x_pix]*u.Jy
                col = Column(np.repeat(-1000.0,len(qt)),name=f"rms_noise_bane{('_%s'%image_depth) if is_suffix else ''}",unit=u.mJy)
                for i in range(len(col)):
                    if isinstance(qt,QTableClass):
                        pos = SkyCoord(qt[i]['ra'],qt[i]['dec'])
                    else:
                        pos = SkyCoord(qt[i]['ra'],qt[i]['dec'],unit=(qt['ra'].unit,qt['dec'].unit))
                    #depth = qt[i]['image_type'] if 'image_type' in qt.colnames else image_depth
                    depth = qt[i]['image_type'] if image_depth is None else image_depth
                    col[i] = get_noise(pos,depth).to(u.mJy).value
                return col
            # add bane rms noise calcs to catalogues
            for module in self.catalogues['modules']:
                for depth in self.catalogues['modules'][module]['tables']:
                    qt = self.catalogues['modules'][module]['tables'][depth]
                    if len(qt) > 0 and 'ra' in qt.colnames and 'dec' in qt.colnames:
                        qt.add_column(create_noise_column(qt,depth))

            # create common catalogue
            def get_translate_func(catalogue_common):
                # get common catalogue fields and info
                def get_common_fields(catalogue_common):
                    common_fields = dict()
                    for field in catalogue_common:
                        fname = list(field.keys())[0]
                        units = eval(re.sub(r"(\w+)",r"u.\1",field[fname]['units'])) if 'units' in field[fname] else None
                        if 'use_if_nan' in field[fname]:
                            use_if_nan = field[fname]['use_if_nan']
                            if not units is None:
                                use_if_nan *= units
                        else:
                            use_if_nan = None
                        for module in field[fname]['modules']:
                            mname = list(module.keys())[0]
                            hname = module[mname]['header'] if not module[mname] is None and 'header' in module[mname] else fname
                            mult_by = module[mname]['mult_by'] if not module[mname] is None and 'mult_by' in module[mname] else None
                            datum = dict()
                            datum['header'] = fname
                            if not units is None:
                                datum['units'] = units
                            if not use_if_nan is None:
                                datum['use_if_nan'] = use_if_nan
                            if not mult_by is None:
                                datum['mult_by'] = mult_by
                            if not mname in common_fields:
                                common_fields[mname] = dict()
                            common_fields[mname][hname] = datum
                    return common_fields
                def get_required_fields(catalogue_common):
                    #return [list(f.keys())[0] for f in catalogue_common]
                    return [list(f.keys())[0] for f in config.get_catalogue_common()]
                common_fields   = get_common_fields(catalogue_common)
                required_fields = get_required_fields(catalogue_common)
                def translate(table,module):
                    meta = common_fields[module]
                    ## debug
                    #if module == 'aegean':            # debug
                    #    print_warning(f" => {meta}")  # debug
                    qt = QTable()
                    qt.add_column(table['id'])
                    for field in meta:
                        qt.add_column(table[field],name=meta[field]['header'])
                        if 'units' in meta[field] and qt[meta[field]['header']].unit != meta[field]['units']:
                            qt[meta[field]['header']] = qt[meta[field]['header']].to(meta[field]['units'])
                        if 'mult_by' in meta[field]:
                            qt[meta[field]['header']] = meta[field]['mult_by']*qt[meta[field]['header']]
                    for field in required_fields:
                        if not field in qt.columns:
                            # TO-DO: 
                            # [1] Need to make more generic, re., should probably handle correct dtype
                            # [2] Also mask is not doing what I'd expected: i.e., I'm not getting cols with '--',
                            #     as in example in https://docs.astropy.org/en/stable/table/masking.html,
                            #        >>> from astropy.table import Table, Column, MaskedColumn
                            #        >>> a = MaskedColumn([1, 2], name='a', mask=[False, True], dtype='i4')
                            #        >>> b = Column([3, 4], name='b', dtype='i8')
                            #        >>> Table([a, b])
                            #        <Table length=2>
                            #          a     b
                            #        int32 int64
                            #        ----- -----
                            #            1     3
                            #           --     4
                            #     So I've put in -1000 for now.
                            col = MaskedColumn([-1000 for _ in range(len(table))],name=field,mask=[False for _ in range(len(table))])
                            qt.add_column(col)
                    # remove nan's
                    ## debug
                    #if module == 'aegean': # debug
                    #    print_warning(" ==> BART!") # debug
                    #    qd = qt # debug
                    #    print_warning(qd[(321.14*u.deg<=qd['ra'])&(qd['ra']<=321.15*u.deg)&(-55.53*u.deg<=qd['dec'])&(qd['dec']<=-55.52*u.deg)]) # debug
                    #    print_warning("[DONE]")
                    for i in range(len(qt)-1,-1,-1):
                        for field in meta:
                            if np.isnan(qt[meta[field]['header']][i]):
                                #qt.remove_row(i)
                                #break
                                if 'use_if_nan' in meta[field]:
                                    ## debug
                                    #if module == 'aegean': # debug
                                    #    tmp = qt[meta[field]['header']][i] # debug
                                    qt[meta[field]['header']][i] = meta[field]['use_if_nan']
                                    ## debug
                                    #if module == 'aegean': # debug
                                    #    print_warning(f"*** ==> FIXED: {tmp} ==> {qt[meta[field]['header']][i]}") # debug
                                else:
                                    ## debug
                                    #if module == 'aegean': # debug
                                    #    print_ok(f"***") # debug
                                    #    print_ok(f" ==>     FIELD: {field}") # debug
                                    #    print_ok(f" ==>      META: {meta[field]}") # debug
                                    #    print_ok(f" ==>    HEADER: {meta[field]['header']}") # debug
                                    #    print_ok(f" ==> qt of {i}: {qt[meta[field]['header']][i]}") # debug
                                    qt.remove_row(i)
                                    break
                    ## debug
                    #if module == 'aegean': # debug
                    #    print_warning(" ==> LISA!") # debug
                    #    qd = qt # debug
                    #    print_warning(qd[(321.14*u.deg<=qd['ra'])&(qd['ra']<=321.15*u.deg)&(-55.53*u.deg<=qd['dec'])&(qd['dec']<=-55.52*u.deg)]) # debug
                    #    print_ok(f"REQUIRED_FIELDS: {required_fields}")
                    #    print_warning("[DONE]")
                    return qt
                return translate
            translate = get_translate_func(config.get_catalogue_common())
            common_catalogue = dict()
            for module in self.catalogues['modules']:
                ## debug
                #if module == 'aegean':  # debug
                #    qd = self.catalogues['modules'][module]['tables']['deep'] # debug
                #    print_warning(qd[(321.14<=qd['ra'])&(qd['ra']<=321.15)&(-55.53<=qd['dec'])&(qd['dec']<=-55.52)]['id','ra','dec','semimajor','semiminor','pa'])
                common_catalogue[module] = {
                    'deep':    translate(self.catalogues['modules'][module]['tables']['deep'],module),
                    'shallow': translate(self.catalogues['modules'][module]['tables']['shallow'],module),
                }

            ## debug
            #print_ok("***<TEST III>***")
            #print(f"cols: {common_catalogue['aegean']['deep'].colnames}")
            #print_warning(f"CATALOGUE_COMMON: {config.get_catalogue_common()}")
            #print_warning("DEEP TEST")
            #qd = common_catalogue['aegean']['deep']
            ##print_ok(qd)
            #print(qd[(321.14*u.deg<=qd['ra'])&(qd['ra']<=321.15*u.deg)&(-55.53*u.deg<=qd['dec'])&(qd['dec']<=-55.52*u.deg)])
            #print_warning("SHALLOW TEST")
            #qs = common_catalogue['aegean']['shallow']
            #print(qs[(321.14*u.deg<=qs['ra'])&(qs['ra']<=321.15*u.deg)&(-55.53*u.deg<=qs['dec'])&(qs['dec']<=-55.52*u.deg)])
            #print_warning("[DONE]")
            #print_ok("***")
            #exit()

            # get the graphics context
            self.graphics_context = config.get_graphics_context()

            # conglomerate external catalogues
            if not catalogues_yml is None:
                extern_config = yml.load(open(catalogues_yml))
                path = re.sub(r"^((.*?/)*).*$",r"\1",catalogues_yml)
                translate = get_translate_func(extern_config['common'])
                for module in extern_config['modules']:
                    mname = list(module.keys())[0]
                    fname = path+module[mname]['fits_file']
                    table = QTable.read(fname)
                    table.add_column([r for r in range(1,len(table)+1)],name='id',index=0)
                    if not 'extern' in self.catalogues:
                        self.catalogues['extern'] = dict()
                    self.catalogues['extern'][mname] = {
                        'file_name': fname,
                        'table': translate(table,mname)
                    }
                    common_catalogue[mname] = {
                        'deep':    self.catalogues['extern'][mname]['table'],
                        'shallow': self.catalogues['extern'][mname]['table']
                    }
                # update graphics context
                self.graphics_context['plot']['colors'].extend(extern_config['graphics']['plot']['colors'])
                self.graphics_context['plot']['labels'].extend(extern_config['graphics']['plot']['labels'])

            # get common catalogue meta data
            common_meta = dict()
            for col in common_catalogue[list(common_catalogue.keys())[0]]['deep'].columns:
                common_meta[col] = {
                    'dtype': common_catalogue[list(common_catalogue.keys())[0]]['deep'][col].dtype,
                    'units': common_catalogue[list(common_catalogue.keys())[0]]['deep'][col].unit
                }

            # create cluster catalogue using common catalogue's [ra,dec] with ellipsoidal-extents
            cluster_parameters = list()
            ref_id_offset = 0
            for module in common_catalogue:
                for image_type in common_catalogue[module]:
                    df = common_catalogue[module][image_type].to_pandas()
                    if len(df) > 0:
                        extents = {
                            'survey': f"{module}_{image_type}",
                            'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                            'cat_nos': df.id.to_numpy(),
                            'max_extent': self.skirt_factor * np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                            'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                            'a_extents': self.skirt_factor * df.extent_semimajor.to_numpy() * u.deg,
                            'b_extents': self.skirt_factor * df.extent_semiminor.to_numpy() * u.deg,
                            't_extents': df.extent_angle.to_numpy() * u.deg
                        }
                        ref_id_offset += len(df)
                        cluster_parameters.append(extents)
            qt = Cluster(cluster_parameters).get_qtable()
            qt['extent_semimajor'] /= self.skirt_factor
            qt['extent_semiminor'] /= self.skirt_factor
            qt.add_column([re.sub(r"^(.*?_)+","",s) for s in qt['survey']],name='image_type',index=3)
            qt['survey'] = [re.sub(r"_$","",re.sub(r"^(.*?_)+.*$",r"\1",s)) for s in qt['survey']]
            qt.rename_column('survey','source_finder')
            qt.add_column(qt['cat_no'],index=2,name='cat_id')
            qt.remove_column('cat_no')

            # trim isolated cluster from extern catalogue
            is_trim_cluster = False
            if is_trim_cluster:
                if 'extern' in self.catalogues:
                    extern_modules = set(self.catalogues['extern'].keys())
                    intern_modules = set(filter(lambda m: not m in extern_modules,set(common_catalogue.keys())))
                    for i in np.flip(np.unique(qt['clump_id'])):
                        clump_i = qt[(qt['clump_id']==i)]
                        if intern_modules.isdisjoint(set(clump_i['source_finder'])):
                            qt.remove_rows(np.array(clump_i['ref_id'])-1)
                    qt['ref_id'] = range(1,len(qt)+1)
                    df = pd.Index(np.unique(qt['clump_id']))
                    qt['clump_id'] = [df.get_loc(i)+1 for i in qt['clump_id']]

            # reconstitute cluster catalogue with rest of common catalogue
            add_fields = list(filter(lambda f: False if f in qt.columns else True,[list(f.keys())[0] for f in config.get_catalogue_common()]))
            add_units  = [common_meta[f]['units'] for f in add_fields]
            add_dtypes = [common_meta[f]['dtype'] for f in add_fields]
            add_cols = QTable(names=add_fields,dtype=add_dtypes)
            def extract(module,image_type,col_id):
                cat = common_catalogue[module][image_type]
                return np.array(cat[(cat['id']==col_id)][add_fields]).tolist()[0]
            for row in qt:
                add_cols.add_row(extract(row['source_finder'],row['image_type'],row['cat_id']))
            for i in range(len(add_fields)):
                if not add_units[i] is None:
                    add_cols[add_fields[i]] *= add_units[i]
            qt = hstack([qt,add_cols])

            # update main catalogue
            self.catalogues['cluster'] = qt

            # clump/cluster decomposition
            print_ok(f"Decomposing clumps...")
            print_ok(f"> clumps={np.max(self.catalogues['cluster']['clump_id'])}, components={np.max(self.catalogues['cluster']['ref_id'])}")
            print_ok(f"> working...")
            qt = self.catalogues['cluster'] # get pointer to the cluster catalogue
            match_idx = 1 # create global match_idx and add 'match_id' column
            qt.add_column(np.array(np.repeat(-1,len(qt))),name="match_id",index=2)
            for clump_i in np.unique(qt['clump_id']):
                clump = qt[(qt['clump_id']==clump_i)] # get ith clump
                cmap = dict() # initialize the match_id to ref_id map
                if len(clump) > 1:
                    matches = dict() # let's do some local matching for our finders/catalogues
                    finders = np.unique([f"{c['source_finder']}_{c['image_type']}" for c in clump])
                    for i in range(len(finders)-1):
                        fi = finders[i] # ith finder; avec, subclump, si, and positions, pi
                        si = clump[((clump['source_finder']==re.sub(r"_\w+$","",fi))&(clump['image_type']==re.sub(r".*?_(\w+)$",r"\1",fi)))]
                        pi = SkyCoord(si['ra'],si['dec'])
                        matches[fi] = {r:list() for r in si['ref_id']}
                        for j in range(i+1,len(finders)):
                            fj = finders[j] # jth (!=ith) finder; avec, subclump, sj, and positions, pj
                            sj = clump[((clump['source_finder']==re.sub(r"_\w+$","",fj))&(clump['image_type']==re.sub(r".*?_(\w+)$",r"\1",fj)))]
                            pj = SkyCoord(sj['ra'],sj['dec'])
                            idx, d2d, d3d = pi.match_to_catalog_sky(pj) # match (i,j)th finders/catalogues

                            # restrict matching to closest distances -- this is the tricky part!
                            cji_idxs = dict() # j-i finder index matcher
                            for cj in np.unique(idx):
                                ci_idx = None
                                ci_d2d = None
                                for ci in range(len(idx)):
                                    if cj==idx[ci]:
                                        if (ci_d2d is None) or (ci_d2d > d2d[ci]):
                                            ci_idx = ci
                                            ci_d2d = d2d[ci]
                                cji_idxs[cj] = ci_idx

                            # tabulate ref_id matches between ith finder and (i+1)th-nth finders
                            for cj in cji_idxs:
                                matches[fi][si[cji_idxs[cj]]['ref_id']].append(sj[cj]['ref_id'])

                    # ok, we are done with this clump, let's map (cmap) our matches -- this is the subtle part!
                    matches = {k:v for c in matches for k,v in matches[c].items()}
                    while len(matches) > 0:
                        cmap[match_idx]=[list(matches.keys())[0]]+matches[list(matches.keys())[0]]
                        for ref_id in cmap[match_idx]:
                            if ref_id in matches:
                                del matches[ref_id]
                        match_idx += 1
                else:
                    cmap[match_idx] = clump['ref_id']
                    match_idx += 1

                # now update our cluster table
                for match_id in cmap:
                    for ref_id in cmap[match_id]:
                        qt['match_id'][(qt['ref_id']==ref_id)]=match_id
            # ok, finally, we have to number the unmatched components
            for i in range(len(qt['ref_id'])):
                if qt['match_id'][i]==-1:
                    qt['match_id'][i] = match_idx
                    match_idx += 1
            print_ok(f"> matches={np.max(self.catalogues['cluster']['match_id'])}")

            # add subclumps
            self.catalogues['cluster'].add_column(
                name = 'subclump_id',
                col = self.__get_subclusters()['subclump_id'],
                index = 2,
            )

            # add rms noise columns
            self.catalogues['cluster'].add_column(create_noise_column(self.catalogues['cluster']),index=12)
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    self.catalogues['extern'][catalogue]['table'].add_column(create_noise_column(self.catalogues['extern'][catalogue]['table'],'deep',True))
                    self.catalogues['extern'][catalogue]['table'].add_column(create_noise_column(self.catalogues['extern'][catalogue]['table'],'shallow',True))

            # add source finder category to cluster table
            source_finder_types = {(list(m.keys()))[0]: m[list(m.keys())[0]]['category'] for m in config.get_catalogue_modules()}
            self.catalogues['cluster'].add_column(
                name = 'source_finder_type',
                col  = [(lambda s: source_finder_types[s] if s in source_finder_types else 'injected')(s) for s in self.catalogues['cluster']['source_finder']],
                index = 6
            )

            # create source finder list
            self.source_finders = np.unique(self.catalogues['cluster']['source_finder'])
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    self.source_finders = self.source_finders[(self.source_finders!=catalogue)]

            # create clump table
            qt = QTable()
            qt.add_column(
                name = 'clump_id',
                col  = np.unique(self.catalogues['cluster']['clump_id'])
            )
            def get_centroid(clump_id):
                qt = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
                ra_disc  = np.sqrt((qt['extent_semimajor']*np.sin(qt['extent_angle']))**2+(qt['extent_semiminor']*np.cos(qt['extent_angle']))**2)
                dec_disc = np.sqrt((qt['extent_semimajor']*np.cos(qt['extent_angle']))**2+(qt['extent_semiminor']*np.sin(qt['extent_angle']))**2)
                ra_min  = qt['ra']  - ra_disc
                ra_max  = qt['ra']  + ra_disc
                dec_min = qt['dec'] - dec_disc
                dec_max = qt['dec'] + dec_disc
                bb_ra_min = min(min(ra_min),min(ra_max))
                bb_ra_max = max(max(ra_min),max(ra_max))
                bb_dec_min = min(min(dec_min),min(dec_max))
                bb_dec_max = max(max(dec_min),max(dec_max))
                ra  = (bb_ra_min+bb_ra_max)/2.0
                dec = (bb_dec_min+bb_dec_max)/2.0
                return SkyCoord(ra,dec)
            centroids = [get_centroid(i) for i in qt['clump_id']]
            qt.add_column(
                name = 'ra',
                col  = [c.ra.to(u.deg) for c in centroids]
            )
            qt.add_column(
                name = 'dec',
                col  = [c.dec.to(u.deg) for c in centroids]
            )
            def get_max_clump_size(clump_id):
                qt = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
                ra_disc  = np.sqrt((qt['extent_semimajor']*np.sin(qt['extent_angle']))**2+(qt['extent_semiminor']*np.cos(qt['extent_angle']))**2)
                dec_disc = np.sqrt((qt['extent_semimajor']*np.cos(qt['extent_angle']))**2+(qt['extent_semiminor']*np.sin(qt['extent_angle']))**2)
                ra_min  = min(qt['ra']  - ra_disc)
                ra_max  = max(qt['ra']  + ra_disc)
                dec_min = min(qt['dec'] - dec_disc)
                dec_max = max(qt['dec'] + dec_disc)
                return max(abs(ra_max-ra_min),abs(dec_max-dec_min))
            qt.add_column(
                name = 'size',
                col = [get_max_clump_size(i).to(u.arcmin) for i in qt['clump_id']]
            )
            qt.add_column(
                name = 'n_clump',
                col = np.zeros(len(qt['clump_id']),dtype=np.int64)
            )
            for module in np.unique(self.catalogues['cluster']['source_finder']):
                qt.add_column(
                    name = f'n_{module}',
                    col = np.zeros(len(qt['clump_id']),dtype=np.int64)
                )
            for i in qt['clump_id']:
                qt_cl = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==i)]
                qt['n_clump'][i-1] = len(qt_cl['clump_id'])
                for module in np.unique(qt_cl['source_finder']):
                    qt[f'n_{module}'][i-1] = len(qt_cl[(qt_cl['source_finder']==module)])
            self.catalogues['clump'] = qt
            
            # update cluster table with residual stats
            is_normalized = True
            u_flux = u.mJy
            u_area = u.arcmin**2 if is_normalized else 1.0
            residual_parameters = {
                'residual_rms': u_flux/(u_area*u.beam),
                'residual_madfm': u_flux/(u_area*u.beam),
                'residual_sumsq': u_flux**2/(u_area*u.beam**2),
            }
            for res_par in reversed(list(residual_parameters.keys())):
                self.catalogues['cluster'].add_column(
                    name = res_par,
                    col  = np.array(np.repeat(-1000.0,len(self.catalogues['cluster'])),dtype=np.float64),
                    index = 14
                )
                self.catalogues['cluster'][res_par] *= residual_parameters[res_par]
            qt = self.catalogues['cluster']
            for i in self.catalogues['clump']['clump_id']:
                qt_cl = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==i)]
                for module in self.source_finders:
                    qt_sf = qt_cl[(qt_cl['source_finder']==module)]
                    if len(qt_sf)>0:
                        for image_type in self.image_types:
                            qt_it = qt_sf[(qt_sf['image_type']==image_type)]
                            if len(qt_it)>0:
                                d = self.catalogues['clump'][(self.catalogues['clump']['clump_id']==i)]
                                rs = self.residuals[module][image_type].get_residual_stats(d['ra'],d['dec'],d['size'])
                                area = d['size'].to(u.arcmin)**2 if is_normalized else 1.0
                                for ref_id in qt_it['ref_id']:
                                    qt['residual_rms'][(qt['ref_id']==ref_id)]   = (rs['rms']*u.Jy).to(u_flux)/(area*u.beam)
                                    qt['residual_madfm'][(qt['ref_id']==ref_id)] = (rs['madfm']*u.Jy).to(u_flux)/(area*u.beam)
                                    qt['residual_sumsq'][(qt['ref_id']==ref_id)] = (rs['sumsq']*u.Jy**2).to(u_flux**2)/(area*u.beam**2)

            # add s/n columns to cluster table
            qt.add_column(
                name = 'sn_bane',
                col = qt['flux_total']/qt['rms_noise_bane'],
                index = 14
            )
            qt.add_column(create_noise_column(self.catalogues['cluster'],'deep',True),index=15)
            qt.add_column(
                name = 'sn_bane_deep',
                col = qt['flux_total']/qt['rms_noise_bane_deep'],
                index = 16
            )
            qt.add_column(create_noise_column(self.catalogues['cluster'],'shallow',True),index=17)
            qt.add_column(
                name = 'sn_bane_shallow',
                col = qt['flux_total']/qt['rms_noise_bane_shallow'],
                index = 18
            )

            # add residual info to clump table
            sf_type_blank_field = max([len(s) for s in np.unique(self.catalogues['cluster']['source_finder_type'])])*' '
            for res_par in residual_parameters:
                for depth in self.image_types:
                    self.catalogues['clump'].add_column(
                        name = f"sf_{res_par.split('_')[1]}_{depth}",
                        col  = np.repeat(' '*max([len(s) for s in self.source_finders]),len(self.catalogues['clump'])),
                    )
                    self.catalogues['clump'].add_column(
                        name = f"{res_par}_{depth}",
                        col  = np.array(np.repeat(-1000.0,len(self.catalogues['clump'])),dtype=np.float64),
                    )
                    self.catalogues['clump'][f'{res_par}_{depth}'] *= residual_parameters[res_par]
            def fetch_clump(i):
                qt = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==i)]
                if 'extern' in self.catalogues:
                    for catalogue in self.catalogues['extern']:
                        qt = qt[(qt['source_finder']!=catalogue)]
                if 'gaussian' in np.unique(qt['source_finder_type']):
                    qt = qt[(qt['source_finder_type']=='gaussian')]
                return qt
            for i in self.catalogues['clump']['clump_id']:
                qt = fetch_clump(i)
                for res_par in residual_parameters:
                    for depth in self.image_types:
                        qt_depth = qt[(qt['image_type']==depth)]
                        if len(qt_depth)>0:
                            is_first = True
                            for j in range(len(qt_depth)):
                                if is_first:
                                    module = qt_depth['source_finder'][j]
                                    sf_type = qt_depth['source_finder_type'][j]
                                    min_par = qt_depth[res_par][j]
                                    is_first = False
                                elif min_par > qt_depth[res_par][j]:
                                    module = qt_depth['source_finder'][j]
                                    sf_type = qt_depth['source_finder_type'][j]
                                    min_par = qt_depth[res_par][j]
                            self.catalogues['clump'][f"sf_{res_par.split('_')[1]}_{depth}"][(self.catalogues['clump']['clump_id']==i)] = module
                            self.catalogues['clump'][f"{res_par}_{depth}"][(self.catalogues['clump']['clump_id']==i)] = min_par

            # do some plotting
            self.plot_colors = {list(m.keys())[0]:m[list(m.keys())[0]] for m in self.graphics_context['plot']['colors']}
            self.plot_labels = {list(m.keys())[0]:m[list(m.keys())[0]] for m in self.graphics_context['plot']['labels']}
            self.plot_collection = self.__collect_plot_info()

            # create cutout files
            self.cutout_files = self.__create_cutouts()

            # tarball it up!
            self.__tarball()
        # done -- phew!
        print_ok(f"[Done]")


    def __plot(self,plot_info,is_show=True,is_save=False):
        def attributize(info):
            if 'xlim' in info:
                plt.xlim(info['xlim'])
            if 'xline' in info:
                plt.axhline(info['xline'],color='black',linestyle=':')
            if 'xticks' in info:
                plt.xticks(info['xticks'])
            if 'ylim' in info:
                plt.ylim(info['ylim'])
            if 'xscale' in info:
                plt.xscale(info['xscale'])
            if 'yscale' in info:
                plt.yscale(info['yscale'])
            if 'grid' in info:
                plt.grid(info['grid'])
            if 'xlabel' in info:
                #plt.xlabel(info['xlabel'],fontsize=16,weight='bold')
                plt.xlabel(info['xlabel'],fontsize=16)
            if 'ylabel' in info:
                #plt.ylabel(info['ylabel'],fontsize=16,weight='bold')
                plt.ylabel(info['ylabel'],fontsize=16)
            if 'title' in info and not info['title'] is None:
                #plt.title(info['title'],fontsize=14,weight='bold')
                plt.title(info['title'],fontsize=14)
            if 'legend' in info:
                #plt.legend(**info['legend'],fontsize=14,prop={'weight':'bold'})
                plt.legend(**info['legend'],fontsize=14)
            if 'curves' in info:
                for curve in info['curves']:
                    kwargs = {
                        'color': curve['c'],
                        'linestyle': curve['ls'],
                    }
                    plt.plot(curve['x'],curve['y'],**kwargs)
            #plt.xticks(fontsize=14,weight='bold')
            plt.xticks(fontsize=14)
            #plt.yticks(fontsize=14,weight='bold')
            plt.yticks(fontsize=14)

            ax = plt.gca()
            # https://e2eml.school/matplotlib_ticks.html
            ax.tick_params(width=1.5)
            ax.tick_params(which='minor',width=1.5)
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.5)
        for info in plot_info:
            fig = plt.figure(figsize=(10,8))
            kwargs = info['data']['kwargs'] if 'kwargs' in info['data'] else {}
            if info['type']=='hist':
                for k in info['data']['n']:
                    plt.hist(info['data']['n'][k],edgecolor=info['data']['color'][k],**kwargs)
            elif info['type']=='step':
                for k in info['data']['y']:
                    #plt.step(info['data']['x'],info['data']['y'][k],color=info['data']['color'][k],where='pre')
                    plt.step(info['data']['x'],info['data']['y'][k],color=info['data']['color'][k],where='mid')
            elif info['type']=='line':
                for k in info['data']['x']:
                    plt.plot(info['data']['x'][k],info['data']['y'][k],color=info['data']['color'][k])
            elif info['type']=='error_bar':
                for k in info['data']['x']:
                    plt.errorbar(info['data']['x'][k],info['data']['y'][k],yerr=info['data']['y_errors'][k],color=info['data']['color'][k],**kwargs)
            elif info['type']=='scatter':
                plt.scatter(info['data']['x'],info['data']['y'],color=info['data']['color'] if 'color' in info['data'] else 'b',s=5)
            attributize(info)
            if is_show:
                plt.show()
            if is_save:
                plt.savefig(info['filename'])
            plt.close(fig)


    def __translate_to_fluxes_column(self,catalogue,is_sn=True,depth=None):
        local_copy = catalogue.copy()
        rms_noise_bane = f"rms_noise_bane{'' if depth is None else ('_%s'%depth)}"
        local_copy.add_column(
            name = 'fluxes',
            col  = (catalogue['flux_total'].to(u.mJy)/(catalogue[rms_noise_bane].to(u.mJy) if is_sn else 1.0)).value
        )
        #local_copy.remove_columns(['flux_total',rms_noise_bane])
        return local_copy


    def __get_flux_bin_info(self,is_sn=True):
        datum = dict()
        for depth in self.image_types:
            # extract flux related info
            cluster = self.catalogues['cluster'][(self.catalogues['cluster']['image_type']==depth)]
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    cluster = cluster[(cluster['source_finder']!=catalogue)]
            cluster = cluster['flux_total','rms_noise_bane']
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    rms_noise_bane = f"rms_noise_bane_{depth}"
                    dcat = self.catalogues['extern'][catalogue]['table']['flux_total',rms_noise_bane]
                    dcat.rename_column(rms_noise_bane,'rms_noise_bane')
                    cluster = vstack([cluster,dcat])

            # extract S/N or fluxes
            cluster = self.__translate_to_fluxes_column(cluster,is_sn)

            # make log-bins
            cluster = cluster[(cluster['fluxes']>0)]
            x_min = min(cluster['fluxes'])
            x_max = max(cluster['fluxes'])
            bins = np.logspace(np.log10(x_min),np.log10(x_max),self.default_bins+1)

            # get S/N or flux bin centers
            fluxes = list()
            for i in range(len(bins)-1):
                fluxes.append((bins[i]+bins[i+1])/2.0)

            # collect info
            datum[depth] = {'fluxes': fluxes, 'bins': bins, 'range': (x_min,x_max)}

        # compilate 
        x_min = min([min(datum[d]['range']) for d in datum])
        x_max = max([max(datum[d]['range']) for d in datum])
        bins = np.logspace(np.log10(x_min),np.log10(x_max),self.default_bins+1)
        fluxes = list()
        for i in range(len(bins)-1):
            fluxes.append((bins[i]+bins[i+1])/2.0)
        datum['global'] = {'fluxes': fluxes, 'bins': bins, 'range': (x_min,x_max)}

        return datum

    def __make_typhon_stats_plots(self):
        plot_info = list()
        x_min = np.min([np.min(self.catalogues['metrics']['statistics'][d]['rms_par']) for d in self.image_types])
        x_max = np.max([np.max(self.catalogues['metrics']['statistics'][d]['rms_par']) for d in self.image_types])
        for depth in self.image_types:
            qt = self.catalogues['metrics']['statistics'][depth]
            source_finders = np.unique(qt['module'])
            x_values = {s:np.unique(qt[(qt['module']==s)]['rms_par']).tolist() for s in source_finders}
            y_values = {s:[np.mean(qt[((qt['module']==s)&(qt['rms_par']==v))]['prd']) for v in x_values[s]] for s in source_finders}
            def miny(s,v,metric='prd'):
                prds = qt[((qt['module']==s)&(qt['rms_par']==v))][metric]
                return np.abs(np.mean(prds)-np.min(prds))
            def maxy(s,v,metric='prd'):
                prds = qt[((qt['module']==s)&(qt['rms_par']==v))][metric]
                return np.abs(np.max(prds)-np.mean(prds))
            y_errors = {s:[[miny(s,v) for v in x_values[s]],[maxy(s,v) for v in x_values[s]]] for s in source_finders}
            plot_info.append({
                'idx': f"hpp{1 if depth=='deep' else 2}",
                'type': 'error_bar',
                'filename': f"typhon_{depth}_optimization_stats.png",
                'html': {
                    'title': f"Typhon {depth.capitalize()} Run Optimization Statistics",
                    'datum': QTable(),
                },
                'data': {
                    'x': x_values,
                    'y': y_values,
                    'y_errors': y_errors,
                    'color': {s:self.plot_colors[s] for s in source_finders},
                    'kwargs': {
                        'fmt': '-o',
                        'capsize': 5,
                        'markersize': 4,
                    },
                },
                #'xlim': (0,x_max),
                #'xlim': (0,1010.0*x_max/1000.0),
                'xlim': (990.0*x_min/1000.0,1010.0*x_max/1000.0),
                'xline': self.percent_real_detections_cut,
                'ylim': (0,105),
                'xlabel': "RMS Parameters",
                'ylabel': 'PRD = ($N_{Image}$ - $N_{Inv.\\;Image}$) x 100 / $N_{Image}$',
                'grid': True,
                'title': f"Typhon: {depth.capitalize()} Percent Real Detections (PRD) Run Statistics",
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in source_finders], 
                    'framealpha': 0.85,
                    'loc': 'upper left',
                },
            })
            #print_warning(f"x_values: {x_values}")
            #print_ok(f"y_values: {y_values}")
            y_values = dict()
            y_errors = dict()
            qt.add_column(
                name = 'prd_cpu_time',
                col = qt['dt_image']+qt['dt_invim']
            )
            metrics = {'prd_cpu_time': 'PRD CPU Time (s)', 'rms': 'RMS', 'madfm': 'MADFM', 'sumsq': 'SumSq'}
            for metric in metrics:
                y_values[metric] = {s:[np.mean(qt[((qt['module']==s)&(qt['rms_par']==v))][metric]) for v in x_values[s]] for s in source_finders}
                y_errors[metric] = {s:[[miny(s,v,metric) for v in x_values[s]],[maxy(s,v,metric) for v in x_values[s]]] for s in source_finders}
                plot_item = {
                    'idx': f"hpp{1 if depth=='deep' else 2}_{metric}",
                    'type': 'error_bar',
                    'filename': f"typhon_{depth}_{metric}_optimization_stats.png",
                    'html': {
                        'title': f"Typhon {depth.capitalize()} Run {metrics[metric]} Optimization Statistics",
                        'datum': QTable(),
                    },
                    'data': {
                        'x': x_values,
                        'y': y_values[metric],
                        'y_errors': y_errors[metric],
                        'color': {s:self.plot_colors[s] for s in source_finders},
                        'kwargs': {
                            'fmt': '-o',
                            'capsize': 5,
                            'markersize': 4,
                        },
                    },
                    #'xlim': (0,x_max),
                    #'xlim': (0,1010.0*x_max/1000.0),
                    'xlim': (990.0*x_min/1000.0,1010.0*x_max/1000.0),
                    #'ylim': (0,105),
                    'xlabel': "RMS Parameters",
                    'ylabel': f"{metrics[metric]}",
                    'grid': True,
                    'title': f"Typhon: {depth.capitalize()} {metrics[metric]} Run Statistics",
                    'legend': {
                        'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in source_finders], 
                        'framealpha': 0.85,
                        'loc': 'upper left',
                    },
                }
                if metric == 'prd_cpu_time':
                    plot_item['yscale'] = 'log'
                plot_info.append(plot_item)
        return plot_info

    def __make_hist_plots(self,is_sn=True):
        # add external catalogue fluxes
        fluxes = dict()
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                datum = self.catalogues['extern'][catalogue]['table']
                datum = datum[(datum['flux_total']>-1000*u.mJy)]
                if is_sn: # use S/N
                    fluxes[catalogue] = {
                        'deep':    (datum['flux_total'].to(u.mJy)/datum['rms_noise_bane_deep'].to(u.mJy)).value,
                        'shallow': (datum['flux_total'].to(u.mJy)/datum['rms_noise_bane_shallow'].to(u.mJy)).value,
                    }
                else:
                    datum = datum['flux_total'].to(u.mJy).value
                    fluxes[catalogue] = {
                        'deep':    datum,
                        'shallow': datum,
                    }

        # filter out extern catalogues from cluster table
        clusters = self.catalogues['cluster'][(self.catalogues['cluster']['flux_total']>-1000*u.mJy)]
        for catalogue in fluxes:
            clusters = clusters[(clusters['source_finder']!=catalogue)]

        # add source finder fluxes
        for module in np.unique(clusters['source_finder']):
            datum = clusters[(clusters['source_finder']==module)]
            deep    = datum[(datum['image_type']=='deep')]
            shallow = datum[(datum['image_type']=='shallow')]
            if is_sn: # use S/N
                fluxes[module] = {
                    'deep':    (deep['flux_total'].to(u.mJy)/deep['rms_noise_bane'].to(u.mJy)).value,
                    'shallow': (shallow['flux_total'].to(u.mJy)/shallow['rms_noise_bane'].to(u.mJy)).value,
                }
            else:
                fluxes[module] = {
                    'deep':    deep['flux_total'].to(u.mJy).value,
                    'shallow': shallow['flux_total'].to(u.mJy).value,
                }
        x_values = list()
        for module in fluxes:
            x_values.extend(fluxes[module]['deep'])
            x_values.extend(fluxes[module]['shallow'])
        x_min = min(x_values)
        x_max = max(x_values)
        del x_values

        # compilate plot info
        plot_info = [{
            'idx': 'dnh1',
            'type': 'hist',
            'filename': f"{'sn' if is_sn else 'flux'}_histogram_deep_image.png",
            'html': {
                'title': f"Deep N <i>vs.</i> {'S/N' if is_sn else 'Flux'} Histogram",
                'datum': QTable(),
            },
            'data': {
                'n': {m:fluxes[m]['deep'] for m in fluxes},
                'color': {m:self.plot_colors[m] for m in fluxes},
                'kwargs': {
                    'histtype': 'stepfilled',
                    'bins': np.logspace(np.log10(x_min),np.log10(x_max),self.default_bins),
                    'facecolor': 'None',
                    'linewidth': 2,
                },
             },
            'xscale': 'log',
            'yscale': 'symlog',
            'xlabel': 'S/N (Deep-Signal/Deep-Noise)' if is_sn else 'Total Flux (mJy)',
            'ylabel': 'N',
            'title': 'Deep Image',
            'legend': {
                'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in fluxes],
                'framealpha': 0.85,
            },
        }]
        plot_info.append(plot_info[0].copy())
        plot_info[1]['idx'] = 'snh2'
        plot_info[1]['filename'] = f"{'sn' if is_sn else 'flux'}_histogram_shallow_image.png" 
        plot_info[1]['html']['title'] = re.sub("Deep","Shallow",plot_info[0]['html']['title'])
        plot_info[1]['data'] = plot_info[0]['data'].copy()
        plot_info[1]['data']['n'] = {m:fluxes[m]['shallow'] for m in fluxes}
        plot_info[1]['title'] = 'Shallow Image'
        plot_info[1]['xlabel'] = 'S/N (Shallow-Signal/Shallow-Noise)' if is_sn else 'Total Flux (mJy)' 

        return plot_info 


    def __binner(self,data,bins):
        counts = np.zeros(len(bins)-1)
        for m_flux in data:
            for i in range(len(bins)-1):
                if bins[i] <= m_flux and (m_flux < bins[i+1] or (i==(len(bins)-2) and m_flux <= bins[i+1])):
                    counts[i] += 1
                    break
        return counts


    #def __get_real_detections(self,is_sn,catalogue,img_depth=None):
    #    # get 'real' detections w.r.t. injected 'catalogue'
    #    def is_reject(subclump):
    #        finders = np.unique(subclump['source_finder'])
    #        if (not catalogue in finders) or len(finders)==1:
    #            return True
    #        return False
    #    match_ids = list()
    #    for i in np.unique(self.catalogues['cluster']['match_id']):
    #        source_finders = np.unique(self.catalogues['cluster'][(self.catalogues['cluster']['match_id']==i)]['source_finder'])
    #        if not (len(source_finders)==1 and catalogue in source_finders):
    #            match_ids.append(i)
    #    real = QTable()
    #    rms_noise_bane = f"rms_noise_bane{'' if img_depth is None else ('_%s'%img_depth)}"
    #    for i in match_ids:
    #        clump = self.catalogues['cluster'][(self.catalogues['cluster']['match_id']==i)]
    #        for depth in self.image_types:
    #            mask = clump['image_type']==depth
    #            if is_reject(clump[mask]):
    #                clump = clump[~mask]
    #        if len(clump) > 0:
    #            clump = clump['ref_id','match_id','source_finder','image_type','flux_total',rms_noise_bane]
    #            real = vstack([real,clump]) if len(real) > 0 else clump
    #    real = real[(real['flux_total']>-1000*u.mJy)]
    #
    #    if is_sn:
    #        real.add_column(name='fluxes',col=(real['flux_total'].to(u.mJy)/real[rms_noise_bane].to(u.mJy)).value)
    #    else:
    #        real.add_column(name='fluxes',col=real['flux_total'].to(u.mJy).value)
    #    real.remove_columns(['flux_total',rms_noise_bane])
    #
    #    return real
    def __get_real_detections(self,is_sn,catalogue,img_depth=None):
        real = QTable()
        rms_noise_bane = f"rms_noise_bane{'' if img_depth is None else ('_%s'%img_depth)}"
        for i in np.unique(self.catalogues['cluster']['match_id']):
            clump = self.catalogues['cluster'][(self.catalogues['cluster']['match_id']==i)]
            source_finders = np.unique(clump['source_finder'])
            if len(source_finders)>1 and catalogue in source_finders:
                clump = clump['ref_id','clump_id','match_id','source_finder','image_type','flux_total',rms_noise_bane]
                real = vstack([real,clump]) if len(real) > 0 else clump
        real = real[(real['flux_total']>-1000*u.mJy)]

        if is_sn:
            real.add_column(name='fluxes',col=(real['flux_total'].to(u.mJy)/real[rms_noise_bane].to(u.mJy)).value)
        else:
            real.add_column(name='fluxes',col=real['flux_total'].to(u.mJy).value)
        real.remove_columns(['flux_total',rms_noise_bane])

        return real


    def __get_flux_bins(self,fluxes):
        fluxes = fluxes[(fluxes>0)]
        x_min = min(fluxes)
        x_max = max(fluxes)
        bins = np.logspace(np.log10(x_min),np.log10(x_max),self.default_bins+1)
        bins[0]  = x_min
        bins[-1] = x_max
        fluxes = list()
        for i in range(len(bins)-1):
            fluxes.append((bins[i]+bins[i+1])/2.0)
        return fluxes,bins


    def __is_close(self, a, b, rel_tol=1e-05, abs_tol=0.0):
        # notes: https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


    def __make_completeness_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"n_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            # get real sources
            real = self.__get_real_detections(is_sn,catalogue)

            # get injected sources
            injected = self.catalogues['extern'][catalogue]['table']
            injected = injected[(injected['flux_total']>-1000*u.mJy)]

            # ok, let's bin stuff...
            match_ids = np.unique(real['match_id'])
            for depth in self.image_types:
                def get_qt(injected_fluxes):
                    fluxes,bins = self.__get_flux_bins(injected_fluxes)
                    qt = QTable()
                    qt.add_column(name='id',col=range(1,len(fluxes)+1))
                    qt.add_column(name='fluxes',col=bins[0:-1])
                    qt.add_column(
                        name='n_injected',
                        col=np.array(self.__binner(injected_fluxes,bins),dtype=np.int64)
                    )
                    for source_finder in self.source_finders:
                        qt.add_column(name=f"n_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                        for match_id in match_ids:
                            subclump = real[((real['match_id']==match_id)&(real['image_type']==depth))]
                            n = len(subclump[(subclump['source_finder']==source_finder)])
                            if n != 0 and len(subclump[(subclump['source_finder']==catalogue)])>0:
                                flux = subclump[(subclump['source_finder']==catalogue)]['fluxes']
                                for i in range(len(bins)-1):
                                    if ((i==0 and self.__is_close(bins[i],flux)) or bins[i] <= flux) and (flux < bins[i+1] or (i==(len(bins)-2) and self.__is_close(flux,bins[i+1]))):
                                        qt[f"n_{source_finder}"][i] += n
                                        break
                        qt.add_column(name=source_finder,col=qt[f"n_{source_finder}"]/np.array(qt['n_injected'],dtype=np.float64))
                    return qt,bins
                injected_fluxes = self.__translate_to_fluxes_column(injected,is_sn,depth)['fluxes']
                qt,bins = get_qt(injected_fluxes)

                # let's rescale!
                idx = get_idx(qt)
                if idx>-1:
                    qt,bins = get_qt(injected_fluxes[(injected_fluxes>bins[idx])])

                #if depth=='deep':
                #    qt.pprint_all()
                #    print(f"bins=[{','.join(['%f' % b for b in bins])}]")
                #    exit()

                # update plot info
                plot_info.append({
                    'idx': 'dcs1' if depth == 'deep' else 'scs2',
                    'type': 'step',
                    'filename': f"{'sn' if is_sn else 'flux'}_completeness_plot_simulated_{depth}_{catalogue.lower()}_image.png",
                    'html': {
                        'title': f"({depth.capitalize()} &cap; Injected)/Injected Completeness <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                        'datum': QTable(),
                    },
                    'data': {
                        'x': qt['fluxes'],
                        'y': {s:qt[s] for s in self.source_finders},
                        'color': {s:self.plot_colors[s] for s in self.source_finders},
                    },
                    'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                    'ylim': (0,1.1),
                    'xscale': 'log',
                    'xlabel': f"S/N ({depth.capitalize()}-Signal/{depth.capitalize()}-Noise)" if is_sn else 'Flux (mJy)',
                    'ylabel': 'Completeness',
                    #'title': f"({depth.capitalize()} \u2229 Injected)/Injected Injected-Catalogue Completeness Plot",
                    'title': f"({depth.capitalize()} \u2229 Injected)/Injected",
                    'legend': {
                        'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders],
                        'framealpha': 0.85,
                        'loc': 'lower right',
                     },
                })
        else:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            cluster = self.__translate_to_fluxes_column(self.catalogues['cluster'],depth='shallow')
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    cluster = cluster[(cluster['source_finder']!=catalogue)]
            def get_qt(cluster_fluxes):
                #fluxes,bins = self.__get_flux_bins(cluster[(cluster['image_type']=='deep')]['fluxes'])
                fluxes,bins = self.__get_flux_bins(cluster_fluxes)
                qt = QTable()
                datum = QTable(names=('bin_no','source_finder','ref_id_deep','clump_size_deep','ref_id_shallow','clump_size_shallow'),dtype=(np.int64,f'S{max([len(s) for s in self.source_finders])}',np.int64,np.int64,np.int64,np.int64))
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = cluster[(cluster['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        #img_types = np.unique(clump['image_type'])
                        if len(pos_deep)==0:
                            continue
                        idx, d2d, d3d = pos_shallow.match_to_catalog_sky(pos_deep)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_deep)):
                            flux = cat_deep['fluxes'][i]
                            for j in range(len(bins)-1):
                                #if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                if ((j==0 and self.__is_close(bins[j],flux)) or bins[j] <= flux) and (flux < bins[j+1] or (j==(len(bins)-2) and self.__is_close(flux,bins[j+1]))):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                        for k in range(len(idx)):
                                            if i==k:
                                                datum.add_row((j,source_finder,cat_deep['ref_id'][i],n_deep,cat_shallow['ref_id'][k],n_shallow))
                                                break
                                    else:
                                        datum.add_row((j,source_finder,cat_deep['ref_id'][i],n_deep,-1,n_shallow))
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt,bins,datum
            cluster_fluxes = cluster[(cluster['image_type']=='deep')]
            qt,bins,datum = get_qt(cluster_fluxes['fluxes'])

            # let's rescale!
            idx = get_idx(qt)
            if idx>-1:
                qt,bins,datum = get_qt(cluster_fluxes[(cluster_fluxes['fluxes']>bins[idx])]['fluxes'])

            #qt.pprint_all()
            #print(f"bins=[{','.join(['%f' % b for b in bins])}]")
            #exit()

            # let's compilate our html data
            datum.add_column(
                name = 'sn_bin_min',
                col = [bins[i] for i in datum['bin_no']],
                index = 1
            )
            datum.add_column(
                name = 'sn_avg',
                col = [(bins[i]+bins[i+1])/2.0 for i in datum['bin_no']],
                index = 2
            )
            datum.add_column(
                name = 'sn_bin_max',
                col = [bins[i+1] for i in datum['bin_no']],
                index = 3
            )
            datum.add_column(
                name = 'completeness',
                col = [qt[s][b] for s,b in zip(datum['source_finder'],datum['bin_no'])],
                index = 4
            )
            datum.add_column(
                name = 'clump_id_deep',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['clump_id'][0] for ref_id in datum['ref_id_deep']],
                index = 7
            )
            datum.add_column(
                name = 'flux_total_deep',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['flux_total'][0] for ref_id in datum['ref_id_deep']],
                index = 9
            )
            datum.add_column(
                name = 'rms_noise_bane_deep',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['rms_noise_bane'][0] for ref_id in datum['ref_id_deep']],
                index = 10
            )
            datum.add_column(
                name = 'sn_deep',
                col = datum['flux_total_deep']/datum['rms_noise_bane_deep'],
                index = 11
            )
            datum.add_column(
                name = 'clump_id_shallow',
                col = [-1 if ref_id < 0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['clump_id'][0] for ref_id in datum['ref_id_shallow']],
                index = 13
            )
            flux_units = self.catalogues['cluster']['flux_total'].unit
            datum.add_column(
                name = 'flux_total_shallow',
                col = [-1000*flux_units if ref_id<0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['flux_total'][0] for ref_id in datum['ref_id_shallow']],
            )
            rms_noise_units = self.catalogues['cluster']['rms_noise_bane'].unit
            datum.add_column(
                name = 'rms_noise_bane_shallow',
                col = [-1000*rms_noise_units if ref_id<0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['rms_noise_bane'][0] for ref_id in datum['ref_id_shallow']],
            )
            datum.add_column(
                name = 'sn_shallow',
                col = [-1000 if n==-1000.0*rms_noise_units else s/n for s,n in zip(datum['flux_total_shallow'],datum['rms_noise_bane_shallow'])],
            )
            datum.sort(['bin_no','source_finder','sn_deep'])

            # compilate plot info
            plot_info.append({
                'idx': 'dscs1',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_completeness_plot_deep_shallow_image.png",
                'html': {
                    'title': f"(Shallow &cap; Deep)/Deep Completeness <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': datum,
                    'fname': f"completeness_deep_shallow_image_plot_data_table.fits",
                 },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                'ylim': (0,1.1),
                'xscale': 'log',
                'xlabel': 'S/N (Deep-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Completeness',
                #'title': '(Shallow \u2229 Deep)/Deep Source Finder Completeness Plots',
                'title': '(Shallow \u2229 Deep)/Deep',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'lower right',
                 },
            })
        return plot_info 


    def __make_reliability_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            # get real sources
            real = self.__get_real_detections(is_sn,catalogue)
            cluster = self.__translate_to_fluxes_column(self.catalogues['cluster']['match_id','source_finder','image_type','flux_total','rms_noise_bane'],is_sn)
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    cluster = cluster[(cluster['source_finder']!=catalogue)]

            # ok, let's bin stuff...
            match_ids = np.unique(cluster['match_id'])
            for depth in self.image_types:
                def get_qt(cluster_fluxes):
                    fluxes,bins = self.__get_flux_bins(cluster_fluxes['fluxes'])
                    qt = QTable()
                    qt.add_column(name='id',col=range(1,len(fluxes)+1))
                    qt.add_column(name='fluxes',col=bins[0:-1])
                    for source_finder in self.source_finders:
                        qt.add_column(
                            name=f"nd_{source_finder}",
                            col  = np.array(self.__binner(cluster_fluxes[(cluster_fluxes['source_finder']==source_finder)]['fluxes'],bins),dtype=np.int64)
                        )
                        qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                        for match_id in match_ids:
                            subclump = real[((real['match_id']==match_id)&(real['image_type']==depth)&(real['source_finder']==source_finder))]
                            n = len(subclump)
                            if n != 0:
                                flux = subclump['fluxes']
                                for i in range(len(bins)-1):
                                    if ((i==0 and self.__is_close(bins[i],flux)) or bins[i] <= flux) and (flux < bins[i+1] or (i==(len(bins)-2) and self.__is_close(flux,bins[i+1]))):
                                        qt[f"ns_{source_finder}"][i] += n
                                        break
                        qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                    return qt,bins
                cluster_fluxes = cluster[(cluster['image_type']==depth)]
                qt,bins = get_qt(cluster_fluxes)

                # let's rescale!
                idx = get_idx(qt)
                if idx>-1:
                    qt,bins = get_qt(cluster_fluxes[(cluster_fluxes['fluxes']>bins[idx])])

                #if depth=='deep':
                #    qt.pprint_all()
                #    print(f"bins=[{','.join(['%f' % b for b in bins])}]")
                #    exit()

                # compilate plot info
                plot_info.append({
                    'idx': 'drs1' if depth=='deep' else 'srs2',
                    'type': 'step',
                    'filename': f"{'sn' if is_sn else 'flux'}_reliability_plot_simulated_{depth}_{catalogue.lower()}_image.png",
                    'html': {
                        'title': f"({depth.capitalize()} &cap; Injected)/{depth.capitalize()} Reliability <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                        'datum': QTable(),
                    },
                    'data': {
                        'x': qt['fluxes'],
                        'y': {s:qt[s] for s in self.source_finders},
                        'color': {s:self.plot_colors[s] for s in self.source_finders},
                    },
                    'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                    'ylim': (0,1.1),
                    'xscale': 'log',
                    'xlabel': f"S/N ({depth.capitalize()}-Signal/{depth.capitalize()}-Noise)" if is_sn else 'Flux (mJy)',
                    'ylabel': 'Reliability',
                    #'title': f"({depth.capitalize()} \u2229 Injected)/{depth.capitalize()} Injected-Catalogue Reliability Plot",
                    'title': f"({depth.capitalize()} \u2229 Injected)/{depth.capitalize()}",
                    'legend': {
                        'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders],
                        'framealpha': 0.85,
                        'loc': 'lower right',
                     },
                })
        else:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            cluster = self.__translate_to_fluxes_column(self.catalogues['cluster'],'shallow')
            if 'extern' in self.catalogues:
                for catalogue in self.catalogues['extern']:
                    cluster = cluster[(cluster['source_finder']!=catalogue)]
            def get_qt(cluster_fluxes):
                #fluxes,bins = self.__get_flux_bins(cluster[(cluster['image_type']=='shallow')]['fluxes'])
                fluxes,bins = self.__get_flux_bins(cluster_fluxes)
                qt = QTable()
                datum = QTable(names=('bin_no','source_finder','ref_id_shallow','clump_size_shallow','ref_id_deep','clump_size_deep'),dtype=(np.int64,f'S{max([len(s) for s in self.source_finders])}',np.int64,np.int64,np.int64,np.int64))
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = cluster[(cluster['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_shallow)==0:
                            continue
                        idx, d2d, d3d = pos_deep.match_to_catalog_sky(pos_shallow)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_shallow)):
                            flux = cat_shallow['fluxes'][i]
                            for j in range(len(bins)-1):
                                #if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                if ((j==0 and self.__is_close(bins[j],flux)) or bins[j] <= flux) and (flux < bins[j+1] or (j==(len(bins)-2) and self.__is_close(flux,bins[j+1]))):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                        for k in range(len(idx)):
                                            if i==k:
                                                datum.add_row((j,source_finder,cat_shallow['ref_id'][i],n_shallow,cat_deep['ref_id'][k],n_deep))
                                                break
                                    else:
                                        datum.add_row((j,source_finder,cat_shallow['ref_id'][i],n_shallow,-1,n_deep))
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt,bins,datum
            cluster_fluxes = cluster[(cluster['image_type']=='shallow')]
            qt,bins,datum = get_qt(cluster_fluxes['fluxes'])

            ## let's rescale!
            idx = get_idx(qt)
            if idx>-1:
                qt,bins,datum = get_qt(cluster_fluxes[(cluster_fluxes['fluxes']>bins[idx])]['fluxes'])

            #qt.pprint_all()
            #print(f"bins=[{','.join(['%f' % b for b in bins])}]")
            #exit()

            # let's compilate our html data
            datum.add_column(
                name = 'sn_bin_min',
                col = [bins[i] for i in datum['bin_no']],
                index = 1
            )
            datum.add_column(
                name = 'sn_avg',
                col = [(bins[i]+bins[i+1])/2.0 for i in datum['bin_no']],
                index = 2
            )
            datum.add_column(
                name = 'sn_bin_max',
                col = [bins[i+1] for i in datum['bin_no']],
                index = 3
            )
            datum.add_column(
                name = 'reliability',
                col = [qt[s][b] for s,b in zip(datum['source_finder'],datum['bin_no'])],
                index = 4
            )
            datum.add_column(
                name = 'clump_id_shallow',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['clump_id'][0] for ref_id in datum['ref_id_shallow']],
                index = 7
            )
            datum.add_column(
                name = 'flux_total_shallow',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['flux_total'][0] for ref_id in datum['ref_id_shallow']],
                index = 9
            )
            datum.add_column(
                name = 'rms_noise_bane_shallow',
                col = [self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['rms_noise_bane'][0] for ref_id in datum['ref_id_shallow']],
                index = 10
            )
            datum.add_column(
                name = 'sn_shallow',
                col = datum['flux_total_shallow']/datum['rms_noise_bane_shallow'],
                index = 11
            )
            datum.add_column(
                name = 'clump_id_deep',
                col = [-1 if ref_id < 0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['clump_id'][0] for ref_id in datum['ref_id_deep']],
                index = 13
            )
            flux_units = self.catalogues['cluster']['flux_total'].unit
            datum.add_column(
                name = 'flux_total_deep',
                col = [-1000*flux_units if ref_id<0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['flux_total'][0] for ref_id in datum['ref_id_deep']],
            )
            rms_noise_units = self.catalogues['cluster']['rms_noise_bane'].unit
            datum.add_column(
                name = 'rms_noise_bane_deep',
                col = [-1000*rms_noise_units if ref_id<0 else self.catalogues['cluster'][(self.catalogues['cluster']['ref_id']==ref_id)]['rms_noise_bane'][0] for ref_id in datum['ref_id_deep']],
            )
            datum.add_column(
                name = 'sn_deep',
                col = [-1000 if n==-1000.0*rms_noise_units else s/n for s,n in zip(datum['flux_total_deep'],datum['rms_noise_bane_deep'])],
            )
            datum.sort(['bin_no','source_finder','sn_shallow'])

            # compilate plot info
            plot_info.append({
                'idx': 'dsrs2',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_reliability_plot_deep_shallow_image.png",
                'html': {
                    'title': f"(Shallow &cap; Deep)/Shallow Reliability <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': datum,
                    'fname': f"reliability_deep_shallow_image_plot_data_table.fits",
                },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                'ylim': (0,1.1),
                'xscale': 'log',
                'xlabel': 'S/N (Shallow-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Reliability',
                #'title': '(Shallow \u2229 Deep)/Shallow Source Finder Reliability Plots',
                'title': '(Shallow \u2229 Deep)/Shallow',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'lower right',
                 },
            })
        return plot_info


    def __make_goodness_of_completeness_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            real = self.__get_real_detections(is_sn,catalogue,'shallow')
            cluster = self.catalogues['cluster']
            for col in ['ra','dec','extent_semimajor','extent_semiminor','extent_angle']:
                real.add_column(
                     name = col,
                     col  = [cluster[(cluster['ref_id']==rid)][col][0] for rid in real['ref_id']]
                )
            def get_qt(real_fluxes):
                #fluxes,bins = self.__get_flux_bins(real[(real['image_type']=='deep')]['fluxes'])
                fluxes,bins = self.__get_flux_bins(real_fluxes)
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = real[(real['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_deep)==0:
                            continue
                        idx, d2d, d3d = pos_shallow.match_to_catalog_sky(pos_deep)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_deep)):
                            flux = cat_deep['fluxes'][i]
                            for j in range(len(bins)-1):
                                if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt,bins
            real_fluxes = real[(real['image_type']=='deep')]
            qt,bins = get_qt(real_fluxes['fluxes'])

            # let's rescale!
            idx = get_idx(qt)
            if idx>-1:
                qt,bins = get_qt(real_fluxes[(real_fluxes['fluxes']>bins[idx])]['fluxes'])

            # compilate plot info
            plot_info.append({
                'idx': 'dsgcs1',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_goodness_of_completeness_plot_deep_shallow_image.png",
                'html': {
                    'title': f"((Shallow&cap;Deep)&cap;Injected)/(Deep&cap;Injected) Goodness of Completeness <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': QTable(),
                 },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                'ylim': (0,1.1),
                'xscale': 'log',
                'xlabel': 'S/N (Deep-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Goodness of Completeness',
                #'title': '((Shallow \u2229 Deep) \u2229 Injected)/(Deep \u2229 Injected) Goodness of Completeness Plots',
                'title': '((Deep \u2229 Injected) \u2229 (Shallow \u2229 Injected)/(Deep \u2229 Injected)',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'lower right',
                 },
            })
            #print(qt)
        return plot_info


    def __make_goodness_of_reliability_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        idx = i-1
                        break
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                return np.int(np.min([idx,idx_max]))

            real = self.__get_real_detections(is_sn,catalogue,'shallow')
            cluster = self.catalogues['cluster']
            for col in ['ra','dec','extent_semimajor','extent_semiminor','extent_angle']:
                real.add_column(
                     name = col,
                     col  = [cluster[(cluster['ref_id']==rid)][col][0] for rid in real['ref_id']]
                )
            def get_qt(real_fluxes):
                #fluxes,bins = self.__get_flux_bins(real[(real['image_type']=='shallow')]['fluxes'])
                fluxes,bins = self.__get_flux_bins(real_fluxes)
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = real[(real['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_shallow)==0:
                            continue
                        idx, d2d, d3d = pos_deep.match_to_catalog_sky(pos_shallow)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_shallow)):
                            flux = cat_shallow['fluxes'][i]
                            for j in range(len(bins)-1):
                                if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt,bins
            real_fluxes = real[(real['image_type']=='shallow')]
            qt,bins = get_qt(real_fluxes['fluxes'])

            # let's rescale!
            idx = get_idx(qt)
            if idx>-1:
                qt,bins = get_qt(real_fluxes[(real_fluxes['fluxes']>bins[idx])]['fluxes'])

            # compilate plot info
            plot_info.append({
                'idx': 'dsgrs2',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_goodness_of_reliability_plot_deep_shallow_image.png",
                'html': {
                    'title': f"((Shallow &cap; Deep) &cap; Injected)/(Shallow &cap; Injected) Goodness of Reliability <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': QTable(),
                },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                'ylim': (0,1.1),
                'xscale': 'log',
                'xlabel': 'S/N (Shallow-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Goodness of Reliability',
                #'title': '((Shallow \u2229 Deep) \u2229 Injected)/(Shallow \u2229 Injected) Goodness of Reliability Plots',
                'title': '((Deep \u2229 Injected) \u2229 (Shallow \u2229 Injected)/(Shallow \u2229 Injected)',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'lower right',
                 },
            })
            #print(qt)
        return plot_info


    def __make_delta_completeness_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                b_idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        b_idx = i-1
                        break
                b_idx = np.int(np.min([b_idx,idx_max]))
                e_idx = -1
                for i in reversed(range(len(qt))):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        e_idx = i+1 if i<len(qt)-1 else -1
                        break
                e_idx = np.int(np.max([e_idx,idx_max+1 if idx_max<len(qt)-1 else len(qt)-1]))
                return [b_idx,e_idx]

            def get_cluster():
                cluster = self.__translate_to_fluxes_column(self.catalogues['cluster'],depth='shallow')
                if 'extern' in self.catalogues:
                    for catalogue in self.catalogues['extern']:
                        cluster = cluster[(cluster['source_finder']!=catalogue)]
                return cluster
            cluster = get_cluster()
            cluster_fluxes = cluster[(cluster['image_type']=='deep')]

            def get_real():
                real = self.__get_real_detections(is_sn,catalogue,'shallow')
                cluster = self.catalogues['cluster']
                for col in ['ra','dec','extent_semimajor','extent_semiminor','extent_angle']:
                    real.add_column(
                         name = col,
                         col  = [cluster[(cluster['ref_id']==rid)][col][0] for rid in real['ref_id']]
                    )
                return real
            real = get_real()
            real_fluxes = real[(real['image_type']=='deep')]

            collective_fluxes = vstack([cluster_fluxes['fluxes'],real_fluxes['fluxes']])
            cfluxes,cbins = self.__get_flux_bins(collective_fluxes['fluxes'])

            def get_dirty(fluxes,bins):
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = cluster[(cluster['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_deep)==0:
                            continue
                        idx, d2d, d3d = pos_shallow.match_to_catalog_sky(pos_deep)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_deep)):
                            flux = cat_deep['fluxes'][i]
                            for j in range(len(bins)-1):
                                if ((j==0 and self.__is_close(bins[j],flux)) or bins[j] <= flux) and (flux < bins[j+1] or (j==(len(bins)-2) and self.__is_close(flux,bins[j+1]))):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt

            def get_clean(fluxes,bins):
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = real[(real['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_deep)==0:
                            continue
                        idx, d2d, d3d = pos_shallow.match_to_catalog_sky(pos_deep)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_deep)):
                            flux = cat_deep['fluxes'][i]
                            for j in range(len(bins)-1):
                                if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt

            def get_delta(fluxes,bins):
                qt_dirty = get_dirty(fluxes,bins)
                qt_clean = get_clean(fluxes,bins)
                b_idx = min(get_idx(qt_dirty)[0],get_idx(qt_clean)[0])
                e_idx = min(get_idx(qt_dirty)[1],get_idx(qt_clean)[1])
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    qt.add_column(
                        name = source_finder,
                        col = qt_clean[source_finder]-qt_dirty[source_finder]
                    )
                return qt,b_idx,e_idx
            qt,b_idx,e_idx = get_delta(cfluxes,cbins)

            if b_idx>-1:
                collective_fluxes=collective_fluxes[(collective_fluxes['fluxes']>cbins[b_idx])]
            if e_idx>-1:
                collective_fluxes=collective_fluxes[(collective_fluxes['fluxes']<cbins[e_idx])]
            if b_idx>-1 or e_idx>-1:
                cfluxes,cbins = self.__get_flux_bins(collective_fluxes['fluxes'])
                qt,b_idx,e_idx = get_delta(cfluxes,cbins)

            # compilate plot info
            plot_info.append({
                'idx': 'deltac',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_residual_completeness_plot_deep_shallow_image.png",
                'html': {
                    'title': f"Delta Completeness <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': QTable(),
                 },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                #'ylim': (0,1.1),
                #'ylim': (-0.12,0.12), # is_paper: uncomment for emu_ext2x2_sim
                'xscale': 'log',
                'xlabel': 'S/N (Deep-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Residual Completeness',
                #'title': 'Residual Completeness Plots',
                'title': None if self.is_paper else 'Residual Completeness Plots',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'upper right',
                 },
            })
        return plot_info


    def __make_delta_reliability_plot(self,is_sn=True,catalogue=None):
        plot_info = list()
        if not catalogue is None:
            def get_idx(qt):
                idx_max = np.int(np.round(9.1*(len(qt)-1)/10.0))
                b_idx = -1
                for i in range(len(qt)):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        b_idx = i-1
                        break
                b_idx = np.int(np.min([b_idx,idx_max]))
                e_idx = -1
                for i in reversed(range(len(qt))):
                    if np.sum([qt[f"ns_{s}"][i] for s in self.source_finders]) > 0:
                        e_idx = i+1 if i<len(qt)-1 else -1
                        break
                e_idx = np.int(np.max([e_idx,idx_max+1 if idx_max<len(qt)-1 else len(qt)-1]))
                return [b_idx,e_idx]

            def get_cluster():
                cluster = self.__translate_to_fluxes_column(self.catalogues['cluster'],'shallow')
                if 'extern' in self.catalogues:
                    for catalogue in self.catalogues['extern']:
                        cluster = cluster[(cluster['source_finder']!=catalogue)]
                return cluster
            cluster = get_cluster()
            cluster_fluxes = cluster[(cluster['image_type']=='deep')]

            def get_real():
                real = self.__get_real_detections(is_sn,catalogue,'shallow')
                cluster = self.catalogues['cluster']
                for col in ['ra','dec','extent_semimajor','extent_semiminor','extent_angle']:
                    real.add_column(
                         name = col,
                         col  = [cluster[(cluster['ref_id']==rid)][col][0] for rid in real['ref_id']]
                    )
                return real
            real = get_real()
            real_fluxes = real[(real['image_type']=='deep')]

            collective_fluxes = vstack([cluster_fluxes['fluxes'],real_fluxes['fluxes']])
            cfluxes,cbins = self.__get_flux_bins(collective_fluxes['fluxes'])

            def get_dirty(fluxes,bins):
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = cluster[(cluster['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_shallow)==0:
                            continue
                        idx, d2d, d3d = pos_deep.match_to_catalog_sky(pos_shallow)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_shallow)):
                            flux = cat_shallow['fluxes'][i]
                            for j in range(len(bins)-1):
                                if ((j==0 and self.__is_close(bins[j],flux)) or bins[j] <= flux) and (flux < bins[j+1] or (j==(len(bins)-2) and self.__is_close(flux,bins[j+1]))):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt

            def get_clean(fluxes,bins):
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    sf_group = real[(real['source_finder']==source_finder)]
                    qt.add_column(name=f"nd_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    qt.add_column(name=f"ns_{source_finder}",col=np.zeros(len(fluxes),dtype=np.int64))
                    ref_id_offset = 0
                    cluster_parameters = list()
                    for depth in self.image_types:
                        df = sf_group[(sf_group['image_type']==depth)].to_pandas()
                        if len(df) > 0:
                            cluster_parameters.append({
                                'survey': depth,
                                'ref_ids': ref_id_offset+(df.index+1).to_numpy(),
                                'cat_nos': df['ref_id'],
                                'max_extent': np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                                'a_extents': df.extent_semimajor.to_numpy() * u.deg,
                                'b_extents': df.extent_semiminor.to_numpy() * u.deg,
                                't_extents': df.extent_angle.to_numpy() * u.deg
                            })
                            ref_id_offset += len(df)
                    ct = Cluster(cluster_parameters).get_qtable()
                    ct.remove_column('ref_id')
                    ct.rename_columns(['cat_no','survey'],['ref_id','image_type'])
                    ct.add_column(
                        name = 'fluxes',
                        col = [sf_group[(sf_group['ref_id']==ref_id)]['fluxes'][0] for ref_id in ct['ref_id']]
                    )
                    for clump_id in np.unique(ct['clump_id']):
                        clump = ct[(ct['clump_id']==clump_id)]
                        cat_deep    = clump[(clump['image_type']=='deep')]
                        cat_shallow = clump[(clump['image_type']=='shallow')]
                        pos_deep    = SkyCoord(cat_deep['ra'],cat_deep['dec'])
                        pos_shallow = SkyCoord(cat_shallow['ra'],cat_shallow['dec'])
                        if len(pos_shallow)==0:
                            continue
                        idx, d2d, d3d = pos_deep.match_to_catalog_sky(pos_shallow)
                        n_deep = len(cat_deep)
                        n_shallow = len(cat_shallow)
                        for i in range(len(cat_shallow)):
                            flux = cat_shallow['fluxes'][i]
                            for j in range(len(bins)-1):
                                if bins[j] <= flux and (flux < bins[j+1] or (j==(len(bins)-2) and flux <= bins[j+1])):
                                    qt[f"nd_{source_finder}"][j] += 1
                                    if i in idx:
                                        qt[f"ns_{source_finder}"][j] += 1
                                    break
                    qt.add_column(name=source_finder,col=qt[f"ns_{source_finder}"]/qt[f"nd_{source_finder}"])
                return qt

            def get_delta(fluxes,bins):
                qt_dirty = get_dirty(fluxes,bins)
                qt_clean = get_clean(fluxes,bins)
                b_idx = min(get_idx(qt_dirty)[0],get_idx(qt_clean)[0])
                e_idx = min(get_idx(qt_dirty)[1],get_idx(qt_clean)[1])
                qt = QTable()
                qt.add_column(name='id',col=range(1,len(fluxes)+1))
                qt.add_column(name='fluxes',col=bins[0:-1])
                for source_finder in self.source_finders:
                    qt.add_column(
                        name = source_finder,
                        col = qt_clean[source_finder]-qt_dirty[source_finder]
                    )
                return qt,b_idx,e_idx
            qt,b_idx,e_idx = get_delta(cfluxes,cbins)

            if b_idx>-1:
                collective_fluxes=collective_fluxes[(collective_fluxes['fluxes']>cbins[b_idx])]
            if e_idx>-1:
                collective_fluxes=collective_fluxes[(collective_fluxes['fluxes']<cbins[e_idx])]
            if b_idx>-1 or e_idx>-1:
                cfluxes,cbins = self.__get_flux_bins(collective_fluxes['fluxes'])
                qt,b_idx,e_idx = get_delta(cfluxes,cbins)

            # compilate plot info
            plot_info.append({
                'idx': 'deltar',
                'type': 'step',
                'filename': f"{'sn' if is_sn else 'flux'}_residual_reliability_plot_deep_shallow_image.png",
                'html': {
                    'title': f"Delta Reliability <i>vs.</i> {'S/N' if is_sn else 'Flux'} Plot",
                    'datum': QTable(),
                 },
                'data': {
                    'x': qt['fluxes'],
                    'y': {s:qt[s] for s in self.source_finders},
                    'color': {s:self.plot_colors[s] for s in self.source_finders},
                },
                'xlim': (qt['fluxes'][0],qt['fluxes'][-1]) if not self.is_paper else self.completeness_and_reliability_sn_plot_limits,
                #'ylim': (0,1.1),
                'xscale': 'log',
                'xlabel': 'S/N (Deep-Signal/Shallow-Noise)' if is_sn else 'Flux (mJy)',
                'ylabel': 'Residual Reliability',
                #'title': 'Residual Reliability Plots',
                'title': None if self.is_paper else 'Residual Reliability Plots',
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in self.source_finders], 
                    'framealpha': 0.85,
                    'loc': 'upper right',
                 },
            })
        return plot_info


    def __make_flux_ratio_plots(self,is_sn=True,catalogue=None):
        #is_sn = False
        plot_info = list()
        def get_complement(child):
            parent = self.catalogues['cluster']
            comp = QTable()
            if not catalogue is None:
                comp.add_column(
                    name = 'ref_id',
                    col = list(set(parent[(parent['source_finder']!=catalogue)]['ref_id'])-set(child[(child['source_finder']!=catalogue)]['ref_id']))
                )
            else:
                comp.add_column(
                    name = 'ref_id',
                    col = list(set(parent['ref_id'])-set(child['ref_id']))
                )
            comp.add_column(
                name = 'clump_id',
                col = [parent[(parent['ref_id']==ref_id)]['clump_id'][0] for ref_id in comp['ref_id']]
            )
            comp.add_column(
                name = 'match_id',
                col = [parent[(parent['ref_id']==ref_id)]['match_id'][0] for ref_id in comp['ref_id']]
            )
            comp.add_column(
                name = 'source_finder',
                col = [parent[(parent['ref_id']==ref_id)]['source_finder'][0] for ref_id in comp['ref_id']]
            )
            comp.add_column(
                name = 'image_type',
                col = [parent[(parent['ref_id']==ref_id)]['image_type'][0] for ref_id in comp['ref_id']]
            )
            if is_sn:
                comp.add_column(
                    name = 'fluxes',
                    col = [parent[(parent['ref_id']==ref_id)]['sn_bane'][0] for ref_id in comp['ref_id']]
                )
            else:
                comp.add_column(
                    name = 'fluxes',
                    col = [parent[(parent['ref_id']==ref_id)]['flux_total'][0].to(u.mJy).value for ref_id in comp['ref_id']]
                )
            if catalogue is None and 'extern' in self.catalogues:
                for module in self.catalogues['extern']:
                    comp = comp[(comp['source_finder']!=module)]
            comp.sort(['ref_id','clump_id','match_id','source_finder','image_type'])
            return comp
        if not catalogue is None:
            # TO-DO: Fix/investigate cases where num injected > 1 for a given match_id
            real = self.__get_real_detections(is_sn,catalogue)
            comp = get_complement(real)
            if is_sn:
                xlims = (1,max(real['fluxes']))
            else:
                xlims = (min(real['fluxes']),max(real['fluxes']))
            x_range = np.logspace(np.log(xlims[0]),np.log(xlims[1]),100)
            ylims = (0,10)

            def get_flux_ratio(ref_id):
                match_id = np.unique(self.catalogues['cluster']['match_id'][(self.catalogues['cluster']['ref_id']==ref_id)])[0]
                depth = np.unique(self.catalogues['cluster']['image_type'][(self.catalogues['cluster']['ref_id']==ref_id)])[0]
                in_flux = self.catalogues['cluster']['flux_total'][((self.catalogues['cluster']['match_id']==match_id)&(self.catalogues['cluster']['source_finder']==catalogue)&(self.catalogues['cluster']['image_type']==depth))][0].to(u.mJy).value
                out_flux = self.catalogues['cluster']['flux_total'][(self.catalogues['cluster']['ref_id']==ref_id)][0].to(u.mJy).value
                return out_flux/in_flux
            real.add_column(
                name = 'flux_ratio',
                col = [get_flux_ratio(rid) for rid in real['ref_id']]
            )
            def get_in_flux(ref_id):
                match_id = np.unique(self.catalogues['cluster']['match_id'][(self.catalogues['cluster']['ref_id']==ref_id)])[0]
                depth = np.unique(self.catalogues['cluster']['image_type'][(self.catalogues['cluster']['ref_id']==ref_id)])[0]
                if is_sn:
                    in_flux = self.catalogues['cluster']['sn_bane'][((self.catalogues['cluster']['match_id']==match_id)&(self.catalogues['cluster']['source_finder']==catalogue)&(self.catalogues['cluster']['image_type']==depth))][0]
                else:
                    in_flux = self.catalogues['cluster']['flux_total'][((self.catalogues['cluster']['match_id']==match_id)&(self.catalogues['cluster']['source_finder']==catalogue)&(self.catalogues['cluster']['image_type']==depth))][0].to(u.mJy).value
                return in_flux
            real.add_column(
                name = 'in_flux',
                col = [get_in_flux(rid) for rid in real['ref_id']]
            )

            #injected_rms = {d:np.std(self.catalogues['cluster'][(self.catalogues['cluster']['source_finder']==catalogue)&(self.catalogues['cluster']['image_type']==d)]['flux_total']).to(u.mJy).value for d in self.image_types}
            #print_warning(injected_rms)
            ## img: simulated/emu_simulated_04.fits
            ## deep:
            ##    Median: 2.0123256945225876 μJy
            ##    Mean:   2.6764882932184264 μJy
            ##    RMS:    21.53401328541804 μJy
            ## shallow:
            ##    Median: 9.558973804951647 μJy
            ##    Mean:   10.092792036823841 μJy
            ##    RMS:    111.52637464380608 μJy
            #injected_rms = {
            #    'deep':     (21.53401328541804*u.uJy).to(u.mJy).value,
            #    'shallow': (111.52637464380608*u.uJy).to(u.mJy).value,
            #}
            #injected_mean = {
            #    'deep':     (2.6764882932184264*u.uJy).to(u.mJy).value,
            #    'shallow': (10.092792036823841*u.uJy).to(u.mJy).value,
            #}
            #print_ok(injected_rms)
            #injected_rms = {d:injected_rms[d]+injected_mean[d] for d in injected_rms}
            #print(injected_rms)
            #print_warning(self.catalogues['metrics']['parameters'])
            #print_ok(self.catalogues['metrics']['rms_box_statistics']['deep'])
            #print(self.catalogues['metrics']['rms_box_statistics']['shallow'])
            injected = dict()
            for depth in self.image_types:
                metrics = self.catalogues['metrics']['rms_box_statistics'][depth]
                injected[depth] = {
                    'mean': metrics[(metrics['is_opt']>0)]['mean'].to(u.mJy).value,
                    'rms':  metrics[(metrics['is_opt']>0)]['rms' ].to(u.mJy).value,
                }
            metrics = self.catalogues['metrics']['parameters']
            sratio = {
                'deep': dict(),
                'shallow': dict(),
            }
            for source_finder in self.source_finders:
                for depth in self.image_types:
                    datum = real[((real['source_finder']==source_finder)&(real['image_type']==depth))]
                    rms = metrics[((metrics['source_finder']==source_finder)&(metrics['image_type']==depth))]['rms_parameter_value'][0]
                    legend = [
                        {'ls': '-' , 'tt': '$1\sigma$ Noise' },
                        {'ls': '--', 'tt': '$3\sigma$ Noise' },
                        {'ls': '-.', 'tt': '$S_{out}$=%.2f$\sigma_{S_{in}}$' % rms},
                        {'ls': ':',  'tt': '$S_{out}$=5$\sigma_{S_{in}}$'},
                    ]
                    if is_sn:
                        for n in [1,2,3]:
                            datum.add_column(
                               name = f'minus_{n}_sigma',
                               col  = 1.0-float(n)/datum['in_flux'],
                            )
                            datum.add_column(
                               name = f'plus_{n}_sigma',
                               col  = 1.0+float(n)/datum['in_flux'],
                            )
                            def get_sratio(ref_id):
                                rp = datum[(datum['ref_id']==ref_id)]
                                flx = rp['in_flux'][0]
                                a = rp[f'minus_{n}_sigma'][0]
                                b = rp[f'plus_{n}_sigma'][0]
                                base = datum[(datum['in_flux']>=flx)]
                                return len(base[((a<=base['flux_ratio'])&(base['flux_ratio']<=b))])/len(base)
                            datum.add_column(
                               name = f'frac_{n}_sigma',
                               col  = [get_sratio(rid) for rid in datum['ref_id']],
                            )
                    #print_warning(datum)
                    #exit()
                    sratio[depth][source_finder]=datum.copy()
                    plot_info.append({
                        'idx': f"frp{source_finder}{'1' if depth=='deep' else 2}",
                        'type': 'scatter',
                        'filename': f"{source_finder}_{depth}_flux_ratio_vs_injected_{catalogue.lower()}_{'sn' if is_sn else 'flux'}_plot.png",
                        'html': {
                            'title': f"{source_finder.capitalize()} {depth.capitalize()} Flux Ratio vs. Injected {catalogue} {'S/N' if is_sn else 'Flux'} Plot",
                            'datum': datum,
                            'fname': f"{source_finder}_{depth}_to_injected_flux_ratio_vs_{'sn' if is_sn else 'flux'}_plot_data_table.fits",
                        },
                        'data': {
                            'x': datum['in_flux'],
                            'y': datum['flux_ratio'],
                            'color': self.plot_colors[source_finder],
                        },
                        'curves': [
                            {'x': x_range, 'y': 1.0+1.0*(1.0 if is_sn else (injected[depth]['mean']+injected[depth]['rms']))/x_range, 'c': 'k', 'ls': '-'},
                            {'x': x_range, 'y': 1.0-1.0*(1.0 if is_sn else (injected[depth]['mean']+injected[depth]['rms']))/x_range, 'c': 'k', 'ls': '-'},
                            {'x': x_range, 'y': 1.0+3.0*(1.0 if is_sn else (injected[depth]['mean']+injected[depth]['rms']))/x_range, 'c': 'k', 'ls': '--'},
                            {'x': x_range, 'y': 1.0-3.0*(1.0 if is_sn else (injected[depth]['mean']+injected[depth]['rms']))/x_range, 'c': 'k', 'ls': '--'},
                            {'x': x_range, 'y': rms*(1.0 if is_sn else injected[depth]['mean']+injected[depth]['rms'])/x_range, 'c': 'k', 'ls': '-.'},
                            {'x': x_range, 'y': 5.0*(1.0 if is_sn else injected[depth]['mean']+injected[depth]['rms'])/x_range, 'c': 'k', 'ls': ':'},
                        ],
                        'xlim': xlims,
                        'xline': 1.0,
                        'ylim': ylims,
                        'xscale': 'log',
                        'xlabel': f"S/N (Injected-Signal/{depth.capitalize()}-Noise)" if is_sn else 'Flux In (mJy)',
                        'ylabel': '$S_{out}/S_{in}$',
                        'title': f"{self.plot_labels[source_finder]} {depth.capitalize()} Flux Ratios Plot",
                        'legend': {
                            'title': f"{self.plot_labels[source_finder]}",
                            'handles': [mlines.Line2D([],[],color='k',label=d['tt'],linestyle=d['ls']) for d in legend],
                            'loc': 'upper right',
                        },
                    })
                    #return plot_info # debug
                    #break # debug
                #break # debug
            #exit() # debug
            for depth in self.image_types:
                fluxes = dict()
                datum = comp[(comp['image_type']==depth)]
                for source_finder in self.source_finders:
                    fluxes[source_finder] = datum[(datum['source_finder']==source_finder)]['fluxes']
                plot_info.append({
                    'idx': f"frfph{1 if depth=='deep' else 2}",
                    'type': 'hist',
                    'filename': f"{'sn' if is_sn else 'flux'}_false_positives_histogram_{depth}_image.png",
                    'html': {
                        'title': f"{depth.capitalize()} Flase-Positives <i>vs.</i> {'S/N' if is_sn else 'Flux'} Histogram",
                        'datum': QTable(), # TO-DO
                        'fname': f"{depth}_injeted_false_positive_flux_ratio_vs_{'sn' if is_sn else 'flux'}_histogram_data_table.fits",
                    },
                    'data': {
                        'n': {m:fluxes[m] for m in fluxes},
                        'color': {m:self.plot_colors[m] for m in fluxes},
                        'kwargs': {
                            'histtype': 'stepfilled',
                            ##'bins': np.logspace(np.log10(np.min(datum['fluxes'].value)),np.log10(np.max(datum['fluxes'].value)),self.default_bins),
                            #'bins': np.logspace(np.log10(np.min(list(datum['fluxes']))),np.log10(np.max(list(datum['fluxes']))),self.default_bins),
                            'bins': np.logspace(np.log10(max(1,np.min(list(datum['fluxes'])))),np.log10(np.max(list(datum['fluxes']))),self.default_bins),
                            'facecolor': 'None',
                            'linewidth': 2,
                        },
                     },
                    'xscale': 'log',
                    #'yscale': 'symlog',
                    'xlabel': f"S/N ({depth.capitalize()}-Signal/{depth.capitalize()}-Noise)" if is_sn else "Total Flux (mJy)",
                    'ylabel': 'False-Positives',
                    'title': f"{depth.capitalize()} False-Positives Histogram",
                    'legend': {
                        'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in fluxes],
                        'framealpha': 0.85,
                    },
                })
            if is_sn:
                modes = ['linear','errors','bins']
                mode = 'linear'
                n=3 # select n-sigme where n = 1,2,3.
                is_diagnostic = True
                for depth in self.image_types:
                    if mode == 'linear':
                        def prt(inp_str):
                            if depth=='deep':
                                print_warning(inp_str)
                            else:
                                print_ok(inp_str)
                        for source_finder in self.source_finders:
                            sratio[depth][source_finder].sort(['in_flux'])
                            sratio[depth][source_finder][f'frac_{n}_sigma'] = 100.0*sratio[depth][source_finder][f'frac_{n}_sigma']
                            if is_diagnostic:
                                prt(f"   ===<{source_finder}::{depth}>===")
                                d=sratio[depth][source_finder][(sratio[depth][source_finder]['in_flux']<=100)]['in_flux','fluxes',f'frac_{n}_sigma']
                                prt(f"S/N ~ 1-10 =>  r = [({d['in_flux'][1]}, {d[f'frac_{n}_sigma'][1]}),({d['in_flux'][-1]}, {d[f'frac_{n}_sigma'][-1]})] %")
                                prt(f"               s = [({100-d['in_flux'][1]}, {100-d[f'frac_{n}_sigma'][1]}), ({100-d['in_flux'][-1]},{100-d[f'frac_{n}_sigma'][-1]})] %")
                                prt(f"               > s_min = {100-min(d[f'frac_{n}_sigma'])}, s_avg = {100-np.mean(d[f'frac_{n}_sigma'])}, s_max = {100-max(d[f'frac_{n}_sigma'])}")
                                prt(f"               > s = {100-np.mean(d[f'frac_{n}_sigma'])}+{max(d[f'frac_{n}_sigma'])-np.mean(d[f'frac_{n}_sigma'])}-{np.mean(d[f'frac_{n}_sigma'])-min(d[f'frac_{n}_sigma'])}")
                        plot_info.append({
                            'idx': f"fr{n}sr{1 if depth=='deep' else 2}",
                            'type': 'line',
                            'filename': f"sn_flux_ratio_{n}sigma_fractions_{depth}_image.png",
                            'html': {
                                'title': f"Flux-Ratio {n}&sigma; Fractions <i>vs.</i> S/N ({depth})",
                                'datum': QTable(),
                            },
                            'data': {
                                'x': {m:sratio[depth][m]['in_flux'] for m in sratio[depth]},
                                'y': {m:sratio[depth][m][f'frac_{n}_sigma'] for m in sratio[depth]},
                                'color': {m:self.plot_colors[m] for m in sratio[depth]},
                                'kwargs': {
                                    'fmt': '-o',
                                    'capsize': 5,
                                    'markersize': 4,
                                },
                            },
                            'xlim': xlims,
                            'xline': 100.0,
                            'ylim': (0,105),
                            'xscale': 'log',
                            'grid': True,
                            'xlabel': f"S/N (Injected-Signal/{depth.capitalize()}-Noise)",
                            'ylabel': f'Flux-Ratio Fraction < {n}$\sigma$ (%)',
                            'title': f"{depth.capitalize()} Fraction of Flux Ratios within {n}$\sigma$",
                            'legend': {
                                'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in sratio[depth]],
                                'framealpha': 0.85,
                                'loc': 'lower left',
                            },
                        })
                    elif mode == 'errors':
                        is_mad = False
                        def madfm(y):
                            md = np.median(y)
                            return np.median(np.abs(md-y))/0.6744888
                        for source_finder in self.source_finders:
                            sratio[depth][source_finder].sort(['in_flux'])
                            sratio[depth][source_finder][f'frac_{n}_sigma'] = 100.0*sratio[depth][source_finder][f'frac_{n}_sigma']
                            sratio[depth][source_finder] = {
                                'data': sratio[depth][source_finder],
                                'x': list(),
                                'y': list(),
                                'dy': list(),
                            }
                            #bins = np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),self.default_bins)
                            bins = np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),7)
                            for i in range(len(bins)-1):
                                data = sratio[depth][source_finder]['data'].copy()
                                if i < len(bins)-2:
                                    data = data[((bins[i]<=data['in_flux'])&(data['in_flux']<bins[i+1]))]
                                else:
                                    data = data[((bins[i]<=data['in_flux'])&(data['in_flux']<=bins[i+1]))]
                                sratio[depth][source_finder]['x'].append((bins[i]+bins[i+1])/2.0)
                                if is_mad:
                                    sratio[depth][source_finder]['y'].append(np.median(data[f'frac_{n}_sigma']))
                                    sratio[depth][source_finder]['dy'].append(madfm(data[f'frac_{n}_sigma']))
                                else:
                                    sratio[depth][source_finder]['y'].append(np.mean(data[f'frac_{n}_sigma']))
                                    sratio[depth][source_finder]['dy'].append(np.std(data[f'frac_{n}_sigma']))
                        y_str = 'Madian$\pm$MADFM' if is_mad else '$\mu\pm\sigma$'
                        plot_info.append({
                            'idx': f"fr{n}sr{1 if depth=='deep' else 2}",
                            'type': 'error_bar',
                            'filename': f"sn_flux_ratio_{n}sigma_fractions_{depth}_image.png",
                            'html': {
                                'title': f"Flux-Ratio {n}&sigma; Fractions <i>vs.</i> S/N ({depth})",
                                'datum': QTable(),
                            },
                            'data': {
                                'x': {m:sratio[depth][m]['x'] for m in sratio[depth]},
                                'y': {m:sratio[depth][m]['y'] for m in sratio[depth]},
                                'y_errors': {m:sratio[depth][m]['dy'] for m in sratio[depth]},
                                'color': {m:self.plot_colors[m] for m in sratio[depth]},
                                'kwargs': {
                                    'fmt': '-o',
                                    'capsize': 5,
                                    'markersize': 4,
                                },
                            },
                            'xlim': xlims,
                            'xline': 100.0,
                            'ylim': (0,105),
                            'xscale': 'log',
                            'grid': True,
                            'xlabel': f"S/N (Injected-Signal/{depth.capitalize()}-Noise)",
                            'ylabel': f'Flux-Ratio Fraction < {n}$\sigma$ (%, {y_str})',
                            'title': f"{depth.capitalize()} Fraction of Flux Ratios within {n}$\sigma$",
                            'legend': {
                                'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in sratio[depth]],
                                'framealpha': 0.85,
                                'loc': 'lower left',
                            },
                        })
                    elif mode == 'bins':
                        bins = np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),25)
                        for source_finder in self.source_finders:
                            sratio[depth][source_finder].sort(['in_flux'])
                            sratio[depth][source_finder][f'frac_{n}_sigma'] = 100.0*sratio[depth][source_finder][f'frac_{n}_sigma']
                            sratio[depth][source_finder] = {
                                'data': sratio[depth][source_finder],
                                'x': list(),
                                'y': list(),
                            }
                            for i in range(len(bins)-1):
                                data = sratio[depth][source_finder]['data'].copy()
                                if i < len(bins)-2:
                                    data = data[((bins[i]<=data['in_flux'])&(data['in_flux']<bins[i+1]))]
                                else:
                                    data = data[((bins[i]<=data['in_flux'])&(data['in_flux']<=bins[i+1]))]
                                sratio[depth][source_finder]['x'].append(bins[i])
                                sratio[depth][source_finder]['x'].append(bins[i+1])
                                sratio[depth][source_finder]['y'].append(np.mean(data[f'frac_{n}_sigma']))
                                sratio[depth][source_finder]['y'].append(np.mean(data[f'frac_{n}_sigma']))
                        plot_info.append({
                            'idx': f"fr{n}sr{1 if depth=='deep' else 2}",
                            'type': 'line',
                            'filename': f"sn_flux_ratio_{n}sigma_fractions_{depth}_image.png",
                            'html': {
                                'title': f"Flux-Ratio {n}&sigma; Fractions <i>vs.</i> S/N ({depth})",
                                'datum': QTable(),
                            },
                            'data': {
                                'x': {m:sratio[depth][m]['x'] for m in sratio[depth]},
                                'y': {m:sratio[depth][m]['y'] for m in sratio[depth]},
                                'color': {m:self.plot_colors[m] for m in sratio[depth]},
                            },
                            'xlim': xlims,
                            'xline': 100.0,
                            'ylim': (0,105),
                            'xscale': 'log',
                            'grid': True,
                            'xlabel': f"S/N (Injected-Signal/{depth.capitalize()}-Noise)",
                            'ylabel': f'Flux-Ratio Fraction < {n}$\sigma$ (%)',
                            'title': f"{depth.capitalize()} Fraction of Flux Ratios within {n}$\sigma$",
                            'legend': {
                                'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in sratio[depth]],
                                'framealpha': 0.85,
                                'loc': 'lower left',
                            },
                        })
            #print_warning(f"x: {plot_info[-2]['data']['x']['aegean']}")
            #print_warning(f"y: {plot_info[-2]['data']['y']['aegean']}")
            #return [plot_info[-2]] # debug
        else:
            def extract_datum(injected,detected):
                qt = detected['ref_id','clump_id','match_id','source_finder','image_type']
                qt.add_column(
                    name = 'fluxes',
                    col = detected['sn_bane_shallow'] if is_sn else detected['flux_total'].to(u.mJy)
                )
                qt.add_column(
                    name = 'flux_ratio',
                    col = detected['flux_total'].to(u.mJy)/injected['flux_total'].to(u.mJy)
                )
                qt.add_column(
                    name = 'in_ref_id',
                    col = injected['ref_id'],
                    index = 1
                )
                qt.add_column(
                    name = 'in_flux',
                    col = injected['sn_bane_shallow'] if is_sn else injected['flux_total'].to(u.mJy)
                )
                return qt
            ref_ids = list()
            real = {s:{'x':list(),'y':list(),'d': QTable()} for s in self.source_finders}
            for match_id in np.unique(self.catalogues['cluster']['match_id']):
                match_set = self.catalogues['cluster'][(self.catalogues['cluster']['match_id']==match_id)]
                for source_finder in self.source_finders:
                    datum = match_set[(match_set['source_finder']==source_finder)]
                    injected = datum[(datum['image_type']=='deep')]
                    detected = datum[(datum['image_type']=='shallow')]
                    if len(injected)==1 and len(detected)==1:
                        real[source_finder]['x'].append((injected['sn_bane_shallow'] if is_sn else injected['flux_total'].to(u.mJy))[0].value)
                        real[source_finder]['y'].append((detected['flux_total'].to(u.mJy)/injected['flux_total'].to(u.mJy))[0].value)
                        if len(real[source_finder]['d'])>0:
                            real[source_finder]['d'] = vstack([real[source_finder]['d'],extract_datum(injected,detected)])
                        else:
                            real[source_finder]['d'] = extract_datum(injected,detected)
                        ref_ids.extend(datum['ref_id'])
            if is_sn:
                for source_finder in self.source_finders:
                    datum = real[source_finder]['d']
                    for n in [1,2,3]:
                        datum.add_column(
                           name = f'minus_{n}_sigma',
                           col  = 1.0-float(n)/datum['in_flux'],
                        )
                        datum.add_column(
                           name = f'plus_{n}_sigma',
                           col  = 1.0+float(n)/datum['in_flux'],
                        )
                        def get_sratio(ref_id):
                            rp = datum[(datum['ref_id']==ref_id)]
                            flx = rp['in_flux'][0]
                            a = rp[f'minus_{n}_sigma'][0]
                            b = rp[f'plus_{n}_sigma'][0]
                            base = datum[(datum['in_flux']>=flx)]
                            return len(base[((a<=base['flux_ratio'])&(base['flux_ratio']<=b))])/len(base)
                        datum.add_column(
                           name = f'frac_{n}_sigma',
                           col  = [get_sratio(rid) for rid in datum['ref_id']],
                        )
            comp = get_complement(QTable([ref_ids],names=('ref_id',)))
            comp = comp[(comp['image_type']=='shallow')]
            xlims = (np.min([np.min(real[s]['x']) for s in self.source_finders]),np.max([np.max(real[s]['x']) for s in self.source_finders]))
            if is_sn:
                xlims = (max(1,xlims[0]),xlims[1])
            x_range = np.logspace(np.log(0.01),np.log(xlims[1]),100)
            #ylims = (0,np.max([np.max(real[s]['y']) for s in self.source_finders]))
            ylims = (0,10)
            injected = dict()
            for depth in self.image_types:
                metrics = self.catalogues['metrics']['rms_box_statistics'][depth]
                injected[depth] = {
                    'mean': metrics[(metrics['is_opt']>0)]['mean'].to(u.mJy).value,
                    'rms':  metrics[(metrics['is_opt']>0)]['rms' ].to(u.mJy).value,
                }
            metrics = self.catalogues['metrics']['parameters']
            for source_finder in self.source_finders:
                rms = metrics[((metrics['source_finder']==source_finder)&(metrics['image_type']=='shallow'))]['rms_parameter_value'][0]
                legend = [
                    {'ls': '-' , 'tt': '$1\sigma$ Noise' },
                    {'ls': '--', 'tt': '$3\sigma$ Noise' },
                    {'ls': '-.', 'tt': '$S_{Shallow}$=%.2f$\sigma_{S_{Deep}}$' % rms},
                    {'ls': ':',  'tt': '$S_{out}$=5$\sigma_{S_{in}}$'},
                ]
                plot_info.append({
                    'idx': f"dsfrp{source_finder}",
                    'type': 'scatter',
                    'filename': f"{source_finder}_shallow_deep_flux_ratio_vs_deep_{'sn' if is_sn else 'flux'}_plot.png",
                    'html': {
                        'title': f"{source_finder.capitalize()} Shallow-Deep Flux Ratio vs. Deep {'S/N' if is_sn else 'Flux'} Plot",
                        #'datum': QTable(), # TO-DO
                        'datum': real[source_finder]['d'],
                        'fname': f"{source_finder}_shallow_to_deep_flux_ratio_vs_{'sn' if is_sn else 'flux'}_plot_data_table.fits",
                    },
                    'data': {
                        'x': real[source_finder]['x'],
                        'y': real[source_finder]['y'],
                        'color': self.plot_colors[source_finder],
                    },
                    'curves': [
                        {'x': x_range, 'y': 1.0+1.0*(1.0 if is_sn else (injected['shallow']['mean']+injected['shallow']['rms']))/x_range, 'c': 'k', 'ls': '-'},
                        {'x': x_range, 'y': 1.0-1.0*(1.0 if is_sn else (injected['shallow']['mean']+injected['shallow']['rms']))/x_range, 'c': 'k', 'ls': '-'},
                        {'x': x_range, 'y': 1.0+3.0*(1.0 if is_sn else (injected['shallow']['mean']+injected['shallow']['rms']))/x_range, 'c': 'k', 'ls': '--'},
                        {'x': x_range, 'y': 1.0-3.0*(1.0 if is_sn else (injected['shallow']['mean']+injected['shallow']['rms']))/x_range, 'c': 'k', 'ls': '--'},
                        {'x': x_range, 'y': rms*(1.0 if is_sn else injected['shallow']['mean']+injected['shallow']['rms'])/x_range, 'c': 'k', 'ls': '-.'},
                        {'x': x_range, 'y': 5.0*(1.0 if is_sn else injected[depth]['mean']+injected[depth]['rms'])/x_range, 'c': 'k', 'ls': ':'},
                    ],
                    'xlim': xlims,
                    'xline': 1.0,
                    'ylim': ylims,
                    'xscale': 'log',
                    'xlabel': f"S/N (Deep-Signal/Shallow-Noise)" if is_sn else 'Deep Flux (mJy)',
                    'ylabel': '$S_{Shallow}/S_{Deep}$',
                    'title': f"{self.plot_labels[source_finder]} Shallow-Deep Flux Ratios Plot",
                    'legend': {
                        'title': f"{self.plot_labels[source_finder]}",
                        'handles': [mlines.Line2D([],[],color='k',label=d['tt'],linestyle=d['ls']) for d in legend],
                        'loc': 'upper right',
                    },
                })
                #break # debug
            fluxes = dict()
            for source_finder in self.source_finders:
                fluxes[source_finder] = comp[(comp['source_finder']==source_finder)]['fluxes']
            plot_info.append({
                'idx': f"dsfrfph",
                'type': 'hist',
                'filename': f"{'sn' if is_sn else 'flux'}_shallow_deep_false_positives_histogram.png",
                'html': {
                    'title': f"Shallow-Deep Flase-Positives <i>vs.</i> {'S/N' if is_sn else 'Flux'} Histogram",
                    'datum': QTable(), # TO-DO
                    'fname': f"shallow_deep_false_positive_flux_ratio_vs_{'sn' if is_sn else 'flux'}_histogram_data_table.fits",
                },
                'data': {
                    'n': {m:fluxes[m] for m in fluxes},
                    'color': {m:self.plot_colors[m] for m in fluxes},
                    'kwargs': {
                        'histtype': 'stepfilled',
                        ##'bins': np.logspace(np.log10(np.min(comp['fluxes'].value)),np.log10(np.max(comp['fluxes'].value)),self.default_bins),
                        #'bins': np.logspace(np.log10(np.min(list(comp['fluxes']))),np.log10(np.max(list(comp['fluxes']))),self.default_bins),
                        'bins': np.logspace(np.log10(max(1,np.min(list(comp['fluxes'])))),np.log10(np.max(list(comp['fluxes']))),self.default_bins),
                        'facecolor': 'None',
                        'linewidth': 2,
                    },
                 },
                'xscale': 'log',
                #'yscale': 'symlog',
                'xlabel': f"S/N (Shallow-Signal/Shallow-Noise)" if is_sn else "Total Flux (mJy)",
                'ylabel': 'False-Positives',
                'title': f"Shallow-Deep False-Positives Histogram",
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in fluxes],
                    'framealpha': 0.85,
                },
            })

        return plot_info


    def __get_subclump_size_distribution_table(self):
        # define some local pars
        qt = self.catalogues['cluster']
        modules = np.unique(qt['source_finder'])
        depths = np.unique(qt['image_type'])
        subclump_ids = np.unique(qt['subclump_id'])

        # get maxium clump size
        max_size = 0
        for module in modules:
            sf = qt[(qt['source_finder']==module)]
            for depth in depths:
                sfd = sf[(sf['image_type']==depth)]
                for subclump_id in subclump_ids:
                    s_size = len(sfd[(sfd['subclump_id']==subclump_id)])
                    if max_size < s_size:
                        max_size = s_size

        # build clump size distribution table
        dst = QTable()
        dst.add_column(
            name = 'counts',
            col = np.arange(1,max_size+1,dtype=np.int64)
        )
        for module in modules:
            sf = qt[(qt['source_finder']==module)]
            for depth in depths:
                col_name = f"{module}_{depth}"
                dst.add_column(
                    name = col_name,
                    col = np.zeros(len(dst),dtype=np.int64)
                )
                sfd = sf[(sf['image_type']==depth)]
                for subclump_id in subclump_ids:
                    s_size = len(sfd[(sfd['subclump_id']==subclump_id)])
                    for i in range(len(dst)):
                        if s_size == dst['counts'][i]:
                            dst[col_name][i] += s_size
                            break
        return dst


    def __make_subclump_hist(self):
        def reduce_aticks(aticks):
            n = 10
            if len(aticks) > n:
                a_min = np.int64(np.floor(min(aticks)))
                a_max = np.int64(np.ceil(max(aticks)))
                step = np.int64(np.ceil((a_max-a_min)/n))
                return [i for i in range(a_min,a_max+1,step)]
            return aticks

        # build clump size distribution table
        dst = self.__get_subclump_size_distribution_table()

        # create plot info 
        x_min = 0.5
        x_max = max(dst['counts'])+0.5
        y_min = 0
        idx = 1
        plot_info = list()
        modules = np.unique(self.catalogues['cluster']['source_finder'])
        for depth in np.unique(self.catalogues['cluster']['image_type']):
            y_max = 12*max([max(dst[f"{m}_{depth}"]) for m in modules])/10
            ords = dict()
            for module in modules:
                freq = dst[f"{module}_{depth}"]
                cnts = dst['counts']
                x = list()
                y = list()
                for i in range(len(freq)-1):
                    x.extend([cnts[i]-0.5, cnts[i]+0.5, cnts[i]+0.5])
                    y.extend([freq[i],     freq[i],     freq[i+1]  ])
                i = len(freq)-1
                x.extend([cnts[i]-0.5, cnts[i]+0.5])
                y.extend([freq[i],     freq[i]])
                ords[module] = {'x': x, 'y': y, 'color': self.plot_colors[module]}
            plot_info.append({
                'idx': f"fc{idx}",
                'type': 'line',
                'filename': f"clump_size_histogram_{depth}_image.png",
                'html': {
                    'title': f"Frequency <i>vs.</i> {depth.capitalize()} Clump Size Historgram",
                    'datum': QTable(),
                },
                'data': {
                    'x': {s:ords[s]['x'] for s in ords},
                    'y': {s:ords[s]['y'] for s in ords},
                    'color': {s:ords[s]['color'] for s in ords},
                },
                'xlim': [x_min,x_max],
                'xticks': reduce_aticks(dst['counts']),
                'ylim': [y_min,12.0*y_max/10.0],
                'yscale': 'symlog',
                'xlabel': 'Clump Size',
                'ylabel': 'Frequency',
                'title': f"{depth.capitalize()} Image Source Finder Clump Distributions",
                'legend': {
                    'handles': [mpatches.Patch(color=self.plot_colors[m],label=self.plot_labels[m]) for m in modules],
                    'framealpha': 0.85,
                    'loc': 'upper right',
                },
            })
            idx += 1
        return plot_info


    def __collect_plot_info(self,is_sn=True):
        plot_info = list()

        # typhon run stats plots
        plot_info.extend(self.__make_typhon_stats_plots())

        # create subclump histograms
        plot_info.extend(self.__make_subclump_hist())

        # histograms
        plot_info.extend(self.__make_hist_plots(is_sn=is_sn))

        # completeness step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_completeness_plot(is_sn=is_sn,catalogue=catalogue))
        plot_info.extend(self.__make_completeness_plot(is_sn=is_sn))

        # reliability step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_reliability_plot(is_sn=is_sn,catalogue=catalogue))
        plot_info.extend(self.__make_reliability_plot(is_sn=is_sn))

        # goodness of completeness step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_goodness_of_completeness_plot(is_sn=is_sn,catalogue=catalogue))

        # goodness reliability step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_goodness_of_reliability_plot(is_sn=is_sn,catalogue=catalogue))

        # delta-completeness step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_delta_completeness_plot(is_sn=is_sn,catalogue=catalogue))

        # delta-reliability step-plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_delta_reliability_plot(is_sn=is_sn,catalogue=catalogue))

        # flux-ratio plots
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_flux_ratio_plots(is_sn=is_sn,catalogue=catalogue))
        plot_info.extend(self.__make_flux_ratio_plots(is_sn=is_sn))

        return plot_info


    def plot_collection_test(self,is_sn=True):
        self.__plot(self.plot_collection)


    def make_typhon_stats_plots(self):
        plot_info = self.__make_typhon_stats_plots()
        self.__plot(plot_info)
        return self


    def make_hist_plots(self,is_sn=True):
        plot_info = self.__make_hist_plots(is_sn)
        self.__plot(plot_info)
        return self


    def make_completeness_plots(self,is_sn=True,is_injected=True):
        plot_info = list()
        if is_injected and 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_completeness_plot(is_sn=is_sn,catalogue=catalogue))
        else:
            plot_info.extend(self.__make_completeness_plot(is_sn=is_sn))
        self.__plot(plot_info)
        return self


    def make_reliability_plots(self,is_sn=True,is_injected=True):
        plot_info = list()
        if is_injected and 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_reliability_plot(is_sn=is_sn,catalogue=catalogue))
        else:
            plot_info.extend(self.__make_reliability_plot(is_sn=is_sn))
        self.__plot(plot_info)
        return self

    def make_goodness_of_completeness_plot(self,is_sn=True):
        plot_info = list()
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_goodness_of_completeness_plot(is_sn=is_sn,catalogue=catalogue))
        self.__plot(plot_info)
        return self


    def make_goodness_of_reliability_plot(self,is_sn=True):
        plot_info = list()
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_goodness_of_reliability_plot(is_sn=is_sn,catalogue=catalogue))
        self.__plot(plot_info)
        return self


    def make_delta_completeness_plot(self,is_sn=True):
        plot_info = list()
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_delta_completeness_plot(is_sn=is_sn,catalogue=catalogue))
        self.__plot(plot_info)
        return self


    def make_delta_reliability_plot(self,is_sn=True):
        plot_info = list()
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_delta_reliability_plot(is_sn=is_sn,catalogue=catalogue))
        self.__plot(plot_info)
        return self


    def make_flux_ratio_plots(self,is_sn=True,is_injected=True):
        plot_info = list()
        if is_injected and 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                plot_info.extend(self.__make_flux_ratio_plots(is_sn=is_sn,catalogue=catalogue))
        else:
            plot_info.extend(self.__make_flux_ratio_plots(is_sn=is_sn))
        self.__plot(plot_info)
        return self


    def print_cluster_table(self,verbose=False):
        print_ok("Cluster Table:")
        if verbose:
            self.catalogues['cluster'].pprint_all()
        else:
            print_ok(self.catalogues['cluster'])
        print_ok(f"[Done]")
        return self


    def print_clump_table(self,verbose=False):
        print_ok("Clump Table:")
        if verbose:
            self.catalogues['clump'].pprint_all()
        else:
            print_ok(self.catalogues['clump'])
        print_ok(f"[Done]")
        return self


    def print_metrics_table(self,verbose=False):
        print_ok("Global Metrics Table:")
        if verbose:
            self.catalogues['metrics']['parameters'].pprint_all()
        else:
            print_ok(self.catalogues['metrics']['parameters'])
        print_ok(f"[Done]")
        return self


    def save_catalogue(self,is_fits=True):
        if is_fits:
            fits_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".hydra.fits",self.deep_image)
            print_ok(f"Writing FITS file: {fits_file}")
            self.catalogues['cluster'].write(fits_file,format='fits',overwrite=True)
        else:
            csv_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".hydra.csv",self.deep_image)
            print_ok(f"Writing CSV file: {csv_file}")
            self.catalogues['cluster'].write(csv_file,format='csv',overwrite=True)
        print_ok(f"[Done]")
        return self


    def save_parameters(self,is_fits=True):
        if is_fits:
            fits_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".hydra.parameters.fits",self.deep_image)
            print_ok(f"Writing Source Finders Input Parameters FITS file: {fits_file}")
            self.catalogues['metrics']['parameters'].write(fits_file,format='fits',overwrite=True)
        else:
            csv_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".hydra.parameters.csv",self.deep_image)
            print_ok(f"Writing Source Finders Input Parameters CSV file: {csv_file}")
            self.catalogues['metrics']['parameters'].write(csv_file,format='csv',overwrite=True)
        print_ok(f"[Done]")
        return self


    def __create_cutouts(self):
        def mk_fname(i,depth,module=None):
            template = {
                'deep': re.sub(r"\.([Ff][Ii][Tt](|[Ss]))$",r".deep.\1",self.deep_image),
                'shallow': self.shallow_image
            }
            m = "" if module is None else f".residual.{module}"
            fname = re.sub(r"^(.*?/)*(.*?)\.(deep|shallow)\.[Ff][Ii][Tt](|[Ss])$",r"\2.hydra.clump_%07d.\3%s.fits" % (i,m),template[depth])
            return f"{cache}/{fname}"
        fits_data = list()
        # RE: TODO's 23 and 24 in cerberus.py
        def kludge(c):
            is_kludge = False
            if is_kludge:
                c['extent_semimajor'][((c['source_finder']=='aegean')|(c['source_finder']=='profound'))] /= 2
                c['extent_semiminor'][((c['source_finder']=='aegean')|(c['source_finder']=='profound'))] /= 2
            return c
        qt = self.catalogues['clump']['clump_id','ra','dec','size'].copy()
        #max_limit = 3.0*u.arcmin
        max_size = 12.0*max(qt['size'])/10.0
        opt_size = (max_size if max_size < self.max_limit else self.max_limit).to(u.arcmin)
        qt.add_column(
            name = 'is_oversized',
            col  = [s > self.max_limit for s in qt['size']]
        )
        qt['size'] = [opt_size if s < opt_size else s for s in qt['size']]
        for i in qt['clump_id']:
            manifest = dict()
            manifest['__main__'] = {d:mk_fname(i,d) for d in self.image_types}
            for module in self.residuals:
                manifest[module] = {d:mk_fname(i,d,module) for d in self.image_types}
            fits_data.append(manifest)

        print_ok(f"Creating Cutouts: OPT_SIZE={opt_size.to(u.arcmin)} (Limiting size: {self.max_limit.to(u.arcmin)})")
        plot_data = list()
        def has_components(clump_id,module,depth):
            if module == '__main__':
                return True
            #df =  self.catalogues['cluster']
            df =  kludge(self.catalogues['cluster'])
            return len(df[((df['clump_id']==clump_id)&(df['source_finder']==module)&(df['image_type']==depth))])>0
        def get_sexadecimal_string(position):
            sexadecimal = "%02d%02d%05.2f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%04.1f" % position.dec.signed_dms)
            return sexadecimal 
        cropper = {
            'deep': HydraCroppingTool(self.deep_image),
            'shallow': HydraCroppingTool(self.shallow_image)
        }
        def get_subclump_box(clump_id,match_id,depth):
            #qt = self.catalogues['cluster']
            qt = kludge(self.catalogues['cluster'])
            qt = qt[((qt['clump_id']==clump_id)&(qt['match_id']==match_id)&(qt['image_type']==depth))]
            ra_disc  = np.sqrt((qt['extent_semimajor']*np.sin(qt['extent_angle']))**2+(qt['extent_semiminor']*np.cos(qt['extent_angle']))**2)
            dec_disc = np.sqrt((qt['extent_semimajor']*np.cos(qt['extent_angle']))**2+(qt['extent_semiminor']*np.sin(qt['extent_angle']))**2)
            ra_min  = qt['ra']  - ra_disc
            ra_max  = qt['ra']  + ra_disc
            dec_min = qt['dec'] - dec_disc
            dec_max = qt['dec'] + dec_disc
            def get_centroid():
                bb_ra_min  = min(min(ra_min),min(ra_max))
                bb_ra_max  = max(max(ra_min),max(ra_max))
                bb_dec_min = min(min(dec_min),min(dec_max))
                bb_dec_max = max(max(dec_min),max(dec_max))
                ra  = (bb_ra_min+bb_ra_max)/2.0
                dec = (bb_dec_min+bb_dec_max)/2.0
                return SkyCoord(ra,dec)
            def get_width():
                wt_ra_min  = min(qt['ra']  - ra_disc)
                wt_ra_max  = max(qt['ra']  + ra_disc)
                return abs(wt_ra_max-wt_ra_min)
            def get_height():
                ht_dec_min = min(qt['dec'] - dec_disc)
                ht_dec_max = max(qt['dec'] + dec_disc)
                return abs(ht_dec_max-ht_dec_min)
            return {'postn': get_centroid(), 'width': get_width(), 'height': get_height()}
        def get_match_ids(clump_id,depth):
            #qt = self.catalogues['cluster']
            qt = kludge(self.catalogues['cluster'])
            qt = qt[((qt['clump_id']==clump_id)&(qt['image_type']==depth))]
            return qt['match_id']
        def get_fits_data_count():
            cnt = 0
            for i in range(len(fits_data)):
                clump_id = qt['clump_id'][i]
                #clump = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
                clump = kludge(self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)])
                for module in fits_data[i]:
                    for depth in fits_data[i][module]:
                        if not has_components(clump_id,module,depth):
                            continue
                        for is_annotated in [False,True]:
                            if is_annotated:
                                is_ok = False
                                datum = clump[(clump['image_type']==depth)]
                                if len(datum)>0:
                                    for source_finder in np.unique(datum['source_finder']):
                                        df = datum[(datum['source_finder']==source_finder)]
                                        if len(df)>0:
                                            is_ok = True
                            else:
                                is_ok = True
                            cnt += 1 if is_ok else 0
            return cnt
        is_use_image = True
        def push(tasks):
            if len(tasks)>0:
                if sys.platform == 'darwin':
                    pool = multiprocessing.Pool()
                else:
                    pool = multiprocessing.Pool(25)
                print_ok(">>      * * *   F L U S H I N G   S T A C K   * * *")
                pool.map(render,tasks.copy())
                print_ok(">>      * * *   < D O N E >   * * *")
                pool.close()
                gc.collect()
        cnt = 1
        tasks = list() 
        fits_files = list()
        fits_data_count = get_fits_data_count()
        for i in range(len(fits_data)):
            clump_id = qt['clump_id'][i]
            ra   = qt['ra'][i]
            dec  = qt['dec'][i]
            size = qt['size'][i].to(u.arcmin)
            #clump = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
            clump = kludge(self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)])
            #compl = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']!=clump_id)]
            compl = kludge(self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']!=clump_id)])
            ra_min  = (ra-size).to(u.deg)
            ra_max  = (ra+size).to(u.deg)
            dec_min = (dec-size).to(u.deg)
            dec_max = (dec+size).to(u.deg)
            compl = compl[((ra_min<compl['ra'])&(compl['ra']<=ra_max)&(dec_min<=compl['dec'])&(compl['dec']<=dec_max))]
            for module in fits_data[i]:
                for depth in fits_data[i][module]:
                    if not has_components(clump_id,module,depth):
                        continue
                    fname = fits_data[i][module][depth]
                    print_ok(f"> Creating[{cnt}/{fits_data_count}]: {fname}")
                    if is_use_image:
                        if module == '__main__':
                            image = cropper[depth].crop_image(ra,dec,size)
                        else:
                            image = self.residuals[module][depth].crop_image(ra,dec,size)
                    else:
                        if module == '__main__':
                            cropper[depth].crop_image(ra,dec,size,fname)
                        else:
                            self.residuals[module][depth].crop_image(ra,dec,size,fname)
                        fits_files.append(fname)
                        image = None
                    datum = clump[(clump['image_type']==depth)]
                    # png plots
                    title = list()
                    title.append(f"[Clump_ID:{clump_id}]")
                    title.append(f"J{get_sexadecimal_string(SkyCoord(ra,dec))}:")
                    if module != '__main__':
                        title.append(f"Residual={self.plot_labels[module]},")
                    title.append(f"Type={depth.capitalize()},")
                    title.append(f"Side={size.value:.2}'")
                    title = " ".join(title)
                    for is_annotated in [False,True]:
                        task = {
                            'plt_title': title,
                            'plt_colors': self.plot_colors,
                            'plt_labels': self.plot_labels,
                            'clump_id': clump_id,
                            'module': module,
                            'depth': depth,
                            'fname': fname,
                            'counter': {'cnt': cnt, 'max': fits_data_count},
                        }
                        if is_annotated:
                            is_ok = False
                            source_finders = list()
                            if len(datum)>0:
                                for source_finder in np.unique(datum['source_finder']):
                                    df = datum[(datum['source_finder']==source_finder)]
                                    if len(df)>0:
                                        source_finders.append(source_finder)
                                        is_ok = True
                            if is_ok:
                                annotation = {
                                    'clump': datum.copy(),
                                    'bpars': list(),
                                }
                                for match_id in get_match_ids(clump_id,depth):
                                    bbox = get_subclump_box(clump_id,match_id,depth) 
                                    annotation['bpars'].append({'bbox': bbox,'match_id':match_id})
                                df = compl[(compl['image_type']==depth)]
                                annotation['compl'] = df.copy()
                                hdls = [mpatches.Patch(color=self.plot_colors[s],label=self.plot_labels[s]) for s in source_finders]
                                hdls.append(mlines.Line2D([],[],color='red',marker='P',linestyle='None',label='Centroid',markerfacecolor='None',markersize=7.5))
                                annotation['positions'] = SkyCoord(ra,dec)
                                annotation['legend_hdls'] = hdls
                                task['annotation'] = annotation
                            pname = re.sub(r"\.fits$",".annotated.png",fname)
                        else:
                            is_ok = True
                            pname = re.sub(r"\.fits$",".png",fname)
                        if is_ok:
                            plot_data.append(pname)
                            task['pname'] = pname
                            if is_use_image:
                                task['image'] = image.copy()
                            tasks.append(task)
                        if (cnt % 1000) == 0:
                            push(tasks)
                            tasks = list()
                        cnt += 1 if is_ok else 0
        push(tasks)
        for fname in fits_files:
            if os.path.isfile(fname):
                os.remove(fname)
        print_ok(f"> [{cnt-1}/{fits_data_count}]")
        print_ok("[Done]")

        return plot_data


    def __html(self,manifest):
        # create html libs list
        html_libs_root_dir = os.path.abspath(get_this_source_file_directory()+"/modules/libs/html")
        archive_root_dir   = re.sub(r"\.tar\.gz$","_dir",re.sub(r"^(.*?/)+","",self.tar_gz_file))
        html_meta = [{
                'kind': 'css',
                'keep': True,
                'system':  f"{html_libs_root_dir}/extern/bootstrap/4.3.1/css/bootstrap.min.css",
                'archive': f"{archive_root_dir}/html/extern/bootstrap/4.3.1/css/bootstrap.min.css",
            },{
                'kind': 'js',
                'keep': True,
                'system':  f"{html_libs_root_dir}/extern/jquery/3.4.1/jquery-3.4.1.min.js",
                'archive': f"{archive_root_dir}/html/extern/jquery/3.4.1/jquery-3.4.1.min.js",
            },{
                'kind': 'js',
                'keep': True,
                'system':  f"{html_libs_root_dir}/extern/bootstrap/4.3.1/js/bootstrap.min.js",
                'archive': f"{archive_root_dir}/html/extern/bootstrap/4.3.1/js/bootstrap.min.js",
        },]

        # create manifest.js html lib contents
        html_file = "manifest.js"
        def get_cutout_stack():
            def get_id(fname):
                return int(re.sub("^.*?(hydra\.clump_(\d+)\.).*$",r"\2",fname))
            def get_field_width(field):
                size = 0
                for datum in manifest['cutouts']:
                    for item in datum[field]:
                        if size < len(item):
                            size = len(item)
                return size
            qt = Table(
                names = ('clump_id','module','depth','is_annotated','cutout'),
                dtype = (
                    np.int64,                         # clump_id
                    f"S{get_field_width('module')}",  # module
                    f"S{get_field_width('depth')}",   # depth
                    np.int64,                         # is_annotated
                    f"S{get_field_width('cutouts')}", # cutout
                )
            )
            for datum in manifest['cutouts']:
                module       = datum['module']
                depth        = datum['depth']
                is_annotated = 1 if datum['is_annotated'] else 0
                for cutout in datum['cutouts']:
                    clump_id = get_id(cutout)
                    qt.add_row((clump_id,module,depth,is_annotated,cutout))
            qt.sort(['clump_id','module','depth','is_annotated'])
            return qt
        def round_sig(x, sig=3):
            try:
                return 0 if x == 0 else np.round(x, sig-np.int64(np.floor(np.log10(np.abs(x))))-1)
            except Exception as e:
                # Notes: I had the following lines,
                #    print_warning(f"ERROR: {e}")
                #    print_warning(f" ==> {x}")
                # and captured two errors, 
                #    hydra.py:1842: RuntimeWarning: overflow encountered in long_scalars
                #      return 0 if x == 0 else np.round(x, sig-np.int64(np.floor(np.log10(np.abs(x))))-1)
                #    ERROR: signed integer is less than minimum
                #     ==> nan
                #    ERROR: signed integer is less than minimum
                #     ==> nan
                # for clump_id=1022, re.,
                #     emu_simulated_2x2.hydra.pybdsf.deep.residual.hydra.tar.gz
                # I suspect this is an underflow error, and so I should really return -1000, but left it
                # a nan (coverting it to NaN for the javascript interpreter.
                #
                # TO-DO: Perhaps?
                #        Investigate further and determine sutible return value. (Would need to capture,
                #        appropriate exception.)
                pass
            return x
        def get_params(clump_id,module,depth):
            qt = self.catalogues['cluster']
            qt = qt[((qt['clump_id']==clump_id)&(qt['source_finder']==module)&(qt['image_type']==depth))]
            params = dict()
            params['n'] = len(qt['source_finder'])
            params['rms']   = round_sig(qt['residual_rms'][0].value)
            params['madfm'] = round_sig(qt['residual_madfm'][0].value)
            params['sumsq'] = round_sig(qt['residual_sumsq'][0].value)
            params['size']  = round_sig(self.catalogues['clump'][(self.catalogues['clump']['clump_id']==clump_id)]['size'].value)
            return params
        str_ref_id        = 'Ref ID'
        str_source_finder = 'Source'
        str_depth         = 'Image'
        str_sn_bane       = 'S/N'
        def get_cluster_table(clump_id,indent=3):
            qt = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
            header = {
                'ref_id': str_ref_id,
                'cat_id': 'Catalogue ID',
                'clump_id': 'Clump ID',
                'subclump_id': 'SubClump ID',
                'match_id': 'Match ID',
                'source_finder': str_source_finder,
                'image_type': str_depth,
                'ra': 'RA&nbsp;(&deg;)',
                'dec': 'Dec&nbsp;(&deg;)',
                'extent_semimajor': 'a<sub>extent</sub>&nbsp;(&Prime;)',
                'extent_semiminor': 'b<sub>extent</sub>&nbsp;(&Prime;)',
                'extent_angle': '&theta;<sub>extent</sub>&nbsp;(&deg;)',
                'flux_total': 'Total&nbsp;Flux<br> mJy',
                'rms_noise_bane': 'RMS&nbsp;Noise<br>mJy',
                'sn_bane': str_sn_bane,
                'flux_peak': 'Peak&nbsp;Flux<br> mJy',
                'residual_rms': 'Residual&nbsp;RMS<br>mJy/(arcmin<sup>2</sup>beam)',
                'residual_madfm': 'Residual&nbsp;MADFM<br>mJy/(arcmin<sup>2</sup>beam)',
                'residual_sumsq': 'Residual&nbsp;SumSq<br>mJy<sup>2</sup>/(arcmin<sup>2</sup>beam)<sup>2</sup>',
            }
            data = Table()
            for h in header:
                data.add_column(qt[h])
            data['extent_semimajor'] = data['extent_semimajor'].to(u.arcsec)
            data['extent_semiminor'] = data['extent_semiminor'].to(u.arcsec)
            data.sort(['match_id','source_finder','image_type','ra','dec'])
            html = list()
            html.append("cluster: {")
            html.append(f"   header: {list(header.values())},")
            html.append("   data: [")
            for i in range(len(data['clump_id'])):
                html.append(re.sub(r" nan,"," NaN,",f"      {list(data[i])},"))
            html.append("   ],")
            html.append("},")
            return [re.sub(r"^ ","",f"{indent*'   '}{h}") for h in html]
        html = list()
        qt = get_cutout_stack()
        html.append("get_cutout_stack = function() {")
        html.append("   var cutout_stack = [")
        for clump_id in np.unique(qt['clump_id']):
            clump = qt[(qt['clump_id']==clump_id)]
            html.append("       {clump_id: %d," % clump_id)
            html.append("        modules: {")
            for module in np.unique(clump['module']):
                datum = clump[(clump['module']==module)]
                html.append("            %s: {" % module)
                for depth in np.unique(datum['depth']):
                    html.append("                %s: {" % depth)
                    d = datum[(datum['depth']==depth)]
                    for c in d:
                        kind = 'annotated' if c['is_annotated'] else 'plain'
                        cname = '"'+c['cutout']+'"'
                        html.append("                    %s: %s," % (kind,cname))
                    if module != "__main__":
                        params = get_params(clump_id,module,depth)
                        html.append("                    params: {")
                        html.append("                       n: %d," % params['n'])
                        html.append("                       rms: %s,"   % re.sub("nan","NaN",re.sub(r"\.$","",re.sub(r"0+$","","%f" % params['rms']))))
                        html.append("                       madfm: %s," % re.sub("nan","NaN",re.sub(r"\.$","",re.sub(r"0+$","","%f" % params['madfm']))))
                        html.append("                       sumsq: %s," % re.sub("nan","NaN",re.sub(r"\.$","",re.sub(r"0+$","","%f" % params['sumsq']))))
                        html.append("                       size: %s,"  % re.sub("nan","NaN",re.sub(r"\.$","",re.sub(r"0+$","","%f" % params['size']))))
                        html.append("                    },")
                    html.append("                },")
                html.append("            },")
            html.append("        },")
            html.extend(get_cluster_table(clump_id))
            html.append("       },")
        html.append("   ]")
        html.append("")
        def fetch_plot_info(at_idx):
            for item in manifest['plot_info']: 
                if item['idx']==at_idx:
                    return item['html']['datum']
            return None
        ds_info = fetch_plot_info('dscs1')
        html.append("   var meta = {")
        html.append("      deep_shallow_completeness: {")
        html.append("         bin_info: {")
        html.append("            range: [0,%d]," % (self.default_bins-1))
        html.append("         },")
        html.append("         datum: [")
        for bin_no in np.unique(ds_info['bin_no']):
            b_cnts = ds_info[(ds_info['bin_no']==bin_no)]
            clumps_str = list()
            for cid in np.unique(b_cnts['clump_id_deep']):
                clp = b_cnts[(b_cnts['clump_id_deep']==cid)]
                lcl = list()
                for rid in clp['ref_id_deep']:
                    s_info = clp[(clp['ref_id_deep']==rid)]
                    lcl.append("{"+f"{s_info['source_finder'][0]}:"+"{"+f"rid:[{rid}{(',%s' % s_info['ref_id_shallow'][0]) if s_info['ref_id_shallow'][0]>-1 else ''}],com:{s_info['completeness'][0]}"+"}}")
                clumps_str.append(f"{cid}:[{','.join(lcl)}]")
            clumps_str = "{"+",".join(clumps_str)+"}"
            html.append("            {bin_no: %d," % bin_no)
            html.append("             sn_min: %s," % f"{b_cnts['sn_bin_min'][0]}")
            html.append("             sn_avg: %s," % f"{b_cnts['sn_avg'][0]}")
            html.append("             sn_max: %s," % f"{b_cnts['sn_bin_max'][0]}")
            html.append("             deep_clump_ids: [%s]," % ",".join(np.unique([f"{clump_id}" for clump_id in b_cnts['clump_id_deep']])))
            html.append("             shallow_clump_ids: [%s]," % ",".join(np.unique([f"{clump_id}" for clump_id in filter(lambda x: x>-1,b_cnts['clump_id_shallow'])]))) 
            html.append("             clumps: %s," % clumps_str)
            html.append("            },")
        html.append("         ],")
        html.append("      },")
        html.append("      deep_shallow_reliability: {")
        html.append("         bin_info: {")
        html.append("            range: [0,%d]," % (self.default_bins-1))
        html.append("         },")
        html.append("         datum: [")
        ds_info = fetch_plot_info('dsrs2')
        for bin_no in np.unique(ds_info['bin_no']):
            b_cnts = ds_info[(ds_info['bin_no']==bin_no)]
            clumps_str = list()
            for cid in np.unique(b_cnts['clump_id_shallow']):
                clp = b_cnts[(b_cnts['clump_id_shallow']==cid)]
                lcl = list()
                for rid in clp['ref_id_shallow']:
                    s_info = clp[(clp['ref_id_shallow']==rid)]
                    lcl.append("{"+f"{s_info['source_finder'][0]}:"+"{"+f"rid:[{rid}{(',%s' % s_info['ref_id_deep'][0]) if s_info['ref_id_deep'][0]>-1 else ''}],rel:{s_info['reliability'][0]}"+"}}")
                clumps_str.append(f"{cid}:[{','.join(lcl)}]")
            clumps_str = "{"+",".join(clumps_str)+"}"
            html.append("            {bin_no: %d," % bin_no)
            html.append("             sn_min: %s," % f"{b_cnts['sn_bin_min'][0]}")
            html.append("             sn_avg: %s," % f"{b_cnts['sn_avg'][0]}")
            html.append("             sn_max: %s," % f"{b_cnts['sn_bin_max'][0]}")
            html.append("             shallow_clump_ids: [%s]," % ",".join(np.unique([f"{clump_id}" for clump_id in b_cnts['clump_id_shallow']])))
            html.append("             deep_clump_ids: [%s]," % ",".join(np.unique([f"{clump_id}" for clump_id in filter(lambda x: x>-1,b_cnts['clump_id_deep'])]))) 
            html.append("             clumps: %s," % clumps_str)
            html.append("            },")
        html.append("         ],")
        html.append("      },")
        clump_sizes = list()
        modules = np.unique(self.catalogues['cluster']['source_finder'])
        for clump_id in np.unique(self.catalogues['cluster']['clump_id']):
            clump = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
            for module in modules:
                sf_clump = clump[(clump['source_finder']==module)]
                for depth in self.image_types:
                    clump_sizes.append(len(sf_clump[(sf_clump['image_type']==depth)]))
        clump_sizes = list(np.unique(np.array(clump_sizes)[(np.array(clump_sizes)>0)]))
        clump_info = list([{'clump_size': cs, 'clump_ids': list()} for cs in clump_sizes])
        for clump_id in np.unique(self.catalogues['cluster']['clump_id']):
            clump = self.catalogues['cluster'][(self.catalogues['cluster']['clump_id']==clump_id)]
            for module in modules:
                sf_clump = clump[(clump['source_finder']==module)]
                for depth in self.image_types:
                    sfd_clump = sf_clump[(sf_clump['image_type']==depth)]
                    n = len(sfd_clump)
                    if n > 0:
                        for i in range(len(clump_info)):
                            if clump_info[i]['clump_size']==n:
                                clump_info[i]['clump_ids'].extend(sf_clump['clump_id'])
                                break
        for i in range(len(clump_info)):
            clump_info[i]['clump_ids'] = list(np.unique(clump_info[i]['clump_ids']))
        html.append("      clumps: [")
        for clump in clump_info:
            html.append("         {clump_size: %d," % clump['clump_size'])
            html.append("          clump_ids: "+f"{clump['clump_ids']},")
            html.append("         },")
        html.append("      ],")
        html.append("   }")
        html.append("")
        html.append("   // ok, let's have '__main__' show up even if there are no annotations, by using the unannoted (plain) image...")
        html.append("   for (var i in cutout_stack) {")
        html.append("      if (!('annotated' in cutout_stack[i]['modules']['__main__']['deep'])) {")
        html.append("         cutout_stack[i]['modules']['__main__']['deep']['annotated']=cutout_stack[i]['modules']['__main__']['deep']['plain']")
        html.append("      }")
        html.append("      if (!('annotated' in cutout_stack[i]['modules']['__main__']['shallow'])) {")
        html.append("         cutout_stack[i]['modules']['__main__']['shallow']['annotated']=cutout_stack[i]['modules']['__main__']['shallow']['plain']")
        html.append("      }")
        html.append("   }")
        html.append("")
        html.append("   // define deep-shallow completeness mode pars")
        html.append("   var dsc_sns = []")
        html.append("   for (var i in meta.deep_shallow_completeness.datum) {")
        html.append("      dsc_sns.push(meta.deep_shallow_completeness.datum[i].sn_avg)")
        html.append("   }")
        html.append("   var dsc_bins = []")
        html.append("   for (var i in meta.deep_shallow_completeness.datum) {")
        html.append("      dsc_bins.push(meta.deep_shallow_completeness.datum[i].bin_no)")
        html.append("   }")
        html.append("")
        html.append("   // define deep-shallow reliability mode pars")
        html.append("   var dsr_sns = []")
        html.append("   for (var i in meta.deep_shallow_reliability.datum) {")
        html.append("      dsr_sns.push(meta.deep_shallow_reliability.datum[i].sn_avg)")
        html.append("   }")
        html.append("   var dsr_bins = []")
        html.append("   for (var i in meta.deep_shallow_reliability.datum) {")
        html.append("      dsr_bins.push(meta.deep_shallow_reliability.datum[i].bin_no)")
        html.append("   }")
        html.append("")
        html.append("   // define clump mode pars")
        html.append("   var clump_sizes = []")
        html.append("   for (var i in meta.clumps) {")
        html.append("      clump_sizes.push(meta.clumps[i].clump_size)")
        html.append("   }")
        html.append("")
        html.append("   var cutout_stack_counter = 0")
        html.append("   var cutout_stack_counter_min = 0")
        html.append("   var cutout_stack_counter_max = cutout_stack.length-1")
        html.append("   var current_stack = cutout_stack")
        html.append("")
        html.append("   return {")
        html.append("      set_mode_clump: function() {")
        html.append("         cutout_stack_counter = 0")
        html.append("         cutout_stack_counter_min = 0")
        html.append("         cutout_stack_counter_max = cutout_stack.length-1")
        html.append("         current_stack = cutout_stack")
        html.append("      },")
        html.append("      get_deep_shallow_completeness_sns: function() {")
        html.append("         return dsc_sns")
        html.append("      },")
        html.append("      get_deep_shallow_completeness_bins: function() {")
        html.append("         return dsc_bins")
        html.append("      },")
        html.append("      set_deep_shallow_completeness_bin: function(bin_no) {")
        html.append("         var clumps = []")
        html.append("         var clump_ids = []")
        html.append("         for (var i in meta.deep_shallow_completeness.datum) {")
        html.append("            if (meta.deep_shallow_completeness.datum[i].bin_no==bin_no) {")
        html.append("               clump_ids.push(...meta.deep_shallow_completeness.datum[i].deep_clump_ids)")
        html.append("               clump_ids.push(...meta.deep_shallow_completeness.datum[i].shallow_clump_ids)")
        html.append("               clump_ids = clump_ids.filter(function(value,index,self) {")
        html.append("                  return self.indexOf(value) === index")
        html.append("               }).sort()")
        html.append("               this.set_mode_clump()")
        html.append("               for (var i in clump_ids) {")
        html.append("                  var clump = this.get_clump(clump_ids[i])")
        html.append("                  clump.bin_no = bin_no")
        html.append("                  clump.mode = 'completeness'")
        html.append("                  clumps.push(clump)")
        html.append("               }")
        html.append("               cutout_stack_counter = 0")
        html.append("               cutout_stack_counter_min = 0")
        html.append("               cutout_stack_counter_max = clumps.length-1")
        html.append("               current_stack = clumps")
        html.append("               break")
        html.append("            }")
        html.append("         }")
        html.append("      },")
        html.append("      get_deep_shallow_reliability_sns: function() {")
        html.append("         return dsr_sns")
        html.append("      },")
        html.append("      get_deep_shallow_reliability_bins: function() {")
        html.append("         return dsr_bins")
        html.append("      },")
        html.append("      set_deep_shallow_reliability_bin: function(bin_no) {")
        html.append("         var clumps = []")
        html.append("         var clump_ids = []")
        html.append("         for (var i in meta.deep_shallow_reliability.datum) {")
        html.append("            if (meta.deep_shallow_reliability.datum[i].bin_no==bin_no) {")
        html.append("               clump_ids.push(...meta.deep_shallow_reliability.datum[i].deep_clump_ids)")
        html.append("               clump_ids.push(...meta.deep_shallow_reliability.datum[i].shallow_clump_ids)")
        html.append("               clump_ids = clump_ids.filter(function(value,index,self) {")
        html.append("                  return self.indexOf(value) === index")
        html.append("               }).sort()")
        html.append("               this.set_mode_clump()")
        html.append("               for (var i in clump_ids) {")
        html.append("                  var clump = this.get_clump(clump_ids[i])")
        html.append("                  clump.bin_no = bin_no")
        html.append("                  clump.mode = 'reliability'")
        html.append("                  clumps.push(clump)")
        html.append("               }")
        html.append("               cutout_stack_counter = 0")
        html.append("               cutout_stack_counter_min = 0")
        html.append("               cutout_stack_counter_max = clumps.length-1")
        html.append("               current_stack = clumps")
        html.append("               break")
        html.append("            }")
        html.append("         }")
        html.append("      },")
        html.append("      get_clump_sizes: function() {")
        html.append("         return clump_sizes")
        html.append("      },")
        html.append("      set_clump_size: function(clump_size) {")
        html.append("         var clumps = []")
        html.append("         var clump_ids = []")
        html.append("         for (var i in meta.clumps) {")
        html.append("            if (meta.clumps[i].clump_size==clump_size) {")
        html.append("               clump_ids.push(...meta.clumps[i].clump_ids)")
        html.append("               clump_ids = clump_ids.filter(function(value,index,self) {")
        html.append("                  return self.indexOf(value) === index")
        html.append("               }).sort()")
        html.append("               this.set_mode_clump()")
        html.append("               for (var i in clump_ids) {")
        html.append("                  var clump = this.get_clump(clump_ids[i])")
        html.append("                  clump.clump_size = clump_size")
        html.append("                  clump.mode = 'clump'")
        html.append("                  clumps.push(clump)")
        html.append("               }")
        html.append("               cutout_stack_counter = 0")
        html.append("               cutout_stack_counter_min = 0")
        html.append("               cutout_stack_counter_max = clumps.length-1")
        html.append("               current_stack = clumps")
        html.append("               break")
        html.append("            }")
        html.append("         }")
        html.append("      },")
        html.append("      get_source_finders: function() {")
        html.append("         var clump = this.get_clump()")
        html.append("         if (clump.bin_no !== undefined && clump.mode !== undefined) {")
        html.append("            var mode   = clump.mode")
        html.append("            var bin_no = clump.bin_no")
        html.append("            var datum = undefined")
        html.append("            if (mode == 'completeness') {")
        html.append("               datum = meta.deep_shallow_completeness.datum")
        html.append("            } else if (mode == 'reliability') {")
        html.append("               datum = meta.deep_shallow_reliability.datum")
        html.append("            }")
        html.append("            if (datum !== undefined) {")
        html.append("               var matches = undefined")
        html.append("               for (var i in datum) {")
        html.append("                   if (datum[i].bin_no == bin_no) {")
        html.append("                      matches = datum[i].clumps[clump.clump_id]")
        html.append("                      break")
        html.append("                   }")
        html.append("               }")
        html.append("               if (matches !== undefined) {")
        html.append("                  var deep = []")
        html.append("                  var shallow = []")
        html.append("                  var ref_ids = []")
        html.append("                  for (var i in matches) {")
        html.append("                     for (var source_finder in matches[i]) {")
        html.append("                        if (mode == 'completeness') {")
        html.append("                           if (!deep.includes(source_finder)) {")
        html.append("                              deep.push(source_finder)")
        html.append("                              if (!shallow.includes(source_finder) && matches[i][source_finder].rid.length>1) {")
        html.append("                                 shallow.push(source_finder)")
        html.append("                              }")
        html.append("                           }")
        html.append("                        } else {")
        html.append("                           if (!shallow.includes(source_finder)) {")
        html.append("                              shallow.push(source_finder)")
        html.append("                              if (!deep.includes(source_finder) && matches[i][source_finder].rid.length>1) {")
        html.append("                                 deep.push(source_finder)")
        html.append("                              }")
        html.append("                           }")
        html.append("                        }")
        html.append("                        var cr_value = (mode=='completeness') ? matches[i][source_finder].com : matches[i][source_finder].rel")
        html.append("                        ref_ids.push([matches[i][source_finder].rid[0],cr_value])")
        html.append("                        if (matches[i][source_finder].rid.length>1) {")
        html.append("                           ref_ids.push([matches[i][source_finder].rid[1],cr_value])")
        html.append("                        }")
        html.append("                     }")
        html.append("                  }")
        html.append("                  return {")
        html.append("                     deep:    (deep.length    > 0) ? deep    : undefined,")
        html.append("                     shallow: (shallow.length > 0) ? shallow : undefined,")
        html.append("                     ref_ids: (ref_ids.length > 0) ? ref_ids : undefined,")
        html.append("                  }")
        html.append("               }")
        html.append("            }")
        html.append("         }")
        html.append("         return undefined")
        html.append("      },")
        html.append("      get_min_clump_id: function() {")
        html.append("         return current_stack[cutout_stack_counter_min].clump_id")
        html.append("      },")
        html.append("      get_current_clump_id: function() {")
        html.append("         return current_stack[cutout_stack_counter].clump_id")
        html.append("      },")
        html.append("      get_max_clump_id: function() {")
        html.append("         return current_stack[cutout_stack_counter_max].clump_id")
        html.append("      },")
        html.append("      get_abs_max_clump_id: function() {")
        html.append("         return cutout_stack[cutout_stack.length-1].clump_id")
        html.append("      },")
        html.append("      get_no_clumps: function() {")
        html.append("         return current_stack.length")
        html.append("      },")
        html.append("      get_next_clump_id: function() {")
        html.append("         if (cutout_stack_counter < cutout_stack_counter_max) {")
        html.append("            cutout_stack_counter += 1")
        html.append("         }")
        html.append("         return current_stack[cutout_stack_counter].clump_id")
        html.append("      },")
        html.append("      get_previous_clump_id: function() {")
        html.append("         if (cutout_stack_counter > cutout_stack_counter_min) {")
        html.append("            cutout_stack_counter -= 1")
        html.append("         }")
        html.append("         return current_stack[cutout_stack_counter].clump_id")
        html.append("      },")
        html.append("      is_first: function() {")
        html.append("         return cutout_stack_counter == cutout_stack_counter_min")
        html.append("      },")
        html.append("      is_last: function() {")
        html.append("         return cutout_stack_counter == cutout_stack_counter_max")
        html.append("      },")
        html.append("      set_clump: function(clump_id) {")
        html.append("         if (clump_id == current_stack[cutout_stack_counter].clump_id) {")
        html.append("            return current_stack[cutout_stack_counter]")
        html.append("         } else {")
        html.append("            for (var i=cutout_stack_counter_min;i<=cutout_stack_counter_max;i++) {")
        html.append("               if (current_stack[i].clump_id==clump_id) {")
        html.append("                  cutout_stack_counter = i")
        html.append("                  return current_stack[cutout_stack_counter]")
        html.append("               }")
        html.append("            }")
        html.append("         }")
        html.append("         return undefined")
        html.append("      },")
        html.append("      get_clump: function(clump_id) {")
        html.append("         if (typeof clump_id !== 'undefined') {")
        html.append("            for (var i in current_stack) {")
        html.append("               if (current_stack[i]['clump_id']==clump_id) {")
        html.append("                  return current_stack[i]")
        html.append("               }")
        html.append("            }")
        html.append("            return undefined")
        html.append("         }")
        html.append("         return current_stack[cutout_stack_counter]")
        html.append("      },")
        html.append("      get_cluster_table: function(clump_id) {")
        html.append("         var clump = ((typeof clump_id !== 'undefined') ? this.get_clump(clump_id) : current_stack[cutout_stack_counter]).cluster")
        html.append("         var idx_ref = clump['header'].findIndex(function(s){return s == '%s'})" % str_ref_id)
        html.append("         var idx_depth = clump['header'].findIndex(function(s){return s == '%s'})" % str_depth)
        html.append("         var table = ''")
        html.append("         table += '<table border=\"1\" cellpadding=\"5px\" class=\"cluster-info-table\">\\n'")
        html.append("         table += '   <thead>\\n'")
        html.append("         table += '      <tr>\\n'")
        html.append("         table += '         <th class=\"mode-title collapse\"></th>\\n'")
        html.append("         for (var i in clump['header']) {")
        html.append("            table += '         <th>'+clump['header'][i]+'</th>\\n'")
        html.append("         }")
        html.append("         table += '      </tr>\\n'")
        html.append("         table += '   </thead>\\n'")
        html.append("         table += '   <tbody>\\n'")
        html.append("         for (var i in clump['data']) {")
        html.append("            table += '      <tr class=\"ref-id-'+clump['data'][i][idx_ref]+' '+clump['data'][i][idx_depth]+'\">\\n'")
        html.append("            row = ''")
        html.append("            row += '<td class=\"mode-data collapse\"></td>'")
        html.append("            for (var j in clump['data'][i]) {")
        html.append("               if (clump['data'][i][j] !== -1000) {")
        html.append("                  row += '<td>'+clump['data'][i][j]+'</td>'")
        html.append("               } else {")
        html.append("                  row += '<td>&mdash;</td>'")
        html.append("               }")
        html.append("            }")
        html.append("            table += '         '+row+'\\n'")
        html.append("            table += '      </tr>\\n'")
        html.append("         }")
        html.append("         table += '   </tbody>\\n'")
        html.append("         table += '</table>\\n'")
        html.append("         return table")
        html.append("      },")
        html.append("   }")
        html.append("}")
        html_meta.append({
            'kind': 'js',
            'keep': False,
            'system': f"{cache}/{html_file}",
            'archive': f"{archive_root_dir}/html/local/js/{html_file}",
            'contents': html,
        })

        # plot table builder
        def fetch_plot(idx):
            for plot in manifest['plot_info']:
                if plot['idx'] == idx:
                    return plot
            return None 
        def layout_plots(indent=0):
            plot_layout = [
                [  'fc1',  'fc2'   ],
                [  'dnh1', 'snh2'  ],
                [  'dcs1', 'scs2'  ],
                [  'drs1', 'srs2'  ],
                ['dsgcs1', 'dsgrs2'],
                ['deltac', 'deltar'],
                [ 'dscs1', 'dsrs2' ],
            ]
            for source_finder in self.source_finders:
                plot_layout.append([f"frp{source_finder}1",f"frp{source_finder}2"])
            plot_layout.append(['frfph1','frfph2'])
            for n in [1,2,3]:
                plot_layout.append([f'fr{n}sr1',f'fr{n}sr2'])
            for i in range(0,len(self.source_finders),2):
                try:
                    plot_layout.append([f"dsfrp{self.source_finders[i]}",f"dsfrp{self.source_finders[i+1]}"])
                except:
                    plot_layout.append([f"dsfrp{self.source_finders[i]}",""])
            plot_layout.append(["dsfrfph",""])
            def fetch_plot(idx):
                for plot in manifest['plot_info']:
                    if plot['idx'] == idx:
                        return plot
                return None 
            rows = list()
            img_width = 750
            sp = 3 * ' '
            for layout in plot_layout:
                row = list() 
                for idx in layout:
                    plot = fetch_plot(idx)
                    if not plot is None:
                        row.append(f"<td><img src='{plot['image']}' width='{img_width}px'/></td>")
                if len(row)>0:
                    rows.append(f"{indent*sp}   <tr>{''.join(row)}</tr>")
            if len(rows)>0:
                table = ""
                table += f"{indent*sp}<table>\n"
                table += "\n".join(rows)+"\n"
                table += f"{indent*sp}</table>"
            else:
                table = None
            return table

        # creates typhon deep/shallow metrics table
        def create_metrics_table(indent=0):
            def translate(txt):
                if txt == 'rms':
                    return 'RMS'
                if txt == 'isl':
                    return 'Island'
                if txt == 'madfm':
                    return 'MADFM'
                if txt == 'sumsq':
                    return 'SumSq'
                return txt.capitalize()
            qt = self.catalogues['metrics']['parameters'].copy()
            qt.rename_columns(qt.colnames,[" ".join([translate(t) for t in n.split("_")]) for n in qt.colnames])
            metrics = qt.to_pandas().to_html(index=False)
            metrics = re.sub(r"<table +","<table style='padding: 5px;text-align: center;' cellpadding='5px' ",metrics)
            for module in self.source_finders:
                metrics = re.sub(f"<td>{module}</td>",f"<td align='left'><b>{self.plot_labels[module]}</b></td>",metrics)
            for depth in self.image_types:
                metrics = re.sub(f"<td>{depth}</td>",f"<td align='left'>{depth.capitalize()}</td>",metrics)
            for param in np.unique(self.catalogues['metrics']['parameters']['rms_parameter_name']):
                metrics = re.sub(f"<td>{param}</td>",f"<td align='left'>{param}</td>",metrics)
            for param in np.unique(self.catalogues['metrics']['parameters']['isl_parameter_name']):
                metrics = re.sub(f"<td>{param}</td>",f"<td align='left'>{param}</td>",metrics)
            if indent > 0:
                sp = 3 * ' '
                metrics = (indent*sp)+("\n"+(indent*sp)).join(metrics.split('\n'))
            return metrics

        # creates mini-table with residual information
        def create_residual_table(is_header=True):
            table = list()
            table.append("<table style='padding: 5px;text-align: center;' cellpadding='5px' border='1'>")
            table.append("   <tr>")
            if is_header:
                table.append("      <th style='font-weight:normal' class='th-mini-1'>N</th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2'>RMS</th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2'>MADFM</th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2'>SumSq</th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2'>Size</th>")
            else:
                table.append("      <th style='font-weight:normal' class='th-mini-1 res-n'></th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2 res-rms'></th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2 res-madfm'></th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2 res-sumsq'></th>")
                table.append("      <th style='font-weight:normal' class='th-mini-2 res-size'></th>")
            table.append("   </tr>")
            table.append("</table>")
            return "".join([re.sub(r"^ +","",t) for t in table])

        # create index.html file contents
        html_file = "index.html"
        html = list()
        html.append("<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01//EN' 'http://www.w3.org/TR/html4/strict.dtd'>")
        html.append("<html lang='en'>")
        html.append("<head>")
        html.append("   <meta charset='utf-8'>")
        html.append("   <meta name='viewport' content='width=device-width, initial-scale=1'>")
        html.append("")
        html.append("   <title>Hydra</title>")
        html.append("")
        for lib in html_meta:
            if lib['kind'] == 'css':
                lib_dir = re.sub(r"^%s/" % re.escape(archive_root_dir),"",lib['archive'])
                html.append(f"   <link rel='stylesheet' href='{lib_dir}'/>")
        html.append("")
        for lib in html_meta:
            if lib['kind'] == 'js':
                lib_dir = re.sub(r"^%s/" % re.escape(archive_root_dir),"",lib['archive'])
                html.append(f"   <script src='{lib_dir}'></script>")
        cutout_img_size = 300
        cutout_col_size = cutout_img_size+15
        html.append("")
        html.append("   <style>")
        html.append("      h1 {")
        html.append("         padding-left: 10px;")
        html.append("         padding-bottom: 10px;")
        html.append("         padding-top: 10px;")
        html.append("      }")
        html.append("      h3 {")
        html.append("         padding-left: 15px;")
        html.append("         padding-bottom: 10px;")
        html.append("         padding-top: 5px;")
        html.append("      }")
        html.append("      h4 {")
        html.append("         padding-left: 15px;")
        html.append("      }")
        html.append("      div.cutout-stack {")
        html.append("         padding-left:  10px;")
        html.append("      }")
        html.append("      .cutout-stack table,")
        html.append("      .cutout-stack th,")
        html.append("      .cutout-stack td {")
        html.append("         padding: 5px ;")
        html.append("         border: 1px solid black;")
        html.append("         text-align: center;")
        html.append("      }")
        res_table_offset = 35 
        html.append("      .cutout-stack td.cutout-viewer {")
        html.append("         min-width:  %dpx;" % (cutout_col_size))
        html.append("         height: %dpx;" % (cutout_col_size+res_table_offset))
        html.append("         min-height: %dpx;" % (cutout_col_size+res_table_offset))
        html.append("      }")
        html.append("      th.th-mini-1 {")
        n_fractional_size = 0.15
        html.append("         min-width:  %dpx;" % (n_fractional_size*cutout_img_size))
        html.append("      }")
        html.append("      th.th-mini-2 {")
        html.append("         min-width:  %dpx;" % ((1.0-n_fractional_size)*cutout_img_size/4.0))
        html.append("      }")
        html.append("      .cutout-stack button {")
        html.append("         width: 100px;")
        html.append("      }")
        html.append("      .cutout-stack button.cutout-stack-mode {")
        html.append("         width: 150px;")
        html.append("      }")
        html.append("      .cutout-stack-slider {")
        html.append("         width: 750px;")
        html.append("      }")
        html.append("   </style>")
        html.append("")
        html.append("   <script>")
        html.append("      $(document).ready(function() {")
        html.append("         var cutout_states = get_cutout_stack()")
        html.append("")
        html.append("         // cutout viewer annotation control")
        html.append("         var is_annotated = false")
        html.append("         var cutout_viewer_annotate = function() {")
        html.append("            if (!is_annotated) {")
        html.append("               $('.cutout-stack .cutout-panel div.annotated').removeClass('collapse')")
        html.append("               $('.cutout-stack .cutout-panel div.plain').addClass('collapse')")
        html.append("               $('.cutout-stack-annotate').text('Annotated')")
        html.append("               $('.cutout-stack-annotate').addClass('active')")
        html.append("               is_annotated = true")
        html.append("            }")
        html.append("         }")
        html.append("         var cutout_viewer_unannotate = function() {")
        html.append("            if (is_annotated) {")
        html.append("               $('.cutout-stack .cutout-panel div.plain').removeClass('collapse')")
        html.append("               $('.cutout-stack .cutout-panel div.annotated').addClass('collapse')")
        html.append("               $('.cutout-stack-annotate').text('Annotate')")
        html.append("               $('.cutout-stack-annotate').removeClass('active')")
        html.append("               is_annotated = false")
        html.append("            }")
        html.append("         }")
        html.append("         var cutout_viewer_annotated = function() {")
        html.append("            if (is_annotated) {")
        html.append("                cutout_viewer_unannotate()")
        html.append("            } else {")
        html.append("                cutout_viewer_annotate()")
        html.append("            }")
        html.append("         }")
        html.append("         $('.cutout-stack-annotate').click(cutout_viewer_annotated)")
        html.append("")
        html.append("         // cutout stack viewer state reset controls")
        html.append("         var clear_cutout_viewer = function() {")
        html.append("            if (is_annotated) {")
        html.append("               $('.cutout-stack .cutout-panel').html(\"<div class='plain collapse'></div><div class='annotated'></div>\")")
        html.append("            } else {")
        html.append("               $('.cutout-stack .cutout-panel').html(\"<div class='plain'></div><div class='annotated collapse'></div>\")")
        html.append("            }")
        html.append("         }")
        html.append("")
        html.append("         // cutout stack button enable/disable controls")
        html.append("         var set_cutout_stack_button_states = function() {")
        html.append("               $('.cutout-stack-previous').prop('disabled',cutout_states.is_first())")
        html.append("               $('.cutout-stack-next').prop('disabled',cutout_states.is_last())")
        html.append("         }")
        html.append("")
        html.append("         // cutout stack counter field control")
        html.append("         var set_cutout_stack_counter = function() {")
        html.append("            $('.cutout-stack-counter').text('Clump ID: '+cutout_states.get_current_clump_id()+'/'+cutout_states.get_abs_max_clump_id())")
        html.append("         }")
        html.append("")
        html.append("         // cutout stack cutout image selector")
        html.append("         var color_bin   = '#87CEFA'")
        html.append("         var color_match = '#FFFACD'")
        html.append("         var update_cutout_viewer = function() {")
        html.append("            var cutout_stack_params = [")
        html.append("               'n',")
        html.append("               'rms',")
        html.append("               'madfm',")
        html.append("               'sumsq',")
        html.append("               'size',")
        html.append("            ]")
        html.append("            var mini_table = \"\"")
        html.append("               + \"<br>\\n\"")
        html.append("               + \"<div style='padding-left:5px'>\\n\"")
        def aesthetics(insert):
            insert = re.sub(r"<th",'"\n%s+ "%s<th' % ((15*" "),(9*" ")),insert)
            insert = re.sub(r"(<tr>|</tr>)",r'"\n%s+ "%s\1' % ((15*" "),(6*" ")),insert)
            insert = re.sub(r"(</table>)",r'"\n%s+ "%s\1' % ((15*" "),(3*" ")),insert)
            insert = re.sub(r"(</(table|tr|th)>)",r'\1\\n',insert)
            insert = re.sub(r"(<(table|tr).*?>)",r'\1\\n',insert)
            return insert
        html.append("               + \"   %s\"" % aesthetics(create_residual_table(False)))
        html.append("               + \"</div>\\n\"")
        html.append("            clear_cutout_viewer()")
        html.append("            set_cutout_stack_counter()")
        html.append("            set_cutout_stack_button_states()")
        html.append("            var clump = cutout_states.get_clump()")
        html.append("            for (var module in clump['modules']) {")
        html.append("               for (var depth in clump['modules'][module]) {")
        html.append("                  for (var complications in clump['modules'][module][depth]) {")
        html_image = "\"<img src='\"+clump['modules'][module][depth][complications]+\"' width='%dpx'>\"+((module=='__main__') ? \"\" : mini_table)" % cutout_img_size
        html.append("                     $('.cutout-stack .cutout-panel.'+module+'.'+depth+' div.'+complications).html(%s)" % html_image)
        html.append("                     if (module != '__main__') {")
        html.append("                        for (var j in cutout_stack_params) {")
        html.append("                           $('.cutout-stack .cutout-panel.'+module+'.'+depth+' th.res-'+cutout_stack_params[j]).text(clump['modules'][module][depth]['params'][cutout_stack_params[j]])")
        html.append("                        }")
        html.append("                     }")
        html.append("                  }")
        html.append("               }")
        html.append("            }")
        html.append("")
        html.append("            $('.cutout-stack div.cluster-table').html(cutout_states.get_cluster_table())")
        html.append("            $('.cutout-viewer.cutout-panel').css('background-color', 'transparent')")
        html.append("            var source_finders = cutout_states.get_source_finders()")
        html.append("            if (cutout_modes_idx!=0 && source_finders !== undefined) {")
        html.append("               var color_deep    = undefined")
        html.append("               var color_shallow = undefined")
        html.append("               switch(cutout_modes_idx) {")
        html.append("                  case 0:")
        html.append("                     break;")
        html.append("                  case 1:")
        html.append("                     $('.cluster-info-table .mode-title').text('Completeness')")
        html.append("                     color_deep    = color_bin") # notes: https://www.w3schools.com/cssref/css_colors.asp
        html.append("                     color_shallow = color_match")
        html.append("                     break;")
        html.append("                  case 2:")
        html.append("                     $('.cluster-info-table .mode-title').text('Reliability')")
        html.append("                     color_shallow = color_bin")
        html.append("                     color_deep    = color_match")
        html.append("                     break;")
        html.append("                  case 3:")
        html.append("                     break;")
        html.append("                  default:")
        html.append("                     break;")
        html.append("               }")
        html.append("               if (source_finders.deep !== undefined) {")
        html.append("                  for (var i in source_finders.deep) {")
        html.append("                     $('.cutout-viewer.cutout-panel.deep.'+source_finders.deep[i]).css('background-color',color_deep)")
        html.append("                  }")
        html.append("               }")
        html.append("               if (source_finders.shallow !== undefined) {")
        html.append("                  for (var i in source_finders.shallow) {")
        html.append("                     $('.cutout-viewer.cutout-panel.shallow.'+source_finders.shallow[i]).css('background-color',color_shallow)")
        html.append("                  }")
        html.append("               }")
        html.append("               $('.cluster-info-table .mode-title,.mode-data').removeClass('collapse')")
        html.append("               if (source_finders.ref_ids !== undefined) {")
        html.append("                  for (var i in source_finders.ref_ids) {")
        html.append("                     var is_deep =  $('.cluster-info-table .ref-id-'+source_finders.ref_ids[i][0]).hasClass('deep')")
        html.append("                     $('.cluster-info-table .ref-id-'+source_finders.ref_ids[i][0]).css('background-color', (is_deep) ? color_deep : color_shallow)")
        html.append("                     $('.cluster-info-table .ref-id-'+source_finders.ref_ids[i][0]+' .mode-data').text(source_finders.ref_ids[i][1])")
        html.append("                  }")
        html.append("               }")
        html.append("            }")
        html.append("         }")
        html.append("")
        html.append("         // cutout stack previous/next controls")
        html.append("         var cutout_stack_previous = function() {")
        html.append("            var clump_id = cutout_states.get_previous_clump_id()")
        html.append("            if (cutout_modes_idx==0) {")
        html.append("               $('.cutout-stack-slider').val(clump_id)")
        html.append("            }")
        html.append("            update_cutout_viewer()")
        html.append("         }")
        html.append("         $('.cutout-stack-previous').click(cutout_stack_previous)")
        html.append("         var cutout_stack_next = function() {")
        html.append("            var clump_id = cutout_states.get_next_clump_id()")
        html.append("            if (cutout_modes_idx==0) {")
        html.append("               $('.cutout-stack-slider').val(clump_id)")
        html.append("            }")
        html.append("            update_cutout_viewer()")
        html.append("         }")
        html.append("         $('.cutout-stack-next').click(cutout_stack_next)")
        html.append("")
        html.append("         // arrow key controls")
        html.append("         var isShitKeyDown = false")
        html.append("         $(document).keydown(function(e) {")
        html.append("            var arrow = {shift: 16,left: 37,up: 38,right: 39,down: 40}")
        html.append("            switch (e.which) {")
        html.append("               case arrow.shift:")
        html.append("                  isShitKeyDown = true")
        html.append("                  break;")
        html.append("               case arrow.right:")
        html.append("                  if (isShitKeyDown) {")
        html.append("                     var idx = $('.cutout-stack-slider').val()")
        html.append("                     var max = $('.cutout-stack-slider').attr('max')")
        html.append("                     if (idx++ < max) {")
        html.append("                        switch(cutout_modes_idx) {")
        html.append("                           case 1:")
        html.append("                              cutout_states.set_deep_shallow_completeness_bin(dsc_bins[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('S/N: '+dsc_sns[idx])")
        html.append("                              break;")
        html.append("                           case 2:")
        html.append("                              cutout_states.set_deep_shallow_reliability_bin(dsr_bins[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('S/N: '+dsr_sns[idx])")
        html.append("                              break;")
        html.append("                           case 3:")
        html.append("                              cutout_states.set_clump_size(clump_sizes[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('Size: '+clump_sizes[idx])")
        html.append("                              break;")
        html.append("                           default:")
        html.append("                              return;")
        html.append("                        }")
        html.append("                        cutout_states.set_clump(cutout_states.get_current_clump_id())")
        html.append("                        $('.cutout-stack-slider').val(idx)")
        html.append("                        $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                        update_cutout_viewer()")
        html.append("                     }")
        html.append("                  } else {")
        html.append("                     cutout_stack_next();")
        html.append("                  }")
        html.append("                  break;")
        html.append("               case arrow.left:")
        html.append("                  if (isShitKeyDown) {")
        html.append("                     var idx = $('.cutout-stack-slider').val()")
        html.append("                     var min = $('.cutout-stack-slider').attr('min')")
        html.append("                     if (min < idx--) {")
        html.append("                        switch(cutout_modes_idx) {")
        html.append("                           case 1:")
        html.append("                              cutout_states.set_deep_shallow_completeness_bin(dsc_bins[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('S/N: '+dsc_sns[idx])")
        html.append("                              break;")
        html.append("                           case 2:")
        html.append("                              cutout_states.set_deep_shallow_reliability_bin(dsr_bins[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('S/N: '+dsr_sns[idx])")
        html.append("                              break;")
        html.append("                           case 3:")
        html.append("                              cutout_states.set_clump_size(clump_sizes[idx])")
        html.append("                              $('.cutout-stack-slider-values-1').text('Size: '+clump_sizes[idx])")
        html.append("                              break;")
        html.append("                           default:")
        html.append("                              return;")
        html.append("                        }")
        html.append("                        cutout_states.set_clump(cutout_states.get_current_clump_id())")
        html.append("                        $('.cutout-stack-slider').val(idx)")
        html.append("                        $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                        update_cutout_viewer()")
        html.append("                     }")
        html.append("                  } else {")
        html.append("                     cutout_stack_previous()")
        html.append("                  }")
        html.append("                  break;")
        html.append("               case arrow.down:")
        html.append("                  if (isShitKeyDown) {")
        html.append("                     rotate_mode_right()")
        html.append("                  } else {")
        html.append("                     cutout_viewer_annotate()")
        html.append("                  }")
        html.append("                  break;")
        html.append("               case arrow.up:")
        html.append("                  if (isShitKeyDown) {")
        html.append("                     rotate_mode_left()")
        html.append("                  } else {")
        html.append("                     cutout_viewer_unannotate()")
        html.append("                  }")
        html.append("                  break;")
        html.append("               default:")
        html.append("                  return;")
        html.append("            }")
        html.append("            e.preventDefault()")
        html.append("         })")
        html.append("         $(document).keyup(function(e) {")
        html.append("            if (e.which == 16) {")
        html.append("               isShitKeyDown = false")
        html.append("            }")
        html.append("         })")
        html.append("")
        html.append("         // clump_id keyboard input selector")
        html.append("         var is_valid_clump_id = function(value) {")
        html.append("            return !isNaN(value)&&cutout_states.get_min_clump_id()<=value&&value<=cutout_states.get_max_clump_id()")
        html.append("         }")
        html.append("         $('.cutout-stack-input').keyup(function(e) {")
        html.append("            if (is_valid_clump_id(this.value)) {")
        html.append("               if (e.key === 'Enter' || e.keyCode === 13) {")
        html.append("                   cutout_states.set_clump(parseInt(this.value,10))")
        html.append("                   update_cutout_viewer()")
        html.append("                   this.value = ''")
        html.append("               }")
        html.append("            } else {")
        html.append("               this.value = this.value.replace(/\w$/,'')")
        html.append("            }")
        html.append("         })")
        html.append("         $('.cutout-stack-input').blur(function(e) {")
        html.append("            if (is_valid_clump_id(this.value)) {")
        html.append("               cutout_states.set_clump(parseInt(this.value,10))")
        html.append("               update_cutout_viewer()")
        html.append("               this.value = ''")
        html.append("            }")
        html.append("         })")
        html.append("")
        html.append("         var dsc_sns     = cutout_states.get_deep_shallow_completeness_sns()")
        html.append("         var dsc_bins    = cutout_states.get_deep_shallow_completeness_bins()")
        html.append("         var dsr_sns     = cutout_states.get_deep_shallow_reliability_sns()")
        html.append("         var dsr_bins    = cutout_states.get_deep_shallow_reliability_bins()")
        html.append("         var clump_sizes = cutout_states.get_clump_sizes()")
        #html.append("         var n_clumps    = cutout_states.get_n_clumps()")
        html.append("         var cutout_modes_idx = 0")
        html.append("         var set_mode = function() {")
        html.append("            switch(cutout_modes_idx) {")
        html.append("               case 0:")
        html.append("                  cutout_states.set_mode_clump()")
        html.append("                  $('.cutout-viewer .deep-title').text('Deep')")
        html.append("                  $('.cutout-viewer .deep-title').css('background-color','transparent')")
        html.append("                  $('.cutout-viewer .shallow-title').text('Shallow')")
        html.append("                  $('.cutout-viewer .shallow-title').css('background-color','transparent')")
        html.append("                  $('.cutout-stack-mode').text('Mode')")
        html.append("                  $('.cutout-stack-mode').removeClass('active')")
        html.append("                  $('.cutout-stack-slider-values-1').addClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-2').addClass('collapse')")
        html.append("                  $('.cutout-stack-input').prop('disabled',false)")
        html.append("                  $('.cutout-stack-slider').attr('min',cutout_states.get_min_clump_id())")
        html.append("                  $('.cutout-stack-slider').attr('value',cutout_states.get_current_clump_id())")
        html.append("                  $('.cutout-stack-slider').attr('max',cutout_states.get_max_clump_id())")
        html.append("                  $('.cutout-stack-slider').val(cutout_states.get_min_clump_id())")
        html.append("                  $('.cutout-stack-slider').addClass('collapse')")
        html.append("                  break;")
        html.append("               case 1:")
        html.append("                  cutout_states.set_deep_shallow_completeness_bin(dsc_bins[0])")
        html.append("                  $('.cutout-viewer .deep-title').html('Deep<br>Bins')")
        html.append("                  $('.cutout-viewer .deep-title').css('background-color',color_bin)")
        html.append("                  $('.cutout-viewer .shallow-title').html('Shallow<br>Matches')")
        html.append("                  $('.cutout-viewer .shallow-title').css('background-color',color_match)")
        html.append("                  $('.cutout-stack-mode').text('Completeness')")
        html.append("                  $('.cutout-stack-mode').addClass('active')")
        html.append("                  $('.cutout-stack-slider-values-1').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-2').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-1').text('S/N: '+dsc_sns[0])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  $('.cutout-stack-input').prop('disabled',true)")
        html.append("                  $('.cutout-stack-slider').attr('min',0)")
        html.append("                  $('.cutout-stack-slider').attr('value',0)")
        html.append("                  $('.cutout-stack-slider').attr('max',dsc_bins.length-1)")
        html.append("                  $('.cutout-stack-slider').val(0)")
        html.append("                  $('.cutout-stack-slider').removeClass('collapse')")
        html.append("                  break;")
        html.append("               case 2:")
        html.append("                  cutout_states.set_deep_shallow_reliability_bin(dsr_bins[0])")
        html.append("                  $('.cutout-viewer .deep-title').html('Deep<br>Matches')")
        html.append("                  $('.cutout-viewer .deep-title').css('background-color',color_match)")
        html.append("                  $('.cutout-viewer .shallow-title').html('Shallow<br>Bins')")
        html.append("                  $('.cutout-viewer .shallow-title').css('background-color',color_bin)")
        html.append("                  $('.cutout-stack-mode').text('Reliability')")
        html.append("                  $('.cutout-stack-mode').addClass('active')")
        html.append("                  $('.cutout-stack-slider-values-1').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-2').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-1').text('S/N: '+dsr_sns[0])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  $('.cutout-stack-input').prop('disabled',true)")
        html.append("                  $('.cutout-stack-slider').attr('min',0)")
        html.append("                  $('.cutout-stack-slider').attr('value',0)")
        html.append("                  $('.cutout-stack-slider').attr('max',dsr_bins.length-1)")
        html.append("                  $('.cutout-stack-slider').val(0)")
        html.append("                  $('.cutout-stack-slider').removeClass('collapse')")
        html.append("                  break;")
        html.append("               case 3:")
        html.append("                  cutout_states.set_clump_size(clump_sizes[0])")
        html.append("                  $('.cutout-viewer .deep-title').text('Deep')")
        html.append("                  $('.cutout-viewer .deep-title').css('background-color','transparent')")
        html.append("                  $('.cutout-viewer .shallow-title').text('Shallow')")
        html.append("                  $('.cutout-viewer .shallow-title').css('background-color','transparent')")
        html.append("                  $('.cutout-stack-mode').text('Clump Size')")
        html.append("                  $('.cutout-stack-mode').addClass('active')")
        html.append("                  $('.cutout-stack-slider-values-1').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-2').removeClass('collapse')")
        html.append("                  $('.cutout-stack-slider-values-1').text('Size: '+clump_sizes[0])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  $('.cutout-stack-input').prop('disabled',true)")
        html.append("                  $('.cutout-stack-slider').attr('min',0)")
        html.append("                  $('.cutout-stack-slider').attr('value',0)")
        html.append("                  $('.cutout-stack-slider').attr('max',clump_sizes.length-1)")
        html.append("                  $('.cutout-stack-slider').val(0)")
        html.append("                  $('.cutout-stack-slider').removeClass('collapse')")
        html.append("                  break;")
        html.append("               default:")
        html.append("                  return;")
        html.append("            }")
        html.append("            update_cutout_viewer()")
        html.append("         }")
        html.append("         var rotate_mode_right = function() {")
        html.append("            cutout_modes_idx = ++cutout_modes_idx % 4")
        html.append("            set_mode()")
        html.append("         }")
        html.append("         var rotate_mode_left = function() {")
        html.append("            cutout_modes_idx = (4+--cutout_modes_idx) % 4")
        html.append("            set_mode()")
        html.append("         }")
        html.append("         $('.cutout-stack-mode').on('click', rotate_mode_right)")
        html.append("")
        html.append("         $('.cutout-stack-slider').on('input change', function() {")
        html.append("            switch(cutout_modes_idx) {")
        html.append("               case 0:")
        html.append("                  cutout_states.set_clump(parseInt(this.value))")
        html.append("                  break;")
        html.append("               case 1:")
        html.append("                  cutout_states.set_deep_shallow_completeness_bin(dsc_bins[this.value])")
        html.append("                  cutout_states.set_clump(cutout_states.get_current_clump_id())")
        html.append("                  $('.cutout-stack-slider-values-1').text('S/N: '+dsc_sns[this.value])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  break;")
        html.append("               case 2:")
        html.append("                  cutout_states.set_deep_shallow_reliability_bin(dsr_bins[this.value])")
        html.append("                  cutout_states.set_clump(cutout_states.get_current_clump_id())")
        html.append("                  $('.cutout-stack-slider-values-1').text('S/N: '+dsr_sns[this.value])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  break;")
        html.append("               case 3:")
        html.append("                  cutout_states.set_clump_size(clump_sizes[this.value])")
        html.append("                  cutout_states.set_clump(cutout_states.get_current_clump_id())")
        html.append("                  $('.cutout-stack-slider-values-1').text('Size: '+clump_sizes[this.value])")
        html.append("                  $('.cutout-stack-slider-values-2').text('Clumps: '+cutout_states.get_no_clumps())")
        html.append("                  break;")
        html.append("               default:")
        html.append("                  return;")
        html.append("            }")
        html.append("            update_cutout_viewer()")
        html.append("         })")
        html.append("")
        html.append("         // initialize cutout stack viewer")
        html.append("         update_cutout_viewer()")
        html.append("      })")
        html.append("   </script>")
        html.append("</head>")
        html.append("<body>")
        html.append("   <div class='cutout-stack' style='padding-left: 15px'>")
        html.append("	   <h1 style='padding-bottom: 0'>Hydra Cutout Viewer</h1>")
        html.append("	   <div style='padding: 0 0 10px 15px;font-size:small;color:red'><b>File: %s</b></div>" % re.sub(r"^(.*?/)+","",self.tar_gz_file))
        html.append("	   <div style='padding: 0 0 10px 15px;font-size:small;color:blue'><b>Units: [RMS]=[MADFM]=mJy/(arcmin<sup>2</sup>*beam), [SumSq]=mJy<sup>2</sup>/(arcmin*beam)<sup>2</sup>, [Size]=arcmin</b></div>")
        html.append("")
        html.append("      <h3>")
        html.append("          <button type='button' class='btn btn-outline-primary cutout-stack-previous' disabled>Previous</button>")
        html.append("          <button type='button' class='btn btn-outline-primary cutout-stack-next'>Next</button>")
        html.append("          <button type='button' class='btn btn-outline-success cutout-stack-annotate'>Annotate</button>")
        html.append("          <input type='text' class='cutout-stack-input' style='vertical-align:middle' placeholder='Input Clump ID' value=''/>")
        html.append("      </h3>")
        html.append("      <h4>")
        html.append("          <table style='width:1200px' class='cutout-stack-controls'>")
        html.append("             <tr>")
        html.append("                <td style='text-align:left;border-style:hidden;width:19%'>")
        html.append("                   <span class='label label-outline-defualt cutout-stack-counter'>Clump ID:</span>")
        html.append("                </td><td style='text-align:left;border-style:hidden;width:27%'>")
        html.append("                   <span style='text-align:left;border-style:hidden' class='label label-outline-defualt collapse cutout-stack-slider-values-1'></span>")
        html.append("                </td><td style='text-align:left;border-style:hidden;width:27%'>")
        html.append("                   <span style='text-align:left;border-style:hidden' class='label label-outline-defualt collapse cutout-stack-slider-values-2'></span>")
        html.append("                </td><td style='text-align:left;border-style:hidden;width:27%'>")
        html.append("                   <span style='text-align:left;border-style:hidden' class='label label-outline-defualt collapse cutout-stack-slider-values-3'></span>")
        html.append("                </td>")
        html.append("             </tr><tr>")
        html.append("                <td style='text-align:right;border-style:hidden'>")
        html.append("                   <button type='button' class='btn btn-outline-success cutout-stack-mode'>Mode</button>")
        html.append("                </td><td colspan=3 style='text-align:left;border-style:hidden'>")
        html.append("                   <input type='range' min='1' max='%s' value='1' class='cutout-stack-slider collapse'>" % max(self.catalogues['cluster']['clump_id']))
        html.append("                </td>")
        html.append("             </tr>")
        html.append("          </table>")
        html.append("      </h4>")
        html.append("")
        html.append("      <table class='cutout-viewer'>")
        html.append("         <tr>")
        html.append("             <th class='cutout-viewer'>Image</th>")
        html.append("             <th class='cutout-viewer __main__'>Cutout</th>")
        def get_color(module):
            color = self.plot_colors[module]
            return color if re.search("#",color) else mcolors.cnames[color]
        for module in self.source_finders:
            html.append("             <th class='cutout-viewer %s'><span style='color:%s'>&#9632;</span> %s Residuals<br><div style='padding-top:5px'>%s</div></th>" % (module,get_color(module),self.plot_labels[module],create_residual_table()))
        html.append("         </tr>")
        for depth in self.image_types:
            html.append("         <tr>")
            html.append("             <th class='cutout-viewer %s-title'>%s</th>" % (depth,depth.capitalize()))
            html.append("             <td class='cutout-viewer cutout-panel __main__ %s' stye='text-align:center' valign='top'></td>" % depth)
            for module in self.source_finders:
                html.append("             <td class='cutout-viewer cutout-panel %s %s' stye='text-align:center' valign='top'></td>" % (module,depth))
            html.append("         </tr>")
        html.append("      </table>")
        html.append("")
        html.append("   <br/>")
        html.append("      <div class='cluster-table'>")
        html.append("      </div>")
        html.append("   </div>")
        html.append("   <br/>")
        html.append("   <hr/>")
        html.append("	<h1 style='padding-bottom: 0'>Hydra: Typhon Run Metrics</h1>")
        html.append("   <table>")
        html.append("      <tr>")
        html.append("         <td><img src='%s' width='750px'/></td>" % fetch_plot('hpp1')['image'])
        html.append("         <td><img src='%s' width='750px'/></td>" % fetch_plot('hpp2')['image'])
        html.append("      </tr>")

        for metric in ['prd_cpu_time','rms','madfm','sumsq']:
            html.append("      <tr>")
            html.append("         <td><img src='%s' width='750px'/></td>" % fetch_plot(f"hpp1_{metric}")['image'])
            html.append("         <td><img src='%s' width='750px'/></td>" % fetch_plot(f"hpp2_{metric}")['image'])
            html.append("      </tr>")

        html.append("   </table>")
        html.append("	<div style='padding: 0 0 10px 15px;font-size:small'><b>Units: [RMS]=[MADFM]=Jy/beam, [SumSq]=(Jy/beam)<sup>2</sup></div>")
        html.append("   <div style='padding-left: 15px'>")
        html.append(create_metrics_table(2))
        html.append("   </div>")
        html.append("   <br/>")
        html.append("   <hr/>")
        html.append("	<h1>Hydra Plots</h1>")
        html.append("   <div style='padding-left: 15px'>")
        html.append(layout_plots(2))
        html.append("   </div>")
        html.append("</body>")
        html.append("</html>")
        html_meta.append({
            'kind': 'index.html',
            'keep': False,
            'system': f"{cache}/{html_file}",
            'archive': f"{archive_root_dir}/{html_file}",
            'contents': html,
        })

        # create html files
        for html in html_meta:
            if not html['keep']:
                with open(html['system'],'w') as fd:
                    fd.write("\n".join(html['contents']))

        # ok, we need to check if external libaries missing and add them
        for html in html_meta.copy():
            if html['keep']:
                fname   = re.sub(r"^(.*?/)+(.*)$",r"\2",html['system'])
                sys_dir = re.sub(r"/+$","",re.sub(r"^((.*?/)+).*$",r"\1",html['system']))
                arc_dir = re.sub(r"/+$","",re.sub(r"^((.*?/)+).*$",r"\1",html['archive']))
                for sfile in glob.glob(f"{sys_dir}/*.min.*"):
                    if not re.search(r"%s$" % re.escape(fname),sfile):
                        sfile = re.sub(r"^(.*?/)+(.*)$",r"\2",sfile)
                        html_meta.append({
                            'kind': 'extra',
                            'keep': True,
                            'system': f"{sys_dir}/{sfile}",
                            'archive': f"{arc_dir}/{sfile}",
                        })

        return html_meta


    def __tarball(self):
        print_ok(f"Archiving Data...")
        html_manifest = dict()

        # create catalogue meta data
        def translate(fname,suffix=None):
            fname = re.sub(r"^(.*?/)+",f"{cache}/",fname)
            suffix = "hydra" + ("" if suffix is None else f".{suffix}")
            fname = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",f".{suffix}.fits",fname)
            return fname
        meta_data = {
            'cluster': {
                'fname': translate(self.deep_image,'cluster_catalogue'),
                'table': self.catalogues['cluster']
            },
            'deep_rms_box_statistics': {
                'fname': translate(self.deep_image,'deep_rms_box_statistics'),
                'table': self.catalogues['metrics']['rms_box_statistics']['deep'],
            },
            'shallow_rms_box_statistics': {
                'fname': translate(self.deep_image,'shallow_rms_box_statistics'),
                'table': self.catalogues['metrics']['rms_box_statistics']['shallow'],
            },
            'deep_statistics': {
                'fname': translate(self.deep_image,'deep_statistics'),
                'table': self.catalogues['metrics']['statistics']['deep'],
            },
            'shallow_statistics': {
                'fname': translate(self.deep_image,'shallow_statistics'),
                'table': self.catalogues['metrics']['statistics']['shallow'],
            },
            'parameters': {
                'fname': translate(self.deep_image,'global_metrics'),
                'table': self.catalogues['metrics']['parameters'],
            },
            'clump': {
                'fname': translate(self.deep_image,'clump_catalogue'),
                'table': self.catalogues['clump'],
            },
        }
        tnames = ['deep_rms_box_statistics','shallow_rms_box_statistics','deep_statistics','shallow_statistics','parameters']
        cnames = list()
        for module in self.catalogues['modules']:
            ## TO-DO: Fix this... :/
            #if module == 'profound':
            #    continue
            for depth in self.catalogues['modules'][module]['tables']:
                meta_data[f"{module}_{depth}"] = {
                    'fname': translate(self.deep_image,f"{module}.{depth}"),
                    'table': self.catalogues['modules'][module]['tables'][depth],
                }
                cnames.append(f"{module}_{depth}")
        if 'extern' in self.catalogues:
            for catalogue in self.catalogues['extern']:
                meta_data[catalogue.lower()] = {
                    'fname': translate(self.deep_image,f"{catalogue.lower()}.external_catalogue"),
                    'table': self.catalogues['extern'][catalogue]['table'],
                }
        for datum in meta_data:
            print_ok(f"> Creating {datum} fits table file: {meta_data[datum]['fname']}") 
            meta_data[datum]['table'].write(meta_data[datum]['fname'],overwrite=True)


        # let's archive stuff...
        print_ok(f"> Creating Archive: {self.tar_gz_file}")
        def strip(fname):
            return re.sub(r"^(.*?/)+","",fname) 
        root_dir = re.sub(r"\.tar\.gz$","_dir",strip(self.tar_gz_file))
        with tarfile.open(self.tar_gz_file,"w:gz") as fd:
            # archive catalogue info
            for datum in meta_data:
                print_ok(f"> Archiving: {meta_data[datum]['fname']}") 
                if datum in tnames:
                    fname = f"{root_dir}/statistics/{strip(meta_data[datum]['fname'])}"
                elif datum in cnames:
                    fname = f"{root_dir}/catalogues/{datum.split('_')[1]}/{strip(meta_data[datum]['fname'])}"
                else:
                    fname = f"{root_dir}/{'injected/' if re.search('.external_catalogue.fits$',meta_data[datum]['fname']) else ''}{strip(meta_data[datum]['fname'])}"
                fd.add(meta_data[datum]['fname'],fname,False)
                os.remove(meta_data[datum]['fname'])

            # archive noise images
            for rms_img in self.rms_imgs:
                fname = f"{cache}/"+re.sub(r"\.[Ff][Ii][Tt]([Ss]|)$",f".hydra.bane.{rms_img}.noise.fits",strip(self.deep_image))
                print_ok(f"> Archiving: {fname}")
                self.rms_imgs[rms_img].write(fname)
                fd.add(fname,f"{root_dir}/noise/{rms_img}/{strip(fname)}")
                os.remove(fname)
            typhon_archive = {
                'deep': self.deep_tar_gz_file,
                'shallow': self.shallow_tar_gz_file,
            }

            # archive .reg and residual files
            work_dir = get_this_source_file_directory()
            def flush_cache(member):
                os.chdir(cache)
                subdir = re.sub(r"^(.*?/).*$",r"\1",member.name)
                if os.path.isdir(subdir):
                    shutil.rmtree(subdir)
                os.chdir(work_dir)
            for depth in typhon_archive:
                reg = re.compile("\.(|shallow\.)(%s)\.reg$" % '|'.join(list(self.catalogues['modules'].keys())))
                res = re.compile("\.(|shallow\.)(%s)\.(residual(|\.model))\.fits$" % '|'.join(list(self.catalogues['modules'].keys())))
                with tarfile.open(typhon_archive[depth],"r") as fd_typhon:
                    for member in fd_typhon.getmembers():
                        fname = strip(member.name)
                        if reg.search(fname): # .reg files
                            fname = f"{cache}/{member.name}"
                            print_ok(f"> Archiving: {fname}")
                            fd_typhon.extract(member,cache)
                            rname = reg.sub(r".hydra.\2.%s.reg" % depth,strip(fname))
                            fd.add(fname,f"{root_dir}/regions/{depth}/{rname}",False)
                            flush_cache(member)
                        elif res.search(fname): # residual files
                            fname = f"{cache}/{member.name}"
                            print_ok(f"> Archiving: {fname}")
                            fd_typhon.extract(member,cache)
                            rname = res.sub(r".hydra.\2.%s.\3.fits" % depth,strip(fname))
                            fd.add(fname,f"{root_dir}/residuals/{depth}/{rname}",False)
                            flush_cache(member)
                    png_prefix = re.sub(r"\.[Ff][Ii][Tt]([Ss]|)$",".hydra",strip(self.deep_image))

            # archive plot files
            html_manifest['plot_info'] = list()
            for plot_info in self.plot_collection:
                info = plot_info.copy()
                fname = f"{cache}/{png_prefix}.{plot_info['filename']}"
                print_ok(f"> Archiving: {fname}")
                info['filename'] = fname
                self.__plot([info],is_show=False,is_save=True)
                aname = f"{root_dir}/plots/{strip(fname)}"
                fd.add(fname,aname,False)
                os.remove(fname)
                html_manifest['plot_info'].append({
                    'idx': plot_info['idx'],
                    'type': plot_info['type'],        
                    'image': re.sub(r"^%s/" % re.escape(root_dir),"",aname),
                    'html': plot_info['html'],
                })
                if len(plot_info['html']['datum'])>0:
                    fname = f"{cache}/{info['html']['fname']}"
                    plot_info['html']['datum'].write(fname,format='fits',overwrite=True)
                    aname = f"{root_dir}/plots/{strip(fname)}"
                    fd.add(fname,aname,False)
                    os.remove(fname)

            # add graphic context file
            fname = f"{cache}/graphics_context.yml"
            with open(fname,'w') as fd_gc:
                yml.dump(self.graphics_context,fd_gc)
            fd.add(fname,f"{root_dir}/plots/{strip(fname)}",False)
            os.remove(fname)

            # add cutouts
            dir_map = dict()
            subdir = f"{root_dir}/cutouts"
            for depth in self.image_types:
                dir_map[f"{subdir}/images/{depth}"] = {
                    're': re.compile(r"\.%s\.png$" % depth),
                    'state': {
                        'module': '__main__',
                        'depth': depth,
                        'is_annotated': False,
                        'cutouts': list(),
                     },
                }
                dir_map[f"{subdir}/images/{depth}/annotated"] = {
                    're': re.compile(r"\.%s\.annotated\.png$" % depth),
                    'state': {
                        'module': '__main__',
                        'depth': depth,
                        'is_annotated': True,
                        'cutouts': list(),
                     },
                }
                for module in self.source_finders:
                    dir_map[f"{subdir}/residuals/{depth}/{module}"] = {
                        're': re.compile(r"\.%s\.residual\.%s\.png$" % (depth,module)),
                        'state': {
                            'module': module,
                            'depth': depth,
                            'is_annotated': False,
                            'cutouts': list(),
                         },
                    }
                    dir_map[f"{subdir}/residuals/{depth}/{module}/annotated"] = {
                        're': re.compile(r"\.%s\.residual\.%s\.annotated\.png$" % (depth,module)),
                        'state': {
                            'module': module,
                            'depth': depth,
                            'is_annotated': True,
                            'cutouts': list(),
                         },
                    }
            for fname in self.cutout_files:
                for dname in dir_map:
                    if dir_map[dname]['re'].search(fname):
                        print_ok(f"> Archiving: {fname}")
                        aname = dname+"/"+re.sub(r"^(.*?/)+","",fname)
                        dir_map[dname]['state']['cutouts'].append(re.sub(r"^%s/" % re.escape(root_dir),"",aname))
                        fd.add(fname,aname,False)
                        os.remove(fname)
                        break

            # update html manifest
            html_manifest['cutouts'] = [dir_map[m]['state'] for m in dir_map]

            # build html
            html_meta = self.__html(html_manifest)
            for html in html_meta:
                print_ok(f"> Arhiving: {html['system']}")
                fd.add(html['system'],html['archive'])
                if not html['keep']:
                    os.remove(html['system'])

            # print contents
            print(f"> Tarball: {self.tar_gz_file}")
            for item in np.sort(fd.getnames()):
                print(f"> o {item}")
            print("> [EOF]")
        print_ok(f"[Done]")
        return self

    def __clusterize(self,source_finder,image_type):
        qt = self.catalogues['cluster']
        qt = qt[((qt['source_finder']==source_finder)&(qt['image_type']==image_type))]
        if len(qt) > 0:
            df = qt.to_pandas()
            cluster_parameters = [{
                'survey': f"{source_finder}_{image_type}",
                'ref_ids': (df.index+1).to_numpy(),
                #'cat_nos': df.id.to_numpy(),
                'cat_nos': df['cat_id'].to_numpy(),
                'max_extent': self.skirt_factor * np.max(np.append(df.extent_semimajor,df.extent_semiminor)) * u.deg,
                'coords': SkyCoord(ra=df.ra,dec=df.dec,frame="icrs",unit="deg"),
                'a_extents': self.skirt_factor * df.extent_semimajor.to_numpy() * u.deg,
                'b_extents': self.skirt_factor * df.extent_semiminor.to_numpy() * u.deg,
                't_extents': df.extent_angle.to_numpy() * u.deg
            }]
            ct = Cluster(cluster_parameters).get_qtable()
            ct['extent_semimajor'] /= self.skirt_factor
            ct['extent_semiminor'] /= self.skirt_factor
            ct['ref_id'] = [qt[(qt['cat_id']==cat_no)]['ref_id'][0] for cat_no in ct['cat_no']]
            #ct['clump_id'] = [qt[(qt['cat_id']==cat_no)]['clump_id'][0] for cat_no in ct['cat_no']]
            ct.add_column(
                    name = 'old_clump_id',
                    col = [qt[(qt['cat_id']==cat_no)]['clump_id'][0] for cat_no in ct['cat_no']],
                    index = 1
            )
            ct.add_column([re.sub(r"^(.*?_)+","",s) for s in ct['survey']],name='image_type',index=4)
            ct['survey'] = [re.sub(r"_$","",re.sub(r"^(.*?_)+.*$",r"\1",s)) for s in ct['survey']]
            ct.rename_column('survey','source_finder')
            #print_ok(ct)
        else:
            ct = QTable()
        return ct

    def __get_subclusters(self):
        qt = self.catalogues['cluster']
        modules = np.unique(qt['source_finder'])
        image_types = np.unique(qt['image_type'])
        ct = QTable()
        for module in modules:
            for image_type in image_types:
                if len(ct)>0:
                    qt = self.__clusterize(module,image_type)
                    if len(qt)>0:
                        qt['clump_id'] += max(ct['clump_id'])
                        ct = vstack([ct,qt])
                else:
                    ct = self.__clusterize(module,image_type)
        #ct.sort(['ref_id','old_clump_id','clump_id','source_finder','image_type'])
        ct.sort(['ref_id'])
        ct.add_column(
            name = 'subclump_id',
            col = np.zeros(len(ct),dtype=np.int64),
            index = 3
        )
        cnt = 1
        found = list()
        for cid in ct['clump_id']:
            if not cid in found:
                ct['subclump_id'][(ct['clump_id']==cid)] = cnt
                found.append(cid)
                cnt += 1
        return ct

    def mk_pybdsf_residual(self):
        print_warning("Hello Kitty!")
        # Modes:
        # Single_Select
        # Single
        # Multi
        # CDS_Single
        # CDS_Multi

        # Mode: Single_Select
        # TO-DO: We're gonna have to cluster on source_finder to get proper stats
        # 1) Let's create a cluster table for each source finder (ref main cluster table)
        # 2) Let's create histograms for these and place on web interface
        #qt = self.__get_subclump_size_distribution_table()
        #print_warning(qt)
        # 3) Then back to here...
        source_finder = 'pybdsf'
        image_type = 'deep'
        img = HydraCroppingTool(self.deep_image)
        qt = self.catalogues['cluster']
        qt = qt[((qt['source_finder']==source_finder)&(qt['image_type']==image_type))]
        catalogue = QTable()
        for subclump_id in np.unique(qt['subclump_id']):
            subclump = qt[(qt['subclump_id']==subclump_id)]
            if len(subclump)==1:
                if len(catalogue)>0:
                    catalogue = vstack([catalogue,subclump])
                else:
                    catalogue = subclump
        residual_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",f".hydra.{source_finder}.{image_type}.residual.fits",self.deep_image)
        model_file    = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",f".hydra.{source_finder}.{image_type}.model.fits",self.deep_image)
        #catalogue.pprint_all()
        img.create_residual_image(catalogue,residual_file,model_file)

        reg_file    = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",f".hydra.{source_finder}.{image_type}.reg",self.deep_image)
        print_ok(f"Creating Region File: {reg_file}")
        reg = list()
        reg.append(f"# Hydra Region File")
        reg.append(f"# Mode: Single_Select")
        reg.append(f"#    o Image File: {self.deep_image}")
        reg.append(f"#    o Residual File: {residual_file}")
        reg.append(f"#    o Model File: {model_file}")
        reg.append(f"#")
        reg.append(f"global color={self.plot_colors[source_finder]}")
        reg.append(f"fk5")
        for row in catalogue:
            reg.append(f"ellipse {row['ra'].to(u.deg).value} {row['dec'].to(u.deg).value} {row['extent_semiminor'].to(u.deg).value} {row['extent_semimajor'].to(u.deg).value} {row['extent_angle'].to(u.deg).value}")
        with open(reg_file,'w') as fd:
            fd.write("\n".join(reg))

        #qt = self.__clusterize(source_finder,image_type)
        #clump_ids = list()
        #for clump_id in qt['clump_id']:
        #    if not clump_id in clump_ids:
        #        clump_ids.append(clump_id)
        ##print_ok(qt)
        ##qt.pprint_all()
        ##print_warning(np.unique(list(map(lambda m: m in clump_ids,qt['clump_id']))))
        ##print_ok(len(np.unique(list(map(lambda m: m in clump_ids,qt['clump_id'])))))
        #print(len(clump_ids))
        ##qt = qt[list(filter(lambda m: m in clump_ids,list(qt['clump_id'])))]
        ##qt = qt[(qt['match_id'] in clump_ids)]
        ##print_warning(qt)
        #qt = self.__get_subclusters()
        #print_ok(qt)
        #qt.pprint_all()
        #self.__make_subclump_hist()

        return self



@click.command()
@click.argument('fits_file',nargs=1)
# TO-DO: This use business needs some work... perhaps it should just reference the FITS_FILE
@click.option('--catalogues-yml',default=None,type=str,help="YML file pointing to catalogoues to match.")
@click.option('--bypass-archived',default=False,is_flag=True,help="Use FITS_FILE_PREFIX.hydra.tar.gz archive.")
@click.option('--use',default=None,type=str,help="Use FITS_FILE_PREFIX.typhon.tar.gz optimization.")
def hydra(
    fits_file,
    catalogues_yml,
    bypass_archived,
    use
):
    """\b
        Performs deep-shallow analaysis.
    """
    print_ok(f"Processing: {fits_file}")
    hydra = Hydra(fits_file,use,catalogues_yml,not bypass_archived) \
        .print_cluster_table() \
        .print_metrics_table() \
        .print_clump_table()

    #hydra.mk_pybdsf_residual()
    #hydra.make_typhon_stats_plots()

    #hydra.make_hist_plots()
    #hydra.make_completeness_plots()
    #hydra.make_completeness_plots(is_injected=False)
    #hydra.make_reliability_plots()
    #hydra.make_reliability_plots(is_injected=False)
    #hydra.plot_collection_test()
    #hydra.make_goodness_of_completeness_plot()
    #hydra.make_goodness_of_reliability_plot()
    #hydra.make_flux_ratio_plots(is_sn=False,is_injected=True)
    #hydra.make_flux_ratio_plots(is_sn=False,is_injected=False)
    #hydra.make_flux_ratio_plots(is_injected=True)
    #hydra.make_flux_ratio_plots(is_injected=False)
    #hydra.make_delta_completeness_plot()
    #hydra.make_delta_reliability_plot()


if __name__ == "__main__":
    hydra()
