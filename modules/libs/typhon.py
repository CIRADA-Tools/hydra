import re
import os
import sys
import gzip
import glob
import time
sys.path.insert(0,"../")
from pathlib import Path
import tarfile
import yaml as yml
from astropy.table import Table
from astropy.table import QTable
import astropy.units as u
from libs.config_tools import Config
from libs.config_tools import cfg_get_build_context
from libs.config_tools import cfg_fetch_frequency
from libs.fits_tools   import sigma_clipper
from libs.fits_tools   import invert_fits_image
from libs.fits_tools   import create_sample
from libs.fits_tools   import FITS
from libs.fits_tools   import load_bane_rms_image
from libs.exceptions   import print_warning
from libs.exceptions   import print_ok
from libs.exceptions   import HydraIOError
from libs.cerberus     import cerberus
from math import ceil


class Typhon:
    def __init__(
        self,
        fits_image_file,
        processing_cache,
        use=None,
        is_diagnostics=False,
        is_fits_catalogue=True,
        is_residual=True,
        is_optimize_rms_box = True,
        percent_real_detections_cut=98.0 # (Defualt: C. Hale, et al, cut)
    ):
        # threshold cut
        self.percent_real_detections_cut = percent_real_detections_cut

        # make sure we have the fits image file
        self.fits_file = fits_image_file
        if not os.path.isfile(self.fits_file):
            print_warning(f"ERROR: File '{fits_file}' does not exist!")
            print_ok("Bye!")
            raise HydraIOError

        # create working cache
        self.cache = processing_cache
        if not os.path.isdir(self.cache):
            print_ok(f"> Creating Cache: {self.cache}")
            os.makedirs(self.cache)
        else:
            print_ok(f"> Using Cache: {self.cache}")

        # define output files
        self.optimization_file = re.sub(r"^(.*?/)*(.*?)\.[Ff][Ii][Tt](|[Ss])$",r"%s/\2.optimizations.yml" % self.cache,self.fits_file)
        self.tar_gz_file = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".typhon.tar.gz",self.fits_file)
        self.local_tar_dir = re.sub(r"^(.*?/)*","",re.sub(r"\.tar\.gz$","_dir",self.tar_gz_file))

        # handle 'use' (!=None) reference optimization case
        self.reference_tar_gz_file = None
        self.reference_fits_file   = None
        if not use is None and use != self.tar_gz_file:
            # set flags (nb: is_bail is set by programmer, depending on choice of behavoir)
            is_reset = False # True if can't get optimization parameters
            is_bail  = True  # policy flag (i.e., what to do if is_reset = True):
                             #     is_bail = True  => Bail
                             #     is_bail = False => Use input fits_image_file

            # check if tar.gz file exists
            self.reference_tar_gz_file = use
            if os.path.isfile(self.reference_tar_gz_file): # check if tar.gz file has optimization file
                has_opt_file = False
                opt_file = re.sub(r"^(.*?/)*(.*?)\.typhon\.tar\.gz$",r"\2.optimizations.yml",self.reference_tar_gz_file)
                with tarfile.open(self.reference_tar_gz_file,'r') as fd:
                    for member in fd.getmembers():
                        fname = re.sub(r"^(.*?/)*","",member.name)
                        if opt_file == fname:
                            has_opt_file = True
                            break
                if not has_opt_file:
                    is_reset = True
            else: # check if corresponding fits file exists, so we can generate an optimzation file
                reference_fits_file_prefix = re.sub(r"\.typhon\.tar\.gz$","",self.reference_tar_gz_file)
                for f in glob.glob(f"{reference_fits_file_prefix}.*"):
                    if re.search(r"^%s\.[Ff][Ii][Tt](|[Ss])$" % re.escape(reference_fits_file_prefix),f):
                        self.reference_fits_file = f
                        break
                if self.reference_fits_file is None:
                    self.reference_tar_gz_file = None
                    is_reset = True

            if is_reset and is_bail:
                print_warning(f"ERROR: Cannot use reference optimization: {self.reference_tar_gz_file}")
                print_ok("Bye!")
                raise HydraIOError

        # set diagnostics flag
        self.is_diagnostics = is_diagnostics

        # set catalogue output format flag
        self.is_fits_catalogue = is_fits_catalogue

        # rms box optimization flag
        self.is_optimize_rms_box = is_optimize_rms_box

        # set create residual and model .fits files flag
        self.is_residual = is_residual


    def optimize(self,use_existing=False):
        print_ok(f"Running Optimizer...")
        tar_gz_file = self.tar_gz_file if self.reference_tar_gz_file is None else self.reference_tar_gz_file
        if use_existing: 
            opt_file = re.sub(r"^(.*?/)*(.*?)\.typhon\.tar\.gz$",r"\2.optimizations.yml",tar_gz_file)
            if os.path.isfile(tar_gz_file): 
                with tarfile.open(tar_gz_file,'r') as fd:
                    for member in fd.getmembers():
                        fname = re.sub(r"^(.*?/)*","",member.name)
                        if opt_file == fname:
                            optimizations = yml.load(fd.extractfile(member).read())
                            for module in optimizations:
                                print_ok(f">> {module}:")
                                for key,value in optimizations[module].items():
                                    print_ok(f">>   {key}: {value}")
                            print_ok("[Done]")
                            return optimizations
            print_ok(f"> Can't find optimization file '{opt_file}' in '{tar_gz_file}'.")
        files = list()
        fits_file = self.fits_file if self.reference_fits_file is None else self.reference_fits_file
        print_ok(f"> Optimizing: {fits_file}")
        fits_sample_file = create_sample(
            input_file  = fits_file,
            output_file = re.sub(r"^(.*?/)*(.*?)\.[Ff][Ii][Tt](|[Ss])$",r"%s/\2.sample.fits" % self.cache,fits_file),
            #size        = (2.0*u.deg,2.0*u.deg)
            size        = (2.5*u.deg,2.5*u.deg)
        )
        files.append(fits_sample_file)
        print_ok(f"> Using: {fits_sample_file}")
        try:
            stats = sigma_clipper(fits_sample_file)
            rms  = stats['rms']*10**6
            mean = stats['mean']*10**6
            median = stats['median']*10**6
            initial_min = stats['initial_min']*10**6
            initial_max = stats['initial_max']*10**6
            final_min = stats['final_min']*10**6
            final_max = stats['final_max']*10**6
            pdfs = [p*100 for p in stats['pdfs']]
            sigmas = [s*10**6 for s in stats['sigmas']]
        except Exception as e:
            print_warning(f"ERROR: {e}")
            print_warning("Whoops! This can't be a fits file.")
            print_ok("Bye!")
            raise HydraIOError
     
        print_ok(f"> Images statistics:")
        print_ok(f">> Median: {median} \u03bcJy")
        print_ok(f">> Mean:   {mean} \u03bcJy")
        print_ok(f">> RMS:    {rms} \u03bcJy")
        print_ok(f">> Convergence(\u03c3's): {sigmas} \u03bcJy")
        print_ok(f">>         ==> PDFs: {pdfs} %")
        print_ok(f">> Initial Range: [{initial_min}, {initial_max}] \u03bcJy")
        print_ok(f">> Final Range:   [{final_min}, {final_max}] \u03bcJy")
    
        n_sigma_minimal = (final_max-mean)/rms
        n_sigma_maximal = (initial_max-mean)/rms
        #n_sigma_max = min(1.5*float(ceil(n_sigma_minimal)),n_sigma_maximal)
        n_sigma_max = float(ceil(n_sigma_minimal))
        #n_sigma_min = 1.5 if 1.5 < n_sigma_max else n_sigma_max/2.0
        n_sigma_min = 1.75 if 1.75 < n_sigma_max else n_sigma_max/2.0
        print_ok(f"> Thresholding:")
        print_ok(f">> Minimal n_\u03C3 = {n_sigma_minimal}")
        print_ok(f">> Maximal n_\u03C3 = {n_sigma_maximal}")
        print_ok(f">> Range n_\u03C3: [{n_sigma_min}, {n_sigma_max}]")
    
        # added inverted image to cache
        fits_inverted_file = invert_fits_image(fits_sample_file,self.cache)
        files.append(fits_inverted_file)
    
        # define container external mount points
        input_path  = os.path.abspath(re.sub(r"^((.*?/)*).*$",r"\1",fits_sample_file))
        output_path = self.cache
        print_ok(f"> Container external mount points:")
        print_ok(f">> o inputs:  {input_path}")
        print_ok(f">> o outputs: {output_path}")
    
        # ok, let's do it!
        cfg = Config()
        sigma_range = [n_sigma_min,n_sigma_max]
        statistics = list()
        optimizations = dict()
        statistics_file = re.sub(r"^(.*?/)*(.*)\.fits$",r"%s/\2.statistics.csv" % self.cache,fits_sample_file)
        files.append(statistics_file)

        rms_box_stats = None
        if self.is_optimize_rms_box:
            print_ok(f"Doing RMS Box optimation...")
            rms_box_stats_file = self.cache+"/"+re.sub("^(.*?/)+","",re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".rms_box_statistcs.csv",fits_file))
            load_bane_rms_image(fits_file,self.cache,None,True)
            print_ok(f">> rms_box_stats_file: {rms_box_stats_file}")
            qt = Table().read(rms_box_stats_file)
            print_ok(f">> Summary:")
            print_ok(qt)
            rms_box_stats = {
                'file': rms_box_stats_file,
                'stats': qt.copy(),
                'box_size': int(qt[(qt['is_opt']!=0)]['box_size'][0]),
                'step_size': int(qt[(qt['is_opt']!=0)]['step_size'][0]),
            }
            files.append(rms_box_stats_file)
            print_ok(f">> Using:")
            print_ok(f">>  o Box_Size: {rms_box_stats['box_size']}")
            print_ok(f">>  o Step_Size: {rms_box_stats['step_size']}")
            print_ok(f"> [Done]")

        for module,constraints in cfg.get_constraints().items():
            #if module != 'pybdsf':
            #    continue
            print_ok(f"> Optimizing: {module}")
            rms_params = rms_param_default = rms_param_name = isl_params = isl_param_default = isl_param_name = None
            for param,value in constraints['parameters'].items():
                if isinstance(value,str) and value == 'rms':
                    rms_params = sigma_range
                    rms_param_name = param
                elif isinstance(value,list):
                    isl_params = [min(value),max(value)]
                    isl_param_name = param
            if 'defaults' in constraints:
                for param,value in constraints['defaults'].items(): 
                    if param == rms_param_name:
                        rms_param_default = value
                    elif param == isl_param_name:
                        isl_param_default = value
            if None in [rms_params, rms_param_default, rms_param_name, isl_params, isl_param_default, isl_param_name]:
                print_warning(f">> Skipping {module}, missing parameter(s).")
                continue
            is_ok = eval(f"lambda {rms_param_name},{isl_param_name}: {constraints['directive']}")
            print_ok(f">> o {rms_param_name}: {rms_params} [default: {rms_param_default}]")
            print_ok(f">> o {isl_param_name}: {isl_params} [default: {isl_param_default}]")
            print_ok(f">> o directive: 'lambda {rms_param_name},{isl_param_name}: {constraints['directive']}'")
    
            # define cerberus lauch function
            def exec_cerberus(fits_file,rms_par,isl_par,is_residual=False):
                #if is_residual:
                #    cerberus([
                #        module,
                #        fits_file,
                #        "--output-dir", output_path,
                #        f"--{re.sub(r'_','',rms_param_name)}".lower(), rms_par,
                #        f"--{re.sub(r'_','',isl_param_name)}".lower(), isl_par,
                #        f"--residual"
                #    ], standalone_mode=False)
                #else:
                #    cerberus([
                #        module,
                #        fits_file,
                #        "--output-dir", output_path,
                #        f"--{re.sub(r'_','',rms_param_name)}".lower(), rms_par,
                #        f"--{re.sub(r'_','',isl_param_name)}".lower(), isl_par,
                #    ], standalone_mode=False)
                cerberus_pars = [
                    module,
                    fits_file,
                    "--output-dir", output_path,
                    f"--{re.sub(r'_','',rms_param_name)}".lower(), rms_par,
                    f"--{re.sub(r'_','',isl_param_name)}".lower(), isl_par,
                ]
                if is_residual:
                    cerberus_pars.append(f"--residual")
                if self.is_optimize_rms_box and cfg.has_rms_box_pars(module):
                    cerberus_pars.extend(["--boxsize", rms_box_stats['box_size']])
                    cerberus_pars.extend(["--stepsize", rms_box_stats['step_size']])
                cerberus(cerberus_pars,standalone_mode=False)
    
            # define module output csv files
            csv_inverted_file = re.sub(r"\.fits$",f".{module}.csv",fits_inverted_file)
            csv_file = re.sub(r"\.inverted(\.%s\.csv)$" % module,r"\1",csv_inverted_file)

            # define residual/model output files
            res_file = re.sub(r"\.inverted\.(%s)\.csv$" % module,r".\1.residual.fits",csv_inverted_file)
            mod_file = re.sub(r"\.inverted\.(%s)\.csv$" % module,r".\1.residual.model.fits",csv_inverted_file)
    
            # define rms iteration parameters
            rms_min = min(rms_params)
            rms_max = max(rms_params)
            ##rms_iters = max(int((rms_max-rms_min)/0.25),4)
            #rms_iters = min(max(int((rms_max-rms_min)/0.25),4),20)
            rms_iters = min(max(int((rms_max-rms_min)/0.125),4),20)
            d_rms = (rms_max-rms_min)/(rms_iters-1)
            print_ok(f">> RMS_ITERS: {rms_iters} over [{rms_min},{rms_max}] @ \u03B4(rms)={d_rms}")
    
            # create constrained isl iteration list func
            isl_min = min(isl_params)
            isl_max = max(isl_params)
            #isl_iters = int((isl_max-isl_min)/0.5)
            isl_iters = min(int((isl_max-isl_min)/0.5),5)
            d_isl = (isl_max-isl_min)/(isl_iters-1)
            print_ok(f">> ISL_ITERS: {isl_iters} over [{isl_min},{isl_max}] @ \u03B4(isl)={d_isl}")
            def get_isl_pars(rms_par):
                isl_pars = list()
                for j in range(isl_iters):
                    isl_par = float(j)*d_isl+isl_min
                    if not is_ok(rms_par,isl_par):
                        isl_par = 0.999*rms_par
                    isl_pars.append(isl_par)
                return sorted(list(set(isl_pars)))
    
            idx = 1
            reg_file = os.path.abspath(re.sub(r"\.fits$",f".{module}.reg",re.sub(r"^(.*?/)*",f"{self.cache}/",fits_sample_file)))
            reg_inverted_file = re.sub(r"(\.%s\.reg)$" % module,r".inverted\1",reg_file)
            is_next_module = False
            rms_par_old = None
            isl_par_old = None
            percent_real_detections_old = None
            percent_real_detections_final = None
            for i in range(rms_iters-1,-1,-1):
                rms_par = float(i)*d_rms+rms_min
                if rms_par_old is None:
                    rms_par_old = rms_par
                for isl_par in get_isl_pars(rms_par):
                    if isl_par_old is None:
                        isl_par_old = isl_par
                    print_ok(f"> Running {module} container: /w")
                    print_ok(f">> o {rms_param_name}={rms_par}")
                    print_ok(f">> o {isl_param_name}={isl_par}")
                    sys.stdout.flush()
                    t0 = time.time()
                    exec_cerberus(fits_sample_file,rms_par,isl_par,True)
                    dt_image = time.time()-t0
                    res_img = FITS(res_file)
                    os.remove(mod_file)
                    os.remove(res_file)
                    sys.stdout.flush()
                    t0 = time.time()
                    exec_cerberus(fits_inverted_file,rms_par,isl_par)
                    dt_invim = time.time()-t0
                    n_normal   = sum(1 for line in open(csv_file))-1
                    n_inverted = sum(1 for line in open(csv_inverted_file))-1
                    percent_real_detections = (n_normal-n_inverted)*100.0/n_normal if n_normal != 0 else 0
                    statistics.append({
                        'module': module,
                        'rms_par': rms_par,
                        'isl_par': isl_par,
                        'file_id': idx,
                        'n_image': n_normal,
                        'n_invim': n_inverted,
                        'prd': percent_real_detections,
                        'rms': res_img.get_rms(),
                        'madfm': res_img.get_madfm(),
                        'sumsq': res_img.get_sumsq(),
                        'dt_image': dt_image,
                        'dt_invim': dt_invim,
                    })
                    del res_img
                    print_ok(f"> Stats:")
                    print_ok(f">> o N_image:          {n_normal} (dt: {dt_image}s)")
                    print_ok(f">> o N_inverted_image: {n_inverted} (dt: {dt_invim}s)")
                    print_ok(f">> % Real Dectections: {percent_real_detections}")
                    print_ok(f"> Renaming:")
                    for fname in [reg_file,reg_inverted_file,csv_file,csv_inverted_file]:
                        idx_fname = re.sub(r"\.%s\.(csv|reg)$" % module,r".%s.%s.\1" % (module,f"{idx:02}"),fname)
                        print_ok(f">> {fname} => {idx_fname}")
                        Path(fname).rename(idx_fname)
                        files.append(idx_fname)
                    idx += 1
                    print_ok(f"> Updating statistics file: {statistics_file}")
                    with open(statistics_file,"w+") as fd:
                        contents = [",".join(statistics[0].keys())]
                        for item in statistics:
                            contents.append(",".join([f"{v}" for v in item.values()]))
                        fd.write("\n".join(contents))
                    #if not is_next_module and percent_real_detections < self.percent_real_detections_cut:
                    if not is_next_module and \
                            not percent_real_detections_old is None and \
                            percent_real_detections_old >= self.percent_real_detections_cut and \
                            percent_real_detections < self.percent_real_detections_cut:
                        is_next_module = True
                        percent_real_detections_final = percent_real_detections
                        if not self.is_diagnostics:
                            break
                    if not is_next_module:
                        isl_par_old = isl_par
                    percent_real_detections_old = percent_real_detections
                if is_next_module and not self.is_diagnostics:
                    break
                if not is_next_module:
                    rms_par_old = rms_par
            if percent_real_detections_final is None:
                max_prd = None
                rms_opt = None
                isl_opt = None
                for datum in statistics:
                    if datum['module'] == module:
                        if max_prd is None:
                            max_prd = datum['prd']
                            rms_opt = datum['rms_par']
                            isl_opt = datum['isl_par']
                        elif max_prd < datum['prd']: 
                            max_prd = datum['prd']
                            rms_opt = datum['rms_par']
                            isl_opt = datum['isl_par']
                if not max_prd is None:
                    rms_par_old = rms_opt
                    isl_par_old = isl_opt
            if self.is_optimize_rms_box and cfg.has_rms_box_pars(module):
                optimizations[module] = {
                    rms_param_name: rms_par_old,
                    isl_param_name: isl_par_old,
                    'boxsize':  rms_box_stats['box_size'],
                    'stepsize': rms_box_stats['step_size'],
                }
            else:
                optimizations[module] = {
                    rms_param_name: rms_par_old,
                    isl_param_name: isl_par_old,
                }
                
            # print module summary
            print_ok("> ***")
            print_ok(f"> Summary statistics: {module}")
            print(f">>  rms   isl   prd")
            for datum in statistics:
                if datum['module'] == module:
                    print_ok(f">> {datum['rms_par']:.3f} {datum['isl_par']:.3f} {datum['prd']}")
            print_ok("> ***")
            print_ok(f"> FILES:")
            for f in files:
                if re.search(module,f):
                    print_ok(f">> o {f}")
            print_ok("> ***")
    
        # create a yml file containing our optimal source finder (module) parameters
        optimization_file = re.sub(r"^(.*?/)*(.*?)\.typhon\.tar\.gz$",r"%s/\2.optimizations.yml" % self.cache,tar_gz_file)
        print_ok(f"> Creating optimization file: {optimization_file}")
        with open(optimization_file,"w") as fd:
            yml.dump(optimizations,fd)
        with open(optimization_file,"r") as fd:
            for line in fd:
                line = re.sub(r"\n","",line)
                print_ok(f">> {line}")
    
        # ok, let's archive our work
        print_ok(f"> Archiving: {tar_gz_file}")
        with tarfile.open(tar_gz_file,"w:gz") as fd:
            # store processing run files under local_tar_dir/optimization/
            for f in files:
                print_ok(f">> o {f}")
                fname = f"{self.local_tar_dir}/optimization/{re.sub(r'^(.*?/)*','',f)}"
                fd.add(f,fname,False)
                os.remove(f)
            # store optimization file under local_tar_dir/
            print_ok(f">> o {optimization_file}")
            fname = f"{self.local_tar_dir}/{re.sub(r'^(.*?/)*','',optimization_file)}"
            fd.add(optimization_file,fname,False)
            os.remove(optimization_file)
        print_ok("> ***")
    
        print_ok(f"[Done]")
        return optimizations


    def get_optimizations(self):
        return self.optimize(use_existing=True)


    def process(self):
        # define some helper functions -- makes code look cleaner/modularized
        def pretty_print(params):
            t_str = ""
            t_str += f"> Typhon is calling Cerberus:\n"
            t_str += f">> cerberus::process(\n"
            t_str += f">>     fits_image_file = '{self.fits_file}',\n"
            t_str += f">>     output_dir = '{self.cache}',\n"
            for module in params:
                for key,value in params[module].items():
                    t_str += f">>     {module}_{re.sub(r'_','',key).lower()} = {value},"+"\n"
            t_str = re.sub(r",\n$","\n",t_str)
            t_str += f">> )"
            print_ok(t_str)
        def build_args(params):
            args = list()
            args.append("process")
            args.append(self.fits_file)
            args.extend(['--output-dir',self.cache])
            for module in params:
                for key,value in params[module].items():
                    args.extend([f"--{module}-{re.sub(r'_','',key).lower()}",value])
            if self.is_fits_catalogue:
                args.append("--fits")
            if self.is_residual:
                args.append("--residual")
            return args
        def get_file_list(params):
            files = list()
            residuals = list()
            prefix = re.sub(r"^(.*?/)*(.*?)\.[Ff][Ii][Tt](|[Ss])$",r"%s/\2" % self.cache,self.fits_file)
            for module in params:
                files.append(f"{prefix}.{module}.{'fits' if self.is_fits_catalogue else 'csv'}")
                files.append(f"{prefix}.{module}.reg")
                if self.is_residual:
                    residuals.append(f"{prefix}.{module}.residual.fits")
                    residuals.append(f"{prefix}.{module}.residual.model.fits")
            files.extend(residuals)
            return files

        # run cerberus
        print_ok(f"Processing: {self.fits_file}")
        optimizations = self.optimize(use_existing=True)
        pretty_print(optimizations)
        args = build_args(optimizations)
        cerberus(args,standalone_mode=False)

        # update/create archive (nb: assumes unqiue root-filenames)
        print_ok(f"> Appending archive: {self.tar_gz_file}")
        files_to_add = get_file_list(optimizations) # optimization output products
        if self.reference_tar_gz_file is None: # update
            # open .tar.gz archive
            with tarfile.open(self.tar_gz_file,'r') as f_zipped:
                # create .tar archive
                tar_file = re.sub(r"\.gz","",self.tar_gz_file)
                with tarfile.open(tar_file,'w') as f_unzipped:
                    # copy evertying from .tar.gz into .tar except old optimization output products
                    replacements = [re.sub(r"^(.*?/)*","",f) for f in files_to_add]
                    for member in f_zipped.getmembers():
                        fname = re.sub(r"^(.*?/)*","",member.name)
                        if not fname in replacements:
                            f_unzipped.addfile(member,f_zipped.extractfile(member.name))
                    # add optimization output products under local_tar_dir/{catalogues|residuals}/ dir
                    for f in files_to_add:
                        print_ok(f">> o {f}")
                        subdir = 'residuals' if re.search(r'\.residual(\.model|)\.fits$',f) else 'catalogues'
                        fname = f"{self.local_tar_dir}/{subdir}/{re.sub(r'^(.*?/)*','',f)}"
                        f_unzipped.add(f,fname,False)
                        os.remove(f)
            # zip the .tar archive overtop of the old .tar.gz archive
            with gzip.open(self.tar_gz_file,'wb') as f_zipped:
                with open(tar_file,'rb') as f_unzipped:
                    f_zipped.writelines(f_unzipped)
            os.remove(tar_file) 
        else: # create
            with tarfile.open(self.tar_gz_file,'w:gz') as fd:
                # add optimization output products under local_tar_dir/{catalogues|residuals}/ dir
                for f in files_to_add:
                    print_ok(f">> o {f}")
                    subdir = 'residuals' if re.search(r'\.residual(\.model|)\.fits$',f) else 'catalogues'
                    fname = f"{self.local_tar_dir}/{subdir}/{re.sub(r'^(.*?/)*','',f)}"
                    fd.add(f,fname,False)
                    os.remove(f)
                # store optimization file under local_tar_dir/ dir
                print_ok(f">> o {self.optimization_file}")
                with open(self.optimization_file,"w") as opt_fd:
                    yml.dump(optimizations,opt_fd)
                fname = f"{self.local_tar_dir}/{re.sub(r'^(.*?/)*','',self.optimization_file)}"
                fd.add(self.optimization_file,fname,False)
                os.remove(self.optimization_file)

        print("[Done]")


    def get_archive_name(self):
        return self.tar_gz_file


    def get_catalogue_names(self):
        modules = set(Config().get_modules())
        cat_suffix = 'fits' if self.is_fits_catalogue else 'csv'
        cat_prefix = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$","",self.fits_file)
        return [re.sub(r"^(.*?/)*","",f"{cat_prefix}.{module}.{cat_suffix}") for module in modules]


    def has_catalogues(self):
        targzf = self.tar_gz_file
        if not targzf is None and os.path.isfile(targzf) and tarfile.is_tarfile(targzf):
            modules = set(Config().get_modules())
            with tarfile.open(targzf,'r') as fd:
                cat_suffix = 'fits' if self.is_fits_catalogue else 'csv'
                catalogues = list(filter(lambda n: re.search(r"/catalogues/.*?\.%s$" % cat_suffix,n),fd.getnames()))
                return modules == set([re.sub(r"^.*?(%s).fits$" % "|".join(modules),r"\1",c) for c in catalogues])
        return False


    def load_catalogues(self):
        modules = {m: {} for m in Config().get_modules()} 
        if not self.has_catalogues():
            self.process()
        with tarfile.open(self.tar_gz_file) as fd:
            cat_suffix = 'fits' if self.is_fits_catalogue else 'csv'
            pat = re.compile(r"^.*?/catalogues/(.*?(\.(%s)\.%s)$)" % ("|".join(list(modules.keys())),cat_suffix))
            for member in fd.getmembers():
                if pat.search(member.name):
                    c_name = pat.sub(r"\1",member.name)
                    module = pat.sub(r"\3",member.name)
                    df = Table.read(fd.extractfile(member),format="fits")
                    modules[module] = {
                        'catalogue': c_name,
                        'astropy_table': df
                    }
        return modules


    def load_optimization_table(self):
        if not self.has_catalogues():
            self.process()
        fname = re.sub(r"\.optimizations\.yml$",".sample.statistics.csv",re.sub(r"^(.*?/)+","",self.optimization_file))
        pat = re.compile(r"%s$" % re.escape(fname))
        with tarfile.open(self.tar_gz_file) as fd:
            for member in fd.getmembers():
                if pat.search(member.name):
                    return Table.read(fd.extractfile(member),format="csv")
        return Table()


    def load_rms_box_optimization_table(self):
        if not self.has_catalogues():
            self.process()
        fname = re.sub(r"\.optimizations\.yml$",".rms_box_statistcs.csv",re.sub(r"^(.*?/)+","",self.optimization_file))
        pat = re.compile(r"%s$" % re.escape(fname))
        with tarfile.open(self.tar_gz_file) as fd:
            for member in fd.getmembers():
                if pat.search(member.name):
                    qt = QTable.read(fd.extractfile(member),format="csv")
                    units = {
                        'box_size': u.pixel,
                        'step_size': u.pixel,
                        'mean': u.Jy,
                        'median': u.Jy,
                        'rms': u.Jy,
                        'madfm': u.Jy,
                        'mad': u.Jy,
                        'sumsq': u.Jy**2,
                    }
                    for col in units:
                        if col in qt.columns:
                            qt[col] = qt[col]*units[col]
                            if units[col] == u.Jy:
                                qt[col] = qt[col].to(u.mJy)
                            elif units[col] == u.Jy**2:
                                qt[col] = qt[col].to(u.mJy**2)
                    return qt
        return Table()

       



