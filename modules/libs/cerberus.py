###########################################################################################################
#                                                                                                         #
#                        W A R N I N G :   A U T O G E N E R A T E D   S C R I P T !                      #
#                                                                                                         #
###########################################################################################################
# Cenerated: 2021-04-17

import re
import os
import sys
sys.path.insert(0,"../")
from libs.fits_tools import FITS
from libs.config_tools import cfg_get_main
from libs.config_tools import cfg_get_build_context
from libs.config_tools import cfg_fetch_frequency
from libs.exceptions import print_warning
from libs.exceptions import print_ok
from libs.exceptions import HydraContainerError
from pathlib import Path

# command line
import click
from libs.click_tools import SpecialHelpOrder
from libs.click_tools import HelpCommandPriority

cmd=HelpCommandPriority()
@click.group(cls=SpecialHelpOrder)
def cerberus():
    """\b
       Runs a collection of source finder modules.

       Cerberus: a multiheaded dog, son of Typhon.
    """
    pass




@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--aegean-seedclip',default=5.0,type=float,help="The clipping value (in σ's) for seeding islands.")
@click.option('--aegean-floodclip',default=4.0,type=float,help="The clipping value (in σ's) for growing islands.")
@click.option('--aegean-boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--aegean-stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--caesar-seedthr',default=5.0,type=float,help="Blob finding threshold (in σ's).")
@click.option('--caesar-mergethr',default=2.6,type=float,help="Blob growth threshold (in σ's)")
@click.option('--profound-skycut',default=2.82,type=float,help="Island threshold (in skyRMS).")
@click.option('--profound-tolerance',default=4.0,type=float,help="Define island separation height.")
@click.option('--pybdsf-threshpix',default=5.0,type=float,help="Island threshold boundary (in σ's).")
@click.option('--pybdsf-threshisl',default=3.0,type=float,help="Source detection threshold (in σ's).")
@click.option('--pybdsf-frequency',default=None,type=float,help="Input frequency (optional: if in fits header).")
@click.option('--pybdsf-boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--pybdsf-stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--selavy-snrcut',default=4.0,type=float,help="Threshold value (in σ's) above µ-noise.")
@click.option('--selavy-growthcut',default=3.0,type=float,help="SNR to grow detections down to.")
@click.option('--selavy-boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--selavy-stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def process(
    fits_image_file,
    output_dir,
    aegean_seedclip,
    aegean_floodclip,
    aegean_boxsize,
    aegean_stepsize,
    caesar_seedthr,
    caesar_mergethr,
    profound_skycut,
    profound_tolerance,
    pybdsf_threshpix,
    pybdsf_threshisl,
    pybdsf_frequency,
    pybdsf_boxsize,
    pybdsf_stepsize,
    selavy_snrcut,
    selavy_growthcut,
    selavy_boxsize,
    selavy_stepsize,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using multiple source finders.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    # pybdsf requires input frequency -- make sure we have it
    fd = FITS(fits_image_file)
    if pybdsf_frequency is None and not fd.has_frequency():
        survey = fd.get_survey()
        pybdsf_frequency = cfg_fetch_frequency(survey)
        if not pybdsf_frequency is None:
            msg = ""
            msg += f"> WARNING: Frequency missing from header: {fits_image_file}\n"
            msg += f">> Using {survey} configuration default:\n"
            msg += ">>   FREQ: {} Hz".format(pybdsf_frequency)
            print_warning(msg)
        else:
            msg = ""
            msg += f"> ERROR: Frequency missing from header: {fits_image_file}\n"
            msg += f">> Require Input:\n"
            msg += ">>    --{}".format(re.sub(r'_','-','pybdsf_frequency'))
            print_warning(msg)
            print_ok("Bye!")
            return
    del fd # free some memory

    cfg = cfg_get_main()
    
    # list the modules being used
    print_ok(f"> Using modules:")
    for module in cfg['modules']:
        key = list(module.keys())[0]
        print_ok(f">   o {module[key]['name']}")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'aegean':
            dcr_cmd += "" if aegean_seedclip is None else " --seedclip {}".format(aegean_seedclip)
            dcr_cmd += "" if aegean_floodclip is None else " --floodclip {}".format(aegean_floodclip)
            dcr_cmd += "" if aegean_boxsize is None else " --box-size {}".format(aegean_boxsize)
            dcr_cmd += "" if aegean_stepsize is None else " --step-size {}".format(aegean_stepsize)
        elif container == 'caesar':
            dcr_cmd += "" if caesar_seedthr is None else " --seedThr {}".format(caesar_seedthr)
            dcr_cmd += "" if caesar_mergethr is None else " --mergeThr {}".format(caesar_mergethr)
        elif container == 'profound':
            dcr_cmd += "" if profound_skycut is None else " --skycut {}".format(profound_skycut)
            dcr_cmd += "" if profound_tolerance is None else " --tolerance {}".format(profound_tolerance)
        elif container == 'pybdsf':
            dcr_cmd += "" if pybdsf_threshpix is None else " --thresh-pix {}".format(pybdsf_threshpix)
            dcr_cmd += "" if pybdsf_threshisl is None else " --thresh-isl {}".format(pybdsf_threshisl)
            dcr_cmd += "" if pybdsf_frequency is None else " --frequency {}".format(pybdsf_frequency)
            dcr_cmd += "" if pybdsf_boxsize is None else " --box-size {}".format(pybdsf_boxsize)
            dcr_cmd += "" if pybdsf_stepsize is None else " --step-size {}".format(pybdsf_stepsize)
        elif container == 'selavy':
            dcr_cmd += "" if selavy_snrcut is None else " --snrCut {}".format(selavy_snrcut)
            dcr_cmd += "" if selavy_growthcut is None else " --growthCut {}".format(selavy_growthcut)
            dcr_cmd += "" if selavy_boxsize is None else " --box-size {}".format(selavy_boxsize)
            dcr_cmd += "" if selavy_stepsize is None else " --step-size {}".format(selavy_stepsize)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")




@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--seedclip',default=5.0,type=float,help="The clipping value (in σ's) for seeding islands.")
@click.option('--floodclip',default=4.0,type=float,help="The clipping value (in σ's) for growing islands.")
@click.option('--boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def aegean(
    fits_image_file,
    output_dir,
    seedclip,
    floodclip,
    boxsize,
    stepsize,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using aegean source finder.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    

    cfg = cfg_get_main()
    
    # print module being used
    print_ok(f"> Using module:")
    print_ok(f">   o aegean")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        if container != 'aegean':
            continue
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'aegean':
            dcr_cmd += "" if seedclip is None else " --seedclip {}".format(seedclip)
            dcr_cmd += "" if floodclip is None else " --floodclip {}".format(floodclip)
            dcr_cmd += "" if boxsize is None else " --box-size {}".format(boxsize)
            dcr_cmd += "" if stepsize is None else " --step-size {}".format(stepsize)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")



@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--seedthr',default=5.0,type=float,help="Blob finding threshold (in σ's).")
@click.option('--mergethr',default=2.6,type=float,help="Blob growth threshold (in σ's)")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def caesar(
    fits_image_file,
    output_dir,
    seedthr,
    mergethr,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using caesar source finder.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    

    cfg = cfg_get_main()
    
    # print module being used
    print_ok(f"> Using module:")
    print_ok(f">   o caesar")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        if container != 'caesar':
            continue
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'caesar':
            dcr_cmd += "" if seedthr is None else " --seedThr {}".format(seedthr)
            dcr_cmd += "" if mergethr is None else " --mergeThr {}".format(mergethr)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")



@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--skycut',default=2.82,type=float,help="Island threshold (in skyRMS).")
@click.option('--tolerance',default=4.0,type=float,help="Define island separation height.")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def profound(
    fits_image_file,
    output_dir,
    skycut,
    tolerance,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using profound source finder.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    

    cfg = cfg_get_main()
    
    # print module being used
    print_ok(f"> Using module:")
    print_ok(f">   o profound")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        if container != 'profound':
            continue
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'profound':
            dcr_cmd += "" if skycut is None else " --skycut {}".format(skycut)
            dcr_cmd += "" if tolerance is None else " --tolerance {}".format(tolerance)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")



@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--threshpix',default=5.0,type=float,help="Island threshold boundary (in σ's).")
@click.option('--threshisl',default=3.0,type=float,help="Source detection threshold (in σ's).")
@click.option('--frequency',default=None,type=float,help="Input frequency (optional: if in fits header).")
@click.option('--boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def pybdsf(
    fits_image_file,
    output_dir,
    threshpix,
    threshisl,
    frequency,
    boxsize,
    stepsize,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using pybdsf source finder.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    # pybdsf requires input frequency -- make sure we have it
    fd = FITS(fits_image_file)
    if frequency is None and not fd.has_frequency():
        survey = fd.get_survey()
        frequency = cfg_fetch_frequency(survey)
        if not frequency is None:
            msg = ""
            msg += f"> WARNING: Frequency missing from header: {fits_image_file}\n"
            msg += f">> Using {survey} configuration default:\n"
            msg += ">>   FREQ: {} Hz".format(frequency)
            print_warning(msg)
        else:
            msg = ""
            msg += f"> ERROR: Frequency missing from header: {fits_image_file}\n"
            msg += f">> Require Input:\n"
            msg += ">>    --{}".format(re.sub(r'_','-','frequency'))
            print_warning(msg)
            print_ok("Bye!")
            return
    del fd # free some memory

    cfg = cfg_get_main()
    
    # print module being used
    print_ok(f"> Using module:")
    print_ok(f">   o pybdsf")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        if container != 'pybdsf':
            continue
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'pybdsf':
            dcr_cmd += "" if threshpix is None else " --thresh-pix {}".format(threshpix)
            dcr_cmd += "" if threshisl is None else " --thresh-isl {}".format(threshisl)
            dcr_cmd += "" if frequency is None else " --frequency {}".format(frequency)
            dcr_cmd += "" if boxsize is None else " --box-size {}".format(boxsize)
            dcr_cmd += "" if stepsize is None else " --step-size {}".format(stepsize)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")



@cerberus.command(help_priority=cmd.prio())
@click.argument('fits_image_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
@click.option('--snrcut',default=4.0,type=float,help="Threshold value (in σ's) above µ-noise.")
@click.option('--growthcut',default=3.0,type=float,help="SNR to grow detections down to.")
@click.option('--boxsize',default=None,type=int,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--stepsize',default=None,type=int,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and module FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def selavy(
    fits_image_file,
    output_dir,
    snrcut,
    growthcut,
    boxsize,
    stepsize,
    fits,
    residual,
    dump
):
    """\b
       Process FITS_IMAGE_FILE using selavy source finder.
    """
    # ok, let's do this!
    print_ok(f"Processing: {fits_image_file}")

    # make sure file exist before we start spawning containers
    if not os.path.isfile(fits_image_file):
        msg = ""
        msg += f"> ERROR: File does not exist: {fits_image_file}"
        print_warning(msg)
        print_ok("Bye!")
        return

    

    cfg = cfg_get_main()
    
    # print module being used
    print_ok(f"> Using module:")
    print_ok(f">   o selavy")
    

    # get file-name and full-paths
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    output_path = Path(output_dir).resolve() if not output_dir is None else input_path

    # loop over source finders and process image
    context = cfg_get_build_context()
    for module in cfg['modules']:
        container = list(module.keys())[0]
        
        if container != 'selavy':
            continue
        args = [re.sub(r" +","",c) for c in context['services'][container]['build']['args']]
        args = {p.split('=')[0]:p.split('=')[1] for p in args}
        # TO-DO: Simplify
        if input_path == output_path:
            dcr_cmd = "docker run --rm -t -v {}:{} {} {} {} {} {}".format(
                input_path,args['input_dir'],   # data in-out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['input_dir'],              # output dir
                input_file                      # image file
            )
        else:
            dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
                input_path, args['input_dir'],  # data in mount point
                output_path,args['output_dir'], # data out mount point
                container,                      # service
                args['input_dir'],              # input dir
                args['processing_dir'],         # processing dir
                args['output_dir'],             # output dir
                input_file                      # image file
            )
        
        if container == 'selavy':
            dcr_cmd += "" if snrcut is None else " --snrCut {}".format(snrcut)
            dcr_cmd += "" if growthcut is None else " --growthCut {}".format(growthcut)
            dcr_cmd += "" if boxsize is None else " --box-size {}".format(boxsize)
            dcr_cmd += "" if stepsize is None else " --step-size {}".format(stepsize)
        dcr_cmd += "" if residual is False else " --residual"
        dcr_cmd += "" if fits is False else " --fits"
        dcr_cmd += "" if dump is False else " --dump"
        print_ok(f"> Launching {module[container]['name']} Container:")
        print_ok(f"> $ {dcr_cmd}")
        if os.system(dcr_cmd) != 0:
            print_warning(f"Container Error: {module[container]['name']}")
            raise HydraContainerError
        # debug
        #if container == 'profound':
        #    os.system(dcr_cmd)

    print_ok(f"[Done]")



if __name__ == "__main__":
    cerberus()