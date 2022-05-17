#!/usr/bin/env Rscript

###############################################################################################
# 
#    * * *   P R O F O U N D   C O N T A I N E R  P R O C E S S I N G   S C R I P T    * * *
#

# initialization
library(magicaxis)
library(FITSio)
library(ProFound)
library(reticulate)
use_python("/usr/bin/python3",required=TRUE)
source_python("python_tools.py")


###############################################################################################
#
#    * * *   C O N F I G U R A T I O N   * * *
#

# profound defaults config
# Note: ProFound is an optical source finder, so skyckut has been
#       adjusted for VLASS, as a starting point. Its true defualt
#       is 1.
profound_cfg <- list(
    skyckut    = 2.82,
    tolerance  = 4.0
)

# catalogue units
catalogue_header <- list(
    segID           = '',
    #groupID         = '',
    uniqueID        = '',
    RAcen           = 'deg',
    Deccen          = 'deg',
    RAmax           = 'deg',
    Decmax          = 'deg',
    xcen            = '',
    ycen            = '',
    xsd             = 'pixel',
    ysd             = 'pixel',
    covxy           = 'pixel',
    corxy           = 'pixel',
    semimaj         = 'deg',
    semimin         = 'deg',
    ang             = 'deg',
    axrat           = '',
    xmax            = '',
    ymax            = '',
    sep             = '',
    flux            = 'Jy',
    flux_err        = 'Jy',
    mag             = '',
    mag_err         = '',
    cenfrac         = '',
    asymm           = '',
    flux_reflect    = '',
    mag_reflect     = '',
    N50             = '',
    N90             = '',
    N100            = '',
    R50             = '',
    R90             = '',
    R100            = '',
    SB_N50          = '',
    SB_N90          = '',
    SB_N100         = '',
    con             = '',
    signif          = '',
    FPlim           = '',
    flux_err_sky    = 'Jy',
    flux_err_skyRMS = 'Jy',
    flux_err_shot   = 'Jy',
    flux_err_cor    = 'Jy',
    cor_seg         = '',
    sky_mean        = '',
    sky_sum         = '',
    skyRMS_mean     = '',
    Nedge           = '',
    Nsky            = '',
    Nobject         = '',
    Nborder         = '',
    Nmask           = '',
    edge_frac       = '',
    edge_excess     = '',
    flag_border     = ''
)


###############################################################################################
#
#    * * *   U T I L I T Y   S C R I P T S   * * *
#

is.symlink <- function(paths) isTRUE(nzchar(Sys.readlink(paths), keepNA=TRUE))

get_pix_to_deg_conversion_factor <- function(tile) {
    fd <- file(tile,'rb')
    hdr <- parseHdr(readFITSheader(fd))
    close(fd)
    cdelt1 <- abs(as.double(hdr[which(hdr=='CDELT1')+1]))
    cdelt2 <- abs(as.double(hdr[which(hdr=='CDELT2')+1]))
    bpa    <- pi*as.double(hdr[which(hdr=='BPA')+1])/180.0
    list(
        fmaj = sqrt((cdelt1*sin(bpa))**2+(cdelt2*cos(bpa))**2),
        fmin = sqrt((cdelt1*cos(bpa))**2+(cdelt2*sin(bpa))**2)
    )
}

# catalogue generator
extract_catalogue <- function(profound_results) {
    # get group ids -- this takes a bit...
    sg <- profound_results$segstats
    get_groupID <- function(segID) {
       gp <- profound_results$group$groupsegID
       get_gid <- function(sid) {
          for (groupID in gp$groupID) {
             if (sid %in% gp[(gp$groupID==groupID),]$segID[[1]])
                return(groupID)
          }
          -1
       }
       ids <- c()
       for (sid in segID) {
          ids <- c(ids,get_gid(sid))
       }
       ids
    }
    #gids <- get_groupID(sg$segID)
    
    # pix convertion routine
    #pix_to_deg <- 1.0
    pix_to_deg <- get_pix_to_deg_conversion_factor(dst)
    
    # create data frame
    cast <- function(datum) {
      if (is.null(datum)) {
        rep('',length(sg$segID))
      } else {
        datum
      }
    }

    get_flux_adu_shift <- function(tile) {
        fd <- file(tile,'rb')
        hdr <- parseHdr(readFITSheader(fd))
        close(fd)
	cdelt1 <- abs(as.double(hdr[which(hdr=='CDELT1')+1]))
	cdelt2 <- abs(as.double(hdr[which(hdr=='CDELT2')+1]))
	bpa    <- pi*as.double(hdr[which(hdr=='BPA')+1])/180.0
	bmaj <- as.double(hdr[which(hdr=='BMAJ')+1])*sqrt((cdelt1*sin(bpa))**2+(cdelt2*cos(bpa))**2)/(cdelt1*cdelt2)
	bmin <- as.double(hdr[which(hdr=='BMIN')+1])*sqrt((cdelt1*cos(bpa))**2+(cdelt2*sin(bpa))**2)/(cdelt1*cdelt2)
        bmaj*bmin*pi/(4.0*log(2.0))
    }
    flux_correction_factor <- get_flux_adu_shift(dst)
    
    # return the catalogue
    data.frame(
        segID            = cast(sg$segID),
        #groupID          = cast(gids),
        uniqueID         = cast(sg$uniqueID),
        RAcen            = cast(sg$RAcen),
        Deccen           = cast(sg$Deccen),
        RAmax            = cast(sg$RAmax),
        Decmax           = cast(sg$Decmax),
        xcen             = cast(sg$xcen),
        ycen             = cast(sg$ycen),
        xsd              = cast(sg$xsd),
        ysd              = cast(sg$ysd),
        covxy            = cast(sg$covxy),
        corxy            = cast(sg$corxy),
	# SEMIMAJ:
        # dfn: Weighted standard deviation along the major axes (i.e.,
        #      the semi-major first moment, so ~2 times this would  be
        #      a typical major axis Kron radius) in units of pix.
	# Convert to typical form in degrees.
        semimaj          = cast(2.0*pix_to_deg$fmaj*sg$semimaj), # 2-sigma
	# SEMIMIN:
        # dfn: Weighted standard deviation along the minor axes (i.e.,
        #      the semi-minor first moment, so ~2 times this would be
        #      a typical minor axis Kron radius) in units of pix 
	# Convert to typical form in degrees.
        semimin          = cast(2.0*pix_to_deg$fmin*sg$semimin), # 2-sigma
	# PA:
        # dfn: The orientation of the semi-major axis in degrees (this
	#      has the convention that 0 = | (vertical), 45 = \\,
	#      90 = - (horizontal), 135  = /, 180 = | (vertical)).
	# Reorient to typical form.
        #ang              = cast(sg$ang-90.0),
        ang              = cast(sg$ang),
        axrat            = cast(sg$axrat),
        xmax             = cast(sg$xmax),
        ymax             = cast(sg$ymax),
        sep              = cast(sg$sep),
        flux             = cast(sg$flux/flux_correction_factor),
        flux_err         = cast(sg$flux_err/flux_correction_factor),
        mag              = cast(sg$mag),
        mag_err          = cast(sg$mag_err),
        cenfrac          = cast(sg$cenfrac),
        asymm            = cast(sg$asymm),
        flux_reflect     = cast(sg$flux_reflect),
        mag_reflect      = cast(sg$mag_reflect),
        N50              = cast(sg$N50),
        N90              = cast(sg$N90),
        N100             = cast(sg$N100),
        R50              = cast(sg$R50),
        R90              = cast(sg$R90),
        R100             = cast(sg$R100),
        SB_N50           = cast(sg$SB_N50),
        SB_N90           = cast(sg$SB_N90),
        SB_N100          = cast(sg$SB_N100),
        con              = cast(sg$con),
        signif           = cast(sg$signif),
        FPlim            = cast(sg$FPlim),
        flux_err_sky     = cast(sg$flux_err_sky),
        flux_err_skyRMS  = cast(sg$flux_err_skyRMS),
        flux_err_shot    = cast(sg$flux_err_shot),
        flux_err_cor     = cast(sg$flux_err_cor),
        cor_seg          = cast(sg$cor_seg),
        sky_mean         = cast(sg$sky_mean),
        sky_sum          = cast(sg$sky_sum),
        skyRMS_mean      = cast(sg$skyRMS_mean),
        Nedge            = cast(sg$Nedge),
        Nsky             = cast(sg$Nsky),
        Nobject          = cast(sg$Nobject),
        Nborder          = cast(sg$Nborder),
        Nmask            = cast(sg$Nmask),
        edge_frac        = cast(sg$edge_frac),
        edge_excess      = cast(sg$edge_excess),
        flag_border      = cast(sg$flag_border),
        stringsAsFactors = FALSE
    )
}


# Rscript profound.R data/ processing/ results/ J165209+212528_s3arcmin_VLASS.fits --fits --skycut 2.2941176470588234 --tolerance 2.2918235294117646
# Rscript profound.R data/ processing/ results/ emu_simulated_04.fits --fits --skycut 3.7222222222222223 --tolerance 3.7185
is_polygons <- FALSE # TO-DO: TRUSE case is too slow -- needs work.
extract_segment_info <- function(fits_image_file,profound_results) { 
    if (is_polygons) {
        if (!is.null(profound_results$segim)) { # guard condition
            # get seg info
            si <- ProFoundSegments(fits_image_file,profound_results$segim)
            si.segments   <- si$get_segments()
            si.boundaries <- si$get_boundaries()
        
            # get stats
            n_easy <- 0
            n_hard <- 0
            for (id in si$get_segment_ids()) {
                if (si$get_boundary(id)$is_hard) {
                    n_hard <- n_hard + 1
                } else {
                    n_easy <- n_easy + 1
                }
            }
        } else {
            si.segments   <- c()
            si.boundaries <- c()
            n_easy <- 0
            n_hard <- 0
        }
        
        # create region file contents
        reg <- c()
        reg <- c(reg,"# ProFound Segments to Polygons Stats:")
        reg <- c(reg,"# ====================================")
        reg <- c(reg,sprintf("# Easy Method: %d",n_easy))
        reg <- c(reg,sprintf("# Hard Method: %d",n_hard))
        reg <- c(reg,sprintf("# Total: %d",n_easy + n_hard))
        reg <- c(reg,"#")
        reg <- c(reg,"global color=red")
        reg <- c(reg,"fk5")
        for (s in si.boundaries) {
            pg <- c("polygon")
            for (i in 1:length(s$ra)) {
                if ( (!is.na(s$ra[i])) && (!is.na(s$dec[i])) ) {
                    pg <- c(pg,s$ra[i],s$dec[i])
                }
            }
            if (length(pg > 3)) {
                pg <- paste(pg,collapse=" ")
                reg <- c(reg,pg)
            }
        }
        
        list(
            segments = si.segments,
            polygons = si.boundaries,
            stats = list(
                n_easy = n_easy,
                n_hard = n_hard
            ),
            ds9 = reg
        )
    } else {
        # create region file contents
	sg <- profound_results$segstats
        pix_to_deg <- get_pix_to_deg_conversion_factor(fits_image_file)
        reg <- c()
        reg <- c(reg,"# ProFound Simple Region File.")
        reg <- c(reg,"#")
        reg <- c(reg,"global color=red")
        reg <- c(reg,"fk5")
	if (length(sg$RAcen)>0) {
	    for (i in 1:length(sg$RAcen)) {
	        reg <- c(reg,paste(
	            'ellipse',
	            sg$RAcen[i],
	            sg$Deccen[i],
	            2.0*pix_to_deg$fmaj*sg$semimaj[i],
	            2.0*pix_to_deg$fmin*sg$semimin[i],
	            sprintf("%fd",sg$ang[i]-90.0),
	            collapse=" "
	        ))
	    }
	}
        
        list(
            segments = list(),
            polygons = list(),
            stats = list(),
            ds9 = reg
        )
    }
}


###############################################################################################
#
#    * * *  C R U D E   C O M A N D - L I N E   I N T E R F A C E   * * *
#

# crude argumenat parser / help utility
args = commandArgs(trailingOnly=TRUE)
flags <- c()
dirs  <- c()
is_skip <- FALSE
for (i in 1:length(args)) {
    if (args[i]  %in% c('--skycut','--tolerance','--fits','--residual','--dump','--help')) {
        flags <- c(flags,args[i])
        if (args[i] == '--skycut' && !is.na(as.numeric(args[i+1]))) {
            profound_cfg['skyckut'] <- as.numeric(args[i+1])
	    is_skip <- TRUE
        }
        if (args[i] == '--tolerance' && !is.na(as.numeric(args[i+1]))) {
            profound_cfg['tolerance'] <- as.numeric(args[i+1])
	    is_skip <- TRUE
        }
    } else if (!is_skip) {
	dirs <- c(dirs,args[i])
    } else {
        is_skip <- FALSE
    }
}
pad <- function(n,text) {
    paste(paste(rep(" ",n),collapse=''),text,collapse='')
}
help <- c(
    "Usage: profound.R [OPTIONS] INPUT_DIR PROCESSING_DIR OUTPUT_DIR FITS_IMAGE_FILE",
    "   ",
    "   ProFound image processing tool.",
    "   ",
    "   inputs:",
    "   ",
    "         INPUT_DIR: location of image_filename.fits",
    "         PROCESSING_DIR: location of scratch directory",
    "         OUTPUT_DIR: location to place results",
    "         FITS_IMAGE_FILE: image_filename.fits (without path)",
    "   ",
    "   outputs:",
    "   ",
    "         OUTPUT_DIR/image_filename.profound.fits",
    "         OUTPUT_DIR/image_filename.profound.reg",
    "   ",
    "   Options:",
    "     --skycut    FLOAT Island threshold (in skyRMS).",
    pad(22,sprintf("[default: %.2f]",profound_cfg['skyckut'])),
    "     --tolerance FLOAT Define island separation height.",
    pad(22,sprintf("(in sykRMS). [default: %.2f]",profound_cfg['tolerance'])),
    "     --fits            Output FITS catalogue. [default: CSV]",
    "     --residual        Output residual and model FITS files.",
    "     --dump            Dump out all processing files.",
    "     --help            Show this message and exit."
)
if ('--help' %in% flags || length(dirs) != 4) {
    cat(help,sep="\n")
    quit(save="no",status=0)
}
input_dir       <- gsub("/+$","",dirs[1])
processing_dir  <- gsub("/+$","",dirs[2])
output_dir      <- gsub("/+$","",dirs[3])
fits_image_file <- dirs[4]
fits     <- if ('--fits'     %in% flags) TRUE else FALSE
residual <- if ('--residual' %in% flags) TRUE else FALSE
dump     <- if ('--dump'     %in% flags) TRUE else FALSE
# -- phew!


###############################################################################################
#
#    * * *   M A I N   S C R I P T    * * *
#

# summarize input paramters
src <- sprintf("%s/%s",input_dir,fits_image_file)
t_str <- c(
    sprintf("> Local Context: %s",src),
    sprintf(">   INPUT_DIR:  %s",input_dir),
    sprintf(">   PROCESSING_DIR: %s",processing_dir),
    sprintf(">   OUTPUT_DIR: %s",output_dir),
    sprintf("> Supported Flags: "),
    sprintf(">    --skycut: %.2f",profound_cfg['skyckut']),
    sprintf(">    --tolerance: %.2f",profound_cfg['tolerance'])
)
cat(t_str,sep="\n")

# check if image file exists
if (!file.exists(src)) {
    e_str <- c(
        sprintf("ERROR: Image file '%s' not found!",src),
	"Bye!"
    )
    cat(e_str,sep="\n")
    quit(save="no",status=0)
}

# link image file from data to processing dir
dst <- sprintf("%s/%s",processing_dir,fits_image_file)
t_str <- c(
    "> Linking image file:",
    sprintf("> $ ln -s %s %s",dst,src)
)
cat(t_str,sep="\n")
# this gaurd condition is usefull for debugging.
#if (!is.symlink(dst) || !file.exists(dst)) {
if (!is.symlink(dst) && !file.exists(dst)) {
    file.symlink(src,dst)
}

# execute profound
t_str <- c(
    "> Running ProFound:",
    ">    profound_results <- profoundProFound(",
             sprintf(">        image      = %s,",dst),
             sprintf(">        skyckut    = %f,",profound_cfg$skyckut),
             sprintf(">        tolerance  = %f,",profound_cfg$tolerance),
             sprintf(">        rotstats   = TRUE,"),
             sprintf(">        boundstats = TRUE,"),
             sprintf(">        nearstats  = TRUE,"),
             sprintf(">        groupstats = TRUE,"),
             sprintf(">        groupby    = 'segim'"),
    ">    )"
)
cat(t_str,sep="\n")
image <- readFITS(dst)
if (length(dim(image$imDat)) > 2) {
    idxs <- ""
    for (i in 1:(length(dim(image$imDat))-2)) {
      idxs <- sprintf("%s,1",idxs)
    }
    image$imDat <- eval(parse(text=sprintf("image$imDat[,%s]",idxs)))
}
# OK, we need to replace all of the 1.13.1 defaults with 1.10.8!
box        <- c(100,100)
boundstats <- TRUE
# TO-DO: Rscript --verbose profound.R data processing results J165209+212528_s3arcmin_VLASS.sample.inverted.fits --skycut 4.0 --tolerance 2.0
time.start <- Sys.time()
profound_results <- tryCatch({
        profoundProFound(
            image         = image,
            segim         = NULL,
            objects       = NULL,
            mask          = NULL,
            skycut        = profound_cfg$skyckut, # default: 1
            pixcut        = 3,
            tolerance     = profound_cfg$tolerance, # default: 4
            ext           = 2,
            reltol        = 0,
            cliptol       = Inf,
            sigma         = 1,
            smooth        = TRUE,
            #SBlim,
            SBdilate      = NULL,
            SBN100        = 100,
            size          = 5,
            shape         = "disc",
            iters         = 6,
            threshold     = 1.05,
            magzero       = 0,
            gain          = NULL,
            pixscale      = 1,
            sky           = NULL,
            skyRMS        = NULL,
            redosegim     = FALSE,
            redosky       = TRUE,
            redoskysize   = 21,
            box           = box, # default: c(100,100)
            grid          = box,
            type          = "bicubic",
            skytype       = "median",
            skyRMStype    = "quanlo",
            roughpedestal = FALSE,
            sigmasel      = 1,
            skypixmin     = prod(box)/2,
            boxadd        = box/2,
            boxiters      = 0,
            iterskyloc    = TRUE,
            deblend       = FALSE,
            df            = 3,
            radtrunc      = 2,
            iterative     = FALSE,
            doclip        = TRUE,
            shiftloc      = FALSE,
            paddim        = TRUE,
            #header,
            verbose       = TRUE,
            plot          = FALSE,
            stats         = TRUE,
            rotstats      = TRUE, # default: FALSE
            boundstats    = boundstats, # default: FALSE,
            nearstats     = boundstats,
            groupstats    = boundstats,
            group         = NULL,
            groupby       = 'segim', # default: 'segim_orig'
            offset        = 1,
            haralickstats = FALSE,
            sortcol       = "segID",
            decreasing    = FALSE,
            lowmemory     = FALSE,
            keepim        = TRUE,
            watershed     = 'ProFound',
            pixelcov      = FALSE,
            deblendtype   = 'fit',
            psf           = NULL,
            fluxweight    = 'sum',
            convtype      = 'brute',
            convmode      = 'extended',
            fluxtype      = 'Raw',
            app_diam      = 1,
            Ndeblendlim   = Inf
        )
    }, warning = function(var) {
        cat(paste("WARNING:",var,"\n"))
	NULL
    }, error = function(err) {
        cat(paste("ERROR:",err,"\n"))
	NULL
    }
)
time.profound <- paste(Sys.time() - time.start)
cat("Profound Processing Time:",time.profound,"\n")
flush.console()

# output profound version
cat("ProFound Version:\n")
profound_results$ProFound.version

# get profound segment polygon outlines
cat("> Extracting Polygon Outlines from ProFound Segments:\n")
img_file <- sprintf("%s/%s",processing_dir,fits_image_file)
time.start = Sys.time()
regions <- extract_segment_info(img_file,profound_results)
time.polygon <- paste(Sys.time() - time.start)
cat("Polygon Extraction Time:",time.polygon,"\n")
flush.console()
if (length(regions$stats)>0) {
    t_str <- c(
        ">> # ProFound Segments to Polygons Stats:",
        ">> # ====================================",
        sprintf(">> # Easy Method: %d",regions$stats$n_easy),
        sprintf(">> # Hard Method: %d",regions$stats$n_hard),
        sprintf(">> # Total: %d",regions$stats$n_easy+regions$stats$n_hard),
        ">> #"
    )
}
cat(t_str,sep="\n")
cat("> [Done]\n")
flush.console()

# write the catalogue
csv_file = sprintf("%s/%s",processing_dir,gsub("\\.[Ff][Ii][Tt]([Ss]|)$",".profound.csv",fits_image_file))
t_str <- c(
    sprintf("> Creating CSV Catalogue: %s",csv_file)
)
cat(t_str,sep="\n")
time.start = Sys.time()
catalogue <- extract_catalogue(profound_results)
time.catalogue <- paste(Sys.time() - time.start)
cat("Catalogue Extraction Time:",time.catalogue,"\n")
flush.console()
write.table(x=catalogue,file=csv_file,sep=",",row.names=FALSE,quote=FALSE)
if (fits) {
    fits_file = sprintf("%s/%s",processing_dir,gsub("\\.[Ff][Ii][Tt]([Ss]|)$",".profound.fits",fits_image_file))
    t_str <- c(
        sprintf("> Converting CSV Catalogue to FITS: %s",fits_file)
    )
    cat(t_str,sep="\n")
    time.start <- Sys.time()
    csv_to_fits(csv_file,catalogue_header,regions$segments,regions$polygons)
    time.csv2fits <- paste(Sys.time()-time.start)
    cat_file = sprintf("%s/%s",output_dir,gsub("\\.[Ff][Ii][Tt]([Ss]|)$",".profound.fits",fits_image_file))
} else {
    cat_file = sprintf("%s/%s",output_dir,gsub("\\.[Ff][Ii][Tt]([Ss]|)$",".profound.csv",fits_image_file))
}
t_str <- c(
    sprintf("> Outputing Catalogue: %s",cat_file),
    sprintf(">> $ cp %s %s",if (fits) fits_file else csv_file,cat_file)
)
cat(t_str,sep="\n")
flush.console()
status <- file.copy(if (fits) fits_file else csv_file,cat_file,overwrite=TRUE)

# create region file
reg_file <- sprintf("%s/%s",output_dir,gsub("\\.[Ff][Ii][Tt]([Ss]|)$",".profound.reg",fits_image_file))
cat("> Creating Region File:",reg_file,"\n")
flush.console()
fd <- file(reg_file)
writeLines(regions$ds9,fd,"\n")
close(fd)
warnings()

# residauls
if (residual) {
    time.start <- Sys.time()
    create_residual_image(profound_results$segim,fits_image_file,processing_dir,output_dir)
    time.residual <- paste(Sys.time()-time.start)
}

# print of time stats info
time.total <-  as.double(time.profound) + 
    as.double(time.polygon) + 
    as.double(time.catalogue) + 
    (if (fits) as.double(time.csv2fits) else 0) + 
    (if (residual) as.double(time.residual) else 0)
times <- c("Time Stats:")
times <- c(times,paste("  Profound:",time.profound,"sec"))
times <- c(times,paste("   Polygon:",time.polygon,"sec"))
times <- c(times,paste(" Catalogue:",time.catalogue,"sec"))
times <- c(times,paste("  CSV2FITS:",if (fits) time.csv2fits else "n/a",if (fits) "sec"))
times <- c(times,paste("  Residaul:",if (residual) time.residual else "n/a",if (residual) "sec"))
times <- c(times,paste("Total:",time.total,"sec"))
cat(times,sep="\n")

cat("[Done]\n")
flush.console()


