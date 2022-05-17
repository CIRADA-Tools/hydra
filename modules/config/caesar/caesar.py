import re
import os
import sys
import glob
import click
from shutil import copy2 as copy
from astropy.table import Table
from astropy.table import QTable
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import astropy.units as u


catalogue_header = {
    'name': None,
    'npix': None,
    'componentId': None,
    'iauName': None,
    'x': None,
    'y': None,
    'x_err': None,
    'y_err': None,
    'x_wcs': u.deg,
    'y_wcs': u.deg,
    'x_wcs_err': u.deg,
    'y_wcs_err': u.deg,
    'nu': u.GHz,
    'Speak': u.Jy/u.beam,
    'Speak_err': u.Jy/u.beam,
    'S': u.Jy,
    'S_err': u.Jy,
    'S_island': u.Jy,
    'S_island_err': u.Jy,
    'beamArea': None,
    'bmaj': None,
    'bmin': None,
    'pa': u.deg,
    'bmaj_err': None,
    'bmin_err': None,
    'pa_err': u.deg,
    'bmaj_wcs': u.arcsec,
    'bmin_wcs': u.arcsec,
    'pa_wcs': u.deg,
    'bmaj_wcs_err': u.arcsec,
    'bmin_wcs_err': u.arcsec,
    'pa_wcs_err': u.deg,
    'bmaj_beam': u.arcsec,
    'bmin_beam': u.arcsec,
    'pa_beam': u.deg,
    'bmaj_deconv_wcs': u.arcsec,
    'bmin_deconv_wcs': u.arcsec,
    'pa_deconv_wcs': u.deg,
    'fitBeamEllipseEccentricityRatio': None,
    'fitBeamEllipseAreaRatio': None,
    'fitBeamEllipseRotAngle': u.deg,
    'bkgSum': u.Jy/u.beam,
    'rmsSum': u.Jy/u.beam,
    'chi2': None,
    'ndf': None,
    'fitQuality': None,
    'fitComponentFlag': None,
    'fitComponentType': None,
}

# root residual filename
global_residaul_filename = "residual.fits"

# caesar defaults config
caesar_cfg = {
    'seedThr': 5.0,
    'mergeThr':  2.6,
}
#caesar_cfg = {
#    ############################################
#    ###    CAESAR CONFIG OPTIONS
#    ############################################
#    'comment_001': "###",
#    'comment_002': "##=====================================================================================================================================================================",
#    'comment_003': "##==                                                                               INPUT                                                                             ==",
#    'comment_004': "##=====================================================================================================================================================================",
#    'inputFile': 'dummy_image_file.fits',                                 # Input image filename (.root#.fits)
#    'inputImage': 'img',                                                  # Input image name in ROOT file
#    'readTileImage': 'false',                                             # Read sub-image (T#F)
#    'tileMinX': 0,                                                        # Min x coords to be read in image
#    'tileMaxX': 0,                                                        # Max x coords to be read in image
#    'tileMinY': 0,                                                        # Min y coords to be read in image
#    'tileMaxY': 0,                                                        # Max y coords to be read in image
#    'comment_005': "###",
#    'comment_006': "###",
#    'comment_007': "##=====================================================================================================================================================================",
#    'comment_008': "##==                                                                             BEAM INFO                                                                           ==",
#    'comment_009': "##=====================================================================================================================================================================",
#    'pixSize': 1.0,                                                       # User-supplied map pixel area in arcsec (pixSize=CDELT, default=1 arcsec)
#    #'beamFWHM': 5,                                                        # User-supplied circular beam FWHM in arcsec (beamFWHM=BMAJ=BMIN, default=6.5 arcsec)
#    'beamFWHM': 6.5,                                                      # User-supplied circular beam FWHM in arcsec (beamFWHM=BMAJ=BMIN, default=6.5 arcsec)
#    'beamBmaj': 10,                                                       # User-supplied elliptical beam bmaj FWHM in arcsec (default=10 arcsec)
#    'beamBmin': 5,                                                        # User-supplied elliptical beam bmin FWHM in arcsec (default=5 arcsec)
#    'beamTheta': 0,                                                       # User-supplied beam theta in deg (default=0)
#    'comment_010': "###",
#    'comment_011': "###",
#    'comment_012': "##=====================================================================================================================================================================",
#    'comment_013': "##==                                                                       DISTRIBUTED PROCESSING                                                                    ==",
#    'comment_014': "##=====================================================================================================================================================================",
#    'nThreads': 1,                                                        # Number of threads used if OPENMP is enabled (-1=all available threads)
#    'splitInTiles': 'false',                                              # Split input image in tiles (default=false)
#    #'tileSizeX': 0,                                                       # Size of tile X (in pixels) to partition the input image
#    'tileSizeX': 1000,                                                    # Size of tile X (in pixels) to partition the input image
#    'tileSizeY': 1000,                                                    # Size of tile Y (in pixels) to partition the input image *NEW*
#    'useTileOverlap': 'false',                                            # Allow for tile overlap
#    'tileStepSizeX': 1,                                                   # Tile step size fraction X to partition the input image (1=no overlap,0.5=half overlap, ...)
#    'tileStepSizeY': 1,                                                   # Tile step size fraction Y to partition the input image (1=no overlap,0.5=half overlap, ...)
#    #'mergeSourcesAtEdge': 'false',                                        # Merge sources found at tile edge by each worker (default=true)
#    'mergeSourcesAtEdge': 'ture',                                         # Merge sources found at tile edge by each worker (default=true)
#    #'mergeSources': 'true',                                               # Merge overlapping sources found by each worker (default=false)
#    'mergeSources': 'false',                                              # Merge overlapping sources found by each worker (default=false)
#    'comment_015': "###",
#    'comment_016': "###",
#    'comment_017': "##=====================================================================================================================================================================",
#    'comment_018': "##==                                                                              LOGGING                                                                            ==",
#    'comment_019': "##=====================================================================================================================================================================",
#    'loggerTarget': 1,                                                    # Logger target (1=CONSOLE, 2=FILE, 3=SYSLOG)
#    'loggerTag': 'logger',                                                # Tag given to the log messages
#    'logLevel': 'INFO',                                                   # Log level threshold (DEBUG>INFO>WARN>ERROR>FATAL)
#    'logFile': 'out.log',                                                 # Log file name (for FILE target only)
#    'appendToLogFile': 'false',                                           # If false a new log file is created, otherwise logs are appended (T#F)
#    'maxLogFileSize': '10MB',                                             # Max size of log file before rotation (e.g. 10MB, 1KB)
#    'maxBackupLogFiles': 2,                                               # Max number of backup files created after file threshold is reached
#    'consoleTarget': 'System.out',                                        # Console target (System.out#System.err)
#    'syslogFacility': 'local6',                                           # Syslog facility used with syslog target
#    'comment_020': "###",
#    'comment_021': "###",
#    'comment_022': "##=====================================================================================================================================================================",
#    'comment_023': "##==                                                                              OUTPUT                                                                             ==",
#    'comment_024': "##=====================================================================================================================================================================",
#    'isInteractiveRun': 'false',                                          # Is interactive run (graph plots enabled) (T#F)
#    'outputFile': 'out-emu_simulated_04.root',                            # Output filename (.root)
#    'outputCatalogFile': 'catalog-emu_simulated_04.dat',                  # Output catalog ascii filename
#    'outputComponentCatalogFile': 'catalog_fitcomp-emu_simulated_04.dat', # Output fitted component catalog ascii filename
#    'ds9RegionFile': 'ds9-emu_simulated_04.reg',                          # DS9 region file (.reg) where to store source catalog
#    'ds9FitRegionFile': 'ds9_fitcomp-emu_simulated_04.reg',               # DS9 region file (.reg) where to store fitted source catalog
#    'ds9RegionFormat': 2,                                                 # DS9 region format (1=ellipse, 2=polygon)
#    'convertDS9RegionsToWCS': 'false',                                    # Convert DS9 regions (contours & ellipses) to WCS (default=false)
#    'ds9WCSType': 0,                                                      # DS9 region WCS output format (0=J2000,1=B1950,2=GAL) (default=0)
#    'inputMapFITSFile': 'input_map.fits',                                 # Output filename where to store input map in FITS format (.fits)
#    'residualMapFITSFile': 'out-emu_simulated_04_res.fits',               # Output filename where to store residual map in FITS format (.fits)
#    'saliencyMapFITSFile': 'saliency_map.fits',                           # Output filename where to store saliency map in FITS format (.fits)
#    'bkgMapFITSFile': 'out-emu_simulated_04_bkg.fits',                    # Output filename where to store bkg map in FITS format (.fits)
#    'noiseMapFITSFile': 'out-emu_simulated_04_rms.fits',                  # Output filename where to store noise map in FITS format (.fits)
#    'significanceMapFITSFile': 'out-emu_simulated_04_significance.fits',  # Output filename where to store significance map in FITS format (.fits)
#    'saveToFile': 'true',                                                 # Save results & maps to output ROOT file (T#F)
#    'saveToCatalogFile': 'true',                                          # Save sources to catalog files (island, fitted components) (T#F)
#    'saveToFITSFile': 'true',                                             # Save results to output FITS file(s) (T#F)
#    'saveDS9Region': 'true',                                              # Save DS9 region files (T#F) (default=T)
#    'saveConfig': 'true',                                                 # Save config options to ROOT file (T#F)
#    'saveSources': 'true',                                                # Save sources to ROOT file (T#F)
#    'saveInputMap': 'false',                                              # Save input map to ROOT file (T#F)
#    'saveBkgMap': 'true',                                                 # Save bkg map to ROOT file (T#F)
#    'saveNoiseMap': 'true',                                               # Save noise map to ROOT file (T#F)
#    'saveResidualMap': 'true',                                            # Save residual map to ROOT file (T#F)
#    'saveSignificanceMap': 'true',                                        # Save significance map to ROOT file (T#F)
#    'saveSaliencyMap': 'false',                                           # Save saliency map to ROOT file (T#F)
#    'saveSegmentedMap': 'false',                                          # Save segmented map computed in extended source search to ROOT file (T#F)
#    'saveEdgenessMap': 'false',                                           # Save edgeness map computed in extended source search to ROOT file (T#F)
#    'saveCurvatureMap': 'false',                                          # Save curvature map to ROOT file (T#F)
#    'comment_025': "###",
#    'comment_026': "###",
#    'comment_027': "##=====================================================================================================================================================================",
#    'comment_028': "##==                                                                           STATS OPTIONS                                                                         ==",
#    'comment_029': "##=====================================================================================================================================================================",
#    'useParallelMedianAlgo': 'true',                                      # Use parallel median algo (based on nth_parallel) (default=true) (T#F)
#    'comment_030': "###",
#    'comment_031': "###",
#    'comment_032': "##=====================================================================================================================================================================",
#    'comment_033': "##==                                                                             BKG OPTIONS                                                                         ==",
#    'comment_034': "##=====================================================================================================================================================================",
#    'useLocalBkg': 'true',                                                # Use local background calculation instead of global bkg (T#F)
#    'localBkgMethod': 1,                                                  # Local background method (1=Grid, 2=Superpixel)
#    'use2ndPassInLocalBkg': 'true',                                       # Use 2nd pass to refine noise calculation in local bkg (T#F)
#    'skipOutliersInLocalBkg': 'false',                                    # Skip outliers (e.g. bright point sources) in local bkg computation (T#F)
#    'bkgEstimator': 2,                                                    # Background estimator (1=Mean,2=Median,3=BiWeight,4=ClippedMedian)
#    'useBeamInfoInBkg': 'true',                                           # Use beam information in bkg box definition (if available) (T#F)
#    'boxSizeX': 10,                                                       # X Size of local background box in #pixels
#    'boxSizeY': 10,                                                       # Y Size of local background box in #pixels
#    'gridSizeX': 0.2,                                                     # X Size of local background grid used for bkg interpolation
#    'gridSizeY': 0.2,                                                     # Y Size of local background grid used for bkg interpolation
#    'sourceBkgBoxBorderSize': 20,                                         # Border size in pixels of box around source bounding box used to estimate bkg for fitting
#    'comment_035': "###",
#    'comment_036': "###",
#    'comment_037': "##=====================================================================================================================================================================",
#    'comment_038': "##==                                                                        FILTERING OPTIONS                                                                        ==",
#    'comment_039': "##=====================================================================================================================================================================",
#    'usePreSmoothing': 'true',                                            # Use a pre-smoothing stage to filter input image for extended source search (T#F)
#    'smoothFilter': 2,                                                    # Smoothing filter (1=gaus,2=guided)
#    'gausFilterKernSize': 5,                                              # Gaussian filter kernel size
#    'gausFilterSigma': 1,                                                 # Gaussian filter sigma
#    'guidedFilterRadius': 12,                                             # Guided filter radius par
#    'guidedFilterColorEps': 0.04,                                         # Guided filter color epsilon parameter
#    'comment_040': "###",
#    'comment_041': "###",
#    'comment_042': "##=====================================================================================================================================================================",
#    'comment_043': "##==                                                                       SOURCE FINDING OPTIONS                                                                    ==",
#    'comment_044': "##=====================================================================================================================================================================",
#    'searchCompactSources': 'true',                                       # Search compact sources (T#F)
#    'minNPix': 5,                                                         # Minimum number of pixel to consider a source
#    'seedThr': 5,                                                         # Seed threshold in flood filling algo for faint sources
#    #'mergeThr': 2.5,                                                      # Merge#aggregation threshold in flood filling algo
#    'mergeThr': 2.6,                                                      # Merge#aggregation threshold in flood filling algo
#    'mergeBelowSeed': 'false',                                            # Aggregate to seed only pixels above merge threshold but below seed threshold (T#F)
#    'searchNegativeExcess': 'false',                                      # Search negative excess together with positive in compact source search
#    #'compactSourceSearchNIters': 2,                                       # Number of iterations to be performed in compact source search (default=10)
#    'compactSourceSearchNIters': 1,                                       # Number of iterations to be performed in compact source search (default=10)
#    'seedThrStep': 0.5,                                                   # Seed threshold decrease step size between iteration (default=1)
#    'comment_045': "###",
#    'comment_046': "###",
#    'comment_047': "##=====================================================================================================================================================================",
#    'comment_048': "##==                                                                    NESTED SOURCE FINDING OPTIONS                                                                ==",
#    'comment_049': "##=====================================================================================================================================================================",
#    'searchNestedSources': 'true', 	                                  # Search for nested sources inside candidate sources (T#F)
#    'blobMaskMethod': 2,                                                  # Blob mask method (1=gaus smooth + laplacian,2=multi-scale LoG)
#    'sourceToBeamAreaThrToSearchNested': 10,                              # Source area#beam thr to add nested sources (e.g. npix>thr*beamArea). NB: thr=0 means always 
#                                                                          # if searchNestedSources is enabled (default=0)
#    'nestedBlobThrFactor': 0,                                             # Threshold (multiple of curvature rms) used for nested blob finding
#    'minNestedMotherDist': 2,                                             # Minimum distance in pixels (in x or y) between nested and parent blob below which nested is
#                                                                          # skipped
#    'maxMatchingPixFraction': 0.5,                                        # Maximum fraction of matching pixels between nested and parent blob above which nested is 
#                                                                          # skipped
#    'nestedBlobPeakZThr': 5,                                              # Nested blob peak significance thr (in curv map) (below thr nested blob is skipped) (default
#                                                                          # ==5 sigmas) 
#    'nestedBlobPeakZMergeThr': 2.5,                                       # Nested blob significance merge thr (in curv map) (default=2.5 sigmas)
#    'nestedBlobMinScale': 1,                                              # Nested blob min search scale (sigma=minscale x beam width) (default=1)
#    #'nestedBlobMaxScale': 2,                                              # Nested blob max search scale (sigma=maxscale x beam width) (default=3)
#    'nestedBlobMaxScale': 3,                                              # Nested blob max search scale (sigma=maxscale x beam width) (default=3)
#    'nestedBlobScaleStep': 1,                                             # Nested blob scale step (scale=minscale + step) (default=1)
#    'nestedBlobKernFactor': 6,                                            # Nested blob curvature#LoG kernel size factor f (kern size=f x sigma) (default=1)
#    'comment_050': "###",
#    'comment_051': "###",
#    'comment_052': "##=====================================================================================================================================================================",
#    'comment_053': "##==                                                                       SOURCE FITTING OPTIONS                                                                    ==",
#    'comment_054': "##=====================================================================================================================================================================",
#    'fitSources': 'true',                                                 # Deblend and fit point-like sources with multi-component gaus fit (T#F)            # NB: Simone's setting...
#    #'fitSources': 'false',                                                # Deblend and fit point-like sources with multi-component gaus fit (T#F)
#    'fitUseThreads': 'false',                                             # Use multithread in source fitting (default=false) (T#F)
#    'fitScaleDataToMax': 'false',                                         # Scale source flux data to max peak flux if true, otherwise scale to mJy units (default=false)
#    'fitMinimizer': 'Minuit2',                                            # Minimizer {Minuit,Minuit2} (default=Minuit) (T#F)
#    'fitMinimizerAlgo': 'minimize',                                       # Minimizer algorithm: {migrad,simplex,scan,minimize,fumili} (default=minimize)
#    'fitPrintLevel': 1,                                                   # Minimizer print level (default=1)
#    'fitStrategy': 2,                                                     # Minimizer strategy (higher means more accurate but more fcn calls) (default=2)
#    #'nBeamsMaxToFit': 100,                                                # Maximum number of beams in compact source for fitting (if above thr fitting not performed)
#    'nBeamsMaxToFit': 20,                                                 # Maximum number of beams in compact source for fitting (if above thr fitting not performed)
#    'fitUseNestedAsComponents': 'false',                                  # If true use nested sources (if any) to estimate fitted components, otherwise estimate blended
#                                                                          #  blobs (default=false)
#    'fitMaxNComponents': 5,                                               # Maximum number of components fitted in a blob
#    'fitWithCentroidLimits': 'true',                                      # Use limits when fitting gaussian centroid (T#F)
#    'fitCentroidLimit': 3,                                                # Source centroid limits in pixels 
#    'fixCentroidInPreFit': 'true',                                        # Fix centroid pars in prefit (T#F)
#    'fitWithBkgLimits': 'true',                                           # Use limits when fitting bkg offset (T#F)
#    'fitWithFixedBkg': 'true',                                            # Fix bkg level parameter in fit (T#F)
#    'fitUseEstimatedBkgLevel': 'true',                                    # Use estimated (avg bkg) as bkg level par in fit (T#F)
#    'fitUseBkgBoxEstimate': 'false',                                      # Use bkg estimated in a box around source (if available) as bkg level par in fit (T#F)
#    'fitBkgLevel': 0,                                                     # Fixed bkg level used when fitWithFixedBkg=true
#    'fitWithAmplLimits': 'true',                                          # Use limits when fitting gaussian amplitude (T#F)
#    'fixAmplInPreFit': 'true',                                            # Fix amplitude par in pre-fit (T#F)
#    #'fitAmplLimit': 0.5,                                                  # Flux amplitude limit around source peak (e.g. Speak*(1+-fitAmplLimit))
#    'fitAmplLimit': 0.3,                                                  # Flux amplitude limit around source peak (e.g. Speak*(1+-fitAmplLimit))
#    'fixSigmaInPreFit': 'false',                                          # Fix sigma in prefit (T#F)
#    'fitWithFixedSigma': 'false',                                         # Fix sigmas in fit (T#F)
#    'fitWithSigmaLimits': 'true',                                         # Use limits when fitting gaussian sigmas (T#F)
#    'fitSigmaLimit': 0.5,                                                 # Gaussian sigma limit around psf or beam (e.g. Bmaj*(1+-fitSigmaLimit))
#    'fitWithFixedTheta': 'false',                                         # Fix gaussian ellipse theta par in fit (T#F)
#    'fixThetaInPreFit': 'false',                                          # Fix theta in prefit (T#F)
#    'fitWithThetaLimits': 'true',                                         # Use limits when fitting gaussian theta par (T#F)
#    #'fitThetaLimit': 360,                                                 # Gaussian theta limit around psf or beam in degrees (e.g. Bpa +- fitThetaLimit)
#    'fitThetaLimit': 5,                                                   # Gaussian theta limit around psf or beam in degrees (e.g. Bpa +- fitThetaLimit)
#    'useFluxZCutInFit': 'false',                                          # Include in fit only source pixels above a given flux significance level (T#F)
#    'fitZCutMin': 2.5,                                                    # Flux significance below which source pixels are not included in the fit
#    'peakMinKernelSize': 3,                                               # Minimum dilation kernel size (in pixels) used to detect peaks (default=3)
#    'peakMaxKernelSize': 7,                                               # Maximum dilation kernel size (in pixels) used to detect peaks (default=7)
#    'peakKernelMultiplicityThr': 1,                                       # Requested peak multiplicity across different dilation kernels (-1=peak found in all given 
#                                                                          # kernels,1=only in one kernel, etc)
#    'peakShiftTolerance': 2,                                              # Shift tolerance (in pixels) used to compare peaks in different dilation kernels (default
#                                                                          # ==1 pixel)
#    'peakZThrMin': 1,                                                     # Minimum peak flux significance (in nsigmas above avg source bkg & noise) below which peak 
#                                                                          # is skipped (default=1)
#    'fitFcnTolerance': 1.e-2,                                             # Fit function minimization tolerance (default=1.e-5)
#    #'fitMaxIters': 1000000,                                               # Fit max number of iterations (default=100000)
#    'fitMaxIters': 100000,                                                # Fit max number of iterations (default=100000)
#    'fitImproveConvergence': 'true',                                      # Try to improve convergence by iterating fit if not converged or converged with pars at 
#                                                                          # limits (default=true)
#    'fitNRetries': 100,                                                   # Number of times fit is repeated (with enlarged limits) if improve convergence flag is 
#                                                                          # enabled (default=1000)
#    #'fitDoFinalMinimizerStep': 'false',                                   # Switch on#off running of final minimizer step after fit convergence with MIGRAD 
#    #                                                                      # (default=true)
#    'fitDoFinalMinimizerStep': 'true ',                                   # Switch on#off running of final minimizer step after fit convergence with MIGRAD 
#                                                                          # (default=true)
#    ###fitFinalMinimizer': 2,                                             # Final minimizer (1=MIGRAD,2=HESS,3=MINOS) (default=2)
#    #'fitChi2RegPar': 1,                                                   # Chi2 regularization par chi2=chi2_signal + regpar*chi2_bkg (default=1)
#    'fitChi2RegPar': 0,                                                   # Chi2 regularization par chi2=chi2_signal + regpar*chi2_bkg (default=1)
#    'fitParBoundIncreaseStepSize': 0.1,                                   # Par bound increase step size (e.g. parmax= parmax_old+(1+nretry)*fitParBoundIncreaseStepSize
#                                                                          # *0.5*#max-min# (default=0.1)
#    'fitRetryWithLessComponents': 'true',                                 # If fit does not converge repeat it iteratively with one component less at each cycle 
#                                                                          # (default=true)
#    'fitApplyRedChi2Cut': 'true',                                         # Apply fit Chi2#NDF cut. Used to set fit quality flag. If Chi2#NDF>cut the good fit cut 
#                                                                          # is not passed (default=true)
#    'fitRedChi2Cut': 5,                                                   # Chi2#NDF cut value (default=5)
#    'fitApplyFitEllipseCuts': 'false',                                    # Apply fit ellipse selection cuts. If not passed fit component is tagged as fake 
#                                                                          # (default=false)
#    'fitEllipseEccentricityRatioMinCut': 0.5,                             # Ellipse eccentricity ratio (fit#beam) min cut value (default=0.5)
#    'fitEllipseEccentricityRatioMaxCut': 1.5,                             # Ellipse eccentricity ratio (fit#beam) max cut value (default=1.5)
#    'fitEllipseAreaRatioMinCut': 0.01,                                    # Ellipse area ratio (fit#beam) min cut value (default=0.01)
#    'fitEllipseAreaRatioMaxCut': 10,                                      # Ellipse area ratio (fit#beam) min cut value (default=10)
#    'fitEllipseRotAngleCut': 45,                                          # Ellipse rot angle diff (#fit-beam#) cut value in degrees (default=45)
#    'comment_055': "###",
#    'comment_056': "###",
#    'comment_057': "##=====================================================================================================================================================================",
#    'comment_058': "##==                                                                      SOURCE SELECTION OPTIONS                                                                   ==",
#    'comment_059': "##=====================================================================================================================================================================",
#    'applySourceSelection': 'true',                                       # Apply selection cuts to sources (T#F)
#    #'useMinBoundingBoxCut': 'true',                                       # Use bounding box cut (T#F)
#    'useMinBoundingBoxCut': 'false',                                      # Use bounding box cut (T#F)
#    'sourceMinBoundingBox': 2,                                            # Minimum bounding box cut (source tagged as bad if below this threshold)
#    #'useCircRatioCut': 'true',                                            # Use circularity ratio cut (T#F)
#    'useCircRatioCut': 'false',                                           # Use circularity ratio cut (T#F)
#    'psCircRatioThr': 0.4,                                                # Circular ratio threshold (source passes point-like cut if above this threshold)
#    #'useElongCut': 'true',                                                # Use elongation cut (T#F)
#    'useElongCut': 'false',                                               # Use elongation cut (T#F)
#    'psElongThr': 0.7,                                                    # Elongation threshold (source passes point-like cut if below this threshold
#    'useEllipseAreaRatioCut': 'false',	                                  # Use Ellipse area ratio cut (T#F)
#    'psEllipseAreaRatioMinThr': 0.6,                                      # Ellipse area ratio min threshold
#    'psEllipseAreaRatioMaxThr': 1.4,                                      # Ellipse area ratio max threshold
#    'useMaxNPixCut': 'false',                                             # Use max npixels cut (T#F)
#    'psMaxNPix': 1000,                                                    # Max number of pixels for point-like sources (source passes point-like cut if below this 
#                                                                          # threshold)
#    'useNBeamsCut': 'true',                                               # Use nBeams cut (T#F)
#    'psNBeamsThr': 10,                                                    # nBeams threshold (sources passes point-like cut if nBeams<thr)
#    'comment_060': "###",
#    'comment_061': "###",
#    'comment_062': "##=====================================================================================================================================================================",
#    'comment_063': "##==                                                                      SOURCE RESIDUAL OPTIONS                                                                    ==",
#    'comment_064': "##=====================================================================================================================================================================",
#    'computeResidualMap': 'false',                                        # Compute compact source residual map (after compact source search) (T#F)
#    #'removeNestedSources': 'false',                                       # Dilate sources nested inside bright sources (T#F)
#    'removeNestedSources': 'true',                                        # Dilate sources nested inside bright sources (T#F)
#    'residualZThr': 5,                                                    # Significance threshold (in sigmas) above which sources of selected type are dilated
#    'residualZHighThr': 10,                                               # Significance threshold (in sigmas) above which sources are always dilated (even if they have 
#                                                                          # nested or different type)
#    'dilateKernelSize': 9,                                                # Size of kernel (odd) to be used in dilation operation
#    'removedSourceType': 2,                                               # Type of bright sources to be dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,
#                                                                          # 3=EXTENDED)
#    'residualModel': 1,                                                   # Model used to replace residual pixel values (1=bkg,2=source median)
#    'residualModelRandomize': 'false',                                    # Randomize pixel values used to replace residual pixels (T#F)
#    'psSubtractionMethod': 1,                                             # Point-source subtraction method (1=dilation, 2=model subtraction (default=1)
#    'residualBkgAroundSource': 'false',                                   # Use bkg around source and not global#local bkg map (default=no)
#    'comment_065': "###",
#    'comment_066': "###",
#    'comment_067': "##=====================================================================================================================================================================",
#    'comment_068': "##==                                                                  EXTENDED SOURCE FINDING OPTIONS                                                                ==",
#    'comment_069': "##=====================================================================================================================================================================",
#    'searchExtendedSources': 'false',                                     # Search extended sources after bright source removal (T#F)
#    'extendedSearchMethod': 4,                                            # Extended source search method (1=WT-thresholding,2=SPSegmentation,3=ActiveContour,4=Saliency 
#                                                                          # thresholding)
#    'useResidualInExtendedSearch': 'true',                                # Use residual image (with selected sources dilated) as input for extended source search
#    'comment_070': "###",
#    'comment_071': "###",
#    'comment_072': "##=====================================================================================================================================================================",
#    'comment_073': "##==                                                            WAVELET TRANSFORM ALGORITHM MAIN OPTIONS                                                             ==",
#    'comment_074': "##=====================================================================================================================================================================",
#    'wtScaleSearchMin': 3,                                                # Minimum Wavelet scale to be used for extended source search
#    'wtScaleSearchMax': 6,                                                # Maximum Wavelet scale to be used for extended source search
#    'comment_075': "###",
#    'comment_076': "###",
#    'comment_077': "##=====================================================================================================================================================================",
#    'comment_078': "##==                                                            ACTIVE CONTOUR ALGORITHMS MAIN OPTIONS                                                               ==",
#    'comment_079': "##=====================================================================================================================================================================",
#    #'acMethod': 2,                                                        # Active contour method (1=Chanvese, 2=LRAC)
#    'acMethod': 1,                                                        # Active contour method (1=Chanvese, 2=LRAC)
#    'acNIters': 1000,                                                     # Maximum number of iterations
#    'acInitLevelSetMethod': 1,                                            # Level set initialization method (1=circle,2=checkerboard,3=saliency)
#    'acInitLevelSetSizePar': 0.1,                                         # Level set size fraction wrt to minimum image size (e.g. circle radius=fraction x image size)
#    'acTolerance': 0.1,                                                   # Tolerance parameter to stop main iteration loop
#    'comment_080': "###",
#    'comment_081': "###",
#    'comment_082': "##=====================================================================================================================================================================",
#    'comment_083': "##==                                                                    CHAN-VESE ALGORITHM OPTIONS                                                                  ==",
#    'comment_084': "##=====================================================================================================================================================================",
#    'cvNItersInner': 5,                                                   # Chan-Vese maximum number of inner iterations
#    'cvNItersReInit': 5,                                                  # Chan-Vese maximum number of iterations performed in re-initialization step
#    'cvTimeStepPar': 0.007,                                               # Chan-Vese time step par
#    'cvWindowSizePar': 1,                                                 # Chan-Vese window size par
#    'cvLambda1Par': 1,                                                    # Chan-Vese lambda1 par
#    'cvLambda2Par': 2,                                                    # Chan-Vese lambda2 par
#    'cvMuPar': 0.5,                                                       # Chan-Vese mu par
#    'cvNuPar': 0,                                                         # Chan-Vese nu par
#    'cvPPar': 1,                                                          # Chan-Vese p par
#    'comment_085': "###",
#    'comment_086': "###",
#    'comment_087': "##=====================================================================================================================================================================",
#    'comment_088': "##==                                                                      LRAC ALGORITHM OPTIONS                                                                     ==",
#    'comment_089': "##=====================================================================================================================================================================",
#    'lracLambdaPar': 0.1,                                                 # Regularization par
#    #'lracRadiusPar': 1,                                                   # Radius of locatization ball par
#    'lracRadiusPar': 10,                                                  # Radius of locatization ball par
#    'lracEpsPar': 0.1,                                                    # Convergence par
#    'comment_090': "###",
#    'comment_091': "###",
#    'comment_092': "##=====================================================================================================================================================================",
#    'comment_093': "##==                                                                     SALIENCY FILTER OPTIONS                                                                     ==",
#    'comment_094': "##=====================================================================================================================================================================",
#    'saliencyUseOptimalThr': 'true',                                      # Use optimal threshold in multiscale saliency thresholding (T#F)
#    #'saliencyThrFactor': 2.8,                                             # Saliency threshold factor for tagging signal regions (thr=<saliency>*factor)
#    'saliencyThrFactor': 1,                                               # Saliency threshold factor for tagging signal regions (thr=<saliency>*factor)
#    'saliencyBkgThrFactor': 1,                                            # Saliency threshold factor for tagging bkg regions (thr=<saliency>*factor)  # NB: Note in online manaul
#    'saliencyImgThrFactor': 1,                                            # Threshold factor to consider a region as significant (thr=<img>*factor)
#    'saliencyResoMin': 20,                                                # Saliency min reso par
#    'saliencyResoMax': 60,                                                # Saliency max reso par
#    'saliencyResoStep': 10,                                               # Saliency reso step par
#    'saliencyUseCurvInDiss': 'false',                                     # Use curvature parameter in dissimilarity estimation (T#F)                  # NB: Note in online manaul
#    'saliencyUseRobustPars': 'false',                                     # Use robust pars in saliency map computation (T#F)
#    'saliencyUseBkgMap': 'false',                                         # Use bkg map in saliency map computation (T#F)
#    'saliencyUseNoiseMap': 'false',                                       # Use noise map in saliency map computation (T#F)
#    'saliencyNNFactor': 1,                                                # Fraction of most similar neighbors used in saliency map computation
#    ###'saliencySpatialRegFactor': 6,                                     # Spatial regularization factor (ruling exp decay in saliency spatial weighting)"
#    'saliencyMultiResoCombThrFactor': 0.7,                                # Fraction of combined salient multi-resolution maps to consider global saliency
#    'saliencyDissExpFalloffPar': 100,                                     # Dissimilarity exponential cutoff parameter (value)
#    'saliencySpatialDistRegPar': 1,                                       # Spatial-color distance regularization par (value, 1=equal weights)
#    'comment_095': "###",
#    'comment_096': "###",
#    'comment_097': "##=====================================================================================================================================================================",
#    'comment_098': "##==                                                                       SP SEGMENTATION OPTIONS                                                                   ==",
#    'comment_099': "##=====================================================================================================================================================================",
#    'spSize': 20,                                                         # Initial superpixel size                                                    # NB: Note in online manaul
#    'spBeta': 1,                                                          # Initial superpixel regularization parameter
#    'spMinArea': 10,                                                      # Initial superpixel min area
#    'spUseLogContrast': 'false',                                          # Use intensity log scale to generate initial superpixel partition (T#F)     # NB: Note in online manaul
#    'comment_100': "###",
#    'comment_101': "###",
#    'comment_102': "##=====================================================================================================================================================================",
#    'comment_103': "##==                                                                       SP MERGING OPTIONS                                                                        ==",
#    'comment_104': "##=====================================================================================================================================================================",
#    'spMergingNSegmentsToStop': 1,                                        # Number of segments below which the segmentation is stopped (default=1)                        # NB: Note in online manaul
#    'spMergingRatio': 0.3,                                                # Fraction of similar segments merged per each hierarchical level (default=0.3)                 # NB: Note in online manaul
#    'spMergingRegPar': 0.5,                                               # Regularization parameters balancing region edge and similarity (default=0.5)                  # NB: Note in online manaul
#    'spMergingMaxDissRatio': 1000,                                        # Mutual segment dissimilarity ratio (R=Diss(ji)#Diss(ij)) above which 1st-order neighbors are
#                                                                          # not merged (default=1000)                                                                     # NB: Note in online manaul
#    'spMergingMaxDissRatio2ndNeighbours': 1.05,                           # Mutual segment dissimilarity ratio (R=Diss(ji)#Diss(ij)) above which 2nd-order neighbors are 
#                                                                          # not merged (default=1.05)                                                                     # NB: Note in online manaul
#    'spMergingDissThreshold': 3,                                          # Absolute dissimilarity ratio (wrt to initial average diss) threshold above which segments are 
#                                                                          # not merged (default=3)                                                                        # NB: Note in online manaul
#    'spMergingEdgeModel': 2,                                              # Superpixel edge model (1=Kirsch,2=Chan-Vese)                                                  # NB: Note in online manaul 
#    'spMergingIncludeSpatialPars': 'false',                               # Include spatial pars in region dissimilarity measure (default=false)                          # NB: Note in online manaul
#    ###
#    ###
#    ###spHierMergingPars': 1  0.3  0.5  3, # Hierarch algo pars: MIN_SEGMENTS, MERGING_RATIO, DIST_REGULARIZATION, DIST_THRESHOLD es 0.25 (value,value,value,value)
#    ###spHierMergingMaxDissRatio': 1000  1.05, # Maximum mutual dissimilarity among regions for merging 1st and 2nd neighbors es. 1.15 (value)
#    ###spMergingEdgeModel': 2,             # Edge model (1=Kirsch,2=Chan-Vese) (value)
#    ###use2ndNeighborsInSPMerging  T,      # Use 2nd-order neighbors in superpixel merging (T#F)
#    ###useCurvatureInSPMerging  T,         # Use curvature params in superpixel merging (T#F)
#    ###useLogContrastInSPGeneration  F,    # Use logarithmic contrast to generate initial partition (T#F)
#    ###usePixelRatioCut  T  0.4,           # Use pixel ratio to tag background regions (if npix#npixtot>cut)  (T#F,value)
#    ###tagSignificantSP  T  1  0.5,        # Tag significant superpixels if fraction of significant subregions is above cut (1=saliency,2=Mahalanobis) (T#F,method,cut)
#    ###
#    ###
#}

def get_config(image,seedthr,mergethr,rms_box={'box_size': None, 'step_size': None}):
    # notes: 
    # [r] https://root.cern.ch/root/htmldoc/guides/users-guide/GettingStarted.html 
    #     root[] .q
    # [2] https://caesar-doc.readthedocs.io/en/latest/tutorials/app_findsource.html
    # [3] https://caesar-doc.readthedocs.io/en/latest/usage/app_options.html#input-options
    #config = [
    #        ['block', 'MAIN'],
    #        ['inputFile',image],
    #        ['block','BKG OPTIONS'],
    #        ['useLocalBkg','true'],
    #        ['bkgEstimator', 2],
    #        ['useBeamInfoInBkg', 'true'],
    #        ['boxSizeX', 20],
    #        ['boxSizeY', 20],
    #        ['gridSizeX', 0.2],
    #        ['gridSizeY', 0.2],
    #        ['block', 'SOURCE FINDING OPTIONS'],
    #        ['searchCompactSources', 'true'],
    #        ['minNPix', 5],
    #        ['seedThr', seedthr],
    #        ['mergeThr', mergethr],
    #        ['compactSourceSearchNIters', 2]
    #        ['seedThrStep', 0.5]
    #]
    contents = list() 
    contents.append("//=======================================")
    contents.append("//==               MAIN                ==")
    contents.append("//=======================================")
    contents.append("inputFile = %s" % strip(image))
    contents.append("###")
    contents.append("//=======================================")
    contents.append("//==            BKG OPTIONS            ==")
    contents.append("//=======================================")
    contents.append("useLocalBkg = true")
    contents.append("bkgEstimator = 2")
    contents.append("useBeamInfoInBkg = true")
    if not rms_box['box_size'] is None and not rms_box['step_size'] is None:
        header = fits.open(image)[0].header
        CDELT = (abs(float(header['CDELT1']))+abs(float(header['CDELT2'])))/2.0
        beam_size = float(header['BMAJ']+header['BMIN'])/(2.0*CDELT)
        box_size = float(rms_box['box_size'])/float(beam_size)
        step_size = float(rms_box['step_size'])/float(rms_box['box_size'])
        contents.append("boxSizeX  = %f" % box_size)
        contents.append("boxSizeY  = %f" % box_size)
        contents.append("gridSizeX  = %f" % step_size)
        contents.append("gridSizeY  = %f" % step_size)
    #contents.append("boxSizeX = 20")
    #contents.append("boxSizeY = 20")
    #contents.append("gridSizeX = 0.2")
    #contents.append("gridSizeY = 0.2")
    contents.append("###")
    contents.append("//=======================================")
    contents.append("//==     SOURCE FINDING OPTIONS        ==")
    contents.append("//=======================================")
    contents.append("searchCompactSources = true")
    #contents.append("minNPix = 5")
    contents.append("seedThr = %f" % seedthr)
    contents.append("mergeThr = %f" % mergethr)
    #contents.append("compactSourceSearchNIters = 2")
    contents.append("compactSourceSearchNIters = 1")
    #contents.append("seedThrStep = 0.5")
    contents.append("###")
    contents.append("//=======================================")
    contents.append("//==   NESTED SOURCE FINDING OPTIONS   ==")
    contents.append("//=======================================")
    contents.append("searchNestedSources = true")
    contents.append("extendedSearchMethod = 4")
    contents.append("blobMaskMethod = 2")
    contents.append("nestedBlobPeakZThr = %f" % seedthr)
    contents.append("nestedBlobPeakZMergeThr = %f" % mergethr)
    contents.append("###")
    contents.append("//=======================================")
    contents.append("//==      SOURCE FITTING OPTIONS       ==")
    contents.append("//=======================================")
    contents.append("fitSources = true")
    contents.append("###")
    contents.append("//=======================================")
    contents.append("//==         RESIDUAL OPTIONS          ==")
    contents.append("//=======================================")
    contents.append("computeResidualMap = true")
    contents.append("removeNestedSources = true")
    contents.append("residualZHighThr = %f" % seedthr)
    contents.append("residualZThr = %f" % mergethr)
    #contents.append("residualZThr = %f" % seedthr)
    contents.append("removedSourceType = -1")
    contents.append("saveResidualMap = true")
    contents.append("residualMapFITSFile = %s" % global_residaul_filename)
    #contents.append("###")
    contents.append("\n") # nb: must have return character after each directive
    return "\n".join(contents)


def get_root_residual_image_extraction_script():
    root_script = list()
    root_script.append('#include <Image.h>')
    root_script.append('')
    root_script.append('void caesar()')
    root_script.append('{')
    root_script.append('    TFile *MyFile = new TFile("output.root","READ");')
    root_script.append('    gFile = MyFile;')
    root_script.append('    gDirectory->ls();')
    root_script.append('    Caesar:Image *img;')
    root_script.append('    MyFile->GetObject("img_residual;1",img);')
    root_script.append('    img->WriteFITS("%s");' % global_residaul_filename)
    root_script.append('}')
    return "\n".join(root_script)+"\n"


def clean_residual_header(image_residual,image_reference):
    hdul_residual  = fits.open(image_residual)
    hdul_reference = fits.open(image_reference)
    data = np.squeeze(hdul_residual[0].data)
    data.shape = hdul_reference[0].data.shape
    header_reference = hdul_reference[0].header.copy()
    fits.PrimaryHDU(data,header=header_reference).writeto(image_residual,overwrite=True)


def strip(fname):
    return re.sub(r"^(.*?/)+","",fname)


def dat_to_csv(dat_file):
    csv_file = re.sub("\.dat$",".csv",dat_file)
    print("> Converting '%s' to '%s'..." % (dat_file,csv_file))
    with open(dat_file,"r") as fd:
        comments = list()
        rows = list()
        for line in fd:
            if re.search(r"^#",line):
                comments.append(line)
            else:
                rows.append(re.sub(r",$","",re.sub("\s+",",",line)))
        header = re.sub(r"\s+",",",re.sub(r"\s+$","",re.sub(r"^#\s*","",comments[-1])))
        header = re.sub(r"\(.*?\)","",header)
        with open(csv_file,"w") as fc:
            fc.write(header+"\n")
            fc.write("\n".join(rows))
    return csv_file

#def create_residual_image(image_filename,catalogue_filename,input_dir,output_dir):
#    def get_wcs(header):
#        header = hdul[0].header
#
#        # minor header fix
#        rotation_matrix_map = {
#            'PC01_01': 'PC1_1',
#            'PC02_01': 'PC2_1',
#            'PC03_01': 'PC3_1',
#            'PC04_01': 'PC4_1',
#            'PC01_02': 'PC1_2',
#            'PC02_02': 'PC2_2',
#            'PC03_02': 'PC3_2',
#            'PC04_02': 'PC4_2',
#            'PC01_03': 'PC1_3',
#            'PC02_03': 'PC2_3',
#            'PC03_03': 'PC3_3',
#            'PC04_03': 'PC4_3',
#            'PC01_04': 'PC1_4',
#            'PC02_04': 'PC2_4',
#            'PC03_04': 'PC3_4',
#            'PC04_04': 'PC4_4',
#        }
#        for key in rotation_matrix_map.keys():
#            if key in header:
#                header.insert(key,(rotation_matrix_map[key],header[key]),after=True)
#                header.remove(key)
#
#        # wcs image header
#        wcs = WCS(header)
#
#        # trim to 2d from nd
#        naxis = wcs.naxis
#        while naxis > 2:
#            wcs = wcs.dropaxis(2)
#            naxis -= 1
#
#        return wcs
#
#    # file definitions
#    image_file     = input_dir+"/"+image_filename
#    catalogue_file = input_dir+"/"+catalogue_filename
#    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".caesar.residual.model.fits",output_dir+"/"+image_filename)
#    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".caesar.residual.fits",output_dir+"/"+image_filename)
#    print("> Creating Residual Image:")
#    print(">>    INPUT_IMAGE: "+image_file)
#    print(">>    INPUT_CATALOGUE: "+catalogue_file)
#
#    # open image and catalogue files
#    hdul = fits.open(image_file)
#    header = hdul[0].header
#    wcs = get_wcs(header)
#    shape = hdul[0].data.shape
#    image = np.squeeze(hdul[0].data)
#    CDELT1 = np.abs(wcs.wcs.cdelt[0])
#    CDELT2 = np.abs(wcs.wcs.cdelt[1])
#    qt = Table.read(catalogue_file)
#    #print(qt) # debug
#
#    # build our model image
#    scale = 8
#    n_x = image.shape[1]-1
#    n_y = image.shape[0]-1
#    model = np.zeros(image.shape)
#    def make_gaussian(ra_0,dec_0,smaj,smin,pa,peak):
#        a = ((np.cos(pa)/smin)**2+(np.sin(pa)/smaj)**2)/2.0
#        b = (np.sin(2.0*pa)/(smin)**2-np.sin(2.0*pa)/(smaj)**2)/2.0
#        c = ((np.sin(pa)/smin)**2+(np.cos(pa)/smaj)**2)/2.0
#        def gaussian(ra,dec):
#            return peak*np.exp(-a*(ra-ra_0)**2-b*(ra-ra_0)*(dec-dec_0)-c*(dec-dec_0)**2)
#        return gaussian
#    for row in qt:
#        # compute (ra,dec) in pixel units
#        ra  = (row['x_wcs']*catalogue_header['x_wcs']).to(u.deg).value
#        dec = (row['y_wcs']*catalogue_header['y_wcs']).to(u.deg).value
#        r_pix = wcs.all_world2pix([[ra,dec]],0,ra_dec_order=True)[0]
#        ra_0  = r_pix[0]
#        dec_0 = r_pix[1]
#
#        # compute semi-major/minor axes in pixel units
#        pa = (row['pa_wcs']*catalogue_header['pa_wcs']).to(u.rad).value
#        a  = (row['bmaj_wcs']*catalogue_header['bmin_wcs']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
#        b  = (row['bmin_wcs']*catalogue_header['bmaj_wcs']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
#        smaj = a*np.sqrt((np.sin(pa)/CDELT1)**2+(np.cos(pa)/CDELT2)**2)
#        smin = b*np.sqrt((np.cos(pa)/CDELT1)**2+(np.sin(pa)/CDELT2)**2)
#
#        # get peak flux
#        pf = row['Speak'] # NB: assumes units of Jy/beam.
#        
#        # create the gaussian function
#        gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,pf)
#
#        # make residuals
#        ra_min  = max(0,int(np.floor(ra_0-scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
#        ra_max  = min(n_x,int(np.ceil(ra_0+scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
#        dec_min = max(0,int(np.floor(dec_0-scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
#        dec_max = min(n_y,int(np.ceil(dec_0+scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
#        i,j = np.mgrid[ra_min:ra_max+1,dec_min:dec_max+1]
#        model[j,i] += gaussian(i,j)
#
#    # create residual image
#    image -= model
#
#    # output model image
#    print(">>    OUTPUT_MODEL: "+model_file)
#    model.shape = shape
#    fits.PrimaryHDU(model,header=header).writeto(model_file,overwrite=True)
#
#    # output residual image
#    print(">>    OUTPUT_IMAGE: "+residual_file)
#    image.shape = shape
#    fits.PrimaryHDU(image,header=header).writeto(residual_file,overwrite=True)
#    print("> [Done]")


def cleanup_reg_file(reg_file):
    if os.path.isfile(reg_file):
        with open(reg_file,"r+") as f:
            file_contents = list()
            for line in f:
                if re.match(r"(ellipse|polygon)",line):
                    # remove the region labels
                    line = re.sub(r"\s*#.*","",line)
                file_contents.append(line)
            f.seek(0)
            f.truncate(0)
            for line in file_contents:
                f.write(line)

@click.command()
@click.argument('input_dir',nargs=1)
@click.argument('processing_dir',nargs=1)
@click.argument('output_dir',nargs=1)
@click.argument('fits_image_file',nargs=1)
@click.option(
    '--seedThr',
    default = caesar_cfg['seedThr'],
    type = float,
    show_default = True,
    help = "Blob finding threshold (in sigmas)."
)
@click.option(
    '--mergeThr',
    default = caesar_cfg['mergeThr'],
    type = float,
    show_default = True,
    help = "Blob growth threshold (in sigmas)."
)
@click.option('--box-size', type=int,default=None,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--step-size',type=int,default=None,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and model FITS files.")
@click.option('--dump',is_flag=True,default=False, help="Dump out all processing files.")
def process(
    input_dir,
    processing_dir,
    output_dir,
    fits_image_file,
    seedthr,
    mergethr,
    box_size,
    step_size,
    fits,
    residual,
    dump
):
    """\b
       Caesar image processing tool.

       inputs:

          \b
          INPUT_DIR: location of image_filename.fits
          PROCESSING_DIR: location of scratch directory
          OUTPUT_DIR: location to place results
          FITS_IMAGE_FILE: image_filename.fits (without path)

       outputs:

          \b
          OUTPUT_DIR/image_filename.caesar.fits
          OUTPUT_DIR/image_filename.caesar.reg
    """
    pwd = os.getcwd()

    # house cleaning
    input_dir      = re.sub(r"/+$","",input_dir)
    processing_dir = re.sub(r"/+$","",processing_dir)
    output_dir     = re.sub(r"/+$","",output_dir)

    src = input_dir+"/"+fits_image_file
    print("Processing: %s" % src)
    print("> Local Context: ")
    print(">   INPUT_DIR:  %s" % input_dir)
    print(">   PROCESSING_DIR: %s" % processing_dir)
    print(">   OUTPUT_DIR: %s" % output_dir)
    print("> Supported Flags: ")
    print(">    --seedThr: %s" % seedthr)
    print(">    --mergeThr: %s" % mergethr)

    # check if image file exists
    if not os.path.isfile(src):
        print("ERROR: Image file '%s' not found!" % src)
        print("Bye!")
        exit()

    # link image file from data to processing dir
    dst = processing_dir+"/"+fits_image_file
    print("> Linking image file:")
    print("> $ ln -s %s %s" % (dst,src))
    # this gaurd condition is usefull for debugging.
    if not (os.path.islink(dst) or os.path.exists(dst)):
        os.symlink(src,dst)

    # setup caesar configurtion
    cfg = get_config(processing_dir+"/"+fits_image_file,seedthr,mergethr,{'box_size': box_size, 'step_size': step_size})
    cfg_file = processing_dir+"/caesar.cfg"
    print("Saving Configuration: %s" % cfg_file)
    print("> "+"\n> ".join(re.sub(r"\n$","",cfg).split("\n")))
    with open(cfg_file,"w") as f:
        f.write(cfg)

    # run caesar
    print("Executing selavy:")
    print("> $ FindSource --config=%s" % cfg_file)
    current_dir = os.getcwd()
    os.chdir(processing_dir)
    os.system("/opt/Software/Sources/caesar-build/FindSource --config=%s" % strip(cfg_file))
    root_script_file = "caesar.C"
    with open(root_script_file,'w') as fd:
        fd.write(get_root_residual_image_extraction_script())
    print("> Extracting residual image from root.")
    print(">> $ root -q %s" % root_script_file)
    os.system("root -q %s" % root_script_file)
    clean_residual_header(
        global_residaul_filename,
        strip(fits_image_file)
    )
    os.chdir(current_dir)

    # convert dat to csv files
    print("Creating .csv files...")
    for fname in glob.glob(processing_dir+"/*.dat"):
        dat_to_csv(fname)

    print("Caesar files: "+processing_dir+"/")
    for f in glob.glob(processing_dir+"/*"):
        f = re.sub(r'^(.*?/)*','',f)
        if f != fits_image_file:
            print("> o "+f)

    # get component catalogue
    print("Creating Catalogue...")
    csv_file = processing_dir+"/catalog_fitcomp.csv"
    if not os.path.isfile(csv_file):
        with open(csv_file,"w") as fd:
            fd.write(",".join(catalogue_header)+"\n")
    if fits:
        cat_file = output_dir+"/"+re.sub(r"\.([Ff][Ii][Tt]([Ss]|))$",".caesar.fits",fits_image_file)
        print("> Mapping '%s' to '%s'..." % (csv_file,cat_file))
        qt = QTable.read(csv_file,format='csv')
        for header in catalogue_header:
            if not catalogue_header[header] is None:
                qt[header] = qt[header]*catalogue_header[header]
        qt.write(cat_file,format='fits',overwrite=True)
    else:
        cat_file = output_dir+"/"+re.sub(r"\.([Ff][Ii][Tt]([Ss]|))$",".caesar.csv",fits_image_file)
        print("> $ cp %s %s" % (csv_file,cat_file))
        copy(csv_file,cat_file)

    # get region file
    print("Getting Region File...")
    reg_in_file = processing_dir+"/ds9_fitcomp.reg"
    reg_out_file = output_dir+"/"+re.sub(r"\.([Ff][Ii][Tt]([Ss]|))$",".caesar.reg",fits_image_file)
    print("> $ cp %s %s" %(reg_in_file,reg_out_file))
    copy(reg_in_file,reg_out_file)
    cleanup_reg_file(reg_out_file)

    # create residual file...
    if residual:
        print("Creating residual and model files.")
        #create_residual_image(fits_image_file,strip(csv_file),processing_dir,output_dir)
        src_residual = processing_dir+"/"+global_residaul_filename
        dst_residual = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".caesar.residual.fits",output_dir+"/"+strip(fits_image_file))
        def make_model():
            # NB: This bad form! That is, we are using fits
            # as a boolean above, so we must re-import...
            # TO-DO: Should fix -- for all modules.
            from astropy.io import fits
            dst_model = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".caesar.residual.model.fits",output_dir+"/"+strip(fits_image_file))
            print("> Creating model file: %s" % dst_model)
            src_image_file = input_dir+"/"+fits_image_file
            hdul_image    = fits.open(src_image_file)
            hdul_residual = fits.open(src_residual)
            image    = np.squeeze(hdul_image[0].data)
            residual = np.squeeze(hdul_residual[0].data)
            model = image - residual
            print(">> Saving...")
            model.shape = hdul_image[0].data.shape
            fits.PrimaryHDU(model,header=hdul_image[0].header).writeto(dst_model,overwrite=True)
        #print("Overwriting: %s" % dst_residual)
        print("> Getting residual file: %s" % dst_residual)
        print(">> $ cp %s %s" % (src_residual,dst_residual))
        os.system("cp %s %s" % (src_residual,dst_residual))
        make_model()

    if dump:
        print("Dumping all files...")
        for in_file in glob.glob(processing_dir+"/*"):
            if strip(in_file) != fits_image_file:
                out_file = output_dir+"/"+strip(in_file)
                print("> $ cp %s %s" % (in_file,out_file))
                copy(in_file,out_file)

    print("[Done]")


if __name__ == "__main__":
    process()




