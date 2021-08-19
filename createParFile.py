import numpy as np
import pylab as pl
import uuid
import pdb
import os

def createParFile(objfile, parfile=None, parfile_dir='/Volumes/Samsung_T5/ESPRESSO/ESPRESSO5/S1D/molecfit/',\
                  obj_dir='/Volumes/Samsung_T5/ESPRESSO/ESPRESSO5/S1D_test/molecfit/',\
                  gdas_prof='', gdas_dir='~/molecfit/share/molecfit/data/profiles/gdas', \
                  outfile='default_result', \
                  listname=None, \
                  outdir='output', trans=True, \
                  wrange_include='/Volumes/Samsung_T5/ESPRESSO/ESPRESSO5/S1D/molecfit/ESPRESSO_include.dat',\
                  wrange_exclude='/Volumes/Samsung_T5/ESPRESSO/ESPRESSO5/S1D/molecfit/ESPRESSO_exclude.dat',\
                  prange_exclude=None, \
                  columns=["WAVELENGTH", "FLUX", "FLUX_ERR", "NULL"], \
                  default_error=0.01, wlgtomicron=0.0001, \
                  vac_air='vac', \
                  plot_type='P', plot_range=True, \
                  ftol=1e-10, xtol=1e-10, \
                  list_molec=['H2O', 'O2'], fit_molec=[True, True,], relcol=[1.0,1.0], \
                  flux_unit=0, fit_back=False, telback=0.1, \
                  fit_cont=True, cont_n=4, cont_const=2500.0, \
                  fit_wlc=True, wlc_n=1, wlc_const=-0.0, \
                  fitresbox=False, kernmode=True, \
                  fit_res_box=False, relres_box=0.0, \
                  fit_res_gauss=True, res_gauss=3.0, \
                  fit_res_lorentz=True, res_lorentz=1.5, \
                  kernfac=30.0, varkern=True, kernel_file=None, \
                  obsdate=None, obsdate_key='MJD-OBS', \
                  utc=None, utc_key='UTC', \
                  telalt=None, telalt_key='ESO TEL1 ALT', \
                  rhum=None, rhum_key='ESO TEL1 AMBI RHUM', \
                  pres=None, pres_key='ESO TEL1 AMBI PRES START', \
                  temp=None, temp_key='ESO TEL1 AMBI TEMP', \
                  m1temp=None, m1temp_key='ESO TEL1 TH M1 TEMP', \
                  geoelev=None, geoelev_key='ESO TEL1 GEOELEV', \
                  longitude=None, longitude_key='ESO TEL1 GEOLON', \
                  latitude=None, latitude_key='ESO TEL1 GEOLAT', \
                  slitw=1.0, slitw_key=None, \
                  pixsc=0.19, pixsc_key=None, \
                  ref_atm="equ.atm", layers=1, emix=5.0, \
                  pwv=-1.0, \
                  clean_mflux=True):


    """
    Script to create a Molecfit parameter file, given various input files
    #vac_air:   vac/air
    #dateKey:   MJD-OBS / DATE_VAL
    #timeKey:   TM-START / TIME_VAL
    #telAltKey: ESO TEL ALT / TELALT_VAL
    #VARFWHM:    0 / 1 
    #FLUXLIM:    -1 or cut value
    
    """
    # create output parfilename
    if parfile==None:
#AS        parfile = user_workdir+outdir+objfile[:objfile.rfind(".fits")]+".molecfit.par"
        parfile = parfile_dir+objfile[:objfile.rfind(".fits")]+".par"
#AS    # create outdir if not present
#AS    if (not os.path.isdir(user_workdir+outdir)):
#AS        os.mkdir(user_workdir+outdir)
    # create pardir if not present
    if (not os.path.isdir(parfile_dir+objfile[:objfile.rfind(".fits")])):
        os.mkdir(parfile_dir+objfile[:objfile.rfind(".fits")])


    # open file handler
    fp = open(parfile, "w")
    # start writing info
    fp.write("\n")
    # files
    fp.write("\nfilename: "+obj_dir+objfile) # object spectrum

    if listname==None:
        fp.write("\nlistname: "+str(listname).lower())# subsequent files to be corrected in list mode
    else:
        fp.write("\nlistname: "+listname)# subsequent files to be corrected in list mode
    fp.write("\ntrans: "+str(int(trans)))# type of spectrum: transmission = True/1, emission=False/0

    # generate colnames
    ncolumns = len(columns)
    if ncolumns>4: raise ValueError('Columns must have no more than 4 elements')
    if ncolumns<4:
        nnone=4-ncolumns
        extra = ["NULL" for x in xrange(nnone)]
        columns.extend(extra)
    if len(columns)!=4: raise ValueError("This should never happen")
    cns=''
    for cn in columns:
        cns=cns+str(cn)+" "
    fp.write("\ncolumns: "+cns)

    fp.write("\ndefault_error: "+str(default_error))
    fp.write("\nwlgtomicron: "+str(wlgtomicron))
    fp.write("\nvac_air: "+vac_air)

    # Exclude and include range files
    if wrange_include==None:
        fp.write("\nwrange_include: "+str(wrange_include).lower()) # none
    else:
        fp.write("\nwrange_include: "+wrange_include) # range to INclude in fit (ASCII or FITS table)
    if wrange_exclude==None:
        fp.write("\nwrange_exclude: "+str(wrange_exclude).lower()) # none
    else:
        fp.write("\nwrange_exclude: "+wrange_exclude) # range to EXclude in fit (")
    if prange_exclude==None:
        fp.write("\nprange_exclude: "+str(prange_exclude).lower()) # none
    else:
        fp.write("\nprange_exclude: "+prange_exclude) # range to EXclude in fit (")


    # results
    fp.write("\noutput_dir: "+outdir+'/')         # output dir
    fp.write("\noutput_name: "+outfile)           # corrected object spectrum
    # plotting
    if plot_type == None: fp.write("\nplot_creation: "+str(plot_type).lower())
    else: fp.write("\nplot_creation: "+str(plot_type))
        
    fp.write("\nplot_range: "+str(int(plot_range)))

    # fit precision
    # fitting
    fp.write("\nftol: "+str(ftol))
    fp.write("\nxtol: "+str(xtol))

    # molec columns: define species present, and which to fit
    fp.write("\nlist_molec: "+" ".join(list_molec))
    fp.write("\nfit_molec: "+" ".join(np.array(np.array(fit_molec,dtype=np.int),dtype=str).tolist()))
    fp.write("\nrelcol: "+" ".join(np.array(np.array(relcol,dtype=np.float),dtype=str).tolist()))

    # background and continuum
    fp.write("\nflux_unit: "+str(int(flux_unit))) # what flux units to use. see manual, 0 most likely
    fp.write("\nfit_back: "+str(int(fit_back))) # fit telescope background?
    fp.write("\ntelback: "+str(telback)) # initial guess for telescope background fit
    fp.write("\nfit_cont: "+str(int(fit_cont))) # fit continuum?
    fp.write("\ncont_n: "+str(cont_n)) # polynomial degree for fit
    fp.write("\ncont_const: "+str(cont_const)) # initial term for continuum fit

    # wavelength solution
    fp.write("\nfit_wlc: "+str(int(fit_wlc))) # fit polynomial for new wavelength solution?
    fp.write("\nwlc_n: "+str(wlc_n)) # order of wavelength fit
    fp.write("\nwlc_const: "+str(wlc_const)) # initial guess for const term

    # resolution fitting
    fp.write("\nfit_res_box: "+str(int(fit_res_box))) # fit resolution by boxcar?
    fp.write("\nrelres_box: "+str(relres_box)) # initial guess of boxcar relative to slitwidth (0-2)
    fp.write("\nkernmode: "+str(int(kernmode))) # Voigt profile approx instead of Gauss and Lorent
    fp.write("\nfit_res_gauss: "+str(int(fit_res_gauss))) # fit resolution with Gaussian?
    fp.write("\nres_gauss: "+str(res_gauss)) # initial FWHM of Gauss in pixels
    fp.write("\nfit_res_lorentz: "+str(int(fit_res_lorentz))) # fit resolution with Lorentian?
    fp.write("\nres_lorentz: "+str(res_lorentz))
    fp.write("\nkernfac: "+str(kernfac)) # size of kernal in FWHM (?!)
    fp.write("\nvarkern: "+str(int(varkern))) # vary kernel size - linear with wavelength
    fp.write("\nkernel_file: "+str(kernel_file).lower()) # instead of parametric kernel, give ascii profile here

    
    # Ambient params

    # make helper function - either use Keyword val or passed value.
    def oneOrOther(thing,key,tag, upper=False):
        if ((thing==None) & (key!=None)):
            # use keyword to get value
            fp.write("\n"+tag)
            fp.write("\n"+tag+"_key: "+key)
        elif ((thing==None) & (key==None)):
           fp.write("\n"+tag)
           fp.write("\n"+tag+"_key")
        elif ((key==None) & (thing!=None)):
            # force use of value, not keyword
            fp.write("\n"+tag+": "+str(thing))
            if (not upper):
                fp.write("\n"+tag+"_key: "+str(key).lower())
            else:
                fp.write("\n"+tag+"_key: "+str(key).upper())
        elif ((key!=None) & (thing!=None)):
            # force use of value, not keyword
            fp.write("\n"+tag+": "+str(thing))
            if (not upper):
                fp.write("\n"+tag+"_key: "+str(key).lower())
            else:
                fp.write("\n"+tag+"_key: "+str(key).upper())
        else:
            raise IOError("Inputs not understood")

    oneOrOther(obsdate,obsdate_key,"obsdate")
    oneOrOther(utc, utc_key, "utc")
    oneOrOther(telalt, telalt_key, "telalt")
    oneOrOther(rhum, rhum_key, "rhum")
    oneOrOther(pres, pres_key, "pres")
    oneOrOther(temp, temp_key, "temp")
    oneOrOther(m1temp, m1temp_key, "m1temp")
    oneOrOther(geoelev, geoelev_key, "geoelev")
    oneOrOther(longitude, longitude_key, "longitude")
    oneOrOther(latitude, latitude_key, "latitude")


    # Instrument params
    oneOrOther(slitw, slitw_key, "slitw", upper=True)
    oneOrOther(pixsc, pixsc_key, "pixsc", upper=True)

    # Atmospheric profiles
    fp.write("\nref_atm: "+ref_atm)
    fp.write("\ngdas_dir: "+gdas_dir)
    fp.write("\ngdas_prof: "+gdas_prof)
    fp.write("\nlayers: "+str(int(layers)))
    fp.write("\nemix: "+str(emix))
    fp.write("\nclean_mflux: "+str(int(trans)))
    fp.write("\npwv: "+str(pwv))
    
    fp.write("\nend\n")
    fp.close()

    return parfile
