import numpy as np
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import pdb

###########################################################################
##
## Makes the input file based on the DSIM standards via the KECK webpage.
## This routine first requires the output from 'determine_priority.py' and
## all subsequent codes in the pipeline. See 'README_mynotes' for details.
##
############################################################################

def input1(maskname,semester,maskPA=0,masknum=1):

    if semester == '2017B':
        hdulist = fits.open(semester+'/COSMOS2015_specz_sdss_priority.fits')
    elif semester == '2018A':
        hdulist = fits.open(semester+'/newfirm_mbs_deepz_sdss_priority.fits')
    else:
        hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority.fits')
    targetdata = hdulist[1].data
    targetINFO = targetdata['target']
    stardata = hdulist[2].data

    if masknum == 30:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
        targetINFO = data_updateinfo['target']

        print((len(targetINFO[np.where(targetINFO==1)])))

    elif masknum == 40:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print((len(targetINFO[np.where(targetINFO==1)])))

    elif masknum == 13:
        hdu_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')
        data_updateinfo = hdu_updateinfo[1].data
        targetINFO = data_updateinfo['target']
        
        print((len(targetINFO[np.where(targetINFO==1)[0]])))

    elif masknum == 14:
        hdu_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')
        data_updateinfo = hdu_updateinfo[1].data
        targetINFO = data_updateinfo['target']
        
        print((len(targetINFO[np.where(targetINFO==1)[0]])))

    elif masknum == 22:
        hdu_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass2.fits')
        data_updateinfo = hdu_updateinfo[1].data
        targetINFO = data_updateinfo['target']
        
        print((len(targetINFO[np.where(targetINFO==1)[0]])))

    elif masknum == 32:
        hdu_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')
        data_updateinfo = hdu_updateinfo[1].data
        targetINFO = data_updateinfo['target']
        
        print((len(targetINFO[np.where(targetINFO==1)[0]])))
        
    else:
        print('all objects are being targeted first time')

    #############################
    ## Read in all necessary data
    #############################
    targetID = targetdata['id']
    targetRA = targetdata['ra']
    targetDEC = targetdata['dec']
    targetFLUX = targetdata['R']
    targetMAG = np.empty(len(targetFLUX))
    for i in range(len(targetFLUX)):
        if targetFLUX[i] <=0:
            targetMAG[i] = -99
        else:
            targetMAG[i] = -2.5*np.log10(targetFLUX[i]) + 25.0
            
    targetPcode = targetdata['Pflag']
    #targetELLIP = targetdata['K_ellip']
    #targetPA = targetdata['K_theta_J2000']
    
    starID = stardata['sID']
    starRA = stardata['ra']
    starDEC = stardata['dec']
    starMAG = stardata['r']
    starPflag = stardata['Pflag']

    if masknum == 30:
        starPflag[np.where(starID == 's15542')[0]] = -2
        print((starPflag[np.where(starID == 's15542')[0]]))
        starPflag[np.where(starID == 's17368')[0]] = -2
        print((starPflag[np.where(starID == 's17368')[0]]))
        starPflag[np.where(starID == 's16768')[0]] = -2
        print((starPflag[np.where(starID == 's16768')[0]]))
        starPflag[np.where(starID == 's11555')[0]] = -2
        print((starPflag[np.where(starID == 's11555')[0]]))

    elif masknum == 14:
        starPflag[np.where(starID == 's13030')[0]] = -1
        print((starPflag[np.where(starID == 's13030')[0]]))

    else:
        starPflag[np.where(starID == 's15542')[0]] = -1
        print((starPflag[np.where(starID == 's15542')[0]]))
        starPflag[np.where(starID == 's17368')[0]] = -1
        print((starPflag[np.where(starID == 's17368')[0]]))
        starPflag[np.where(starID == 's16768')[0]] = -1
        print((starPflag[np.where(starID == 's16768')[0]]))
        starPflag[np.where(starID == 's11555')[0]] = -1
        print((starPflag[np.where(starID == 's11555')[0]]))
        starPflag[np.where(starID == 's13030')[0]] = -2
        print((starPflag[np.where(starID == 's13030')[0]]))

    ###############################
    ## Define and fill final arrays
    ###############################
    finalID = np.array([], dtype = int)
    finalRA = np.array([])
    finalDEC = np.array([])
    finalMAG = np.array([])
    finalPcode = np.array([], dtype = int)
    finalPA = np.array([])
    
    finalEQ = np.empty(len(finalID))
    finalEQ[:] = 2000.0
    finalPBand = np.chararray(len(finalID))
    finalPBand = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0


    ##################################
    ## Begin with the target selection
    ##################################
    print('target selection...')

    #choose this for each pass, depending on what must be on the mask
    if masknum == 30:
        pflagcut = ((targetPcode > 5000) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 40:
        pflagcut = ((targetPcode > 5000) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 13:
        pflagcut = ((targetPcode > 5000) & (targetINFO != 1)) | ((targetPcode == 100) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 14:
        pflagcut = (((targetPcode > 5000) & (targetINFO != 1)) | ((targetPcode == 100) & (targetINFO != 1)))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 15:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 22:
        pflagcut = ((targetPcode > 5000) & (targetINFO != 1)) | ((targetPcode == 100) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 23:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 32:
        pflagcut = ((targetPcode > 5000) & (targetINFO != 1)) | ((targetPcode == 100) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 33:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))
        
    else:
        pflagcut = ((targetPcode == 10000)) #| (targetPcode == 100))

    tID = targetID[pflagcut]
    tRA = targetRA[pflagcut]
    tDEC = targetDEC[pflagcut]
    tMAG = targetMAG[pflagcut]
    tPcode = targetPcode[pflagcut]
    #tELLIP = targetELLIP[pflagcut]
    #tPA = targetPA[pflagcut]

    objPA = np.empty(len(tMAG),dtype = '<f8')
    if masknum == 30:
        objPA[:] = 40.0
    elif masknum == 10:
        objPA[:] = -78.5
    elif masknum == 20:
        objPA[:] = 2.5
    elif masknum == 22:
        objPA[:] = maskPA+5.0
    elif masknum == 40:
        objPA[:] = -3
    elif masknum == 13:
        objPA[:] = -78.5
    elif masknum == 12:
        objPA[:] = -78.5
    elif masknum == 14:
        objPA[:] = maskPA-5.0
    elif masknum == 32:
        objPA[:] = maskPA-5.0
    else:
        objPA[:] = maskPA+5.0
        #print('wrong PA...')

    finalID = np.append(finalID, tID)
    finalRA = np.append(finalRA, tRA)
    finalDEC = np.append(finalDEC, tDEC)
    finalMAG = np.append(finalMAG, tMAG)
    finalPcode = np.append(finalPcode, tPcode)
    finalPA = np.append(finalPA, objPA)

    print('done with targets')

    print(finalID)

    
    ##################################################
    ## Move to guide and alignment star selection
    ##################################################

    print('alignment stars')

    aligncut = starPflag == -2
    alignID = starID[aligncut]
    alignRA = starRA[aligncut]
    alignDEC = starDEC[aligncut]
    alignMAG = starMAG[aligncut]
    alignPcode = starPflag[aligncut]
    alignPA = np.empty(len(alignID))
    alignPA[:]=0

    finalID = np.append(finalID,alignID)
    finalRA = np.append(finalRA,alignRA)
    finalDEC = np.append(finalDEC,alignDEC)
    finalMAG = np.append(finalMAG,alignMAG)
    finalPcode = np.append(finalPcode,alignPcode)
    finalPA = np.append(finalPA,alignPA)

    print('guide stars...')

    guidecut = starPflag == -1
    guideID = starID[guidecut]
    guideRA = starRA[guidecut]
    guideDEC = starDEC[guidecut]
    guideMAG = starMAG[guidecut]
    guidePcode = starPflag[guidecut]
    guidePA = np.empty(len(guideID))
    guidePA[:]=0

    finalID = np.append(finalID,guideID)
    finalRA = np.append(finalRA,guideRA)
    finalDEC = np.append(finalDEC,guideDEC)
    finalMAG = np.append(finalMAG,guideMAG)
    finalMAG = np.round(finalMAG,decimals=2)
    finalPcode = np.append(finalPcode,guidePcode)
    finalPcode = np.round(finalPcode,decimals=0)
    finalPA = np.append(finalPA,guidePA)
    finalPA = np.round(finalPA,decimals=2)

    print('done with stars')

    ##############################################################
    ## Change the RA and DEC format to match DSIM requirements
    ##############################################################

    print('changing coords...')

    c = SkyCoord(finalRA*u.degree, finalDEC*u.degree)
    outputRA = np.chararray(len(finalID),itemsize = 20)
    outputDEC = np.chararray(len(finalID),itemsize = 20)

    h_ra = c.ra.hms[0]
    m_ra = c.ra.hms[1]
    s_ra = c.ra.hms[2]
    d_dec = c.dec.dms[0]
    m_dec = c.dec.dms[1]
    s_dec = c.dec.dms[2]

    for i in range(len(finalID)):

        if s_ra[i] < 10:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':0'+np.str(np.round(s_ra[i],decimals=3))
        else:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':'+np.str(np.round(s_ra[i],decimals=3))

        if np.sign(s_dec[i]) == 1:
            if s_dec[i] < 10:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

        else: 
            if s_dec[i] < 10:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

    print('done')

    ##############################################################
    ## Add the extra columns which are the same for every object
    ##############################################################
    finalEQ = np.empty(len(finalID), dtype = int)
    finalEQ[:] = 2000
    finalPBand = np.chararray(len(finalID))
    finalPBand[:] = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0

    print((len(finalSEL[np.where(finalSEL == 1)[0]])))
    print((finalID[np.where(finalSEL == 1)[0]]))

    print((len(finalID)))
    print((len(finalRA)))
    print((len(finalDEC)))
    print((len(finalEQ)))
    print((len(finalMAG)))
    print((len(finalPBand)))
    print((len(finalPcode)))
    print((len(finalSAM)))
    print((len(finalSEL)))
    print((len(finalPA)))
    

    ##############################################################
    ## Combine all the columns into one table and save output
    ##############################################################

    print('saving...')

    t = Table([finalID,outputRA,outputDEC,finalEQ,finalMAG,finalPBand,finalPcode,finalSAM,finalSEL,finalPA], names=('ID', 'RA', 'DEC','Equinox','Magnitude','PassBand','Pcode','Sample','Select','SlitPA'), meta={'2018A': 'Keck/DEIMOS mask'})

    t.write(semester+'/maskdesign/'+maskname+'.input', format='ascii', overwrite=True)
    
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

def input2(maskname, semester, maskPA, masknum):

    if semester == '2018A':
        hdulist = fits.open('2018A/newfirm_mbs_deepz_sdss_priority.fits')
    else:
        hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority.fits')
        
    targetdata = hdulist[1].data
    targetINFO = targetdata['target']
    stardata = hdulist[2].data

    if masknum == 10:
        hdum = fits.open('masks_round2/mask10a.fits')
    elif masknum == 12:
        hdum = fits.open('maskdesign/'+semester+'/m12a.fits')
    elif masknum == 13:
        hdum = fits.open('maskdesign/'+semester+'/m13a.fits')
    elif masknum == 14:
        hdum = fits.open('maskdesign/'+semester+'/m14a.fits')
    elif masknum == 20:
        hdum = fits.open('masks_round2/mask20a.fits')
    elif masknum == 22:
        hdum = fits.open('maskdesign/'+semester+'/m22a.fits')
    elif masknum == 30:
        hdum = fits.open('masks_round2/mask30a.fits')
    elif masknum == 32:
        hdum = fits.open('maskdesign/'+semester+'/m32a.fits')
    elif masknum == 40:
        hdum = fits.open('masks_round2/mask40a.fits')
    elif masknum == 70:
        hdum = fits.open(semester+'/maskdesign/m17B70a.fits')
    elif masknum == 15:
        hdum = fits.open(semester+'/maskdesign/m18A15a.fits')
    elif masknum == 23:
        hdum = fits.open(semester+'/maskdesign/m18A23a.fits')
    elif masknum == 33:
        hdum = fits.open(semester+'/maskdesign/m18A33a.fits')
    else:
        print('wrong...try again!')

    mdata = hdum[1].data
    mID = mdata['OBJECT']

    if masknum == 30:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 40:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 13:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 14:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 22:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 32:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))
        
    else:
        print('all objects are being targeted first time')

    #############################
    ## Read in all necessary data
    #############################
    targetID = targetdata['id']
    targetRA = targetdata['ra']
    targetDEC = targetdata['dec']
    targetFLUX = targetdata['R']
    targetMAG = np.empty(len(targetFLUX))
    for i in range(len(targetFLUX)):
        if targetFLUX[i] <=0:
            targetMAG[i] = -99
        else:
            targetMAG[i] = -2.5*np.log10(targetFLUX[i]) + 25.0
            
    targetPcode = targetdata['Pflag']
    targetELLIP = targetdata['K_ellip']
    targetPA = targetdata['K_theta_J2000']
    
    starID = stardata['sID']
    starRA = stardata['ra']
    starDEC = stardata['dec']
    starMAG = stardata['r']
    starPflag = stardata['Pflag']

    if masknum == 30:
        starPflag[np.where(starID == 's15542')[0]] = -2
        print(starPflag[np.where(starID == 's15542')[0]])
        starPflag[np.where(starID == 's17368')[0]] = -2
        print(starPflag[np.where(starID == 's17368')[0]])
        starPflag[np.where(starID == 's16768')[0]] = -2
        print(starPflag[np.where(starID == 's16768')[0]])
        starPflag[np.where(starID == 's11555')[0]] = -2
        print(starPflag[np.where(starID == 's11555')[0]])

    elif masknum == 14:
        starPflag[np.where(starID == 's13030')[0]] = -1
        print((starPflag[np.where(starID == 's13030')[0]]))

    else:
        starPflag[np.where(starID == 's15542')[0]] = -1
        print((starPflag[np.where(starID == 's15542')[0]]))
        starPflag[np.where(starID == 's17368')[0]] = -1
        print((starPflag[np.where(starID == 's17368')[0]]))
        starPflag[np.where(starID == 's16768')[0]] = -1
        print((starPflag[np.where(starID == 's16768')[0]]))
        starPflag[np.where(starID == 's11555')[0]] = -1
        print((starPflag[np.where(starID == 's11555')[0]]))
        starPflag[np.where(starID == 's13030')[0]] = -2
        print((starPflag[np.where(starID == 's13030')[0]]))

    ###############################
    ## Define and fill final arrays
    ###############################
    finalID = np.array([], dtype = int)
    finalRA = np.array([])
    finalDEC = np.array([])
    finalMAG = np.array([])
    finalPcode = np.array([], dtype = int)
    finalPA = np.array([])
    
    finalEQ = np.empty(len(finalID))
    finalEQ[:] = 2000.0
    finalPBand = np.chararray(len(finalID))
    finalPBand = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0


    ##################################
    ## Begin with the target selection
    ##################################
    print('target selection...')

    #choose this for each pass, depending on what must be on the mask
    if masknum == 30:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 40:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    #if masknum == 12:
        #pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        #print(targetINFO[pflagcut])
        #print(np.max(targetINFO[pflagcut]))
        #print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 13:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 14:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 22:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 32:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) | (targetPcode == 100)

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 70:
        pflagcut = ((targetPcode == 10000) | (targetPcode == 5000))

    elif masknum == 23:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)) | ((targetPcode == 100)) | ((targetPcode == 50)))

        print((targetINFO[pflagcut]))
        print(targetPcode[pflagcut])
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 15:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)) | ((targetPcode == 100)) | ((targetPcode == 50)))

        print((targetINFO[pflagcut]))
        print(targetPcode[pflagcut])
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 33:
        pflagcut = (((targetPcode == 10000)) | ((targetPcode == 5000)) | ((targetPcode == 100)) | ((targetPcode == 50)))

        print((targetINFO[pflagcut]))
        print(targetPcode[pflagcut])
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))
        
    else:
        pflagcut = targetPcode > 50

    tID = targetID[pflagcut]
    tRA = targetRA[pflagcut]
    tDEC = targetDEC[pflagcut]
    tMAG = targetMAG[pflagcut]
    tPcode = targetPcode[pflagcut]
    tELLIP = targetELLIP[pflagcut]
    tPA = targetPA[pflagcut]

    objPA = np.empty(len(tPA),dtype = '<f8')
    if masknum == 30:
        objPA[:] = 40.0

    elif masknum == 10:
        objPA[:] = -78.5

    elif masknum == 12:
        objPA[:] = -78.5

    elif masknum == 13:
        objPA[:] = -78.5

    elif masknum == 14:
        objPA[:] = maskPA-5.0

    elif masknum == 20:
        objPA[:] = 2.5

    elif masknum == 22:
        objPA[:] = maskPA+5.0

    elif masknum == 32:
        objPA[:] = maskPA-5.0

    elif masknum == 33:
        objPA[:] = maskPA-5.0
        
    elif masknum == 40:
        objPA[:] = -3
    else:
        objPA[:] = maskPA+5.0
        #print('wrong PA...')

    finalID = np.append(finalID, tID)
    finalRA = np.append(finalRA, tRA)
    finalDEC = np.append(finalDEC, tDEC)
    finalMAG = np.append(finalMAG, tMAG)
    finalPcode = np.append(finalPcode, tPcode)
    finalPA = np.append(finalPA, objPA)

    print('done with targets')
    
    ##################################################
    ## Move to guide and alignment star selection
    ##################################################

    print('alignment stars')

    aligncut = starPflag == -2
    alignID = starID[aligncut]
    alignRA = starRA[aligncut]
    alignDEC = starDEC[aligncut]
    alignMAG = starMAG[aligncut]
    alignPcode = starPflag[aligncut]
    alignPA = np.empty(len(alignID))
    alignPA[:]=0

    finalID = np.append(finalID,alignID)
    finalRA = np.append(finalRA,alignRA)
    finalDEC = np.append(finalDEC,alignDEC)
    finalMAG = np.append(finalMAG,alignMAG)
    finalPcode = np.append(finalPcode,alignPcode)
    finalPA = np.append(finalPA,alignPA)

    print('guide stars...')

    guidecut = starPflag == -1
    guideID = starID[guidecut]
    guideRA = starRA[guidecut]
    guideDEC = starDEC[guidecut]
    guideMAG = starMAG[guidecut]
    guidePcode = starPflag[guidecut]
    guidePA = np.empty(len(guideID))
    guidePA[:]=0

    finalID = np.append(finalID,guideID)
    finalRA = np.append(finalRA,guideRA)
    finalDEC = np.append(finalDEC,guideDEC)
    finalMAG = np.append(finalMAG,guideMAG)
    finalMAG = np.round(finalMAG,decimals=2)
    finalPcode = np.append(finalPcode,guidePcode)
    finalPcode = np.round(finalPcode,decimals=0)
    finalPA = np.append(finalPA,guidePA)
    finalPA = np.round(finalPA,decimals=2)

    print('done with stars')

    ##############################################################
    ## Change the RA and DEC format to match DSIM requirements
    ##############################################################

    print('changing coords...')

    c = SkyCoord(finalRA*u.degree, finalDEC*u.degree)
    outputRA = np.chararray(len(finalID),itemsize = 20)
    outputDEC = np.chararray(len(finalID),itemsize = 20)

    h_ra = c.ra.hms[0]
    m_ra = c.ra.hms[1]
    s_ra = c.ra.hms[2]
    d_dec = c.dec.dms[0]
    m_dec = c.dec.dms[1]
    s_dec = c.dec.dms[2]

    for i in range(len(finalID)):

        if s_ra[i] < 10:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':0'+np.str(np.round(s_ra[i],decimals=3))
        else:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':'+np.str(np.round(s_ra[i],decimals=3))

        if np.sign(s_dec[i]) == 1:
            if s_dec[i] < 10:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

        else: 
            if s_dec[i] < 10:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

    print('done')

    ##############################################################
    ## Add the extra columns which are the same for every object
    ##############################################################
    finalEQ = np.empty(len(finalID), dtype = int)
    finalEQ[:] = 2000
    finalPBand = np.chararray(len(finalID))
    finalPBand[:] = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0

    for i in range(len(mID)):

        cut = finalID == mID[i]
        finalSEL[cut] = 1

    print(len(finalSEL[np.where(finalSEL == 1)[0]]))
    print(finalID[np.where(finalSEL == 1)[0]])

    print(len(finalID))
    print(len(finalRA))
    print(len(finalDEC))
    print(len(finalEQ))
    print(len(finalMAG))
    print(len(finalPBand))
    print(len(finalPcode))
    print(np.unique(finalPcode))
    print(len(finalSAM))
    print(len(finalSEL))
    print(len(finalPA))
    

    ##############################################################
    ## Combine all the columns into one table and save output
    ##############################################################

    print('saving...')

    t = Table([finalID,outputRA,outputDEC,finalEQ,finalMAG,finalPBand,finalPcode,finalSAM,finalSEL,finalPA], names=('ID', 'RA', 'DEC','Equinox','Magnitude','PassBand','Pcode','Sample','Select','SlitPA'), meta={'2018A': 'Keck/DEIMOS'})

    t.write(semester+'/maskdesign/'+maskname+'.input', format='ascii',overwrite=True)


############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

def input3(maskname, semester, maskPA, masknum):

    if semester == '2018A':
        hdulist = fits.open(semester+'/newfirm_mbs_deepz_sdss_priority.fits')
    else:
        hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority.fits')
        
    targetdata = hdulist[1].data
    targetINFO = targetdata['target']
    stardata = hdulist[2].data

    if masknum == 10:
        hdum = fits.open('masks_round2/mask10b.fits')
    elif masknum == 12:
        hdum = fits.open('maskdesign/'+semester+'/m12b.fits')
    elif masknum == 13:
        hdum = fits.open('maskdesign/'+semester+'/m13b.fits')
    elif masknum == 14:
        hdum = fits.open('maskdesign/'+semester+'/m14b.fits')
    elif masknum == 20:
        hdum = fits.open('masks_round2/mask20b.fits')
    elif masknum == 22:
        hdum = fits.open('maskdesign/'+semester+'/m22b.fits')
    elif masknum == 30:
        hdum = fits.open('masks_round2/mask30b.fits')
    elif masknum == 32:
        hdum = fits.open('maskdesign/'+semester+'/m32b.fits')
    elif masknum == 40:
        hdum = fits.open('masks_round2/mask40b.fits')
    elif masknum == 70:
        hdum = fits.open(semester+'/maskdesign/m17B70b.fits')
    elif masknum == 15:
        hdum = fits.open(semester+'/maskdesign/m18A15b.fits')
    elif masknum == 23:
        hdum = fits.open(semester+'/maskdesign/m18A23b.fits')
    elif masknum == 33:
        hdum = fits.open(semester+'/maskdesign/m18A33b.fits')
    else:
        print('wrong...try again!')

    mdata = hdum[1].data
    mID = mdata['OBJECT']

    if masknum == 30:
        data_updateinfo = fits.open('newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 40:
        data_updateinfo = fits.open('newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 13:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 14:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 22:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    elif masknum == 32:
        data_updateinfo = fits.open('maskdesign/'+semester+'/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
        targetINFO = data_updateinfo['target']

        print(len(targetINFO[np.where(targetINFO==1)]))

    else:
        print('all objects are being targeted first time')

    #############################
    ## Read in all necessary data
    #############################
    targetID = targetdata['id']
    targetRA = targetdata['ra']
    targetDEC = targetdata['dec']
    targetFLUX = targetdata['R']
    targetZ = targetdata['z_peak']
    targetMAG = np.empty(len(targetFLUX))
    for i in range(len(targetFLUX)):
        if targetFLUX[i] <=0:
            targetMAG[i] = -99
        else:
            targetMAG[i] = -2.5*np.log10(targetFLUX[i]) + 25.0
            
    targetPcode = targetdata['Pflag']
    targetELLIP = targetdata['K_ellip']
    targetPA = targetdata['K_theta_J2000']
    
    starID = stardata['sID']
    starRA = stardata['ra']
    starDEC = stardata['dec']
    starMAG = stardata['r']
    starPflag = stardata['Pflag']

    if masknum == 30:
        starPflag[np.where(starID == 's15542')[0]] = -2
        print(starPflag[np.where(starID == 's15542')[0]])
        starPflag[np.where(starID == 's17368')[0]] = -2
        print(starPflag[np.where(starID == 's17368')[0]])
        starPflag[np.where(starID == 's16768')[0]] = -2
        print(starPflag[np.where(starID == 's16768')[0]])
        starPflag[np.where(starID == 's11555')[0]] = -2
        print(starPflag[np.where(starID == 's11555')[0]])

    elif masknum == 14:
        starPflag[np.where(starID == 's13030')[0]] = -1
        print((starPflag[np.where(starID == 's13030')[0]]))

    else:
        starPflag[np.where(starID == 's15542')[0]] = -1
        print((starPflag[np.where(starID == 's15542')[0]]))
        starPflag[np.where(starID == 's17368')[0]] = -1
        print((starPflag[np.where(starID == 's17368')[0]]))
        starPflag[np.where(starID == 's16768')[0]] = -1
        print((starPflag[np.where(starID == 's16768')[0]]))
        starPflag[np.where(starID == 's11555')[0]] = -1
        print((starPflag[np.where(starID == 's11555')[0]]))
        starPflag[np.where(starID == 's13030')[0]] = -2
        print((starPflag[np.where(starID == 's13030')[0]]))

    ###############################
    ## Define and fill final arrays
    ###############################
    finalID = np.array([], dtype = int)
    finalRA = np.array([])
    finalDEC = np.array([])
    finalMAG = np.array([])
    finalPcode = np.array([], dtype = int)
    finalPA = np.array([])
    
    finalEQ = np.empty(len(finalID))
    finalEQ[:] = 2000.0
    finalPBand = np.chararray(len(finalID))
    finalPBand = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0


    ##################################
    ## Begin with the target selection
    ##################################
    print('target selection...')

    #choose this for each pass, depending on what must be on the mask
    if masknum == 30:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 40:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    #elif masknum == 12:
        #pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        #print(targetINFO[pflagcut])
        #print(np.max(targetINFO[pflagcut]))
        #print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 13:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ ((targetPcode == 100) & (targetINFO != 1)) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 14:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ ((targetPcode == 100)) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 15:
        pflagcut = (((((targetPcode == 10000)) | ((targetPcode == 5000)) | ((targetPcode == 100)) | ((targetPcode == 50)))) | (((targetPcode == 1000)| (targetPcode == 10) | (targetPcode == 5)) & (targetINFO != 1)) | ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5) & (targetINFO != 1)) | ((targetPcode == 1) & (targetMAG < 25.5) & (targetMAG > 18.) & (targetINFO != 1) & (((targetZ > 1.3) & (targetZ < 2.3)) | (targetZ < 0.5))))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 23:
        pflagcut = (((targetPcode >= 5) & (targetINFO != 1)) | ((targetPcode == 100)) | ((targetPcode == 50))) | ((targetPcode == 1) & (targetMAG < 25.5) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) | (targetZ < 0.3))) | ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5) & (targetINFO != 1))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 22:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 32:
        pflagcut = ((targetPcode > 500) & (targetINFO != 1)) ^ (targetPcode == 100) ^ (targetPcode == 10) ^ (targetPcode == 5) ^ ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) ^ (targetZ < 0.3)))

        print(targetINFO[pflagcut])
        print(np.max(targetINFO[pflagcut]))
        print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    elif masknum == 33:
        pflagcut = (((targetPcode >= 5) & (targetINFO != 1)) | ((targetPcode == 100)) | ((targetPcode == 50)) | ((targetPcode == 1) & (targetMAG < 25.5) & (targetMAG > 18.) & (targetINFO != 1) & ((targetZ > 1.3) | (targetZ < 0.3))) | ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)  & (targetINFO != 1)))

        print((targetINFO[pflagcut]))
        print((np.max(targetINFO[pflagcut])))
        print((len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]])))

    elif masknum == 70:
        pflagcut = ((targetPcode == 10000) ^ (targetPcode == 5000) ^ (targetPcode == 1000))
        
    else:
        pflagcut = (targetPcode >= 5) | ((targetPcode == 1) & (targetZ < 1.3) & (targetZ > 0.3) & (targetMAG > 25.5)) ^ ((targetPcode == 1) & (targetMAG < 24.) & (targetMAG > 18.) & ((targetZ > 1.3) ^ (targetZ < 0.3)))


    tID = targetID[pflagcut]
    tRA = targetRA[pflagcut]
    tDEC = targetDEC[pflagcut]
    tMAG = targetMAG[pflagcut]
    tPcode = targetPcode[pflagcut]
    tELLIP = targetELLIP[pflagcut]
    tPA = targetPA[pflagcut]

    objPA = np.empty(len(tPA),dtype = '<f8')
    if masknum == 30:
        objPA[:] = 40.0

    elif masknum == 10:
        objPA[:] = -78.5

    elif masknum == 12:
        objPA[:] = -78.5

    elif masknum == 13:
        objPA[:] = -78.5

    elif masknum == 14:
        objPA[:] = maskPA-5.0

    elif masknum == 20:
        objPA[:] = 2.5

    elif masknum == 22:
        objPA[:] = maskPA+5.0

    elif masknum == 32:
        objPA[:] = maskPA-5.0

    elif masknum == 33:
        objPA[:] = maskPA-5.0
        
    elif masknum == 40:
        objPA[:] = -3
    else:
        objPA[:] = maskPA+5.0
        #print('wrong PA...')
    
    finalID = np.append(finalID, tID)
    finalRA = np.append(finalRA, tRA)
    finalDEC = np.append(finalDEC, tDEC)
    finalMAG = np.append(finalMAG, tMAG)
    finalPcode = np.append(finalPcode, tPcode)
    finalPA = np.append(finalPA, objPA)

    print('done with targets')

    
    ##################################################
    ## Move to guide and alignment star selection
    ##################################################

    print('alignment stars')

    aligncut = starPflag == -2
    alignID = starID[aligncut]
    alignRA = starRA[aligncut]
    alignDEC = starDEC[aligncut]
    alignMAG = starMAG[aligncut]
    alignPcode = starPflag[aligncut]
    alignPA = np.empty(len(alignID))
    alignPA[:]=0

    finalID = np.append(finalID,alignID)
    finalRA = np.append(finalRA,alignRA)
    finalDEC = np.append(finalDEC,alignDEC)
    finalMAG = np.append(finalMAG,alignMAG)
    finalPcode = np.append(finalPcode,alignPcode)
    finalPA = np.append(finalPA,alignPA)

    print('guide stars...')

    guidecut = starPflag == -1
    guideID = starID[guidecut]
    guideRA = starRA[guidecut]
    guideDEC = starDEC[guidecut]
    guideMAG = starMAG[guidecut]
    guidePcode = starPflag[guidecut]
    guidePA = np.empty(len(guideID))
    guidePA[:]=0

    finalID = np.append(finalID,guideID)
    finalRA = np.append(finalRA,guideRA)
    finalDEC = np.append(finalDEC,guideDEC)
    finalMAG = np.append(finalMAG,guideMAG)
    finalMAG = np.round(finalMAG,decimals=2)
    finalPcode = np.append(finalPcode,guidePcode)
    finalPcode = np.round(finalPcode,decimals=0)
    finalPA = np.append(finalPA,guidePA)
    finalPA = np.round(finalPA,decimals=2)

    print('done with stars')

    ##############################################################
    ## Change the RA and DEC format to match DSIM requirements
    ##############################################################

    print('changing coords...')

    c = SkyCoord(finalRA*u.degree, finalDEC*u.degree)
    outputRA = np.chararray(len(finalID),itemsize = 20)
    outputDEC = np.chararray(len(finalID),itemsize = 20)

    h_ra = c.ra.hms[0]
    m_ra = c.ra.hms[1]
    s_ra = c.ra.hms[2]
    d_dec = c.dec.dms[0]
    m_dec = c.dec.dms[1]
    s_dec = c.dec.dms[2]

    for i in range(len(finalID)):

        if s_ra[i] < 10:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':0'+np.str(np.round(s_ra[i],decimals=3))
        else:
            outputRA[i] = np.str(np.int(h_ra[i]))+':'+np.str(np.int(m_ra[i]))+':'+np.str(np.round(s_ra[i],decimals=3))

        if np.sign(s_dec[i]) == 1:
            if s_dec[i] < 10:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '+'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

        else: 
            if s_dec[i] < 10:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':0'+np.str(np.round(s_dec[i],decimals=2))
            else:
                outputDEC[i] = '-'+np.str(np.int(d_dec[i]))+':'+np.str(np.int(m_dec[i]))+':'+np.str(np.round(s_dec[i],decimals=2))

    print('done')

    ##############################################################
    ## Add the extra columns which are the same for every object
    ##############################################################
    finalEQ = np.empty(len(finalID), dtype = int)
    finalEQ[:] = 2000
    finalPBand = np.chararray(len(finalID))
    finalPBand[:] = 'r'
    finalSAM = np.empty(len(finalID), dtype = int)
    finalSAM[:] = 1
    finalSEL = np.empty(len(finalID), dtype = int)
    finalSEL[:] = 0
    print(mID)
    print(finalID)
    

    for i in range(len(mID)):
        print(mID[i])
        cut = finalID == mID[i]
        print(np.unique(cut))
        finalSEL[cut] = 1

    print(len(finalSEL[np.where(finalSEL == 1)[0]]))
    print(finalID[np.where(finalSEL == 1)[0]])

    print(len(finalID))
    print(len(finalRA))
    print(len(finalDEC))
    print(len(finalEQ))
    print(len(finalMAG))
    print(len(finalPBand))
    print(len(finalPcode))
    print(len(finalSAM))
    print(len(finalSEL))
    print(len(finalPA))
    

    ##############################################################
    ## Combine all the columns into one table and save output
    ##############################################################

    print('saving...')

    t = Table([finalID,outputRA,outputDEC,finalEQ,finalMAG,finalPBand,finalPcode,finalSAM,finalSEL,finalPA], names=('ID', 'RA', 'DEC','Equinox','Magnitude','PassBand','Pcode','Sample','Select','SlitPA'), meta={'2018a': 'deimos mask'})

    t.write(semester+'/maskdesign/'+maskname+'.input', format='ascii',overwrite=True)
