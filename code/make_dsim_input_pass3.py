import numpy as np
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)

###########################################################################
##
## Makes the input file based on the DSIM standards via the KECK webpage.
## This routine first requires the output from 'determine_priority.py' and
## all subsequent codes in the pipeline. See 'README_mynotes' for details.
##
############################################################################

def input1(maskname, maskPA, passnum, inputnum):

    hdulist = fits.open('maskdesign/newfirm_mbs_deepz_sdss_priority_pass2.fits')
    targetdata = hdulist[1].data
    stardata = hdulist[2].data

    data_updateinfo = fits.open('maskdesign/newfirm_mbs_deepz_sdss_priority_pass3.fits')[1].data
    targetINFO = data_updateinfo['target']

    len(targetINFO[np.where(targetINFO==1)])

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
    #targetINFO = targetdata['target']
    
    starID = stardata['sID']
    starRA = stardata['ra']
    starDEC = stardata['dec']
    starMAG = stardata['r']
    starPflag = stardata['Pflag']

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
    pflagcut = ((targetPcode == 10000) & (targetINFO != 1)) ^ ((targetPcode == 100) & (targetINFO != 1)) ^ ((targetPcode == 1000) & (targetINFO != 1)) ^ ((targetPcode == 10) & (targetINFO != 1))
    print(targetINFO[pflagcut])
    print(np.max(targetINFO[pflagcut]))
    print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    tID = targetID[pflagcut]
    tRA = targetRA[pflagcut]
    tDEC = targetDEC[pflagcut]
    tMAG = targetMAG[pflagcut]
    tPcode = targetPcode[pflagcut]
    tELLIP = targetELLIP[pflagcut]
    tPA = targetPA[pflagcut]

    objPA = np.empty(len(tPA),dtype = '<f8')
    
    for i in range(len(tELLIP)):

        if (tELLIP[i] > 0.3):

            PAdiff = tPA[i] - maskPA

            if np.abs(PAdiff) < 30:

                objPA[i] = tPA[i]

            else:

                if np.sign(PAdiff) == 1:
                    objPA[i] = maskPA + 30
                else:
                    objPA[i] = maskPA - 30

        else:

            objPA[i] = maskPA + 5.0

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

    t = Table([finalID,outputRA,outputDEC,finalEQ,finalMAG,finalPBand,finalPcode,finalSAM,finalSEL,finalPA], names=('ID', 'RA', 'DEC','Equinox','Magnitude','PassBand','Pcode','Sample','Select','SlitPA'), meta={'2016a': 'deimos mask'})

    t.write('dsim_input/'+maskname+'.input', format='ascii')

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

def input2(maskname, maskPA, passnum, inputnum):

    hdulist = fits.open('maskdesign/newfirm_mbs_deepz_sdss_priority_pass'+np.str(passnum)+'.fits')#_mag'+np.str(magnify)+'EX.fits')
    targetdata = hdulist[1].data
    stardata = hdulist[2].data

    data_updateinfo = fits.open('maskdesign/newfirm_mbs_deepz_sdss_priority_pass2.fits')[1].data
    targetINFO = data_updateinfo['target']

    if inputnum == 1:
        hdum = fits.open('maskdesign/mask3a.fits')
    elif inputnum == 2:
        hdum = fits.open('maskdesign/mask4a.fits')
    else:
        print('wrong...try again!')

    mdata = hdum[1].data
    mID = mdata['OBJECT']

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
    pflagcut = ((targetPcode == 10000) & (targetINFO != 1)) ^ (targetPcode == 100) ^ ((targetPcode == 1000) & (targetINFO != 1)) ^ (targetPcode == 10)
    print(targetINFO[pflagcut])
    print(np.max(targetINFO[pflagcut]))
    print(len(targetINFO[pflagcut][np.where(targetINFO[pflagcut]==1)[0]]))

    tID = targetID[pflagcut]
    tRA = targetRA[pflagcut]
    tDEC = targetDEC[pflagcut]
    tMAG = targetMAG[pflagcut]
    tPcode = targetPcode[pflagcut]
    tELLIP = targetELLIP[pflagcut]
    tPA = targetPA[pflagcut]

    objPA = np.empty(len(tPA),dtype = '<f8')
    
    for i in range(len(tELLIP)):

        if (tELLIP[i] > 0.3):

            PAdiff = tPA[i] - maskPA

            if np.abs(PAdiff) < 30:

                objPA[i] = tPA[i]

            else:

                if np.sign(PAdiff) == 1:
                    objPA[i] = maskPA + 30
                else:
                    objPA[i] = maskPA - 30

        else:

            objPA[i] = maskPA + 5.0

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

    t = Table([finalID,outputRA,outputDEC,finalEQ,finalMAG,finalPBand,finalPcode,finalSAM,finalSEL,finalPA], names=('ID', 'RA', 'DEC','Equinox','Magnitude','PassBand','Pcode','Sample','Select','SlitPA'), meta={'2016a': 'deimos mask'})

    t.write('dsim_input/'+maskname+'.input', format='ascii')
