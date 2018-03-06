import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as Coord


def newcat(masklist):
    """This will take the 'newfirm_mbs_deepz_sdss_offset.fits' file and update
    the spfz column with new redshift values from recent observations

    This will follow very closely the catalog_combine.specmatch() function

    This file can then be sent into the 'determine_priority.py' routine
    in order to assign priorities and generate the input for dsimulator
    One main difference is that we will not be matching via astrometric threshold,
    we will be using the newfirm_mbs survey ID's.

    All priorities should already be placed at 1 for any possible science target

    I need to generalize this using better directory structure and modules
    masklist is a list of strings that are the mask names
    example: 'm1016A.2016-05-05'
    """
    
    newID = np.array([])
    newZ = np.array([])
    newQ = np.array([])

    for jj in range(len(masklist)):
        hdu = fits.open('/Users/spf/zspec/zspec.spf.'+masklist[jj]+'.fits')
        data = hdu[1].data

        newID = np.append(newID,data.OBJNAME)
        newZ = np.append(newZ,data.Z)
        newQ = np.append(newQ,data.ZQUALITY)

    return newID, newZ, newQ

################################################################################
## the newcat function needs to be properly called
################################################################################
def newobs(masklist, instrument, semester):
    """
    """
    
    # the new catalog data
    newID, newZ, newQ = newcat(masklist)
    print(len(newID))
    qualitycut = np.where(newQ >= 3)[0]
    qualityID = newID[qualitycut]
    qualityZ = newZ[qualitycut]
    failurecut = np.where(newQ == 2)[0]
    failureID = newID[failurecut]
    failureZ = newZ[failurecut]
    #starcut = np.where(newQ==-1)[0]
    #starID = newID[starcut]
    #starZ = newZ[starcut]
    
    #read in fits data for matching
    if semester == '2017A':
        nfhdu = fits.open('/Users/spf/keck_obs/deimos/red1_catalog_2016A.fits')
    else:
        nfhdu = fits.open('/Users/spf/keck_obs/'+instrument+'/2017A/newfirm_mbs_deepz_sdss_offset.fits')
    
    nfdata = nfhdu[1].data

    nfID_byte = []
    for ii in range(len(nfdata.id)):
        nfID_byte = np.append(nfID_byte,bytes(np.str(nfdata.id[ii]), encoding="UTF-8"))

    nf_id = nfdata.id
    nf_spfz = nfdata.spfz
    nf_target = nfdata.target

    for jj in range(len(qualityID)):
        
        newcut = np.where(qualityID[jj].strip()==nfID_byte)[0]
        nf_spfz[newcut] = qualityZ[jj]

    for kk in range(len(newID)):

        newcut2 = np.where(newID[kk].strip()==nfID_byte)[0]
        nf_target[newcut2] += 1        

    #for kk in range(len(starID)):
        
        #newcut = np.where(starID[kk].strip()==nfID_byte)[0]
        #nf_spfz[newcut] = -99
        
    nfhdu.writeto('/Users/spf/keck_analysis/deimos/red1_catalog_'+semester+'masktracking.fits')
    print(nf_spfz)
    print(nfhdu[1].data.target)
    return nfhdu#, nf_spfz
    

