import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as Coord


def check(ra,dec,z,rproj=1000.0,dz=0.05):

    projtest = len(ra) != len(dec)
    ztest = len(ra) != len(z)

    output = np.empty(len(ra))
    output[:]=0

    if (projtest^ztest):
        print('coordinate objects need to be the same size')
        print('Try Again...')

        return output

    else:

    ###########################################################################
    ## Determine whether each object is nearby a group
    ###########################################################################

        groupdata = Table.read('/Users/spf/research/deep/groupcat_deimos16a.txt', format = 'ascii')
        groupRA = np.array(groupdata['RA'])
        groupDEC = np.array(groupdata['DEC'])
        groupZ = np.array(groupdata['Z'])
        groupPriority = np.array(groupdata['PRIORITY'])

        for ii in range(len(ra)):

            radiff = ra[ii] - groupRA
            decdiff = dec[ii] - groupDEC
            zdiff = z[ii] - groupZ

            convert = cosmo.arcsec_per_kpc_comoving(groupZ)
            projcut = (convert*rproj)*(u.kpc/u.arcsec)/3600.
            #print(projcut)

            ratest = np.abs(radiff) < projcut
            dectest = np.abs(decdiff) < projcut
            ztest = np.abs(zdiff) < dz


            if (len(radiff[ratest & dectest & ztest])>0):
                output[ii] = 1

            else:
                output[ii] = output[ii]

        return output

            






            
