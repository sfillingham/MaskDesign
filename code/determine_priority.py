import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.coordinates import SkyCoord
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import pdb

###########################################################################
##
## This will select galaxies in the NEWFIRM_MBS above a given magnitude in
## specified band, see below. initial_priority will take all targets and
## assign priorities to all objects. dsim_priority will take the latest output
## of DSIM and the latest priority file and lower the priority of objects which
## have already been targeted on previous masks.
##
## maglim, band corresponds to the magnitude and broad band filter which the
## cut will be based on.
##
## rproj and dz allow for a general selection around the group centers
##
## The final output will be a file that is ready for DSIMULATOR input routine
## 'make_dsim_input.py'
##
## All priorities should already be placed at 1 for any possible science target
##
##
############################################################################
def pass1(rproj,dz):

    groupdata = Table.read('/Users/Sean/research/deep/groupcat_deimos16a.txt', format = 'ascii')
    groupRA = np.array(groupdata['RA'])
    groupDEC = np.array(groupdata['DEC'])
    groupZ = np.array(groupdata['Z'])
    groupPriority = np.array(groupdata['PRIORITY'])
    
    hdulist = fits.open('newfirm_mbs_deepz_sdss_offset.fits')
    nf = hdulist[1].data
    sdss = hdulist[2].data

    nf_pflag = nf['Pflag']
    nf_starflag = nf['star_flag']

    starcut = nf_starflag == 0

    nf_id = nf['id']
    print len(nf_id)
    nf_pflag = nf['Pflag']
    nf_ra = nf['ra']
    nf_dec = nf['dec']
    nf_zpeak = nf['z_peak']
    nf_zspec = nf['z_spec']
    nf_deep = nf['deepz']
    nf_nstar = nf['Near_star']
    nf_rflux = nf['R']
    nf_rmag = -2.5*np.log10(nf_rflux)+25.0
    nf_kflux = nf['K']
    nf_kmag = -2.5*np.log10(nf_kflux)+25.0

    for i in range(len(nf_id)):
        if ((nf_kmag[i] > 18) & (nf_rflux[i] <= 0)):
            nf_rmag[i] = 99
            nf_rflux[i] = 10**(-29.6)
        else:
            nf_rmag[i] = nf_rmag[i]

    weird_rcut = (nf_rflux > 0)

    nf_pflag[np.where(weird_rcut)[0]] == -50

    sciencecut = starcut & weird_rcut

    #nf_id = nf_id[weird_rcut]
    #nf_pflag = nf_pflag[weird_rcut]
    #nf_ra = nf_ra[weird_rcut]
    #nf_dec = nf_dec[weird_rcut]
    #nf_zpeak = nf_zpeak[weird_rcut]
    #nf_zspec = nf_zspec[weird_rcut]
    #nf_deep = nf_deep[weird_rcut]
    #nf_nstar = nf_nstar[weird_rcut]
    #nf_rflux = nf_rflux[weird_rcut]
    #nf_rmag = nf_rmag[weird_rcut]
    #nf_kmag = nf_kmag[weird_rcut]

    sdss_ra = sdss['ra']
    sdss_dec = sdss['dec']
    sdss_pflag = sdss['Pflag']


    ###########################################################################
    ## Ensure the alignment and guide stars are uniquely identified...hopefully
    ###########################################################################

    #sdss_catalog1 = SkyCoord(ra=sdss_ra*u.degree,dec=sdss_dec*u.degree)
    #sdss_catalog2 = SkyCoord(ra=sdss_ra*u.degree,dec=sdss_dec*u.degree)

    #idx, d2d, d3d = sdss_catalog1.match_coordinates_sky(sdss_catalog2,nthneighbor=2)

    #nthcut = d2d > 0.0002778*u.degree

    ###########################################################################
    ## Determine whether each object is nearby a group
    ###########################################################################

    Gz = groupZ#[lowcut]
    Gra = groupRA#[lowcut]
    Gdec = groupDEC#[lowcut]

    convert = cosmo.arcsec_per_kpc_comoving(Gz)
    print convert
    projcut = (convert*rproj)*(u.kpc/u.arcsec)/3600.
    print projcut

    groupflag = np.empty(len(nf_ra))
    groupflag[:]=0

    for j in range(len(Gz)):

        ra_diff = nf_ra - Gra[j]
        dec_diff = nf_dec - Gdec[j]
        z_diff = nf_zpeak - Gz[j]

        ratest = np.abs(ra_diff) < projcut[j]
        dectest = np.abs(dec_diff) < projcut[j]
        ztest = np.abs(z_diff) < dz

        groupflag[np.where(ratest & dectest & ztest)[0]] = 1
    
    #pdb.set_trace()
        
    ###########################################################################
    ## Tier 6: Setup during initial catalog combine by setting everything
    ## to Pflag = 1
    ##############################
    ## Tier 5: This group is for really faint objects near groups.
    ## not been targeted by DEEP, set at a Pflag = 5
    ###########################################################################
    grouptest = groupflag == 1
    deeptest = (nf_deep < 3) & (nf_deep != -1)
    magtest = nf_rmag > 25.5

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 5

    print np.max(nf_pflag)
    print len(nf_pflag[np.where(nf_pflag==5)[0]])

    ###########################################################################
    ## Tier 4: This group is, for intermediate faint objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, set at a Pflag = 10
    ###########################################################################
    #galtest = nf_starflag[i] < 1
    deeptest = (nf_deep < 3) & (nf_deep != -1)
    ztest = (nf_zpeak < 1.3) & (nf_zpeak > 0.3)
    magtest = (nf_rmag < 25.5) & (nf_rmag > 24)

    nf_pflag[np.where(ztest & deeptest & magtest & sciencecut)[0]] = 10

    print np.max(nf_pflag)
    print len(nf_pflag[np.where(nf_pflag==10)[0]])
    
    ####################################################################
    ## Tier 3: This group is for intermediate faint objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 100
    ####################################################################
    grouptest = groupflag == 1
    deeptest = (nf_deep < 3) & (nf_deep != -1)
    magtest = (nf_rmag < 25.5) & (nf_rmag > 24)

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 100

    print np.max(nf_pflag)
    print len(nf_pflag[np.where(nf_pflag==100)[0]])

    ###########################################################################
    ## Tier 2: This group is, for relatively bright objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, at a Pflag = 1000
    ###########################################################################
    #galtest = nf_starflag[i] < 1
    deeptest = (nf_deep < 3) & (nf_deep != -1)
    ztest = (nf_zpeak < 1.3) & (nf_zpeak > 0.3)
    magtest = (nf_rmag < 24)

    nf_pflag[np.where(ztest & deeptest & magtest & sciencecut)[0]] = 1000

    print np.max(nf_pflag)
    print len(nf_pflag[np.where(nf_pflag==1000)[0]])
        
    ####################################################################
    ## Tier 1: This group is for relatively bright objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 10000
    ####################################################################
    grouptest = groupflag == 1
    deeptest = (nf_deep < 3) & (nf_deep != -1)
    magtest = (nf_rmag < 24)

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 10000

    print np.max(nf_pflag)
    print len(nf_pflag[np.where(nf_pflag==10000)[0]])

    #######################################################################
    ## Write the new hdu to a file that will be used to generate the DSIM
    ## input files.
    #######################################################################

    #hdulist.writeto('maskdesign/newfirm_mbs_deepz_sdss_priority_pass1.fits')

    return hdulist


#######################################################################
#######################################################################
##
## This is the Pass 2 priority function.
## It will take the DSIM output from Pass 1 and reassign priorities
## based on the science and what has already been targeted.
##
#######################################################################
#######################################################################


def pass2(rproj,dz):

    #Group center data
    groupdata = Table.read('/Users/Sean/research/deep/groupcat_deimos16a.txt', format = 'ascii')
    groupRA = np.array(groupdata['RA'])
    groupDEC = np.array(groupdata['DEC'])
    groupZ = np.array(groupdata['Z'])
    groupPriority = np.array(groupdata['PRIORITY'])

    #target data
    hdulist = fits.open('newfirm_mbs_deepz_sdss_offset.fits')
    nf = hdulist[1].data
    sdss = hdulist[2].data

    nf_starflag = nf['star_flag']
    starcut = nf_starflag == 0

    nf_id = nf['id'][starcut]
    nf_pflag = nf['Pflag'][starcut]
    nf_ra = nf['ra'][starcut]
    nf_dec = nf['dec'][starcut]
    nf_zpeak = nf['z_peak'][starcut]
    nf_zspec = nf['z_spec'][starcut]
    nf_deep = nf['deepz'][starcut]
    nf_nstar = nf['Near_star'][starcut]
    nf_rflux = nf['R'][starcut]
    nf_rmag = -2.5*np.log10(nf_rflux)+25.0
    nf_kflux = nf['K'][starcut]
    nf_kmag = -2.5*np.log10(nf_kflux)+25.0
    nf_targ = nf['target'][starcut]

    sdss_ra = sdss['ra'][starcut]
    sdss_dec = sdss['dec'][starcut]
    sdss_pflag = sdss['Pflag'][starcut]

    #Pass 1 outputs, should be 2 mask files in the FITS format
    hdu1 = fits.open('maskdesign/mask1.fits')
    hdu2 = fits.open('maskdesign/mask2.fits')
    
    m1 = hdu1[1].data
    m2 = hdu2[1].data

    m1_ID = m1['ObjectId']
    m2_ID = m2['ObjectId']

    print len(m1_ID)
    print len(m2_ID)

    for i in range(len(m1_ID)):

        targcut = m1_ID[i] == nf_id
        nf_targ[targcut] = 1

    for i in range(len(m2_ID)):

        targcut = m2_ID[i] == nf_id

        if (nf_targ[targcut]==1):
            nf_targ[targcut] = 3
        else:
            nf_targ[targcut] = 2
        

    ###########################################################################
    ## Tier 6: Setup during initial catalog combine by setting everything
    ## to Pflag = 1
    ##############################
    ## Tier 5: This group is for really faint objects near groups.
    ## not been targeted by DEEP, set at a Pflag = 5
    ###########################################################################

    #lowcut = groupPriority == 2
    Gz = groupZ#[lowcut]
    Gra = groupRA#[lowcut]
    Gdec = groupDEC#[lowcut]

    convert = cosmo.arcsec_per_kpc_comoving(Gz)
    print convert
    projcut = (convert*rproj)*(u.kpc/u.arcsec)/3600.
    print projcut

    for j in range(len(Gz)):

        ra_diff = nf_ra - Gra[j]
        dec_diff = nf_dec - Gdec[j]
        z_diff = nf_zpeak - Gz[j]

        for i in range(len(ra_diff)):

            projtest = (np.abs(ra_diff[i]) < projcut[j]) and (np.abs(dec_diff[i]) < projcut[j])
            ztest = np.abs(z_diff[i]) < dz
            deeptest = (nf_deep[i] < 3) & (nf_deep[i] != -1)
            magtest = nf_rmag[i] > 25.5

            finalcut = projtest & ztest & deeptest & magtest

            if finalcut:

                nf_pflag[i] = 5

            else:

                nf_pflag[i] = nf_pflag[i]

    ###########################################################################
    ## Tier 4: This group is, for intermediate faint objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, set at a Pflag = 10
    ###########################################################################

    for i in range(len(nf_id)):

        #galtest = nf_starflag[i] < 1
        deeptest = (nf_deep[i] < 3) & (nf_deep[i] != -1)
        ztest = (nf_zpeak[i] < 1.3) & (nf_zpeak[i] > 0.3)
        magtest = (nf_rmag[i] < 25.5) & (nf_rmag[i] > 24)

        if (deeptest & ztest & magtest):

            nf_pflag[i] = 10

        else:

            nf_pflag[i] = nf_pflag[i]

    
    ####################################################################
    ## Tier 3: This group is for intermediate faint objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 100
    ####################################################################

    #lowcut = groupPriority == 2
    Gz = groupZ#[lowcut]
    Gra = groupRA#[lowcut]
    Gdec = groupDEC#[lowcut]

    convert = cosmo.arcsec_per_kpc_comoving(Gz)
    print convert
    projcut = (convert*rproj)*(u.kpc/u.arcsec)/3600.
    print projcut

    for j in range(len(Gz)):

        ra_diff = nf_ra - Gra[j]
        dec_diff = nf_dec - Gdec[j]
        z_diff = nf_zpeak - Gz[j]

        for i in range(len(ra_diff)):

            projtest = (np.abs(ra_diff[i]) < projcut[j]) and (np.abs(dec_diff[i]) < projcut[j])
            ztest = np.abs(z_diff[i]) < dz
            deeptest = (nf_deep[i] < 3) & (nf_deep[i] != -1)
            magtest = (nf_rmag[i] < 25.5) & (nf_rmag[i] > 24)

            finalcut = projtest & ztest & deeptest & magtest

            if finalcut:

                nf_pflag[i] = 100

            else:

                nf_pflag[i] = nf_pflag[i]

    ###########################################################################
    ## Tier 2: This group is, for relatively bright objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, at a Pflag = 1000
    ###########################################################################

    for i in range(len(nf_id)):

        #galtest = nf_starflag[i] < 1
        deeptest = (nf_deep[i] < 3) & (nf_deep[i] != -1)
        ztest = (nf_zpeak[i] < 1.3) & (nf_zpeak[i] > 0.3)
        magtest = (nf_rmag[i] < 24)
        targetcut = nf_targ[i] == 0

        if (deeptest & ztest & magtest & targetcut):

            nf_pflag[i] = 1000

        else:

            nf_pflag[i] = nf_pflag[i]
        
    ####################################################################
    ## Tier 1: This group is for relatively bright objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 10000
    ####################################################################

    #lowcut = groupPriority == 2
    Gz = groupZ#[lowcut]
    Gra = groupRA#[lowcut]
    Gdec = groupDEC#[lowcut]

    convert = cosmo.arcsec_per_kpc_comoving(Gz)
    print convert
    projcut = (convert*rproj)*(u.kpc/u.arcsec)/3600.
    print projcut

    for j in range(len(Gz)):

        ra_diff = nf_ra - Gra[j]
        dec_diff = nf_dec - Gdec[j]
        z_diff = nf_zpeak - Gz[j]

        for i in range(len(ra_diff)):

            projtest = (np.abs(ra_diff[i]) < projcut[j]) and (np.abs(dec_diff[i]) < projcut[j])
            ztest = np.abs(z_diff[i]) < dz
            deeptest = (nf_deep[i] < 3) & (nf_deep[i] != -1)
            magtest = (nf_rmag[i] < 24)
            targetcut = nf_targ[i] == 0

            finalcut = projtest & ztest & deeptest & magtest & targetcut

            if finalcut:

                nf_pflag[i] = 10000

            else:

                nf_pflag[i] = nf_pflag[i]

    #######################################################################
    ## Write the new hdu to a file that will be used to generate the DSIM
    ## input files.
    #######################################################################

    hdulist.writeto('maskdesign/newfirm_mbs_deepz_sdss_priority_pass2.fits')
    
