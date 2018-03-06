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
def assign(semester,rproj=1000.0,dv=1000.0):
    """ This will determine the 3D distance a galaxy is from a group center and assign a priority value based on Rband mag
    
    Parameters
    ----------
    semester : string
        The observational semester that this mask will be observed.
    rproj : float
        The maximum projected distance, think 2D distance in RA and DEC, that two objects can be apart and still possibly a group member.
    dv : float
        The maximum differential velocity between the group center and the galaxy in question for group membership.

    Returns
    -------
    Nothing
    Writes the output table to a fits file in the 'semester' directory, newfirm_mbs_deepz_sdss_priority.fits
    
    """

    if semester == '2017B':
        grouphdu = fits.open('/Users/spf/geec/geec_groupdata.fits')
        groupdata = grouphdu[1].data
        groupRA = groupdata.RAJ2000
        groupDEC = groupdata.DEJ2000
        groupZ = groupdata.z
        groupPriority = np.empty(len(groupZ))
        groupPriority[:] = 1

        zcosmoshdu = fits.open('/Users/spf/zcosmos/zcosmos_groupcat.fits')
        zcosmosdata = zcosmoshdu[1].data
        cosmoscut = (zcosmosdata.RAJ2000<150.)&(zcosmosdata.RAJ2000>149.9)&(zcosmosdata.DEJ2000<2.42)&(zcosmosdata.DEJ2000>2.15)&(zcosmosdata.Ng>10)
        zcosmosdata = zcosmosdata[cosmoscut]
        print(len(zcosmosdata))
        zcosmosRA = zcosmosdata.RAJ2000
        zcosmosDEC = zcosmosdata.DEJ2000
        zcosmosZ = zcosmosdata.__z_
        zcosmosPriority = np.empty(len(zcosmosZ))
        zcosmosPriority[:] = 2
        
    
    else:
        groupdata = Table.read('/Users/spf/research/deep/groupcat_deimos16a.txt', format = 'ascii')
        groupRA = np.array(groupdata['RA'])
        groupDEC = np.array(groupdata['DEC'])
        groupZ = np.array(groupdata['Z'])
        groupPriority = np.array(groupdata['PRIORITY'])

    
    if semester == '2017B':
        hdulist = fits.open('/Users/spf/keck_obs/deimos/2017B/newfirm_mbs_deepz_sdss_offset.fits')
        nf = hdulist[1].data
        sdss = hdulist[2].data

    elif semester == '2018A':
        hdulist = fits.open('/Users/spf/keck_obs/deimos/red1_catalog_2017A.fits')
        nf = hdulist[1].data
        sdss = hdulist[2].data

    else:
        hdulist = fits.open('/Users/spf/keck_obs/deimos/'+semester+'/newfirm_mbs_deepz_sdss_offset_newobs.fits')
        nf = hdulist[1].data
        sdss = hdulist[2].data

    nf_pflag = nf['Pflag']
    nf_starflag = nf['star_flag']

    starcut = nf_starflag == 0

    nf_id = nf['id']
    print(len(nf_id))
    nf_pflag = nf['Pflag']
    nf_ra = nf['ra']
    nf_dec = nf['dec']
    nf_zpeak = nf['z_peak']
    #nf_zspec = nf['z_spec']
    #nf_deep = nf['zcosmos']
    nf_deep = nf['deepz']
    
    if semester == '2017B':
        nf_rmag = nf.r
        nf_kmag = nf.Ks
        weird_rcut = (nf_rmag == -99)
    else:
        nf_rflux = nf['R']
        nf_rmag = -2.5*np.log10(nf_rflux)+25.0
        nf_rmag[np.isnan(nf_rmag)] = 999
        nf_kflux = nf['K']
        nf_kmag = -2.5*np.log10(nf_kflux)+25.0
        nf_nstar = nf['Near_star']

        for i in range(len(nf_id)):
            if ((nf_kmag[i] > 18) & (nf_rflux[i] <= 0)):
                nf_rmag[i] = 99
                nf_rflux[i] = 10**(-29.6)
            else:
                nf_rmag[i] = nf_rmag[i]

        weird_rcut = (nf_rflux > 0)
        
    nf_spfz = nf['spfz']
    nf_pflag[np.where(weird_rcut)[0]] == -50

    sciencecut = starcut #& weird_rcut

    sdss_ra = sdss['ra']
    sdss_dec = sdss['dec']
    sdss_pflag = sdss['Pflag']


    ###########################################################################
    ## Ensure the alignment and guide stars are uniquely identified...
    ###########################################################################

    #sdss_catalog1 = SkyCoord(ra=sdss_ra*u.degree,dec=sdss_dec*u.degree)
    #sdss_catalog2 = SkyCoord(ra=sdss_ra*u.degree,dec=sdss_dec*u.degree)

    #idx, d2d, d3d = sdss_catalog1.match_coordinates_sky(sdss_catalog2,nthneighbor=2)

    #nthcut = d2d > 0.0002778*u.degree

    ###########################################################################
    ## Group Membership estimate
    ###########################################################################

    zcosmosflag = np.empty(len(nf_ra))
    zcosmosflag[:]=0
    groupflag = np.empty(len(nf_ra))
    groupflag[:]=0
        
    if semester == '2017A':
        
        ###########################################################################
        ## Determine whether each object is nearby a geec group
        ###########################################################################

        Gz = groupZ#[lowcut]
        Gra = groupRA#[lowcut]
        Gdec = groupDEC#[lowcut]

        convert = cosmo.arcsec_per_kpc_comoving(Gz)
        print(convert)
        projcut = (convert * rproj) * (u.kpc / u.arcsec) / 3600.
        print(projcut)

        for j in range(len(Gz)):

            dz = dv * (1 + Gz[j]) / (3.e5)
            
            ra_diff = nf_ra - Gra[j]
            dec_diff = nf_dec - Gdec[j]
            z_diff = nf_zpeak - Gz[j]

            ratest = np.abs(ra_diff) < projcut[j]
            dectest = np.abs(dec_diff) < projcut[j]
            ztest = np.abs(z_diff) < dz

            groupflag[np.where(ratest & dectest & ztest)[0]] = 1


        ###########################################################################
        ## Determine whether each object is nearby a zcosmos group
        ###########################################################################
        Gz = zcosmosZ#[lowcut]
        Gra = zcosmosRA#[lowcut]
        Gdec = zcosmosDEC#[lowcut]

        convert = cosmo.arcsec_per_kpc_comoving(Gz)
        print(convert)
        projcut = (convert * rproj) * (u.kpc / u.arcsec) / 3600.
        print(projcut)

        for j in range(len(Gz)):

            dz = dv * (1 + Gz[j]) / (3.e5)
            
            ra_diff = nf_ra - Gra[j]
            dec_diff = nf_dec - Gdec[j]
            z_diff = nf_zpeak - Gz[j]

            ratest = np.abs(ra_diff) < projcut[j]
            dectest = np.abs(dec_diff) < projcut[j]
            ztest = np.abs(z_diff) < dz

            zcosmosflag[np.where(ratest & dectest & ztest)[0]] = 1

    else:

        ###########################################################################
        ## Determine whether each object is nearby a group from catalog
        ###########################################################################

        Gz = groupZ
        Gra = groupRA
        Gdec = groupDEC
        

        convert = cosmo.arcsec_per_kpc_comoving(Gz)
        print(convert)
        projcut = (convert * rproj) * (u.kpc / u.arcsec) / 3600.
        print(projcut)

        for j in range(len(Gz)):

            dz = dv * (1 + Gz[j]) / (3.e5)
            
            ra_diff = nf_ra - Gra[j]
            dec_diff = nf_dec - Gdec[j]
            z_diff = nf_zpeak - Gz[j]

            ratest = np.abs(ra_diff) < projcut[j]
            dectest = np.abs(dec_diff) < projcut[j]
            ztest = np.abs(z_diff) < dz

            groupflag[np.where(ratest & dectest & ztest)[0]] = groupPriority[j]

    ###########################################################################
    ## Tier 6: Setup during initial catalog combine by setting everything
    ## to Pflag = 1
    ## This step should be done when the file is initialized in the
    ## catalog_combine routine
    ##############################
    ## Tier 5: This group is for really faint objects near groups.
    ## not been targeted by DEEP, set at a Pflag = 5
    ###########################################################################
    grouptest = groupflag == 1
    grouptest2 = groupflag == 2
    zcosmostest = zcosmosflag == 1
    deeptest = (nf_deep < 0) & (nf_spfz < 0)
    magtest = nf_rmag > 25.5

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 5
    nf_pflag[np.where(grouptest2 & deeptest & magtest & sciencecut)[0]] = 5
    nf_pflag[np.where(zcosmostest & deeptest & magtest & sciencecut)[0]] = 5

    print(np.unique(nf_pflag))
    print('very faint group members')
    print(len(nf_pflag[np.where(nf_pflag==5)[0]]))

    ###########################################################################
    ## Tier 4: This group is, for intermediate faint objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, set at a Pflag = 10
    ###########################################################################
    deeptest = (nf_deep < 0) & (nf_spfz < 0)
    ztest = (nf_zpeak < 1.3) & (nf_zpeak > 0.3)
    magtest = (nf_rmag < 25.5) & (nf_rmag > 24)

    nf_pflag[np.where(ztest & deeptest & magtest & sciencecut)[0]] = 10

    print(np.unique(nf_pflag))
    print(len(nf_pflag[np.where(nf_pflag==10)[0]]))
    
    ####################################################################
    ## Tier 3: This group is for intermediate faint objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 100
    ####################################################################
    grouptest = groupflag == 1
    grouptest2 = groupflag == 2
    zcosmostest = zcosmosflag == 1
    deeptest = (nf_deep < 0) & (nf_spfz < 0)
    magtest = (nf_rmag < 25.5) & (nf_rmag > 24)
    print(np.unique(nf_rmag))

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 100
    nf_pflag[np.where(grouptest2 & deeptest & magtest & sciencecut)[0]] = 50
    nf_pflag[np.where(zcosmostest & deeptest & magtest & sciencecut)[0]] = 50

    print(np.unique(nf_pflag))
    print('geec faint members')
    print(len(nf_pflag[np.where(nf_pflag==100)[0]]))
    print('zcosmos faint members')
    print(len(nf_pflag[np.where(nf_pflag==50)[0]]))

    ###########################################################################
    ## Tier 2: This group is, for relatively bright objects that have
    ## not been targeted by DEEP, and at 0.3 < z < 1.3, at a Pflag = 1000
    ###########################################################################
    deeptest = (nf_deep < 0) & (nf_spfz < 0)
    ztest = (nf_zpeak < 1.3) & (nf_zpeak > 0.3)
    magtest = (nf_rmag < 24.2)

    nf_pflag[np.where(ztest & deeptest & magtest & sciencecut)[0]] = 1000

    print(np.unique(nf_pflag))
    print('bright field')
    print(len(nf_pflag[np.where(nf_pflag==1000)[0]]))
        
    ####################################################################
    ## Tier 1: This group is for relatively bright objects near groups,
    ## that have not been targeted by DEEP, at a Pflag = 10000
    ####################################################################
    grouptest = groupflag == 1
    grouptest2 = groupflag == 2
    zcosmostest = zcosmosflag == 1
    deeptest = (nf_deep < 0) & (nf_spfz < 0)
    magtest = (nf_rmag < 24.2)

    nf_pflag[np.where(grouptest & deeptest & magtest & sciencecut)[0]] = 10000
    nf_pflag[np.where(grouptest2 & deeptest & magtest & sciencecut)[0]] = 5000
    nf_pflag[np.where(zcosmostest & deeptest & magtest & sciencecut)[0]] = 5000

    print(np.unique(nf_pflag))
    print('group1 bright members')
    print(len(nf_pflag[np.where(nf_pflag==10000)[0]]))
    print('group2 bright members')
    print(len(nf_pflag[np.where(nf_pflag==5000)[0]]))
    print('group1 faint members')
    print(len(nf_pflag[np.where(nf_pflag==100)[0]]))
    print('group2 faint members')
    print(len(nf_pflag[np.where(nf_pflag==50)[0]]))
    print(np.sum(grouptest))
    print(np.sum(grouptest & (nf_spfz > 0)))
    print(np.sum(grouptest2))
    print(np.sum(grouptest2 & (nf_spfz > 0)))


    #######################################################################
    ## Write the new hdu to a file that will be used to generate the DSIM
    ## input files.
    #######################################################################

    if semester == '2017B':
        hdulist.writeto(semester+'/COSMOS2015_specz_sdss_priority.fits', overwrite=True)
    else:
        hdulist.writeto(semester+'/newfirm_mbs_deepz_sdss_priority.fits', overwrite=True)

    #return hdulist


