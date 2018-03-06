import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.coordinates as Coord

################################################################################
## Specify the threshold in degrees. This is the maximum distance that any two
## matched objects can lie in different catalog.
################################################################################
def specmatch(instrument,semester,field,threshold='default'):

    if threshold == 'default':

        limit = 0.0002778

    else:

        limit = threshold

    #read in fits data for matching
    if field == 'egs':
        deephdulist = fits.open('/Users/spf/research/deep/zcat.alldeep.2015jul27.fits')
        specdata = deephdulist[1].data
        spec_ra = specdata['RA']
        spec_dec = specdata['DEC']
        spec_zspec = specdata['Z']
        spec_ID = specdata['OBJNO']
        zqual = specdata['ZQUALITY']
        
    elif field == 'cosmos':
        specdata = ascii.read('/Users/spf/cosmos/cosmos_candels_public.dat')
        spec_ra = specdata['RA']
        spec_dec = specdata['DEC']
        spec_zspec = specdata['z_spec']
        zqual = specdata['z_quality']
        zqual[np.where(zqual == 'A')[0]] = 4.0
        zqual[np.where(zqual == 'C')[0]] = 1.0
        
    else:
        print('choose a real field')

    if field == 'cosmos':
        nfhdulist = fits.open('/Users/spf/keck_obs/'+instrument+'/'+semester+'/COSMOS2015_catalog.fits')
    else:
        nfhdulist = fits.open('/Users/spf/keck_obs/'+instrument+'/'+semester+'/newfirm_mbs_catalog.fits')
    nfdata = nfhdulist[1].data
    nf_ra = nfdata['ra']
    nf_dec = nfdata['dec']
    if field == 'egs':
        nf_zspec = nfdata['deepz']
    elif field == 'cosmos':
        nf_zspec = nfdata['zcosmos']
    else:
        print('choose a real field')
    nf_id = nfdata['ID']

    deepID_out = np.empty(len(nf_ra))
    nmbsID = np.empty(len(nf_ra))

    speccatalog = SkyCoord(ra=spec_ra*u.degree, dec = spec_dec*u.degree)
    nfcatalog = SkyCoord(ra=nf_ra*u.degree, dec = nf_dec*u.degree)

    idx,d2d,d3d = nfcatalog.match_to_catalog_sky(speccatalog)

    d2d = d2d/u.degree

    #use the threshold to determine matching objects and insert specz accordingly
    for i in range(len(idx)):
        location = idx[i]
        zq = zqual[location].astype(np.float)
        testgood = (d2d[i] < limit) and (zq > 2.8)
        testbad = (d2d[i] < limit) and (zq < 2.8) and (zq > 0)
        teststar = (d2d[i] < limit) and (zq == -1)
        #print test
        

        if testgood == True:
            nf_zspec[i] = zqual[location]
            #deepID_out[i] = deep_ID[location]
            #nmbsID[i] = nf_id[i]
        elif testbad == True:
            nf_zspec[i] = zqual[location]
            #deepID_out[i] = deep_ID[location]
            #nmbsID[i] = nf_id[location]
        elif teststar == True:
            nf_zspec[i] = zqual[location]
            #deepID_out[i] = deep_ID[location]
            #nmbsID[i] = nf_id[location]
        else:
            nf_zspec[i] = -99
            #deepID_out[i] = -99
            #nmbsID[i] = nf_id[location]

    if field == 'cosmos':
        nfhdulist.writeto(semester+'/COSMOS2015_specz_catalog.fits')
    else:
        nfhdulist.writeto(semester+'/newfirm_mbs_specz_catalog.fits')

    #return deepID_out, nmbsID

################################################################################################
## Check for and remove any duplicate stars or galaxies
## from the SDSS catalog.
################################################################################################

def doops(input_ra,input_dec):
    print(len(input_ra))

    c = SkyCoord(ra=input_ra*u.degree,dec=input_dec*u.degree)
    h_ra = c.ra.hms[0]
    m_ra = c.ra.hms[1]
    s_ra = c.ra.hms[2]
    d_dec = c.dec.dms[0]
    m_dec = c.dec.dms[1]
    s_dec = c.dec.dms[2]

    ra_check1 = np.round(s_ra,decimals=2)
    RAvalue1, RAindex1 = np.unique(ra_check1, return_index=True)
    #print RAvalue1
    #print RAindex1
    dec_check1 = np.round(s_dec[RAindex1], decimals=1)
    DECvalue1, DECindex1 = np.unique(dec_check1, return_index=True)
    #print len(DECvalue1)


    dec_check2 = np.round(s_dec, decimals=1)
    DECvalue2, DECindex2 = np.unique(dec_check, return_index=True)
    #print len(DECvalue2)
    ra_check2 = np.round(s_ra[DECindex2],decimals=2)
    RAvalue2, RAindex2 = np.unique(ra_check2, return_index=True)
    #print len(RAvalue2)

    combineindex = np.append(DECindex1,RAindex2)

    uniqueindex, uniqueindex_index = np.unique(combineindex,return_index=True)

    
    #print DECvalue

    #idxc,idxc2,d2d,d3d = c.search_around_sky(c,0.0002778*u.degree)
    #print idxc[:5]
    #print c[idxc][:5]
    #print c[:5]
    #print idxc2[:5]
    #print c[idxc2][:5]
    
    #idx,d2d,d3d = Coord.match_coordinates_sky(c,c,nthneighbor = 2)
    #doublecut = d2d < 0.00027778*u.degree
    #singlecut = d2d >= 0.00027778*u.degree
    #singleidx = idx[singlecut]
    #print len(singleidx)
    #c2 = c[doublecut]
    #idx2,d2d2,d3d2 = c2.match_to_catalog_sky(c2)
    #print idx2
    #print len(nidx)
    
    #value, index = np.unique(idxc, return_index=True)
    #print value #These are the indices to keep and apply to the input coords
    #print index
    #print len(value)

    #new_c = c[value]
    #print len(new_c)

    #nidx,nd2d,nd3d = Coord.match_coordinates_sky(new_c,new_c,nthneighbor = 2)
    #doublecut2 = nd2d < 0.00027778*u.degree
    #singlecut2 = nd2d >= 0.00027778*u.degree
    #singleidx2 = nidx[singlecut2]
    #print len(singleidx2)
    #nidx2 = nidx[doublecut2]
    #print nidx2
    #print len(nidx2)

    #value2, index2 = np.unique(nidx2, return_index=True)
    #print value2 #These are the indices to keep and apply to the input coords
    #print index
    #print len(value2)
    
    #uniqueindex = np.append(singleidx,value)

    #print len(uniqueindex)
    
    return uniqueindex

################################################################################################
## Add the SDSS objects to the master catalog, BLAH...I should add more of a description
################################################################################################

def sdss_match(instrument, semester, field):

    
    if field == 'egs':
        sdsshdulist = fits.open('/Users/spf/research/sdss/egs/NMBS_DR7_sfillingham.fit')
        nfhdulist = fits.open('/Users/spf/keck_obs/'+instrument+'/'+semester+'/newfirm_mbs_specz_catalog.fits')
        
    elif field == 'cosmos':
        sdsshdulist = fits.open('/Users/spf/research/sdss/cosmos/COSMOS_r22_sfillingham.fit')
        nfhdulist = fits.open('/Users/spf/keck_obs/'+instrument+'/'+semester+'/COSMOS2015_specz_catalog.fits')
        
    else:
        print('choose a real field')
        
    sdssdata = sdsshdulist[1].data
    Rtest = sdssdata['r']
    star_flagtest = sdssdata['type_r']
    mode = sdssdata['mode']

    ###############################################
    ## Establish cuts for stars and galaxies
    ###############################################
    modecut = mode == 1
    Rscicut = (Rtest < 22) & (Rtest > 18)
    Raligncut = (Rtest < 18) & (Rtest > 16)
    Rguidecut = (Rtest < 16) & (Rtest > 10)

    galcut = star_flagtest == 3
    starcut = star_flagtest == 6

    scicut = Rscicut & modecut #& galcut
    aligncut = Raligncut & starcut & modecut
    guidecut = Rguidecut & starcut & modecut

    starcut = aligncut ^ guidecut
    finalcut = starcut ^ scicut 

    ID = sdssdata['ObjID'][finalcut]
    ra = sdssdata['ra'][finalcut]
    dec = sdssdata['dec'][finalcut]
    z = sdssdata['z'][finalcut]
    I = sdssdata['u'][finalcut]
    R = sdssdata['r'][finalcut]
    G = sdssdata['g'][finalcut]
    U = sdssdata['u'][finalcut]
    star_flag = sdssdata['type_r'][finalcut]

    #remove any duplicate objects
    #uniqueindex = doops(ra,dec)
    #ID = ID[uniqueindex]
    #ra = ra[uniqueindex]
    #dec = dec[uniqueindex]
    #z = z[uniqueindex]
    #I = I[uniqueindex]
    #R = R[uniqueindex]
    #G = G[uniqueindex]
    #U = U[uniqueindex]
    #star_flag = star_flag[uniqueindex]
    #print len(ID)

    pflag = np.empty(len(ID))
    #sloanID = np.empty(len(ID))
    sloanID = np.chararray(len(ID),itemsize = 10)

    for i in range(len(ID)):
        sloanID[i] = 's'+np.str(10000+i)
        #sloanID[i] = 50000+i
        
        Rscicut = (R[i] < 22.) and (R[i] > 18.)
        Raligncut = (R[i] < 18.) and (R[i] > 16.)
        Rguidecut = (R[i] < 16.) and (R[i] > 10.)

        galcut = star_flag[i] == 3
        starcut = star_flag[i] == 6

        scicut = Rscicut and galcut
        starscicut = Rscicut and starcut
        aligncut = Raligncut and starcut
        guidecut = Rguidecut and starcut

        if scicut:
            pflag[i] = 1

        elif starscicut:
            pflag[i] = -3

        elif aligncut:
            pflag[i] = -2

        elif guidecut:
            pflag[i] = -1

        else:
            pflag[i] = 0

        
    Ks = np.empty(len(ID))
    Ks[:] = -1
    K = np.empty(len(ID))
    K[:] = -1
    H = np.empty(len(ID))
    H[:] = -1
    H2 = np.empty(len(ID))
    H2[:] = -1
    H1 = np.empty(len(ID))
    H1[:] = -1
    J = np.empty(len(ID))
    J[:] = -1
    J3 = np.empty(len(ID))
    J3[:] = -1
    J2 = np.empty(len(ID))
    J2[:] = -1
    J1 = np.empty(len(ID))
    J1[:] = -1
    mstar = np.empty(len(ID))
    mstar[:] = -1
    redspec = np.empty(len(ID))
    redspec[:] = -1
    redpeak = np.empty(len(ID))
    redpeak[:] = -1
    K_ellip = np.empty(len(ID))
    K_ellip[:] = -1
    K_theta = np.empty(len(ID))
    K_theta[:] = -1
    deepz = np.empty(len(ID))
    deepz[:] = -1
    N_star = np.empty(len(ID))
    N_star[:] = -1

    col1 = fits.Column(name='id',format='K',array=ID)
    col2 = fits.Column(name='ra',format='D',array=ra)
    col3 = fits.Column(name='dec',format='D',array=dec)
    col4 = fits.Column(name='Ks',format='E',array=Ks)
    col5 = fits.Column(name='K',format='E',array=K)
    col6 = fits.Column(name='H',format='E',array=H)
    col7 = fits.Column(name='H2',format='E',array=H2)
    col8 = fits.Column(name='H1',format='E',array=H1)
    col9 = fits.Column(name='J',format='E',array=J)
    col10 = fits.Column(name='J3',format='E',array=J3)
    col11 = fits.Column(name='J2',format='E',array=J2)
    col12 = fits.Column(name='J1',format='E',array=J1)
    col13 = fits.Column(name='z',format='E',array=z)
    col14 = fits.Column(name='I',format='E',array=I)
    col15 = fits.Column(name='R',format='E',array=R)
    col16 = fits.Column(name='G',format='E',array=G)
    col17 = fits.Column(name='U',format='E',array=U)
    col18 = fits.Column(name='mstar',format='E',array=mstar)
    col19 = fits.Column(name='z_spec',format='E',array=redspec)
    col20 = fits.Column(name='z_peak',format='E',array=redpeak)
    col21 = fits.Column(name='type_r',format='I',array=star_flag)
    col22 = fits.Column(name='K_ellip',format='E',array=K_ellip)
    col23 = fits.Column(name='K_theta_J2000',format='E',array=K_theta)
    col24 = fits.Column(name='Near_Star',format='I',array=N_star)
    col25 = fits.Column(name='Pflag',format='I',array=pflag)
    col26 = fits.Column(name='deepz',format='E',array=deepz)
    col27 = fits.Column(name='sID',format='A10',array=sloanID)

    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27])

    tbhdu = fits.BinTableHDU.from_columns(cols)

    nfhdulist.append(tbhdu)

    if field == 'cosmos':
        nfhdulist.writeto(semester+'/COSMOS2015_specz_sdss_catalog.fits')
    else:
        nfhdulist.writeto(semester+'/newfirm_mbs_specz_sdss_catalog.fits')

    
