import numpy as np
import astropy.io.fits as fits
from astropy.table import Table

def get_catalog(semester,field):
    ### This routine pulls the relevant columns from a photometric catalog needed for DEIMOS mask design ###

    ##############################################################################################################
    ## Read in the relevant photometric data files for a given field flag ##
    ##############################################################################################################
    if field == 'egs':
        cat_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.v5.1.cat', format = 'ascii')
        mstar_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout', format = 'ascii')
        z_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.redshifts/aegis-n2.deblend.v5.1.zout', format = 'ascii')
        rf1_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.153-154.rf', format = 'ascii')
        rf2_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.153-155.rf', format = 'ascii')
        rf3_data = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.155-161.rf', format = 'ascii')

        cat_paramlist = np.array(['id','x','y','ra','dec','Ks','K','H','H2','H1','J','J3','J2','J1','z','I','R','G','U','z_spec','star_flag'])
        mstar_paramlist = np.array(['col1','col7']) #'col1' = 'id', 'col7' = 'lmass'
        z_paramlist = np.array(['id','z_spec','z_peak'])

    elif field == 'cosmos':

        cat_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.v5.1.cat', format = 'ascii')
        mstar_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.sps/cosmos-1.bc03.del.deblend.v5.1.fout', format = 'ascii')
        z_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.redshifts/cosmos-1.deblend.v5.1.zout', format = 'ascii')
        rf1_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.rfcolors/cosmos-1.deblend.v5.1.153-154.rf', format = 'ascii')
        rf2_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.rfcolors/cosmos-1.deblend.v5.1.153-155.rf', format = 'ascii')
        rf3_data = Table.read('/Users/spf/newfirmMBS/cosmos/cosmos-1.deblend.rfcolors/cosmos-1.deblend.v5.1.155-161.rf', format = 'ascii')

        cat_paramlist = np.array(['id','x','y','ra','dec','Ks','K','H','H2','H1','J','J3','J2','J1','z','I','R','G','U','z_spec','star_flag'])
        mstar_paramlist = np.array(['col1','col7']) #'col1' = 'id', 'col7' = 'lmass'
        z_paramlist = np.array(['id','z_spec','z_peak'])

    elif field == 'COSMOS2015':
        cat_hdu = fits.open('/Users/spf/cosmos/COSMOS2015/COSMOS2015_Laigle+_v1.1.fits.gz')
        cat_data = cat_hdu[1].data
    else:
        print('choose a real field')

    ##############################################################################################################
    ## Grab the necessary data from the photometric catalogs in order to design the DEIMOS masks ##
    ##############################################################################################################
    if field == 'COSMOS2015':
        
        ID = cat_data.NUMBER
        x = cat_data.X_IMAGE
        y = cat_data.Y_IMAGE
        ra = cat_data.ALPHA_J2000
        dec = cat_data.DELTA_J2000
        Ks = cat_data.Ks_FLUX_APER2
        r = cat_data.r_FLUX_APER2
        H = cat_data.H_FLUX_APER2
        J = cat_data.J_FLUX_APER2
        z = cat_data.zp_FLUX_APER2
        U = cat_data.u_FLUX_APER2
        mstar = cat_data.MASS_BEST
        cosmosflag = cat_data.FLAG_COSMOS
        MU = cat_data.MU
        MV = cat_data.MV
        MJ = cat_data.MJ
        z_phot = cat_data.ZPDF
        star_flag = cat_data.TYPE

        pflag = np.empty(len(ID))
        for i in range(len(pflag)):
            if star_flag[i] == 0:
                pflag[i] = 1
            elif star_flag[i] == 1:
                pflag[i] = -3
            else:
                pflag[i] = -4

        sID = np.empty(len(ID))
        sID[:] = -1
        targ = np.empty(len(ID))
        targ[:] = 0
        spfz = np.empty(len(ID))
        spfz[:] = -1

    else:
        #data from the general catalog
        ID = np.array(cat_data['id'])
        x = np.array(cat_data['x'])
        y = np.array(cat_data['y'])
        ra = np.array(cat_data['ra'])
        dec = np.array(cat_data['dec'])
        Ks = np.array(cat_data['Ks'])
        K = np.array(cat_data['K'])
        H = np.array(cat_data['H'])
        H2 = np.array(cat_data['H2'])
        H1 = np.array(cat_data['H1'])
        J = np.array(cat_data['J'])
        J3 = np.array(cat_data['J3'])
        J2 = np.array(cat_data['J2'])
        J1 = np.array(cat_data['J1'])
        z = np.array(cat_data['z'])
        I = np.array(cat_data['I'])
        R = np.array(cat_data['R'])
        G = np.array(cat_data['G'])
        U = np.array(cat_data['U'])
        z_spec = np.array(cat_data['z_spec'])
        star_flag = np.array(cat_data['star_flag'])
        K_ellip = np.array(cat_data['K_ellip'])
        K_theta = np.array(cat_data['K_theta_J2000'])
        N_star = np.array(cat_data['Near_Star'])
    
        pflag = np.empty(len(ID))
        for i in range(len(pflag)):
            if cat_data['star_flag'][i] < 1:
                pflag[i] = 1
            else:
                pflag[i] = -3

        sID = np.empty(len(ID))
        sID[:] = -1
        targ = np.empty(len(ID))
        targ[:] = 0
        spfz = np.empty(len(ID))
        spfz[:] = -1

        #data from the SPS catalog
        msID = np.array(mstar_data['col1'])
        mstar = np.array(mstar_data['col7'])

        #data from the photometric spec catalog
        redID = np.array(z_data['id'])
        redspec = np.array(z_data['z_spec'])
        redpeak = np.array(z_data['z_peak'])

    if field == 'egs':
        deepz = np.empty(len(ID))
        deepz[:] = -1
    elif field == 'cosmos':
        zcosmos = np.empty(len(ID))
        zcosmos[:] = -1
    elif field == 'COSMOS2015':
        zcosmos = np.empty(len(ID))
        zcosmos[:] = -1
    else:
        print('choose a real field')


    #######################################################
    ## Put the above data into a single FITS file
    #######################################################

    if field == 'COSMOS2015':
        col1 = fits.Column(name='id',format='K',array=ID)
        col2 = fits.Column(name='x_pixel',format='E',array=x)
        col3 = fits.Column(name='y_pixel',format='E',array=y)
        col4 = fits.Column(name='ra',format='D',array=ra)
        col5 = fits.Column(name='dec',format='D',array=dec)
        col6 = fits.Column(name='Ks',format='E',array=Ks)
        col7 = fits.Column(name='r',format='E',array=r)
        col8 = fits.Column(name='H',format='E',array=H)
        col9 = fits.Column(name='J',format='E',array=J)
        col10 = fits.Column(name='z',format='E',array=z)
        col11 = fits.Column(name='U',format='E',array=U)
        col12 = fits.Column(name='mstar',format='E',array=mstar)
        col13 = fits.Column(name='COSMOS_FLAG',format='E',array=cosmosflag)
        col14 = fits.Column(name='MU',format='E',array=MU)
        col15 = fits.Column(name='MV',format='E',array=MV)
        col16 = fits.Column(name='MJ',format='E',array=MJ)
        col17 = fits.Column(name='z_peak',format='E',array=z_phot)
        col18 = fits.Column(name='star_flag',format='I',array=star_flag)
        col19 = fits.Column(name='Pflag',format='I',array=pflag)
        col20 = fits.Column(name='zcosmos',format='E',array=zcosmos)
        col21 = fits.Column(name='sID',format='E',array=sID)
        col22 = fits.Column(name='target',format='I',array=targ)
        col23 = fits.Column(name='spfz',format='E',array=spfz)

        cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,
                             col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,
                             col21,col22,col23])

        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto('/Users/spf/keck_obs/deimos/'+semester+'/COSMOS2015_catalog.fits')

    else:

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
        col21 = fits.Column(name='star_flag',format='I',array=star_flag)
        col22 = fits.Column(name='K_ellip',format='E',array=K_ellip)
        col23 = fits.Column(name='K_theta_J2000',format='E',array=K_theta)
        col24 = fits.Column(name='Near_Star',format='I',array=N_star)
        col25 = fits.Column(name='Pflag',format='I',array=pflag)
        if field == 'egs':
            col26 = fits.Column(name='deepz',format='E',array=deepz)
        elif field == 'cosmos':
            col26 = fits.Column(name='zcosmos',format='E',array=zcosmos)
        else:
            print('choose a real field')
        col27 = fits.Column(name='sID',format='E',array=sID)
        col28 = fits.Column(name='target',format='I',array=targ)
        col29 = fits.Column(name='spfz',format='E',array=spfz)

        cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,
                             col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,
                             col21,col22,col23,col24,col25,col26,col27,col28,col29])

        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto('/Users/spf/keck_obs/deimos/'+semester+'/newfirm_mbs_catalog.fits')
    
