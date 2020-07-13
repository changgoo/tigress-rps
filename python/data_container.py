import xarray as xr
import numpy as np
import pandas as pd
from pyathena import hst_reader,get_params
from pyathena import set_units
from pyathena.utils import compare_files

import os
import astropy.constants as ac
import astropy.units as au

unit=set_units(muH=1.4271)
Myr=unit['time'].to('Myr').value

def get_snr(sn,taxis,tbin='auto',snth=100.):
    sntime=sn['time']*Myr
    time=np.array(taxis)*Myr
    if tbin is 'auto':
        tbin=0.
        dtbin=0.1
        snrmean=0.
        while (tbin < 40) & (snrmean < snth):
            tbin += dtbin
            idx=np.less(sntime[np.newaxis,:],time[:,np.newaxis]) & \
            np.greater(sntime[np.newaxis,:],(time[:,np.newaxis]-tbin))
            snr=idx.sum(axis=1)
            snrmean=snr.mean()

        snr = snr/tbin
    else:
        idx=np.less(sntime[np.newaxis,:],time[:,np.newaxis]) & \
        np.greater(sntime[np.newaxis,:],(time[:,np.newaxis]-tbin))
        snr=idx.sum(axis=1)/tbin
    snr=xr.DataArray(snr,coords=[taxis])
    return tbin,time,snr

def get_flux_norms(Sigma_snr,ESN=1.e51,mstar=95.5,Mej=10,ZSN=0.2,vcool=200):
    momunit=(1.0*(au.erg/(au.km/au.s))/au.kpc**2/au.yr).to('Msun*km/(s*kpc^2*yr)').value
    eunit=(1.0*au.erg/au.kpc**2/au.yr).to('erg/(yr*kpc^2)').value
    flux_norms={'mass':Sigma_snr*mstar,
                'mom':ESN/vcool*Sigma_snr*momunit/2.,
                'mom_kin':ESN/vcool*Sigma_snr*momunit/2.,
                'mom_th':ESN/vcool*Sigma_snr*momunit/2.,
                'energy':ESN*Sigma_snr*eunit,
                'metal':Mej*ZSN*Sigma_snr,
                'metal_sn':Mej*ZSN*Sigma_snr,
               }
    return flux_norms

#--------------------------------------------------------------------------------------------------
# Flux data class
#--------------------------------------------------------------------------------------------------
class data_container(object):
    '''
        A data container for RPS analysis
    '''
    def __init__(self,pid,pdir=None,base='../',create=True,read_only=False,
                 fluxes=['totflux','netflux','origflux','outflux','influx','totinflux']):
        if pdir is None: pdir=pid + '/'
        self.base=base
        self.pid=pid
        self.pdir=pdir
        self.fluxes=fluxes

        snfname='{}{}hst/{}.sn'.format(base,pdir,pid)
        if os.path.isfile(snfname): self.sn=hst_reader(snfname)
        else: self.sn=pd.read_pickle(snfname+'.p')

        hstfname='{}{}hst/{}.hst'.format(base,pdir,pid)
        if os.path.isfile(hstfname): self.hst=hst_reader(hstfname)
        else: self.hst=pd.read_pickle(hstfname+'.p')

        hstzpfname='{}{}hst/{}.hst_zp.p'.format(base,pdir,pid)
        self.hstzp=pd.read_pickle(hstzpfname)

        parfname='{}{}{}.par'.format(base,pdir,pid)
        self.par=get_params(parfname)
        if 'Omega' in self.par: self.torb=2*np.pi/self.par['Omega']

        self.par=get_params(parfname)
        if 'Omega' in self.par: self.torb=2*np.pi/self.par['Omega']

        if 'd_icm' in self.par:
            self.nicm=self.par['d_icm']
            self.vicm=self.par['v_icm']
            self.csicm=self.par['cs_icm']
            self.Zicm=self.par['Z_icm']

        if read_only: return

        self.set_units()
        self.initial_processing()
        datadir='{}{}/flux/'.format(self.base,self.pid)
        fluxfile='{}{}_totflux.nc'.format(datadir,self.pid)
        if create:# or (not compare_files('./data_container.py',fluxfile)):
            dset,dset_icm,dset_ism=self.merge_zprof()
            #if self.Precon: self.calculate_Pterms(dset)
            self.calculate_flux(dset)
            self.calculate_norm()
            self.create_flux_data()
        else:
            self.read_flux_data()
            self.calculate_norm()
            self.add_ZISM()
            self.calculate_loading()
            self.calculate_vout()
        self.colors={'cold':'C0','unstable':'C6','warm':'C2',
                     'int':'#a6d608','hot':'C3','whole':'C7','cool':'C0','cooint':'#32cd32',
                     'i':'#a6d608','h':'C3'
                    }

    def load_vpdf(self):
        base=self.base
        pid=self.pid
        if not hasattr(self,'vpdf'):
            fname='{}{}/vz_pdf/{}_tavg.nc'.format(base,pid,pid)
            if os.path.isfile(fname): self.vpdf=xr.open_dataset(fname)
            else: print('{} is not found'.format(fname))
            self.vpdf['mvpdf']=self.vpdf.v*self.vpdf['mpdf']
            self.vpdf['mvpdf_H']=self.vpdf.v*self.vpdf['mpdf_H']
            self.vpdf['mvpdf_2H']=self.vpdf.v*self.vpdf['mpdf_2H']
        if not hasattr(self,'vstd'):
            fname='{}{}/vz_pdf/{}_tstd.nc'.format(base,pid,pid)
            if os.path.isfile(fname): self.vstd=xr.open_dataset(fname)
            else: print('{} is not found'.format(fname))
            self.vstd['mvpdf']=np.abs(self.vstd.v)*self.vstd['mpdf']
            self.vstd['mvpdf_H']=np.abs(self.vstd.v)*self.vstd['mpdf_H']
            self.vstd['mvpdf_2H']=np.abs(self.vstd.v)*self.vstd['mpdf_2H']
        if not hasattr(self,'vBpdf'):
            fname='{}{}/vz_pdf/{}_tavg.vB.nc'.format(base,pid,pid)
            if os.path.isfile(fname): self.vBpdf=xr.open_dataset(fname)
            else: print('{} is not found'.format(fname))
        if not hasattr(self,'vBstd'):
            fname='{}{}/vz_pdf/{}_tstd.vB.nc'.format(base,pid,pid)
            if os.path.isfile(fname): self.vBstd=xr.open_dataset(fname)
            else: print('{} is not found'.format(fname))

    def add_ZISM(self,z0=50.):
        for flux in self.fluxes:
            fl=getattr(self,flux)
            flmid=fl.sel(zaxis=slice(-z0,z0))
            fl['Z']=(fl['metal']/fl['mass'])
            fl['ZISM']=(flmid['metal']/flmid['mass']).mean(dim='zaxis').sel(phase='cool')
        return

    def calculate_loading(self):
        for flux in self.fluxes:
            fl=getattr(self,flux)
            for var in ['mass','mom','energy','metal']:
                fl['{}_loading'.format(var)]=(fl[var]/fl['{}_norm'.format(var)]).where(fl['{}_norm'.format(var)]>0.)
            fl['mom_kin_loading'.format(var)]=(fl['mom_kin']/fl['mom_norm']).where(fl['mom_norm']>0.)
            fl['fmass_sn_cum']=((fl['Z']-0.02)/(self.ZSN-0.02))
            fl['fmass_sn_cum']=fl['fmass_sn_cum'].clip(0,None).fillna(0)
            fl['fmetal_sn_cum']=self.ZSN/fl['Z']*fl['fmass_sn_cum']
            fl['mass_sn_cl']=fl['fmass_sn_cum']*fl['mass']
            fl['metal_sn_cl']=fl['metal_sn']
            fl['mass_sn_cl_loading']=(fl['mass_sn_cl']/fl['mass_norm']).where(fl['mass_norm']>0.)
            fl['metal_sn_cl_loading']=(fl['metal_sn_cl']/fl['metal_norm']).where(fl['metal_norm']>0.)

            fl['fmass_sn']=((fl['Z']-fl['ZISM'])/(self.ZSN-fl['ZISM']))
            fl['fmass_sn']=fl['fmass_sn'].clip(0,None).fillna(0)
            fl['fmetal_sn']=self.ZSN/fl['Z']*fl['fmass_sn']
            fl['mass_sn']=fl['fmass_sn']*fl['mass']
            fl['metal_sn']=self.ZSN*fl['mass_sn']
            fl['mass_sn_loading']=(fl['mass_sn']/fl['mass_norm']).where(fl['mass_norm']>0.)
            fl['metal_sn_loading']=(fl['metal_sn']/fl['metal_norm']).where(fl['metal_norm']>0.)
            fl['enrichment']=(fl['Z']/fl['ZISM']).where(fl['ZISM']>0.)

    def calculate_vout(self):
        units=self.flux_units
        for flux in self.fluxes:
            fl=getattr(self,flux)
            fl['vout']=(fl['mass']/fl['d'])
            fl['vout_flux']=(fl['mom_kin']/fl['mass'])
            pot_energy=fl['energy_eg']/fl['mass']*units['mass']/units['energy']
            specific_energy=fl['energy']/fl['mass']*units['mass']/units['energy']
            fl['vB']=np.sqrt(2.0*specific_energy.where(specific_energy>=0))
            fl['vesc']=np.sqrt(2.0*pot_energy.where(pot_energy>=0))
            fl['Z']=(fl['metal']/fl['mass'])

    def create_flux_data(self):
        dc=self
        datadir='{}{}/flux/'.format(self.base,self.pid)
        if not os.path.isdir(datadir): os.mkdir(datadir)

        print('Creating flux files...')
        for flux in self.fluxes:
            fl = getattr(dc,flux)
            for field in dc.flux_norms:
                fl['{}_norm'.format(field)]=dc.flux_norms[field]
            f='{}{}_{}.nc'.format(datadir,self.pid,flux)
            fl.to_netcdf(f,mode='w')
            fl.close()

    def read_flux_data(self):
        datadir='{}{}/flux/'.format(self.base,self.pid)
        for flux in self.fluxes:
            f='{}{}_{}.nc'.format(datadir,self.pid,flux)
            with xr.open_dataset(f) as ds:
                setattr(self,flux,ds)
        self.taxis=getattr(self,flux).taxis
        self.H['taxis']=self.taxis
        self.H_2p['taxis']=self.taxis
        self.calculate_norm()
        for flux in self.fluxes:
            fl = getattr(self,flux)
            for field in self.flux_norms:
                fl['{}_norm'.format(field)]=self.flux_norms[field]
        print('Reading flux files')

    def merge_zprof(self):
        print('Merging zprof phases...')
        base=self.base
        pdir=self.pdir
        pid=self.pid
        dset=None
        for i,ph in enumerate(['phase1','phase2','phase3','phase4','phase5']):
            fname='{}{}zprof_merged/{}.phase{}.zprof.nc'.format(base,pdir,pid,i+1)
            with xr.open_dataarray(fname) as da:
               ds=da.to_dataset(dim='fields')
               ds=ds.expand_dims('phase')
               ds=ds.assign_coords(phase=[ph])
               if dset is None:
                   dset = ds
               else:
                   dset=xr.concat([dset,ds],dim='phase')

        print('Merging zprof-icm phases...')
        dset_icm=None
        for i,ph in enumerate(['phase1','phase2','phase3','phase4','phase5']):
            fname='{}{}zprof_merged/{}.phase{}-icm.zprof.nc'.format(base,pdir,pid,i+1)
            with xr.open_dataarray(fname) as da:
               ds=da.to_dataset(dim='fields')
               ds=ds.expand_dims('phase')
               ds=ds.assign_coords(phase=[ph])
               if dset_icm is None:
                   dset_icm = ds
               else:
                   dset_icm=xr.concat([dset_icm,ds],dim='phase')

        dset=self.redefine_phase(dset)
        dset_icm=self.redefine_phase(dset_icm).assign_coords(taxis=dset.taxis)

        dset_ism = dset - dset_icm
        dset_ism = dset_ism.where(dset_ism['d']>0).fillna(0.)

        return dset,dset_icm,dset_ism

    def redefine_phase(self,data_orig):
        data_cool=data_orig.sel(phase=['phase1','phase2','phase3']).sum(dim='phase').expand_dims('phase')
        data_cool.coords['phase']=['cool']
        data_sum=data_orig.sum(dim='phase').expand_dims('phase')
        data_sum.coords['phase']=['all']
        zpdset=xr.concat([data_orig,data_cool,data_sum],dim='phase')
        zpdset.coords['phase']=['cold','unstable','warm','int','hot','cool','whole']

        return zpdset

    def set_units(self):
        from pyathena import set_units
        import astropy.constants as c
        import astropy.units as u

        unit=set_units(muH=1.4271)

        unit['G_code'] = (c.G/unit['length']**3)*unit['mass']*unit['time']**2

        self.to_Myr=unit['time'].to('Myr').value
        self.to_Msun=unit['mass'].to('Msun').value
        self.unit=unit
        self.ZSN=0.2

    def initial_processing(self):
        par=self.par
        Lx=self.par['x1max']-self.par['x1min']
        Ly=self.par['x2max']-self.par['x2min']
        area=Lx*Ly

        unit=self.unit
        flux_units={'d':(unit['density']).to('Msun*s/(kpc^2*yr*km)').value,
                    'mass':(unit['density']*unit['velocity']).to('Msun/(kpc^2*yr)').value,
                    'mom':(unit['density']*unit['velocity']**2).to('Msun*km/(s*kpc^2*yr)').value,
                    'energy':(unit['density']*unit['velocity']**3).to('erg/(kpc^2*yr)').value,
                    'metal':(unit['density']*unit['velocity']).to('Msun/(kpc^2*yr)').value,
                   }

        self.flux_units=flux_units
        self.area=area
        self.H=self.hstzp['H'].to_xarray()
        self.H_2p=self.hstzp['H_2p'].to_xarray()

    def calculate_norm(self):
        area = self.area
        tbin0,time,snr=get_snr(self.sn,self.taxis,tbin='auto',snth=10.)
        self.snr=snr
        self.Sigma_snr=snr/area
        self.flux_norms=get_flux_norms(self.Sigma_snr)

    def calculate_Pterms(self,data_orig):
        fM=data_orig['pFzd']
        fp=data_orig['pFzM3']
        fEth=data_orig['pFzP']

        fpP=0.4*fEth*fM/fp

        data_orig['pP']=fpP

        fM=data_orig['mFzd']
        fp=data_orig['mFzM3']
        fEth=data_orig['mFzP']

        fmP=0.4*fEth*fM/fp

        data_orig['mP']=fmP

    def calculate_flux(self,zpdset):
        area = self.area
        zhat=zpdset.zaxis/abs(zpdset.zaxis)
        upper=zpdset.where(zpdset.zaxis > 0,drop=True)
        lower=zpdset.where(zpdset.zaxis < 0,drop=True)

        flux=xr.Dataset()
        outflux=xr.Dataset()
        outflux_u=xr.Dataset()
        outflux_l=xr.Dataset()
        influx=xr.Dataset()
        influx_u=xr.Dataset()
        influx_l=xr.Dataset()
        flux['d']=zpdset['pd']*zhat-zpdset['md']*zhat
        flux['mass']=(zpdset['pFzd']+zpdset['mFzd'])
        flux['mom_kin']=zpdset['pFzM3']*zhat-zpdset['mFzM3']*zhat
        flux['mom_th']=zpdset['pP']*zhat-zpdset['mP']*zhat
        flux['mom']=flux['mom_kin']+flux['mom_th']
        flux['energy_eg']=(zpdset['pFzEge']+zpdset['mFzEge'])
        flux['energy_sg']=(zpdset['pFzEgsg']+zpdset['mFzEgsg'])
        flux['energy_mag']=(zpdset['pSzEm1']+zpdset['pSzEm2']+zpdset['pSzvB1']+zpdset['pSzvB2'])
        flux['energy_mag']+=(zpdset['mSzEm1']+zpdset['mSzEm2']+zpdset['mSzvB1']+zpdset['mSzvB2'])
        flux['energy_kx']=(zpdset['pFzE1']+zpdset['mFzE1'])
        flux['energy_ky']=(zpdset['pFzE2']+zpdset['mFzE2'])
        flux['energy_kz']=(zpdset['pFzE3']+zpdset['mFzE3'])
        flux['energy_k']=flux['energy_kx']+flux['energy_ky']+flux['energy_kz']
        flux['energy_th']=(zpdset['pFzP']+zpdset['mFzP'])
        flux['energy_z']=flux['energy_kz']+flux['energy_th']
        flux['energy']=flux['energy_k']+flux['energy_th']
        flux['metal']=(zpdset['pFzs1']+zpdset['mFzs1'])
        flux['metal_sn']=(zpdset['pFzs2']+zpdset['mFzs2']+zpdset['pFzs3']+zpdset['mFzs3'])

        outflux_u['d']=upper['pd']
        outflux_l['d']=lower['md'].assign_coords(zaxis=-lower.zaxis)
        outflux_u['mass']=upper['pFzd']
        outflux_l['mass']=-lower['mFzd'].assign_coords(zaxis=-lower.zaxis)

        outflux_u['mom_kin']=upper['pFzM3']
        outflux_l['mom_kin']=lower['mFzM3'].assign_coords(zaxis=-lower.zaxis)
        outflux_u['mom_th']=upper['pP']
        outflux_l['mom_th']=lower['mP'].assign_coords(zaxis=-lower.zaxis)
        outflux_u['mom']=outflux_u['mom_kin']+outflux_u['mom_th']
        outflux_l['mom']=outflux_l['mom_kin']+outflux_l['mom_th']

        outflux_u['energy_eg']=(upper['pFzEge'])
        outflux_l['energy_eg']=-(lower['mFzEge']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_sg']=(upper['pFzEgsg'])
        outflux_l['energy_sg']=-(lower['mFzEgsg']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_mag']=(upper['pSzEm1']+upper['pSzEm2']+upper['pSzvB1']+upper['pSzvB2'])
        outflux_l['energy_mag']=-(lower['mSzEm1']+lower['mSzEm2']+lower['mSzvB1']+lower['mSzvB2']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_kx']=(upper['pFzE1'])
        outflux_l['energy_kx']=-(lower['mFzE1']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_ky']=(upper['pFzE2'])
        outflux_l['energy_ky']=-(lower['mFzE2']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_kz']=(upper['pFzE3'])
        outflux_l['energy_kz']=-(lower['mFzE3']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_th']=(upper['pFzP'])
        outflux_l['energy_th']=-(lower['mFzP']).assign_coords(zaxis=-lower.zaxis)
        outflux_u['energy_k']=outflux_u['energy_kx']+outflux_u['energy_ky']+outflux_u['energy_kz']
        outflux_l['energy_k']=outflux_l['energy_kx']+outflux_l['energy_ky']+outflux_l['energy_kz']
        outflux_u['energy_z']=outflux_u['energy_kz']+outflux_u['energy_th']
        outflux_l['energy_z']=outflux_l['energy_kz']+outflux_l['energy_th']
        outflux_u['energy']=outflux_u['energy_k']+outflux_u['energy_th']
        outflux_l['energy']=outflux_l['energy_k']+outflux_l['energy_th']

        outflux_u['metal']=upper['pFzs1']
        outflux_l['metal']=-lower['mFzs1'].assign_coords(zaxis=-lower.zaxis)
        outflux_u['metal_sn']=upper['pFzs2']+upper['pFzs3']
        outflux_l['metal_sn']=-(lower['mFzs2']+lower['mFzs3']).assign_coords(zaxis=-lower.zaxis)

        influx_u['d']=upper['md']
        influx_l['d']=lower['pd'].assign_coords(zaxis=-lower.zaxis)
        influx_u['mass']=upper['mFzd']
        influx_l['mass']=-lower['pFzd'].assign_coords(zaxis=-lower.zaxis)

        influx_u['mom_kin']=upper['mFzM3']
        influx_l['mom_kin']=lower['pFzM3'].assign_coords(zaxis=-lower.zaxis)
        influx_u['mom_th']=upper['mP']
        influx_l['mom_th']=lower['pP'].assign_coords(zaxis=-lower.zaxis)
        influx_u['mom']=influx_u['mom_kin']+influx_u['mom_th']
        influx_l['mom']=influx_l['mom_kin']+influx_l['mom_th']

        influx_u['energy_eg']=(upper['mFzEge'])
        influx_l['energy_eg']=-(lower['pFzEge']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_sg']=(upper['mFzEgsg'])
        influx_l['energy_sg']=-(lower['pFzEgsg']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_mag']=(upper['mSzEm1']+upper['mSzEm2']+upper['mSzvB1']+upper['mSzvB2'])
        influx_l['energy_mag']=-(lower['pSzEm1']+lower['pSzEm2']+lower['pSzvB1']+lower['pSzvB2']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_kx']=(upper['mFzE1'])
        influx_l['energy_kx']=-(lower['pFzE1']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_ky']=(upper['mFzE2'])
        influx_l['energy_ky']=-(lower['pFzE2']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_kz']=(upper['mFzE3'])
        influx_l['energy_kz']=-(lower['pFzE3']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_th']=(upper['mFzP'])
        influx_l['energy_th']=-(lower['pFzP']).assign_coords(zaxis=-lower.zaxis)
        influx_u['energy_k']=influx_u['energy_kx']+influx_u['energy_ky']+influx_u['energy_kz']
        influx_l['energy_k']=influx_l['energy_kx']+influx_l['energy_ky']+influx_l['energy_kz']
        influx_u['energy_z']=influx_u['energy_kz']+influx_u['energy_th']
        influx_l['energy_z']=influx_l['energy_kz']+influx_l['energy_th']
        influx_u['energy']=influx_u['energy_k']+influx_u['energy_th']
        influx_l['energy']=influx_l['energy_k']+influx_l['energy_th']

        influx_u['metal']=upper['mFzs1']
        influx_l['metal']=-lower['pFzs1'].assign_coords(zaxis=-lower.zaxis)
        influx_u['metal_sn']=upper['mFzs2']+upper['mFzs3']
        influx_l['metal_sn']=-(lower['pFzs2']+lower['pFzs3']).assign_coords(zaxis=-lower.zaxis)

        flux_units=self.flux_units
        for k in flux:
            unit_key = k.split('_')[0]
            if (k.startswith('mom')) or (k=='d'):
                flux[k] *= flux_units[unit_key]/area
            else:
                flux[k] *= flux_units[unit_key]*zhat/area
            outflux_u[k] *= flux_units[unit_key]/area
            outflux_l[k] *= flux_units[unit_key]/area
            influx_u[k] *= flux_units[unit_key]/area
            influx_l[k] *= flux_units[unit_key]/area
        netflux=flux.assign_coords(zabs=abs(flux.zaxis)).groupby('zabs').sum(dim='zaxis')
        netflux=netflux.rename(zabs='zaxis')

        outflux=xr.concat([outflux_l.assign_coords(zaxis=-outflux_l.zaxis),outflux_u],dim='zaxis')
        influx=xr.concat([influx_l.assign_coords(zaxis=-influx_l.zaxis),influx_u],dim='zaxis')

        # flux at full domain
        self.flux=flux
        self.outflux=outflux
        self.influx=influx
        self.origflux=flux*zhat
        self.origflux['fV']=zpdset['A']/self.area
        self.taxis=self.origflux.taxis
        self.zaxis=self.origflux.zaxis

        # flux at half domain
        self.netflux=netflux
        self.totflux=outflux_u+outflux_l
        self.totinflux=influx_u+influx_l

#--------------------------------------------------------------------------------------------------
# Joint PDF data class
#--------------------------------------------------------------------------------------------------

class pdf_data(object):
    def __init__(self,pid,base='../data/',dc=None):
        self.pid=pid
        self.base=base
        if dc is None: dc=data_container(pid,base=base,create=False,fluxes=['totflux'])
        self.sfr=dc.hstzp['sfr10_hst'].to_xarray().sel(taxis=slice(0.5*dc.torb,1.5*dc.torb)).mean().data
        self.ZISM=dc.totflux['ZISM'].sel(taxis=slice(0.5*dc.torb,1.5*dc.torb)).mean().data
        self.zlist=['500','1000','H','2H']
        self.Nx=dc.par['Nx1']
        self.Ny=dc.par['Nx2']
        self.Nt=len(dc.totflux['mass'].sel(taxis=slice(dc.torb*0.5,dc.torb*1.5)).taxis)
        units=dict()
        units['massflux']=(unit['density']*unit['velocity']).to('Msun/(kpc^2*yr)').value
        units['metalflux']=(unit['density']*unit['velocity']).to('Msun/(kpc^2*yr)').value
        units['energyflux']=(unit['density']*unit['velocity']**3).to('erg/(kpc^2*yr)').value
        units['momflux']=(unit['density']*unit['velocity']**2).to('(Msun*km)/(kpc^2*yr*s)').value
        self.units = units

    def Z_model(self,vB):

        yZcool=10.**(np.max([np.log10(1.15*(self.sfr)**0.05),np.log10(1.)],axis=0))
        vmax=3.2e3
        ZSN=0.2
        expo=1.7

        self.yZcool=yZcool

        return (vB/vmax)**expo*(ZSN-self.ZISM)+self.ZISM*yZcool

    def vB_projections(self,z,dvbin=0.02,wf=None):
        pdf=self.pdfdict[z]
        dbin=pdf.attrs['dbin']
        vbin=np.log10(pdf['vBz'])
        vB_bins=np.arange(0,4,dvbin)
        if wf is None:
            pdf1d_vB=pdf.groupby_bins(vbin,vB_bins).sum()*dbin
            cumpdf_vB=pdf1d_vB.to_array().data[:,::-1].cumsum(axis=1)[:,::-1]*dbin
            cumpdf_vB=xr.DataArray(cumpdf_vB,coords=pdf1d_vB.to_array().coords).to_dataset('variable')
        else:
            wpdf=pdf[wf].groupby_bins(vbin,vB_bins).sum()
            pdf1d_vB=(pdf*pdf[wf]).groupby_bins(vbin,vB_bins).sum()/wpdf
            cumpdf_vB=(pdf1d_vB.to_array().data[:,::-1].cumsum(axis=1)/(np.arange(len(vB_bins)-1)+1)[np.newaxis,:])[:,::-1]
            cumpdf_vB=xr.DataArray(cumpdf_vB,coords=pdf1d_vB.to_array().coords).to_dataset('variable')
        pdf1d_vB.coords['vBz_bins'] = [v.left for v in pdf1d_vB['vBz_bins'].values]
        cumpdf_vB.coords['vBz_bins'] = [v.left for v in cumpdf_vB['vBz_bins'].values]

        return pdf1d_vB,cumpdf_vB

    def load_all(self):
        self.pdfdict=dict()

        for z0 in self.zlist:
            pdf,pdfe=self._pdf_load(z0)
            for k in pdfe:
                if k in pdf:
                    pass
                else:
                    pdf[k]=pdfe[k]
                    pdf.attrs[k]=pdfe.attrs[k]
            self._pdf_add_terms(pdf)
            self.pdfdict[z0]=pdf

    def load_z(self,z0):
        pdf,pdfe=self._pdf_load(z0)
        pdf['momflux']=pdfe['momflux']
        pdf.attrs['momflux']=pdfe.attrs['momflux']
        pdf.attrs['sfr']=self.sfr
        pdf.attrs['ZISM']=self.ZISM
        pdf.attrs['NxNyNt']=self.Nx*self.Ny*self.Nt
        for f,u in self.units.items():
            pdf.attrs[f+'_unit']=u
        vout=10.**pdf.vout
        cs=10.**pdf.cs
        pdf['vBz'] = np.sqrt(5*(cs)**2+(vout)**2)

        return pdf

    def _pdf_add_terms(self,pdf):
        # total energy flux
        # energyflux + poyntingflux
        pdf.attrs['etotflux'] = pdf.attrs['energyflux']+pdf.attrs['poyntingflux']
        pdf['etotflux'] = pdf['energyflux']*pdf.attrs['energyflux']+pdf['poyntingflux']*pdf.attrs['poyntingflux']
        pdf['etotflux'] = pdf['etotflux']/pdf.attrs['etotflux']

        # metalicity
        pdf.attrs['Z_MF'] = pdf.attrs['metalflux']/pdf.attrs['massflux']
        pdf.attrs['Z_V'] = pdf.attrs['Z']/pdf.attrs['volume']
        pdf.attrs['yZ_V'] = pdf.attrs['yZ']/pdf.attrs['volume']

        pdf['Z_MF'] = pdf['metalflux']/pdf['massflux']*pdf.attrs['Z_MF']
        pdf['Z_V']  = pdf['Z']/pdf['volume']*pdf.attrs['Z_V']
        pdf['yZ_V'] = pdf['yZ']/pdf['volume']*pdf.attrs['yZ_V']
        pdf['yZ_MF'] = pdf['Z_MF']/self.ZISM

        # ratio
        pdf.attrs['ratio-mommag-momflux'] = np.abs(pdf.attrs['momflux_mag'])/np.abs(pdf.attrs['momflux'])
        pdf.attrs['ratio-poynting-etotflux'] = pdf.attrs['poyntingflux']/pdf.attrs['etotflux']
        pdf.attrs['ratio-ekinz-ekin'] = pdf.attrs['energyflux_kin_z']/pdf.attrs['energyflux_kin']
        pdf.attrs['ratio-ekinh-ekin'] = 1-pdf.attrs['energyflux_kin_z']/pdf.attrs['energyflux_kin']
        pdf.attrs['ratio-eth-ekin'] = pdf.attrs['energyflux_th']/pdf.attrs['energyflux_kin']

        pdf['ratio-mommag-momflux'] = np.abs(pdf['momflux_mag'])/np.abs(pdf['momflux'])*pdf.attrs['ratio-mommag-momflux']
        pdf['ratio-poynting-etotflux'] = pdf['poyntingflux']/pdf['etotflux']*pdf.attrs['ratio-poynting-etotflux']
        pdf['ratio-ekinz-ekin'] = pdf['energyflux_kin_z']/pdf['energyflux_kin']*pdf.attrs['ratio-ekinz-ekin']
        pdf['ratio-ekinh-ekin'] = 1-pdf['energyflux_kin_z']/pdf['energyflux_kin']*pdf.attrs['ratio-ekinz-ekin']
        pdf['ratio-eth-ekin'] = pdf['energyflux_th']/pdf['energyflux_kin']*pdf.attrs['ratio-eth-ekin']

        # vB
        vout=10.**pdf.vout
        cs=10.**pdf.cs
        pdf.attrs['vB'] = np.sqrt(2.0*pdf.attrs['energyflux']/pdf.attrs['massflux'])
        pdf['vB'] = np.sqrt(pdf['energyflux']/pdf['massflux'])*pdf.attrs['vB']
        pdf['vBz'] = np.sqrt(5*(cs)**2+(vout)**2)
        pdf['ratio-vBz-vB'] = (pdf['vBz']/pdf['vB'])**2

        # yZ_model
        pdf['Z_model_vB'] = self.Z_model(pdf['vB'])
        pdf['Z_model_vBz'] = self.Z_model(pdf['vBz'])
        pdf['ratio-Z-vB'] = pdf['Z_model_vB']/pdf['Z_MF']
        pdf['ratio-Z-vBz'] = pdf['Z_model_vBz']/pdf['Z_MF']

        # derived fluxes
        bias=0.1*np.log10(pdf['vBz'])+0.6
        pdf['metalflux_from_massflux'] = pdf['massflux']*pdf['Z_model_vBz']/pdf.attrs['Z_MF']
        pdf['momflux_from_massflux'] = pdf['massflux']*(vout+cs**2/vout)*(pdf.attrs['massflux']/pdf.attrs['momflux'])
        pdf['etotflux_from_massflux'] = pdf['massflux']*0.5*pdf['vBz']**2*(pdf.attrs['massflux']/pdf.attrs['etotflux'])/bias
        pdf['energyflux_from_massflux'] = pdf['massflux']*0.5*pdf['vBz']**2*(pdf.attrs['massflux']/pdf.attrs['energyflux'])/bias
        pdf['ratio-metalflux'] = pdf['metalflux_from_massflux']/pdf['metalflux']
        pdf['ratio-momflux'] = pdf['momflux_from_massflux']/pdf['momflux']
        pdf['ratio-energyflux'] = pdf['energyflux_from_massflux']/pdf['energyflux']
        pdf['ratio-etotflux'] = pdf['energyflux_from_massflux']/pdf['etotflux']

        pdf.attrs['ratio-metalflux'] = ((pdf['ratio-metalflux']*pdf['massflux']).sum()/pdf['massflux'].sum()).data
        pdf.attrs['ratio-momflux'] = ((pdf['ratio-momflux']*pdf['massflux']).sum()/pdf['massflux'].sum()).data
        pdf.attrs['ratio-energyflux'] = ((pdf['ratio-energyflux']*pdf['massflux']).sum()/pdf['massflux'].sum()).data
        pdf.attrs['ratio-etotflux'] = ((pdf['ratio-etotflux']*pdf['massflux']).sum()/pdf['massflux'].sum()).data

    def _pdf_load(self,z):
        base=self.base
        pid=self.pid
        with xr.open_dataset('{base}{pid}/vz_pdf/{pid}.pdf-out.{z}.nc'.format(**dict(base=base,pid=pid,z=z))) as ds:
            pdfout = ds.load()
        with xr.open_dataset('{base}{pid}/vz_pdf/{pid}.pdf-out-extra.{z}.nc'.format(**dict(base=base,pid=pid,z=z))) as ds:
            pdfextra = ds.load()
        return pdfout,pdfextra
