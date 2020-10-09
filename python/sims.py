from startup import *

import pyathena as pa
units=pa.set_units(muH=1.4271)
import gravity

class data_container(object):
    """Basic data container
    """
    def __init__(self,pid,base='../data/'):
        self.pid=pid
        self.pdir=pid+'/'
        self.base=base
        self._set_filename()
        self.set_units()
        self.load_all()
        self.grav = gravity.gravpot(self.par,self.par['x3max'])

    def _set_filename(self):
        """set file name with assumed path
        """
        pid=self.pid
        pdir=self.pdir
        base=self.base
        self.files=dict(par='{}{}/{}.par'.format(base,pdir,pid),
                        hst='{}{}/hst/{}.hst.p'.format(base,pdir,pid),
                        sn='{}{}/hst/{}.sn.p'.format(base,pdir,pid),
                        hstzp='{}{}/hst/{}.hst_zp.p'.format(base,pdir,pid),
                        zpall=glob.glob('{}{}/zprof_merged/{}.phase?.zprof.nc'.format(base,pdir,pid)),
                        zpicm=glob.glob('{}{}/zprof_merged/{}.phase?-icm.zprof.nc'.format(base,pdir,pid)),
                       )

    def load_all(self):
        """read in history and parameter files
        """
        self._load_hst()
        self._load_par()

    def _load_par(self):
        self.par = pa.get_params(self.files['par'])
        self.dz = self.par['x3max']/self.par['Nx3']*2
        self.Lx = self.par['x1max']-self.par['x1min']
        self.Ly = self.par['x2max']-self.par['x2min']
        self.Lz = self.par['x3max']-self.par['x3min']
        self.area = self.Lx*self.Ly

    def _load_hst(self):
        self.hst = pd.read_pickle(self.files['hst'])
        self.sn = pd.read_pickle(self.files['sn'])
        self.hstzp = pd.read_pickle(self.files['hstzp'])

    def load_zprof(self):
        """read in z-profiles
        """
        print('Merging zprof phases...')
        zpall = self.files['zpall']
        zpall.sort()
        plist=['phase1','phase2','phase3','phase4','phase5']
        dset = None
        for f,ph in zip(zpall,plist):
            with xr.open_dataarray(f) as da:
                ds=da.to_dataset(dim='fields')
                ds=ds.expand_dims('phase')
                ds=ds.assign_coords(phase=[ph])
                if dset is None:
                    dset = ds
                else:
                    dset=xr.concat([dset,ds],dim='phase')

        print('Merging zprof-icm phases...')
        zpicm = self.files['zpicm']
        zpicm.sort()
        dset_icm=None

        for f,ph in zip(zpicm,plist):
            with xr.open_dataarray(f) as da:
                ds=da.to_dataset(dim='fields')
                ds=ds.expand_dims('phase')
                ds=ds.assign_coords(phase=[ph])
                if dset_icm is None:
                    dset_icm = ds
                else:
                    dset_icm=xr.concat([dset_icm,ds],dim='phase')

        # rename phase
        dset = self._redefine_phase(dset)

        # separate species
        if dset_icm is not None:
            dset_icm = self._redefine_phase(dset_icm).assign_coords(taxis=dset.taxis)

            dset_ism = dset - dset_icm
            dset_ism = dset_ism.where(dset_ism['d']>0).fillna(0.)
            self.isicm = True
        else:
            dset_ism = dset
            dset_icm = dset - dset_ism
            self.isicm = False

        return dset,dset_icm,dset_ism

    def set_units(self):
        """set units and a few conversion constants
        """
        units = pa.set_units(muH=1.4271)
        self.to_pok=(units['pressure']/ac.k_B).cgs.value
        self.to_Myr=units['time'].to('Myr').value
        self.to_flux=(units['density']*units['velocity']).to('Msun/(kpc^2*yr)').value
        self.to_Eflux=(units['density']*units['velocity']**3).to('erg/(kpc^2*yr)').value
        self.to_gacc=((units['pressure']/units['length'])/units['density']).cgs.value

    def _redefine_phase(self,data_orig):
        """rename phase
        """
        data_cool=data_orig.sel(phase=['phase1','phase2','phase3']).sum(dim='phase').expand_dims('phase')
        data_cool.coords['phase']=['cool']
        data_cu=data_orig.sel(phase=['phase1','phase2']).sum(dim='phase').expand_dims('phase')
        data_cu.coords['phase']=['cu']
        data_sum=data_orig.sum(dim='phase').expand_dims('phase')
        data_sum.coords['phase']=['all']
        zpdset=xr.concat([data_orig,data_cool,data_cu,data_sum],dim='phase')
        zpdset.coords['phase']=['cold','unstable','warm','int','hot','cool','cu','whole']

        return zpdset

    def calc_icm_params(self):
        """Calculate parameters related to the ICM
        """
        if 'd_icm' in self.par:
            p = self.par
            picm = dict()
            picm['nicm'] = p['d_icm']
            picm['vicm'] = p['v_icm']
            picm['Zicm'] = p['Z_icm']*0.02
            picm['Mdot'] = (p['d_icm']*p['v_icm']*units['density']*units['velocity']).to('Msun/(kpc^2*yr)')
            picm['picm_ram'] = (p['d_icm']*p['v_icm']**2*units['density']*units['velocity']**2/ac.k_B).cgs
            picm['picm_th'] = (p['d_icm']*p['cs_icm']**2*units['density']*units['velocity']**2/ac.k_B).cgs
            picm['picm'] = picm['picm_ram'] + picm['picm_th']
            picm['Sigma_gas0'] = self.hstzp['surf'].iloc[0] * ac.M_sun/ac.pc**2
            picm['Sigma_star'] = p['SurfS'] * ac.M_sun/ac.pc**2
            picm['Wext'] = (2*np.pi*ac.G*picm['Sigma_gas0']*picm['Sigma_star']/ac.k_B).cgs
            picm['ratio'] = picm['picm']/picm['Wext']
            picm['massflux'] = p['d_icm']*p['v_icm']
            picm['metalflux'] = p['d_icm']*p['v_icm']*p['Z_icm']*0.02
            picm['momflux'] = p['d_icm']*(p['v_icm']**2+p['cs_icm']**2)
            picm['energyflux'] = p['d_icm']*p['v_icm']*(p['v_icm']**2+5*p['cs_icm']**2)/2
            self.picm = picm

    def plot_grav(self):
        g=self.grav
        plt.plot(g.z,g.gz,label='total')
        plt.plot(g.z,g.gz_star(g.z),label='star')
        plt.plot(g.z,g.gz_dm(g.z),label='DM')

class allmodels(object):
    """Load all data containers
    """
    def __init__(self,models):
        self.models=models
        self.load_sims()

    def load_sims(self):
        self.dc = dict()
        for m in self.models:
            dc = data_container(self.models[m])
            dc.calc_icm_params()
            tot,icm,ism=dc.load_zprof()
            dc.tot = tot
            dc.icm = icm
            dc.ism = ism
            if m.startswith('P'):
                name = 'ICM-{}'.format(m)
            else:
                name = 'ISM-only'
            dc.name = name
            self.add_fields(dc)
            self.dc[m] = dc

    def print_table(self):
        for m in self.models:
            dc = self.dc[m]
            if hasattr(dc,'picm'):
                p = dc.picm
                print(m,p['nicm']/1.e-4,p['vicm']/1.e3,p['Mdot']/1.e-3,p['picm']/1.e4,p['ratio'])
        print(p['Wext'])

    def add_fields(self,dc):
        for sp in ['tot','icm','ism']:
            zp = getattr(dc,sp)
            zsgn = zp.zaxis/np.abs(zp.zaxis)
            zp['vz'] = zp['M3']/zp['d']
            zp['vout'] = zp['M3']*zsgn/zp['d']
            zp['sicm'] = zp['s4']/zp['d']
            zp['Z'] = zp['s1']/zp['d']
            zp['massflux'] = zp['pFzd']+zp['mFzd']
            zp['Pram'] = zp['pFzM3']+zp['mFzM3']
            zp['pP']=0.4*zp['pFzP']*(zp['pd']/zp['pFzd'])
            zp['mP']=0.4*zp['mFzP']*(zp['md']/zp['mFzd'])
            zp['Pth'] = zp['pP']+zp['mP']
            zp['Pimag'] = zp['PB1']+zp['PB2']-zp['PB3']
            zp['Ptot'] = zp['Pram']+zp['P']+zp['Pimag']
            zp['momflux'] = zp['Pram']+zp['P']#+zp['Pimag']
            zp['energy_eg']=(zp['pFzEge']+zp['mFzEge'])
            zp['energy_sg']=(zp['pFzEgsg']+zp['mFzEgsg'])
            zp['energy_mag']=(zp['pSzEm1']+zp['pSzEm2']+zp['pSzvB1']+zp['pSzvB2'])
            zp['energy_mag']+=(zp['mSzEm1']+zp['mSzEm2']+zp['mSzvB1']+zp['mSzvB2'])
            zp['energy_kx']=(zp['pFzE1']+zp['mFzE1'])
            zp['energy_ky']=(zp['pFzE2']+zp['mFzE2'])
            zp['energy_kz']=(zp['pFzE3']+zp['mFzE3'])
            zp['energy_k']=zp['energy_kx']+zp['energy_ky']+zp['energy_kz']
            zp['energy_th']=(zp['pFzP']+zp['mFzP'])
            zp['energy_z']=zp['energy_kz']+zp['energy_th']
            zp['energyflux']=zp['energy_k']+zp['energy_th']
            zp['metalflux']=(zp['pFzs1']+zp['mFzs1'])
            zp['metalflux_sn']=(zp['pFzs2']+zp['mFzs2']+zp['pFzs3']+zp['mFzs3'])
            zp['vB'] = np.sqrt(2.0*zp['energyflux']/zp['massflux'])

class viz_container(data_container):
    """Data container for visualization
    """
    def __init__(self, pid, base='../data/'):
        data_container.__init__(self,pid,base)

        self.slcfiles = glob.glob('{}{}/slice/{}.*.p'.format(base,pid,pid))
        self.prjfiles = glob.glob('{}{}/surf/{}.*.surf.p'.format(base,pid,pid))
        self.slcfiles.sort()
        self.prjfiles.sort()

    def load_slc(self,i,verbose=True):
        try:
            f = self.slcfiles[i]
        except IndexError:
            print('[Error] idx {} > max idx {}'.format(i,len(self.slcfiles)))
            return

        slc = pd.read_pickle(f)
        if verbose:
            print('reading slices at t = {:.2f}'.format(slc['time']))

        return slc


    def load_prj(self,i,verbose=True):
        try:
            f = self.prjfiles[i]
        except IndexError:
            print('[Error] idx {} > max idx {}'.format(i,len(self.prjfiles)))
            return

        prj = pd.read_pickle(f)
        if verbose:
            print('reading slices at t = {:.2f}'.format(prj['time']))

        return prj

    def get_GG_condition(self):
        """Calculate the critical surface density at which P_ICM=W_GG
        """
        self.calc_icm_params()
        picm = self.picm
        surf_crit = (picm['picm']*ac.k_B/(2*np.pi*ac.G*picm['Sigma_star'])).to('Msun/pc^2')

        return surf_crit
