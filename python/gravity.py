import astropy.constants as c
import astropy.units as u
import numpy as np

class gravpot(object):
    def __init__(self,param,zmax):
        self.param = param
        self.surf_s = param['SurfS']*c.M_sun/c.pc**2
        self.rho_dm = param['rhodm']*c.M_sun/c.pc**3
        self.z0 = param['zstar']*c.pc
        self.R0 = param['R0']*c.pc
        self.z = np.arange(0,zmax,1)*u.pc
        self.phi = self.phi_ext(self.z)
        self.gz = self.gz_ext(self.z)

    def gz_ext(self,z):
        return self.gz_star(z)+self.gz_dm(z)

    def phi_ext(self,z):
        return self.phi_star(z)+self.phi_dm(z)

    def gz_star(self,z):
        gz=-2*np.pi*c.G*self.surf_s*z/np.sqrt(z**2+self.z0**2)
        return gz.to('cm/s**2')

    def gz_dm(self,z):
        gz=-4*np.pi*c.G*self.rho_dm*z/(1+(z/self.R0)**2)
        return gz.to('cm/s**2')

    def phi_star(self,z):
        phi=2*np.pi*c.G*(self.surf_s*(np.sqrt(z**2+self.z0**2)-self.z0))
        return phi.to('km**2/s**2')

    def phi_dm(self,z):
        phi=2*np.pi*c.G*self.rho_dm*self.R0**2*np.log(1+z**2/self.R0**2)
        return phi.to('km**2/s**2')

    def vesc(self,z):
        return np.sqrt(2*(self.phi_ext(z).to('km^2/s^2')))

    def phi_sg(self,z,surfg,zg):
        phi=np.log(np.cosh((z/(zg*u.pc)).cgs.value))*2*np.pi*c.G*surfg*c.M_sun/c.pc**2*zg*u.pc
        return phi.to('km**2/s**2')
