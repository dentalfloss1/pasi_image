import h5py
import numpy as np
f = h5py.File('OVRO_LWA-CST_er-3p5_sig-0p002.h5','r')
f.keys()
f['X-pol_Efields']['Etheta(Mag)']
f['X-pol_Efields']
f['X-pol_Efields'].keys()
f['X-pol_Efields']['Ephi(Mag)']
freq =  f['Freq(Hz)']
phi = f['phi_pts']
theta = f['theta_pts']
# do x
xpol = f['X-pol_Efields']
xpolephi  = np.array(xpol['Ephi(Mag)'])
xpolephimag  = np.array(xpol['Ephi(Mag)'])
xpolephiphase  = np.array(xpol['Ephi(Phase)'])
xpolethetamag = np.array(xpol['Etheta(Mag)'])
xpolethetaphase = np.array(xpol['Etheta(Phase)'])
xphi = xpolephimag * np.exp(1j*np.deg2rad(xpolephiphase))
xtheta = xpolethetamag * np.exp(1j*np.deg2rad(xpolethetaphase))
stokesIX = np.real(np.conj(xphi)*xphi + np.conj(xtheta)*xtheta)
stokesQX = np.real(np.conj(xphi)*xphi - np.conj(xtheta)*xtheta)
stokesUX = 2*np.real(xphi*np.conj(xtheta))
stokesVX = -2*np.imag(np.conj(xphi)*xtheta)

# now y
ypol = f['Y-pol_Efields']
ypolephi  = np.array(ypol['Ephi(Mag)'])
ypolephimag  = np.array(ypol['Ephi(Mag)'])
ypolephiphase  = np.array(ypol['Ephi(Phase)'])
ypolethetamag = np.array(ypol['Etheta(Mag)'])
ypolethetaphase = np.array(ypol['Etheta(Phase)'])
yphi = ypolephimag * np.exp(1j*np.deg2rad(ypolephiphase))
ytheta = ypolethetamag * np.exp(1j*np.deg2rad(ypolethetaphase))
stokesIY = np.real(np.conj(yphi)*yphi + np.conj(ytheta)*ytheta)
stokesQY = np.real(np.conj(yphi)*yphi - np.conj(ytheta)*ytheta)
stokesUY = 2*np.real(yphi*np.conj(ytheta))
stokesVY = -2*np.imag(np.conj(yphi)*ytheta)

np.savez('CSTbeam.npz',[stokesIX,stokesQX,stokesUX,stokesVX,stokesIY,stokesQY,stokesUY,stokesVY])

