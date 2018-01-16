from plotting import *
import numpy as np
import scipy.linalg as la
import astropy.coordinates as coord
import astropy.units as u

satnames=['Fornax', 'Sculptor', 'Carina', 'Draco', 'LeoI', 'UrsaMinor', 'Sextans', 'LeoII', 'Sagittarius', 'SMC', 'LMC']
obsprops=np.loadtxt('MW_dwarfs_6D.txt', skiprows=1, usecols=(4,5,9,10,11,12,13,14,15,16))
obspos=obsprops[:,0]
obsposerr=obsprops[:,1]
obsvel=obsprops[:,2]
obsvelerr=obsprops[:,3]
obsj=obsprops[:,4]
obsjerr=obsprops[:,5]
gcnames=np.genfromtxt('allgc_spacemotion.txt', usecols=0, dtype='str')[2:]
gcs=np.loadtxt('allgc_spacemotion.txt', usecols=(1,2,3,4,5,6), skiprows=2)
xs=gcs[:,0]
ys=gcs[:,1]
zs=gcs[:,2]
def vec_mag(x,y,z):
    mag = la.norm([x,y,z])
    return mag
gcsr= map(vec_mag, xs,ys,zs) 
gcsv=vec_mag(gcs[:,3], gcs[:,4], gcs[:,5])
gcsj=[la.norm(np.cross([x,y,z], [vx, vy,vz])) for x, y, z, vx, vy, vz in zip(gcs[:,0], gcs[:,1], gcs[:,2], gcs[:,3], gcs[:,4], gcs[:,5])] 

plt.figure()
for i in range(len(satnames[:-1])):
    plt.errorbar(obspos[i], obsj[i], xerr=obsposerr[i], yerr=obsjerr[i], fmt='o', mec='None', label='%s'%satnames[i], ms=12, elinewidth=2)
plt.errorbar(obspos[-1], obsj[-1], xerr=obsposerr[-1], yerr=obsjerr[-1], fmt='s', c='k', mfc='black', mec='None', label='%s'%satnames[-1], ms=12, elinewidth=2)
plt.plot(gcsr, gcsj, 'o', ms=10, mec='k', mfc='none', label='Sohn GCs')
plt.xlabel(r'$\rm r^{obs}\; [kpc]$', fontsize=14)
plt.ylabel(r'$\rm j^{obs}\; [kpc \; km \; s^{-1}]$', fontsize=14)
plt.ylim(50., 5e4)
plt.yscale("log")
plt.xscale("log")
plt.legend(ncol=2, fontsize=9)
plt.savefig('test_plot_1.pdf')


