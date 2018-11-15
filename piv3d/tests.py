

import numpy as np
import matplotlib.pyplot as plt


filename='/media/stephane/DATA/Experimental_data/Turbulence3d/20181010/raw_flowfield.npy'

a = np.load(filename)

n = np.prod(a.shape[1:])
N = 1000
nanratio = [np.sum(np.isnan(a[i,...].flatten()))/n for i in range(N)]

#truncate to keep only the 60000
N = 60000
b = a[:N,...]
(N,Nx,Ny,Nd) = b.shape
Nz = 40
Nt = int(N/Nz)
b = b.reshape((Nt,Nz,Nx,Ny,2))

print(b.shape)

umean = np.nanmean(b,axis=0)
urms = np.nanstd(b,axis=0)

#display

for i in range(21,0,-1):
	fig,axs=plt.subplots(2,1,figsize=(12,8)); axs=axs.flatten()
	axs[0].imshow(umean[0:20,i,:,1])
	axs[1].imshow(urms[0:20:,i,:,1])
	plt.show()
	plt.pause(0.2)

for i in range(20):
	fig,axs=plt.subplots(2,1,figsize=(12,8)); axs=axs.flatten()
	axs[0].imshow(umean[i,...,1])
	axs[1].imshow(urms[i,...,1])
	plt.show()
	plt.pause(0.2)
