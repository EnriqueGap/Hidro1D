import matplotlib.pyplot as plt
import numpy as np
import imageio
import os

for i in range(10):
	a=str(i)
	rho=np.loadtxt('rho-0'+a+'.dat')
	vx=np.loadtxt('vx-0'+a+'.dat')
	vy=np.loadtxt('vy-0'+a+'.dat')
	energy=np.loadtxt('pressure-0'+a+'.dat')

	plt.subplots(figsize=(12,24))
	plt.subplot(4,1,1)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Density a t='+a)
	plt.imshow(rho, cmap='viridis', vmin=0, vmax=1)
	plt.colorbar(shrink=0.6)

	plt.subplot(4,1,2)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Vx a t='+a)
	plt.imshow(vx, cmap='inferno', vmin=0, vmax=1)
	plt.colorbar(shrink=0.6)

	plt.subplot(4,1,3)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Vy a t='+a)
	plt.imshow(vy, cmap='plasma', vmin=0, vmax=1)
	plt.colorbar(shrink=0.6)

	plt.subplot(4,1,4)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Pressure a t='+a)
	plt.imshow(energy, cmap='turbo', vmin=0, vmax=1)
	plt.colorbar(shrink=0.6)

	plt.savefig('im'+a+'.png')

with imageio.get_writer("test.gif", mode='I') as writer:
	for i in range(10):
		writer.append_data(imageio.imread('im'+str(i)+'.png'))
os.system('eog test.gif')
