import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import imageio
import os

for i in range(10):
	a=str(i)
	data=ascii.read('hd-0'+a+'.dat')
	plt.subplots(figsize=(12,18))
	plt.subplot(3,1,1)
	plt.xlabel('x')
	plt.ylabel(r'$\rho$')
	#plt.xlim((0,500))
	#plt.ylim((-0.1,1.1))
	plt.grid()
	plt.title('Densidad a t='+a)
	plt.plot(data['col1'],data['col2'])
	plt.subplot(3,1,2)
	plt.xlabel('x')
	plt.ylabel('u')
	#plt.xlim((0,500))
	#plt.ylim((-0.1,1.1))
	plt.grid()
	plt.title('Velocidad a t='+a)
	plt.plot(data['col1'],data['col3'])
	plt.subplot(3,1,3)
	plt.xlabel('x')
	plt.ylabel('P')
	#plt.xlim((0,500))
	#plt.ylim((-0.1,1.1))
	plt.grid()
	plt.title('Presión a t='+a)
	plt.plot(data['col1'],data['col4'])

	plt.savefig('im'+a+'.png')

images = []
for i in range(10):
        images.append(imageio.imread('im'+str(i)+'.png'))
imageio.mimsave('mygif.gif', images, fps=3)

os.system('eog mygif.gif')
