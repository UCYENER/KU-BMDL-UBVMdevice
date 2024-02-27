import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from os import system
from sys import exit


def importData(filename):
	print(f'$$$ Importing {filename}')
	data = np.genfromtxt(filename, delimiter=None, skip_header=81)
	# 3 cols: x(time in s), y(voltage in V), Uncertaintyn
	time = np.asarray(data[:,0])
	volt = np.asarray(data[:,1])
	volt = volt - np.mean(volt)
	return time,volt

def calcPressure(volt):
	volt_pkpk = np.max(volt)-np.min(volt)
	sensitivity = 447*1E-9 # [V/Pa]	(447 mV/MPa) 
	rho = 998 # kg/m3
	c_water = 1474 # m/s
	pressure = volt_pkpk / sensitivity # [Pa]
	return pressure

def Main():
	fileprefix = 'data\\20230510_freq_sweep_of_transducer3f'
	pressureArr = np.zeros(101)
	f = np.linspace(1,3,101)
	for i in range(101):
		_, volt = importData(f'{fileprefix}{i:03}.txt')
		pressure = calcPressure(volt)
		pressureArr[i] = pressure 

	pressureArr /= np.max(pressureArr) # normalize
	pressureArr = 10*np.log10(pressureArr / np.max(pressureArr))


	idx = np.argmin(np.abs(pressureArr-np.max(pressureArr)))
	print(f'Max @ {f[idx]}')
	
	# idx1 = np.argmin(np.abs(pressureArr[:51] + 3))
	# idx2 = np.argmin(np.abs(pressureArr[51:] + 3))
	# BW = 100*2* np.abs(f[idx1]-f[idx2+51]) / (np.abs(f[idx1]+f[idx2+51]))
	# print(f'-3dB BW = {BW}%')
	# print(f'-3dB center frequency = {(np.abs(f[idx1]+f[idx2+51]))/2} MHz')

	

	plt.figure()
	plt.rc('font', size=14)
	plt.title('Frequency Sweep')
	plt.plot(f, pressureArr)
	plt.xlabel('Frequency [MHz]')
	plt.ylabel('Pressure [dB]')


	f = f.reshape(f.size,1)
	pressureArr = pressureArr.reshape(pressureArr.size,1)
	data = np.concatenate((f*1e6,pressureArr), axis=1)
	np.savetxt("freq_sweep.csv", data, delimiter=",", header='Frequency [MHz],Pressure [dB]')


if __name__ == '__main__':
	system('cls')
	Main()
	plt.show()