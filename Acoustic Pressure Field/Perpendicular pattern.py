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

def importMapData(filename):
	print(f'$$$ Importing {filename}')
	data = np.genfromtxt(filename, delimiter=None, skip_header=4)
	# (n+1)x(m+1) array, normalized to 1
	x = np.asarray(data[0,1:])
	y = np.asarray(data[1:,0])
	x,y = np.meshgrid(x,y)
	data = np.asarray(data[1:,1:])
	return x,y,data

def calcPressure(volt):
	volt = volt - np.mean(volt)
	volt_pkpk = np.max(volt)-np.min(volt)
	sensitivity = 447*1E-9 # [V/Pa]	(447 mV/MPa) 
	rho = 998 # kg/m3
	c_water = 1474 # m/s
	pressure = volt_pkpk / sensitivity # [Pa]
	return pressure


def Main():
	fileprefix = 'data\\'
	mapfile = fileprefix + 'tx_with_ml_2d_scan_at_focus3_UMSmap.txt'

	x,y,data = importMapData(mapfile)
	data = data/np.max(data)

	plt.figure(figsize=(10,7))
	plt.rc('font', size=14)
	plt.title('NORMALIZED VOLTAGE SQUARED INTEGRAL')
	plt.pcolormesh(x,y,data, cmap=mpl.cm.jet)# Greys, inferno etc. (matplotlib.cm)
	plt.colorbar(label='Magnitude')
	plt.xlabel('X [mm]')
	plt.ylabel('Y [mm]')


	
	#######################################


	interp_x = np.linspace(x[0,0],x[0,-1], x.shape[1]*10)
	interp_y = np.linspace(y[0,0],y[-1,0], y.shape[0]*10)
	interp_x,interp_y = np.meshgrid(interp_x,interp_y)

	data = np.transpose(data)

	interp_data = interpn((x[0,:],y[:,0]), data, (interp_x,interp_y))
	
	plt.figure()
	plt.rc('font', size=14)
	plt.gca().set_aspect('equal')
	plt.title('INTERPOLATED DISTRIBUTION')
	plt.pcolormesh(interp_x, interp_y, interp_data, cmap=mpl.cm.jet)
	plt.colorbar(label='Magnitude')
	plt.xlabel('X [mm]')
	plt.ylabel('Y [mm]')

	##########################################################
	##########################################################
	##########################################################
	##########################################################

	fileprefix = 'data\\tx_with_ml_2d_scan_at_focus3'
	pressureMap = np.zeros(data.shape)
	pressureMap = np.transpose(pressureMap)
	for i in range(20):
		for j in range(20):		
			fullfilename = f'{fileprefix}Y{i:03}X{j:03}.txt'
			_, volt = importData(fullfilename)
			pressure = calcPressure(volt)
			pressureMap[i,j] = pressure 

	pressureMap = np.transpose(pressureMap)
	pressureMap = interpn((x[0,:],y[:,0]), pressureMap, (interp_x,interp_y), method='linear',bounds_error=False,fill_value=None)
	pressureMap /= np.max(pressureMap) # normalize
	pressureMap = 10*np.log10(pressureMap / np.max(pressureMap))
	
	datatofile = np.zeros((pressureMap.shape[0]+1 , pressureMap.shape[1]+1))
	datatofile[0,1:] = interp_x[0,:]
	datatofile[1:,0] = interp_y[:,0]
	datatofile[1:,1:] = pressureMap

	np.savetxt("PressureMap_perpendicular.csv", datatofile, delimiter=",", header='1st row: X, 1st col: Y')
	#######################################


	plt.figure()
	plt.rc('font', size=14)
	plt.title('Pressure Distribution')
	plt.pcolormesh(interp_x, interp_y, pressureMap, cmap=mpl.cm.jet)
	plt.colorbar(label='Acoustic Pressure [dB]')
	plt.xlabel('X [mm]')
	plt.ylabel('Y [mm]')

	

if __name__ == '__main__':
	system('cls')
	Main()
	plt.show()