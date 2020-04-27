'''This code plots the CV measurement figures from files located in the 'data' folder'''

__author__ = "Shirajum Monira"

import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import numpy as np
import math 
import glob

def straight_lines(x,a,b): #defines straight line 
	first = a*x + b
	return first

file = sorted(glob.glob('*.txt')) #returns an unsorted list of all files with .txt extension
#guess values for the straight line fitting
guess1 = (0.00015,0.012)
guess2 = (0.0006,0)
guess3 = (-0.000183,0.011)
guess4 = (-0.000275,0)
#dictionary with keys as index value of file list. The keys contain guess values of particular files and voltage values to calculate range of each straight line
file_dictionary = {
	
	0:(guess1,guess2,25,18),
	2:(guess1,guess2,25,18),
	4:(guess3,guess4,-40,-30),
}

for index in [0,2,4]:

	data_needle = np.loadtxt(file[index],skiprows=1) #data from closed needle files 
	data_open_needle = np.loadtxt(file[index+1],skiprows=1) #data from open needle files

	source_voltage = data_needle[:,0] #reads column source voltage 
	#subtracts baseline from main data and calculates the true capacitances in pF
	true_capacitance1 = (data_needle[:,1]-data_open_needle[:,1])*1e12
	true_capacitance2 = (data_needle[:,3]-data_open_needle[:,3])*1e12
	true_capacitance3 = (data_needle[:,5]-data_open_needle[:,5])*1e12

	#print(true_capacitance1[0],true_capacitance2[0],true_capacitance3[0])
	print('You are now looking at file -- ',file[index])
	
	plt.figure()
	#plotting 1/c^2 vs source voltage
	plt.plot(source_voltage,1/(true_capacitance1)**2)
	plt.plot(source_voltage,1/(true_capacitance2)**2)
	plt.plot(source_voltage,1/(true_capacitance3)**2)
	plt.suptitle('Detector True Capacitance vs. Bias Voltage')
	plt.xlabel('Bias Voltage (V)')
	plt.ylabel('Capacitance$^{-2}$ in (pF$^{-2}$)')

	#creates mask with voltage values in absolute forms to fix range of each straight line
	v1_index = np.argmax(np.abs(source_voltage) > np.abs(file_dictionary[index][2]))
	v2_index = np.argmax(np.abs(source_voltage) > np.abs(file_dictionary[index][3]))
	#using mask to generate range for line fitting later
	new_source_voltage1 = source_voltage[v1_index:]
	new_source_voltage2 = source_voltage[:v2_index]
	new_source_voltage3 = source_voltage[:v1_index]
	new_capacitance1 = true_capacitance1[v1_index:]
	new_capacitance2 = true_capacitance1[:v2_index]

	#fits first straight line
	g1,cov = curve_fit(straight_lines,new_source_voltage1,1/(new_capacitance1)**2,file_dictionary[index][0])
	plt.plot(source_voltage,straight_lines(source_voltage,g1[0],g1[1]))
	#fits second straight line
	g2,cov = curve_fit(straight_lines,new_source_voltage2,1/(new_capacitance2)**2,file_dictionary[index][1])
	plt.plot(new_source_voltage3,straight_lines(new_source_voltage3,g2[0],g2[1]))

	#calculating depletion voltage by equating the straight lines at the intersection point
	depletion_voltage = (g2[1]-g1[1])/(g1[0]-g2[0])
	plt.axvline(depletion_voltage) #drops a vertical line from intersection to the X-axis
	print('The full depletion voltage calculated is =',depletion_voltage)
	plt.savefig('figure3.pdf', bbox_inches = 'tight')

plt.show()
