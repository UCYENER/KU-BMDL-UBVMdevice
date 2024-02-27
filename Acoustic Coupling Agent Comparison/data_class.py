import numpy as np 
import matplotlib.pyplot as plt
import signal_processing as sp
import os
import pandas as pd


class Data():

    def __init__(self, data_folder, file_namings, number_of_measurements):
        self.data_folder = data_folder
        self.file_namings = file_namings
        self.measurement_number = number_of_measurements
       



    def process(self):
        
        echo_amplitudes = np.empty((self.measurement_number, len(self.file_namings)))
        
        
        for i, naming in enumerate(self.file_namings):
            for mea in range(self.measurement_number):
                filename = f"{self.data_folder}{naming}{str(mea+1)}.txt"
                time_s, volt_V = self.read_raw_data(filename)
                sampling_rate = (time_s.size-1) / (time_s[-1] - time_s[0])
                
                filtered_volt_V = sp.apply_bandpass_filter(volt_V, [1e6,3e6], 4, sampling_rate)
                time_s, filtered_volt_V = self.section_signal(time_s, filtered_volt_V, [20e-6, 45e-6])

                Vpp = np.abs( np.max(filtered_volt_V) - np.min(filtered_volt_V) )
                echo_amplitudes[mea,i] = Vpp
                
        measurement_indexes = np.arange(1,11)
        plt.figure(figsize=(8,6))
        plt.rc("font", size=14, family="Arial")
        for i in range(len(self.file_namings)):
            plt.scatter(measurement_indexes, 1e3*echo_amplitudes[:,i], 75, label=self.file_namings[i])
        plt.xlabel("Measurements"); plt.ylabel("Echo Peak Amplitude [mV]")
        plt.legend()
        plt.savefig("each measurement result.PNG", dpi=300, bbox_inches="tight")
        plt.savefig("each measurement result.SVG", dpi=300, bbox_inches="tight")
        
        exportDict = {}
        for i,agent in enumerate(self.file_namings):
            exportDict[f"{str(agent)} [V]"] = echo_amplitudes[:,i]
        df = pd.DataFrame(exportDict, index=measurement_indexes)
        df.to_csv("each measurement result.csv")
        
        
        
        means = np.mean(echo_amplitudes, axis=0)
        stds = np.std(echo_amplitudes, axis=0)
        errors = stds / np.sqrt(self.measurement_number)
        
        plt.figure(figsize=(8,6))
        plt.rc("font", size=14, family="Arial")
        x_axis = ["Without\nCoupling Agent", "Liquid\nUS Gel", "Solid\nUS Gel", "Silbione\nGel", "Overcured\nSilbione Gel"]
        # plt.title("Acoustic Coupling Agent Comparison")
        plt.errorbar(x_axis, 1e3*means, yerr=1e3*errors, fmt="o", lw=1.5, color="red", ecolor="black", capsize=10, 
                     capthick=2, barsabove=True)
        plt.xlabel(""); plt.ylabel("Echo Peak Amplitude [mV]")
        plt.savefig("With error bars.PNG", dpi=300, bbox_inches="tight")
        plt.savefig("With error bars.SVG", dpi=300, bbox_inches="tight")
        
        

        exportDict = {"Mean [V]": means, "Standard Dev. [V]": stds, "Error [V]": errors}
        df = pd.DataFrame(exportDict, index=self.file_namings)
        df.to_csv("with error bars.csv")
        




        
        
        
    def section_signal(self, time_s, filtered_voltage_V, limits):
        idxlow = np.argmin( np.abs( time_s - limits[0] ) )
        idxhigh = np.argmin( np.abs( time_s - limits[1] ) ) 
        return time_s[idxlow:idxhigh], filtered_voltage_V[idxlow:idxhigh]
        
    
    def read_raw_data(self, filename):
        data = np.genfromtxt(filename, delimiter=None, skip_header=81)
        time = np.asarray(data[:,0])
        voltage = np.asarray(data[:,1])
        voltage = np.nan_to_num(voltage, copy=True, nan=0.0, posinf=np.nanmax(voltage), neginf=np.nanmin(voltage))
        return time, voltage-np.mean(voltage)
    

    def show_plots(self):
        plt.show()
