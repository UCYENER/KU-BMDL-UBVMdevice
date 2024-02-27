import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

import os
from sys import exit



class Data():



    def __init__(self, 
                 csv_file_prefix: str, 
                 file_number_range: range,
                 fft_time_window_us: list):
        
        self.csv_file_prefix = csv_file_prefix
        self.file_number_range = file_number_range
        self.fft_time_window_us = fft_time_window_us
        
        self.export_folder = "00_python_export\\"
        if not os.path.exists(self.export_folder):
            os.mkdir(self.export_folder)

        self.averaged_data_export_path = self.export_folder + "Averaged Signal.csv"

        if os.path.isfile(self.averaged_data_export_path):
            self.read_from_exported_file()
        else:
            self.apply_averaging_and_save_averaged()

        self.full_time_us = self.averaged_time_us.copy()
        self.full_averaged_voltage_V = self.averaged_voltage_V.copy()

        self.averaged_voltage_V -= np.mean(self.averaged_voltage_V)

        index1 = np.argmin(np.abs(self.averaged_time_us - self.fft_time_window_us[0]))
        index2 = np.argmin(np.abs(self.averaged_time_us - self.fft_time_window_us[1]))

        self.averaged_time_us = self.averaged_time_us[index1:index2]
        self.averaged_voltage_V = self.averaged_voltage_V[index1:index2]



    



    def apply_averaging_and_save_averaged(self):
        print("Reading all waveform data files...")
        for i, filenum in enumerate(self.file_number_range):
            data = np.genfromtxt(f"{self.csv_file_prefix}{filenum:0>2}.csv", delimiter=',', skip_header=3)
            if 0==i:
                self.averaged_time_us = np.asarray(data[:,0])
                self.all_voltages_V = np.empty((data.shape[0], len(self.file_number_range)))
            self.all_voltages_V[:,i] = np.asarray(data[:,1])-np.mean(np.asarray(data[:,1]))
        print("Applying averaging...")
        self.averaged_voltage_V = np.mean(self.all_voltages_V, axis=1)
        print("Exporting the averaged data...")
        export_data = np.concatenate((self.averaged_time_us.reshape(self.averaged_time_us.size,1), self.averaged_voltage_V.reshape(self.averaged_voltage_V.size,1)), axis=1)
        np.savetxt(self.averaged_data_export_path, export_data, delimiter=",", header="Time us, Voltage V")







    def read_from_exported_file(self):
        print("Reading from the exported data file...")
        data = np.genfromtxt(self.averaged_data_export_path, delimiter=",", skip_header=2)
        self.averaged_time_us = np.asarray(data[:,0])
        self.averaged_voltage_V = np.asarray(data[:,1])






    def plot_signal(self, arrays: np.array, title: str, axis_labels: str, save:bool = False):
        print(f"Plotting {title}...")
        fig1, ax1 = plt.subplots()
        ax1.set_title(title)
        ax1.set_xlabel(axis_labels[0])
        ax1.set_ylabel(axis_labels[1])
        ax1.plot(arrays[0], arrays[1], lw=1)
        if True==save:
            fig1.savefig(f"{self.export_folder}\\{title}.png",bbox_inches='tight',dpi=300)
        plt.show()

    





    def apply_fft_on_interval(self):
        print("Taking FFT of the specified time interval...")
        self.sampling_frequency_Hz = (self.averaged_time_us.size - 1) / (self.averaged_time_us[-1]-self.averaged_time_us[0]) * 1e6
        length = self.averaged_voltage_V.size
        fft_window = np.hanning(int(length))
        s2 = np.sum(fft_window**2)
        BIN = 1/(length/self.sampling_frequency_Hz)
        D_fft = np.zeros(int(length))
        D_fft = np.sqrt(  np.divide( (abs(np.fft.fft(self.averaged_voltage_V*fft_window))**2)*2 ,  self.sampling_frequency_Hz*s2  ) )
        D_fft_half = D_fft[0:(D_fft.size//2)]
        self.frequency = np.transpose(BIN*np.arange(0,D_fft_half.size))
        self.fft_amp = D_fft_half

        export_data = np.concatenate((self.frequency.reshape(self.frequency.size,1), self.fft_amp.reshape(self.fft_amp.size,1)),axis=1)
        np.savetxt(f"{self.export_folder}FFT result.csv", export_data, delimiter=",", header="Frequency Hz, Amplitude V")






    def smoothen_fft_result(self):
        print("Smoothening the FFT results...")
        cs = CubicSpline(self.frequency,self.fft_amp, bc_type='natural')
        self.smooth_frequency = np.linspace(self.frequency[0],self.frequency[-1],self.frequency.size*5)
        self.smooth_fft_amp = cs(self.smooth_frequency)
        ref = np.max(self.smooth_fft_amp)
        self.smooth_fft_amp = 10 * np.log10(self.smooth_fft_amp / ref)
        self.smooth_fft_amp = self.smooth_fft_amp - np.max(self.smooth_fft_amp)

        export_data = np.concatenate((self.smooth_frequency.reshape(self.smooth_frequency.size,1), self.smooth_fft_amp.reshape(self.smooth_fft_amp.size,1)),axis=1)
        np.savetxt(f"{self.export_folder}Smooth FFT result.csv", export_data, delimiter=",", header="Frequency Hz, Amplitude dB")





    def limit_fft_result(self, limits: list):
        print("Limiting the frequency range of the FFT results...")
        idx1 = np.argmin(np.abs(self.frequency - limits[0]))
        idx2 = np.argmin(np.abs(self.frequency - limits[1]))
        self.frequency = self.frequency[idx1:idx2]
        self.fft_amp = self.fft_amp[idx1:idx2]






    def calculate_bandwidths(self):
        print("Calculating the -3dB and -6dB bandwidths...")
        left_half = self.smooth_fft_amp[:self.smooth_fft_amp.size//2]
        right_half = self.smooth_fft_amp[self.smooth_fft_amp.size//2:]
        
        index1 = np.argmin(np.abs(left_half + 6))
        index2 = np.argmin(np.abs(right_half + 6)) + self.smooth_fft_amp.size//2
        self.indexes3dB = [index1, index2]
        freq1 = self.smooth_frequency[index1]
        freq2 = self.smooth_frequency[index2]
        self.center_frequency_6dB = np.average([freq1,freq2])
        self.BW_6dB_percentage = (freq2-freq1)/self.center_frequency_6dB * 100

        index1 = np.argmin(np.abs(left_half + 3))
        index2 = np.argmin(np.abs(right_half + 3)) + self.smooth_fft_amp.size//2
        self.indexes6dB = [index1, index2]
        freq1 = self.smooth_frequency[index1]
        freq2 = self.smooth_frequency[index2]
        self.center_frequency_3dB = np.average([freq1,freq2])
        self.BW_3dB_percentage = (freq2-freq1)/self.center_frequency_3dB * 100


    