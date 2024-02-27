import numpy as np 
import matplotlib.pyplot as plt
import signal_processing as sp
import os
import pandas as pd


class Data:

    def __init__(self, data_folder, file_namings, number_of_measurements):
        self.data_folder = data_folder
        self.file_namings = file_namings
        self.measurement_number = number_of_measurements
       

    def process_and_compare_averaged(self):
        averaged_waveforms = {}
        for i, naming in enumerate(self.file_namings):
            print(f"Reading {naming} files...")
            for mea in range(self.measurement_number):
                filename = f"{self.data_folder}{naming}\\{naming}_{mea+1:0>2}.csv"
                time_s, volt_V = self.read_raw_data(filename)
                if mea==0:
                    all_waveforms = np.empty((volt_V.size, self.measurement_number))
                all_waveforms[:,mea] = volt_V
                
            averaged_waveform = np.mean(all_waveforms, 1)
            
            sampling_rate = (time_s.size-1) / (time_s[-1] - time_s[0])
            filtered_volt_V = sp.apply_bandpass_filter(averaged_waveform, [1e6,3e6], 3, sampling_rate)
            time_s, filtered_volt_V = self.section_signal(time_s, filtered_volt_V, [100e-6, 135e-6])
            
            averaged_waveforms[naming] = np.concatenate((time_s.reshape(time_s.size,1), filtered_volt_V.reshape(filtered_volt_V.size,1)), axis=1)
            

        before_waveform = {"time (us)": 1e6*averaged_waveforms[self.file_namings[0]][:,0], "voltage (mV)": 1e3*averaged_waveforms[self.file_namings[0]][:,1]}
        before_waveform = pd.DataFrame(data=before_waveform)
        before_waveform.to_csv("Output\\before_waveform.csv")

        after_waveform = {"time (us)": 1e6*averaged_waveforms[self.file_namings[1]][:,0], "voltage (mV)": 1e3*averaged_waveforms[self.file_namings[1]][:,1]}
        after_waveform = pd.DataFrame(data=after_waveform)
        after_waveform.to_csv("Output\\after_waveform.csv")



        plt.figure(figsize=(8,6))
        plt.rc("font", size=14, family="Arial")
        # plt.title("Repetitive Bending Test")
        plt.plot(1e6*averaged_waveforms[self.file_namings[0]][:,0], 1e3*averaged_waveforms[self.file_namings[0]][:,1], lw=1.5, color="orange", alpha=0.9, label="Before")
        plt.plot(1e6*averaged_waveforms[self.file_namings[1]][:,0], 1e3*averaged_waveforms[self.file_namings[1]][:,1], lw=1.5, color="blue", alpha=0.5, label="After")
        plt.xlabel("Time, us"); plt.ylabel("Amplitude, mV")
        plt.legend()
        plt.savefig("averaged waveforms.PNG", dpi=300, bbox_inches="tight")
        plt.savefig("averaged waveforms.SVG", dpi=300, bbox_inches="tight")
        
        
        f_0, out_0 = self.apply_fft_on_interval(averaged_waveforms[self.file_namings[0]][:,1], sampling_rate)
        f_1, out_1 = self.apply_fft_on_interval(averaged_waveforms[self.file_namings[1]][:,1], sampling_rate)
        
        idx1 = np.argmin(np.abs(f_0 - 1e6))
        idx2 = np.argmin(np.abs(f_0 - 3e6))
        
        f_0 = f_0[idx1:idx2]; f_1 = f_1[idx1:idx2]
        out_0 = out_0[idx1:idx2]; out_1 = out_1[idx1:idx2]


        before_fft = {"frequency (Hz)": f_0, "voltage (mV)": out_0}
        before_fft = pd.DataFrame(data=before_fft)
        before_fft.to_csv("Output\\before_fft.csv")

        after_fft = {"frequency (Hz)": f_1, "voltage (mV)": out_1}
        after_fft = pd.DataFrame(data=after_fft)
        after_fft.to_csv("Output\\after_fft.csv")
        
        plt.figure(figsize=(8,6))
        plt.rc("font", size=14, family="Arial")
        # plt.title("Repetitive Bending Test FFT")
        plt.semilogy(1e-6*f_0, out_0*1e3, lw=2.5, color="red", label="Before")
        plt.semilogy(1e-6*f_1, out_1*1e3, lw=2.5, color="blue", label="After")
        plt.xlabel("Frequency, MHz"); plt.ylabel("Amplitude, mV")
        plt.legend()
        plt.savefig("repetitive bending waveform fft.PNG", dpi=300, bbox_inches="tight")
        plt.savefig("repetitive bending waveform fft.SVG", dpi=300, bbox_inches="tight")
        


    def apply_fft_on_interval(self, volt_V, fs_Hz):
        length = volt_V.size
        fft_window = np.hanning(int(length))
        s2 = np.sum(fft_window**2)
        BIN = 1/(length/fs_Hz)
        D_fft = np.zeros(int(length))
        D_fft = np.sqrt(  np.divide( (abs(np.fft.fft(volt_V*fft_window))**2)*2 ,  fs_Hz*s2  ) )
        D_fft_half = D_fft[0:(D_fft.size//2)]
        frequency = np.transpose(BIN*np.arange(0,D_fft_half.size))
        fft_amp = D_fft_half
        return frequency, fft_amp
    


    def process_with_error_bars(self):             
        echo_amplitudes = np.empty((self.measurement_number, len(self.file_namings)))
        for i, naming in enumerate(self.file_namings):
            print(f"Reading {naming} files...")
            for mea in range(self.measurement_number):
                filename = f"{self.data_folder}{naming}\\{naming}_{mea+1:0>2}.csv"
                time_s, volt_V = self.read_raw_data(filename)
                sampling_rate = (time_s.size-1) / (time_s[-1] - time_s[0])
                filtered_volt_V = sp.apply_bandpass_filter(volt_V, [1e6,3e6], 3, sampling_rate)
                time_s, filtered_volt_V = self.section_signal(time_s, filtered_volt_V, [100e-6, 140e-6])
                Vpp = np.abs( np.max(filtered_volt_V) - np.min(filtered_volt_V) )
                echo_amplitudes[mea,i] = Vpp

        means = np.mean(echo_amplitudes, axis=0)
        print(f"Percent reduction in echo amplitude: {np.abs(means[-1]-means[0]) / means[0] * 100:.1f} % ")
        standard_devs = np.std(echo_amplitudes, axis=0)
        errors = standard_devs / np.sqrt(self.measurement_number)

        plt.figure(figsize=(8,6))
        plt.rc("font", size=14, family="Arial")
        # plt.title("Repetitive Bending Test (n=25)")
        plt.bar(["Before", "After"], means*1e3, color="red", alpha=0.4, width=0.2)
        plt.errorbar(["Before", "After"], means*1e3, yerr=errors*1e3, fmt="*", color="black", ecolor="black", capsize=10, capthick=1, barsabove=True)
        plt.xlabel(""); plt.ylabel("Echo Peak Amplitude [mV]")
        plt.savefig("With error bars.PNG", dpi=300, bbox_inches="tight")
        plt.savefig("With error bars.SVG", dpi=300, bbox_inches="tight")


        indiv_echo_amp_datapoints = {"Voltage (mV) - before": echo_amplitudes[:,0]*1e3, "Voltage (mV) - after": echo_amplitudes[:,1]*1e3}
        indiv_echo_amp_datapoints = pd.DataFrame(data=indiv_echo_amp_datapoints)
        indiv_echo_amp_datapoints.to_csv("Output\\indiv_echo_amp_datapoints.csv")


        comparison = {"Mean Voltages (mV)": means*1e3, "Error (mV)": errors*1e3}
        comparison = pd.DataFrame(data=comparison, index=["before", "after"])
        comparison.to_csv("Output\\comparison.csv")




    def section_signal(self, time_s, filtered_voltage_V, limits):
        idxlow = np.argmin( np.abs( time_s - limits[0] ) )
        idxhigh = np.argmin( np.abs( time_s - limits[1] ) ) 
        return time_s[idxlow:idxhigh], filtered_voltage_V[idxlow:idxhigh]
        
    
    def read_raw_data(self, filename):
        data = np.genfromtxt(filename, delimiter=",", skip_header=3, encoding="utf8")
        time = np.asarray(data[:,0])/1e6 # us -> s 
        voltage = np.asarray(data[:,1])/1e3 # mV -> V
        voltage = np.nan_to_num(voltage, copy=True, nan=0.0, posinf=np.nanmax(voltage), neginf=np.nanmin(voltage))
        return time, voltage-np.mean(voltage)
    

    def show_plots(self):
        plt.show()
