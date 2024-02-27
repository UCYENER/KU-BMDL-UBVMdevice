import numpy as np 
import matplotlib.pyplot as plt


class Data():

    def __init__(self, data_folder, file_namings, number_of_measurements, sample_thickness_mm):
        self.data_folder = data_folder
        self.file_namings = file_namings
        self.measurement_number = number_of_measurements
        self.sample_thickness_mm = sample_thickness_mm
        
        self.get_sampling_rate(file_namings[0])
       



    def process(self):
        
        averaged_volt_V = {}
        FFT_results = {}

        for case_name in self.file_namings:
            for mea_num in range(self.measurement_number):
                full_filename = f"{self.data_folder}{case_name}_{str(mea_num+1)}.txt"
                time_s, volt_V = self.read_raw_data(full_filename)
                if mea_num == 0:
                    all_mea_volt_V = np.empty((volt_V.size, self.measurement_number))
                    
                all_mea_volt_V[:, mea_num] = volt_V
            
            averaged_volt_V[case_name] = np.mean(all_mea_volt_V, axis=1)
            freq_Hz, out_V = self.apply_fft_on_interval(averaged_volt_V[case_name])
            FFT_results[case_name] = np.concatenate((freq_Hz.reshape(freq_Hz.size,1), out_V.reshape(out_V.size,1)), axis=1)
        
        
        without_f = FFT_results[self.file_namings[0]][:,0]
        idx1 = np.argmin(np.abs(without_f - 0.5e6)) 
        idx2 = np.argmin(np.abs(without_f - 5e6))
        without_out = FFT_results[self.file_namings[0]][:,1]

        without_f = without_f[idx1:idx2]
        without_out = without_out[idx1:idx2]

        with_f = FFT_results[self.file_namings[1]][:,0]
        idx1 = np.argmin(np.abs(with_f - 0.5e6)) 
        idx2 = np.argmin(np.abs(with_f - 5e6))
        with_out = FFT_results[self.file_namings[1]][:,1]

        with_f = with_f[idx1:idx2]
        with_out = with_out[idx1:idx2]

        plt.figure()
        plt.title("FFT of the Received Waveforms")
        plt.semilogy(without_f, without_out, label=self.file_namings[0])
        plt.semilogy(with_f, with_out, label=self.file_namings[1])
        plt.legend()
        # plt.savefig("FFT plot.PNG", dpi=300, bbox_inches="tight")
        # plt.savefig("FFT plot.SVG", dpi=300, bbox_inches="tight")
        
        
        # alpha = 20 * (1/h) * log10(A_water / A_sample)
        alpha = 20 * (1/self.sample_thickness_mm) * np.log10(FFT_results[self.file_namings[0]][:,1] / FFT_results[self.file_namings[1]][:,1])
        frequency = FFT_results[self.file_namings[0]][:,0]
        
        # clip the frequency range
        idx1 = np.argmin(np.abs(frequency - 0.5e6)) 
        idx2 = np.argmin(np.abs(frequency - 5e6))
        frequency = frequency[idx1:idx2]
        alpha = alpha[idx1:idx2]
        
        plt.figure()
        plt.title("Amplitude Attenuation")
        plt.plot(frequency*1e-6, alpha, lw=2, color="orange")
        plt.xlabel("Frequency, MHz"); plt.ylabel("Attenuation, dB/mm")
        # plt.savefig("amplitude attenuation.PNG", dpi=300, bbox_inches="tight")
        # plt.savefig("amplitude attenuation.SVG", dpi=300, bbox_inches="tight")
        
        
        
        idx_2MHz = np.argmin(np.abs(frequency - 2.05e6))
        
        AVERAGE_ATTENUATION = np.mean(alpha[idx_2MHz-2:idx_2MHz+2])
        print(frequency[idx_2MHz-2], frequency[idx_2MHz+2])
        
        print(f"\nAmplitude attenuation = {AVERAGE_ATTENUATION:.2f} dB/mm @2MHz")
        print("Calcualated by taking average of the following frequencies [MHz]:")
        [print(f"{i:.3f}",end="  ") for i in frequency[idx_2MHz-2:idx_2MHz+2]/1e6]
        
            



    
    
    
    def get_sampling_rate(self, filename):
        time_s, _ = self.read_raw_data(f"{self.data_folder}{filename}_1.txt")
        self.sampling_rate_Hz = (time_s.size-1) / (time_s[-1] - time_s[0])
        print(f"Sampling rate: {self.sampling_rate_Hz/1e6:.1f} MSa/s")
        print(f"Sampled time interval: {(time_s[-1]-time_s[0])*1e6:.1f} us")
        print(f"FFT resolution: {self.sampling_rate_Hz/time_s.size/1e3:.3f} kHz")
        
        
        
        
        
    def apply_fft_on_interval(self, volt_V):
        length = volt_V.size
        fft_window = np.hanning(int(length))
        s2 = np.sum(fft_window**2)
        BIN = 1/(length/self.sampling_rate_Hz)
        D_fft = np.zeros(int(length))
        D_fft = np.sqrt(  np.divide( (abs(np.fft.fft(volt_V*fft_window))**2)*2 ,  self.sampling_rate_Hz*s2  ) )
        D_fft_half = D_fft[0:(D_fft.size//2)]
        frequency = np.transpose(BIN*np.arange(0,D_fft_half.size))
        fft_amp = D_fft_half
        return frequency, fft_amp
    
    
    
    
    def read_raw_data(self, filename):
        data = np.genfromtxt(filename, delimiter=None, skip_header=81)
        time = np.asarray(data[:,0])
        voltage = np.asarray(data[:,1])
        voltage = np.nan_to_num(voltage, copy=True, nan=0.0, posinf=np.nanmax(voltage), neginf=np.nanmin(voltage))
        return time, voltage-np.mean(voltage)
    



    def show_plots(self):
        plt.show()
