import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.integrate import simpson, cumulative_trapezoid
import matplotlib as mpl
import os

class Data1D():
    """
    This class is desiged to handle the 1D scan UMS scan (x-scan or y-scan) 
    The scanned axis is labeled x.
    """


    def __init__(self,
                 PA_map_output_name:str,
                 PA_raw_output_prefix_name:str,
                 initial_distance:float,
                 center_frequency_Hz:float,
                 pulse_repetition_frequency_Hz:float,
                 scanned_axis_letter:str,
                 export_folder:str = ""):
        """
        PA_map_output_name: the path to the UMSmap file of the scan
        PA_raw_output_prefix_name: the path to the individual volt vs time data files, excluding the "Y000X000.txt" part
        center_frequency_Hz: center frequency used in the test, in Herz
        pulse_repetition_frequency_Hz: pulse repetition frequency used in the test, in Herz 
        export_folder: the path that the figures and exported data files will be saved"""
        
        print("Data class is initialized.\n")
        self.PA_map_output_name = PA_map_output_name
        self.PA_raw_output_prefix_name = PA_raw_output_prefix_name 
        
        self.export_folder = export_folder
        if not os.path.exists(self.export_folder):
            os.mkdir(self.export_folder)

        self.initial_distance = initial_distance
        self.center_frequency_Hz = center_frequency_Hz
        self.center_frequency_MHz = self.center_frequency_Hz/1e6
        self.PRF_Hz =  pulse_repetition_frequency_Hz
        
        
        self.scanned_axis_letter = str(scanned_axis_letter)
        if str(scanned_axis_letter).upper() not in ["X", "Y"]:
            print("The scanned axis letter is not recognized. Try using 'x' or 'y'.\n")
            exit()    

        self.read_PA_map_output_file()
        self.x = self.x + self.initial_distance
        
        self.plot_for_starting_distance_calculation()






    def read_PA_map_output_file(self):
        """Reads the UMSmap file of the scan."""
        print("Reading the PA map file output.")
        data = np.genfromtxt(self.PA_map_output_name, delimiter=None, skip_header=81)
        self.x = np.asarray(data[:,0])       
        self.x_count = self.x.shape[0]   
        print("Done...\n")



    

    def plot_for_starting_distance_calculation(self):

        full_filename = f"{self.PA_raw_output_prefix_name}X{0:0>3}.txt"
        data = np.genfromtxt(full_filename, delimiter=None, skip_header=81)
        time = np.asarray(data[:,0])
        volt = np.asarray(data[:,1])

        fig1, ax1 = plt.subplots()
        ax1.set_title("Signal when distance is minimum")
        ax1.set_xlabel("Time, s")
        ax1.set_ylabel("Voltage, V")
        ax1.plot(time,volt, lw=2)



    def calculate_intensities_and_MI(self, time, volt, x_index):
        """Calculates the intensities and some other parameters from a given set of time and voltage array.
        https://www.brl.uiuc.edu/Downloads/sakai/SakaiChapter2.pdf
        https://scholarworks.sjsu.edu/cgi/viewcontent.cgi?article=4780&context=etd_theses"""

        # sensitivity of the hydrophone is 447 mV/MPa
        # ... = 447*1e-3 V/MPa = 447*1e-9 V/Pa
        HYDROPHONE_SENSITIVITY = 447*1E-9 # [V/Pa]
        WATER_DENSITY = 998 # water density [kg/m3]
        WAVE_SPEED = 1480 # speed of sound in water [m/s]
        ATTENUATION_CONST = 0.3 # acoustic attenuation in water [dB/cm/MHz]
        DISTANCE = self.x[x_index]/ 10 # between the hydrophone and transducer (cm)
        
        # pulse intensity integral, I(PA), W*s/cm2  
        PII_Wscm2 = cumulative_trapezoid(volt**2, time, initial=0) / (WATER_DENSITY*WAVE_SPEED*HYDROPHONE_SENSITIVITY**2) 
        PII_Wscm2 *= 1e-4 # here, PII_Wscm2 is an array with a same size as x

        t10_index = np.argmin(np.abs(PII_Wscm2 - 0.1*np.max(PII_Wscm2)))
        t90_index = np.argmin(np.abs(PII_Wscm2 - 0.9*np.max(PII_Wscm2)))
        PD = 1.25 * np.abs( (time[t90_index]-time[t10_index]) ) # pulse duration in seconds
        
        # PII_Wscm2 = float(PII_Wscm2[-1]) # the last index of PII_Wscm2 is the actual integral
        # OR
        PII_Wscm2 = PII_Wscm2[t90_index]-PII_Wscm2[t10_index]
        
        PII_derated = PII_Wscm2 * np.exp(-0.115*ATTENUATION_CONST*self.center_frequency_MHz*DISTANCE) # W*s/cm2

        I_spta_Wcm2 = PII_derated * self.PRF_Hz  # Spatial peak temporal average intensity, W/cm2
        I_sppa_Wcm2 = PII_derated / PD # Spatial peak pulse average intensity, W/cm2
        
        pr_Pa = np.abs(np.min(volt))/HYDROPHONE_SENSITIVITY # Peak rarefactional pressure (most negative pressure)
        pr_MPa = pr_Pa * 1e-6
        pr_derated_MPa = pr_MPa * np.exp(-0.115*ATTENUATION_CONST*self.center_frequency_MHz*DISTANCE) # derated peak rarefactional pressure, Pa
        
        return I_spta_Wcm2, I_sppa_Wcm2, PII_Wscm2, pr_derated_MPa






    def perform_acoustic_calculations(self):
        """Used to calculate the intensities and mechanical index of the whole scan."""
        print("Calculating the intensities and mechanical index...")

        self.mapped_I_spta_Wcm2 = np.zeros(self.x_count)
        self.mapped_I_sppa_Wcm2 = self.mapped_I_spta_Wcm2.copy()
        self.mapped_PII_Wscm2 = self.mapped_I_spta_Wcm2.copy()
        self.mapped_pr_derated_MPa = self.mapped_I_spta_Wcm2.copy()
        self.pr = self.mapped_I_spta_Wcm2.copy()

        for x_index in range(self.x_count):
            full_filename = f"{self.PA_raw_output_prefix_name}{self.scanned_axis_letter}{x_index:0>3}.txt"
            if x_index % 20 == 0:
                print(f"__ file {x_index}/{self.x_count} processed.")
            
            data = np.genfromtxt(full_filename, delimiter=None, skip_header=81)
            time = np.asarray(data[:,0])
            volt = np.asarray(data[:,1])
            volt -= np.mean(volt)
            
            I_spta_Wcm2, I_sppa_Wcm2, PII_Wscm2, pr_derated_MPa = self.calculate_intensities_and_MI(time, volt, x_index)
            
            
            self.mapped_I_spta_Wcm2[x_index] = I_spta_Wcm2
            self.mapped_I_sppa_Wcm2[x_index] = I_sppa_Wcm2
            self.mapped_PII_Wscm2[x_index] = PII_Wscm2
            self.mapped_pr_derated_MPa[x_index] = pr_derated_MPa

        PII_max_index = np.unravel_index(self.mapped_PII_Wscm2.argmax(), self.mapped_PII_Wscm2.shape)
        pr_derated__ = self.mapped_pr_derated_MPa[ PII_max_index]


        # The refernce MHz here, but it results in huge mechanical indexes
        self.mechanical_index = pr_derated__ / np.sqrt(self.center_frequency_MHz)
        
        
        
        self.max_I_spta_Wcm2 = np.max(self.mapped_I_spta_Wcm2)
        self.max_I_sppa_Wcm2 = np.max(self.mapped_I_sppa_Wcm2)
        print("Done...\n")






    def interpolate_intensity_maps(self, K:int):
        """Interpolates the 2D x,y,data arrays."""
        self.interpolated_x = np.linspace(self.x[0],self.x[-1], int(self.x_count*K))
    
        interpolate1 = interp1d(self.x, self.mapped_I_spta_Wcm2)
        self.interpolated_mapped_I_spta_Wcm2 =  interpolate1(self.interpolated_x)
        
        interpolate2 = interp1d(self.x, self.mapped_I_sppa_Wcm2)
        self.interpolated_mapped_I_sppa_Wcm2 = interpolate2(self.interpolated_x)





    def export_calculated_parameters(self):
        """Exports the calculated parameters to external files in csv format."""

        print("Exporting calculated parameters.")
        export = np.zeros((self.x.size, 2))
        export[:,0] = self.x
        export[:,1] = self.mapped_I_spta_Wcm2
        np.savetxt(f"{self.export_folder}\\calculated_Ispta.csv", export, delimiter=",", header="x mm, Ispta W/cm2")

        export = np.zeros((self.x.size, 2))
        export[:,0] = self.x
        export[:,1] = self.mapped_I_sppa_Wcm2
        np.savetxt(f"{self.export_folder}\\calculated_Isppa.csv", export, delimiter=",", header="Ispta W/cm2 (1st row: x and 1st column: y) - (0_0 is a placeholder)")
        
        export = np.zeros((self.interpolated_x.size, 2))
        export[:,0] = self.interpolated_x
        export[:,1] = self.interpolated_mapped_I_spta_Wcm2
        np.savetxt(f"{self.export_folder}\\calculated_Ispta_interpolated.csv", export, delimiter=",", header="Ispta W/cm2 (1st row: x and 1st column: y) - (0_0 is a placeholder)")
        
        export = np.zeros((self.interpolated_x.size, 2))
        export[:,0] = self.interpolated_x
        export[:,1] = self.interpolated_mapped_I_sppa_Wcm2
        np.savetxt(f"{self.export_folder}\\calculated_Isppa_interpolated.csv", export, delimiter=",", header="Isppa W/cm2 (1st row: x and 1st column: y) - (0_0 is a placeholder)")
  
        print("Done...\n")






    def plot_calculated_parameters(self,    
                                   matrices_to_be_plotted:tuple,
                                   title:str,
                                   axis_labels:tuple,
                                   colorbar_label:str,
                                   colormap=mpl.cm.jet,
                                   save=False):
        """Used to plot a 2D colormap plot.
        matrices_to_be_plotted: 3 element tuple, containing x(meshgrid), y(meshgrid), 2D array, respectively."""
        print(f"Plotting <<{title}>>.")

        fig1, ax1 = plt.subplots()
        plt.rc("font", size=14, family="Arial")
        ax1.set_title(title)
        ax1.set_xlabel(axis_labels[0])
        ax1.set_ylabel(axis_labels[1])
        ax1.plot(matrices_to_be_plotted[0], matrices_to_be_plotted[1], lw=2, color="#ed7e77", label=colorbar_label)
        ax1.legend(loc=1)
        if save:
            fig1.savefig(f"{self.export_folder}\\{title}.png",bbox_inches='tight',dpi=300)
            fig1.savefig(f"{self.export_folder}\\{title}.svg",bbox_inches='tight',dpi=300)
        print("Done...\n")


    def show_plots(self):
        plt.show()
    


