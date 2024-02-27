import numpy as np 
import matplotlib.pyplot as plt
import signal_processing as sp
import os, json
import random

class Data():

    def __init__(self, settings_json:str, shape:str):
        
        with open(settings_json, 'r') as json_file:
            settings = json.load(json_file)
        
        self.export_folder = settings["general"]["export_folder"] # folder to export files and images
        if (not os.path.isdir(self.export_folder)): os.mkdir(self.export_folder) # create the folder if it's not there
        self.center_frequency_Hz = settings["general"]["center_frequency_Hz"]
        self.bandpass_cutoffs = settings["general"]["bandpass_cutoffs"]
        self.bandpass_order = settings["general"]["bandpass_order"]
        self.are_limits_set_correctly = settings["general"]["are_limits_set_correctly"]
        self.x_index_range_for_fitting = settings["general"]["x_index_range_for_fitting"]
        self.y_index_range_for_fitting = settings["general"]["y_index_range_for_fitting"]
        self.fitting_alg = settings["general"]["fitting_alg"]
        
        
        self.PA_map_output_name = settings[shape]["map_file"] # UMSmap filename
        self.PA_raw_output_prefix_name = settings[shape]["raw_file"] # Raw signal filename, excluding Y000X000.txt part
        self.anterior_limits = settings[shape]["anterior_limits"]
        self.posterior_limits = settings[shape]["posterior_limits"]
        
        self.object_shape = shape
        self.get_scan_parameters()
        
        if self.x_index_range_for_fitting[1] >= self.x_count: # index out of range
            self.x_index_range_for_fitting[1] = self.x_count
        if self.y_index_range_for_fitting[1] >= self.y_count: # index out of range
            self.y_index_range_for_fitting[1] = self.y_count
        

        self.get_sampling_rate()
        



    def __str__(self):
        return f"<<{self.PA_map_output_name}>>\n({self.x_count}x{self.y_count}) scan."





    def get_scan_parameters(self):
        data = np.genfromtxt(self.PA_map_output_name, delimiter=None, skip_header=4)
        self.x = np.asarray(data[0,1:])
        self.y = np.asarray(data[1:,0])
        self.x_count = self.x.size
        self.y_count = self.y.size
        




    def get_sampling_rate(self):
        center_file = f"{self.PA_raw_output_prefix_name}Y000X000.txt"
        data = np.genfromtxt(center_file, delimiter=None, skip_header=81)
        time_s = np.asarray(data[:,0]); volt_V = np.asarray(data[:,1])
        self.sampling_rate = (time_s.size-1) / (time_s[-1]-time_s[0])
        
        if not self.are_limits_set_correctly:
            print("Print a raw waveform to set the anterior and posterior limits.")
            print("Then, set <are_limits_set_correctly> key to true to run the program.\n")
            exit()

    



    def process(self):
        
        self.filtered_waveforms_us_mV = {}
        self.enveloped_waveforms_us_mV = {}
        self.peak_times_us = {}
        self.peak_voltages_mV = {}
        self.anterior_coordinates = np.empty((self.x_count*self.y_count,3))
        self.posterior_coordinates = self.anterior_coordinates.copy()
        
        
        for x_index in range(self.x_count):
            for y_index in range(self.y_count):
                full_filename = f"{self.PA_raw_output_prefix_name}Y{y_index:0>3}X{x_index:0>3}.txt"
                data = np.genfromtxt(full_filename, delimiter=None, skip_header=81)
                time_s = np.asarray(data[:,0])
                volt_V = np.asarray(data[:,1])
                volt_V -= np.mean(volt_V)

                filtered_volt_V = sp.apply_bandpass_filter(volt_V, self.bandpass_cutoffs, self.bandpass_order, self.sampling_rate)
                amplified_volt_V = sp.amplify_signal(filtered_volt_V, gain=3)
                enveloped_volt_V = sp.get_envelope(amplified_volt_V)
                                
                # if x_index == 2 and y_index == 2:
                #     self.plot_processed_data_plots_and_exit(time_s, volt_V, filtered_volt_V, amplified_volt_V, enveloped_volt_V)
                
                peak_times, peak_voltages = sp.detect_echo_peaks(time_s, enveloped_volt_V, self.anterior_limits, self.posterior_limits)
                self.filtered_waveforms_us_mV[f"x{x_index}_y{y_index}"] = (time_s*1e6, filtered_volt_V*1e3)
                self.enveloped_waveforms_us_mV[f"x{x_index}_y{y_index}"] = (time_s*1e6, enveloped_volt_V*1e3)
                self.peak_times_us[f"x{x_index}_y{y_index}"] = peak_times*1e6
                self.peak_voltages_mV[f"x{x_index}_y{y_index}"] = peak_voltages*1e3
                self.anterior_coordinates[x_index*self.y_count+y_index,:] = np.array([self.x[x_index], self.y[y_index],sp.tof_to_distance_mm(peak_times[0]*1e6)])
                self.posterior_coordinates[x_index*self.y_count+y_index,:] = np.array([self.x[x_index], self.y[y_index],sp.tof_to_distance_mm(peak_times[1]*1e6)])

            
        self.filter_the_selected_measurements()
        
        # # effect of number of coordinates on the volume estimation accuracy
        # sample_count = 16
        # num_of_coords = self.selected_anterior_coordinates.shape[0]
        # selected_indexes = random.sample(range(num_of_coords), sample_count)
        # self.selected_anterior_coordinates = self.selected_anterior_coordinates[list(selected_indexes)]
        # self.selected_posterior_coordinates = self.selected_posterior_coordinates[list(selected_indexes)]
        # print(selected_indexes)
        
        self.selected_coordinates = np.concatenate((self.selected_anterior_coordinates, self.selected_posterior_coordinates), axis=0)
        
        
        
        if self.fitting_alg.lower() == "spherical":
            self.calculated_volume, self.radius, self.center = sp.calculate_the_sphere_volume(self.selected_coordinates)
            self.visualize_selected_coordinates_sphere()
            print(f"Volume: {self.calculated_volume:.2f} mL\n{self.radius:.2f} mm radius and center at {self.center}")
        elif self.fitting_alg.lower() == "ellipsoid":
            self.calculated_volume, self.params, self.center = sp.calculate_the_ellipsoid_volume(self.selected_coordinates)
            self.visualize_selected_coordinates_ellipsoid()
            print(f"Volume: {self.calculated_volume:.2f} mL\n{self.params} mm semi-axis lengths\nCenter at {self.center} mm")
 
 
 
 
 
 
    def filter_the_selected_measurements(self):
        selected_length = np.abs(self.x_index_range_for_fitting[1] - self.x_index_range_for_fitting[0]) * np.abs(self.y_index_range_for_fitting[1] - self.y_index_range_for_fitting[0])
        counter = 0
        self.selected_anterior_coordinates = np.empty(( selected_length, 3))
        self.selected_posterior_coordinates = self.selected_anterior_coordinates.copy()
        
        for x_index in range(self.x_index_range_for_fitting[0], self.x_index_range_for_fitting[1]):
            for y_index in range(self.y_index_range_for_fitting[0], self.y_index_range_for_fitting[1]):
                self.selected_anterior_coordinates[counter,:] = np.array([self.x[x_index], self.y[y_index], sp.tof_to_distance_mm(self.peak_times_us[f"x{x_index}_y{y_index}"][0])])
                self.selected_posterior_coordinates[counter,:] = np.array([self.x[x_index], self.y[y_index], sp.tof_to_distance_mm(self.peak_times_us[f"x{x_index}_y{y_index}"][1])])
                counter += 1
 
 
 
 
 
    def visualize_selected_coordinates_sphere(self):
        """Scatters anterior and posterior coordinates in 3D space."""
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d', aspect="equal")
        ax.set_title(f"{self.object_shape} flask, {self.fitting_alg} fitting - {self.calculated_volume:.1f}mL")
        ax.set_xlabel("X, mm"); ax.set_ylabel("Y, mm"); ax.set_zlabel("Z, mm")
        ax.scatter3D(self.selected_anterior_coordinates[:,0], self.selected_anterior_coordinates[:,1], self.selected_anterior_coordinates[:,2], marker="^", color="red", label="Anterior")
        ax.scatter3D(self.selected_posterior_coordinates[:,0], self.selected_posterior_coordinates[:,1], self.selected_posterior_coordinates[:,2], marker="o", color="blue", label="Posterior")
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = self.radius * np.outer(np.cos(u), np.sin(v)) + self.center[0]
        y = self.radius * np.outer(np.sin(u), np.sin(v)) + self.center[1]
        z = self.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + self.center[2]
        ax.plot_surface(x, y, z, color='pink', alpha=0.4, label="Fitted sphere")
        
        limits = np.array([ ax.get_xlim3d(),
                           ax.get_ylim3d(),
                           ax.get_zlim3d() ])
        origin = np.mean(limits, axis=1)
        max_span = np.max(np.diff(limits, axis=1))
        new_limits = np.array( [    [origin[0]-max_span/2 , origin[0]+max_span/2],
                                    [origin[1]-max_span/2 , origin[1]+max_span/2],
                                    [origin[2]-max_span/2 , origin[2]+max_span/2] ] )
        ax.set_xlim3d(new_limits[0])
        ax.set_ylim3d(new_limits[1])
        ax.set_zlim3d(new_limits[2])
        
        fig.savefig(f"{self.export_folder}\\{self.object_shape}_spherical_fit.png", bbox_inches="tight")

        
        
    def visualize_selected_coordinates_ellipsoid(self):
        """Scatters anterior and posterior coordinates in 3D space."""
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d', aspect="equal")
        ax.set_title(f"{self.object_shape} flask, {self.fitting_alg} fitting - {self.calculated_volume:.1f}mL")
        ax.set_xlabel("X, mm"); ax.set_ylabel("Y, mm"); ax.set_zlabel("Z, mm")
        ax.scatter3D(self.selected_anterior_coordinates[:,0], self.selected_anterior_coordinates[:,1], self.selected_anterior_coordinates[:,2], marker="^", color="red", label="Anterior")
        ax.scatter3D(self.selected_posterior_coordinates[:,0], self.selected_posterior_coordinates[:,1], self.selected_posterior_coordinates[:,2], marker="o", color="blue", label="Posterior")
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = self.params[0] * np.outer(np.cos(u), np.sin(v)) + self.center[0]
        y = self.params[1] * np.outer(np.sin(u), np.sin(v)) + self.center[1]
        z = self.params[2] * np.outer(np.ones(np.size(u)), np.cos(v)) + self.center[2]
        ax.plot_surface(x, y, z, color='pink', alpha=0.4, label="Fitted ellipsoid")
        
        limits = np.array([ ax.get_xlim3d(),
                           ax.get_ylim3d(),
                           ax.get_zlim3d() ])
        origin = np.mean(limits, axis=1)
        max_span = np.max(np.diff(limits, axis=1))
        new_limits = np.array( [    [origin[0]-max_span/2 , origin[0]+max_span/2],
                                    [origin[1]-max_span/2 , origin[1]+max_span/2],
                                    [origin[2]-max_span/2 , origin[2]+max_span/2] ] )
        ax.set_xlim3d(new_limits[0])
        ax.set_ylim3d(new_limits[1])
        ax.set_zlim3d(new_limits[2])
        
        fig.savefig(f"{self.export_folder}\\{self.object_shape}_ellipsoid_fit.png", bbox_inches="tight")

    
        

    def export_calculated_parameters(self, name, array2export, header):
        """Exports the calculated parameters to external files in csv format."""
        np.savetxt(f"{self.export_folder}\\{name}.csv", array2export, delimiter=",", header=header)









    def plot_processed_data_plots_and_exit(self, time_s, volt_V, filtered_volt_V, amplified_volt_V, enveloped_volt_V):
        plt.figure()
        plt.plot(time_s*1e6, volt_V, color="black", label="Raw")
        plt.plot(time_s*1e6, filtered_volt_V, color="blue", label="Filt")
        plt.legend()
                
        plt.figure()
        plt.plot(time_s*1e6, volt_V, color="black", label="Raw")
        plt.plot(time_s*1e6, amplified_volt_V, color="red", label="Filt+Ampl")
        plt.legend()

        plt.figure()
        plt.plot(time_s*1e6, enveloped_volt_V, color="green", label="Enve")
        plt.legend()
        
        plt.show()
        exit()



    def show_plots(self):
        plt.show()
