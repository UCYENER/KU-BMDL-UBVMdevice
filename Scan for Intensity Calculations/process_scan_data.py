from data_class import Data1D
from os import system
from sys import exit 


if __name__ == "__main__":
    system("cls")

    # axial scan (1D)
    mapfile_path = "Axial Scan\\alignment_UMSProfile.txt" # UMSmap file path in the data folder
    raw_file_prefix = "Axial Scan\\alignment" # raw file path, exclude the X and Y index parts and file format
    export_folder = "00_python_output_axial_scan" # the FOLDER you want to export into


    pulse_period_Seconds = 2.5 # the period of the pulses in seconds
    pulse_repetition_frequency_Hz = 1/pulse_period_Seconds
    center_frequency_Hz = 2e6 # frequency used in the experiment
    scanned_axis_letter = "x"
    initial_distance = 3.3e-6 * 1500 * 1e3 # mm - UPDATE THIS AFTER EXAMINING THE FIRST DATA FILE
    
    data = Data1D(mapfile_path, raw_file_prefix, 
                  initial_distance, center_frequency_Hz, 
                  pulse_repetition_frequency_Hz, 
                  scanned_axis_letter, export_folder)

    

    data.perform_acoustic_calculations()
    data.interpolate_intensity_maps(10) # the argument multiples the number of X and Y coordinate points
        
    print(f"Max Ispta = {data.max_I_spta_Wcm2:.2E} W/cm2")
    print(f"Max Isppa = {data.max_I_sppa_Wcm2:.2E} W/cm2")
    print(f"Mechanical Index = {data.mechanical_index:.2E}\n")
    

    data.plot_calculated_parameters((data.x, data.mapped_I_spta_Wcm2), "I_spta", ("X [mm]", "Intensity"), "I_spta W/cm2", save=True)
    data.plot_calculated_parameters((data.x, data.mapped_I_sppa_Wcm2), "I_sppa", ("X [mm]", "Intensity"), "I_sppa W/cm2", save=True)

    # data.plot_calculated_parameters((data.interpolated_x, data.interpolated_mapped_I_spta_Wcm2), "Interpolated I_spta", ("X [mm]", "Y [mm]"), "I_spta W/cm2", save=True)
    # data.plot_calculated_parameters((data.interpolated_x, data.interpolated_mapped_I_sppa_Wcm2), "Interpolated I_sppa", ("X [mm]", "Y [mm]"), "I_sppa W/cm2", save=True)


    data.export_calculated_parameters()
    
    data.show_plots()
    
    