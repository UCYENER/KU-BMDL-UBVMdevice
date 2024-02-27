from os import system
from sys import exit
from data_class import Data
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 

""" 
Experiment: "in-vitro tests for volume calculation algorithm validation"
Here, the aim is to take data at as many points as possible using a commercial transducer and detect
the anterior and posterior walls of the differently shaped flasks, then use various combinations of
those points/coordinates with our algorithm and show that using only 4 transducers is enough for our
application. Additionally, ellipsoid fitting can be done on the points, since we will have more than
8 coordinates to use during the fitting. By doing this, we can have the opportunity to compare the
performances of spherical and ellipsoid fitting.
"""


def plot_every_waveform_1by1(data, shape):
    for x_index in range(data.x_count):
        for y_index in range(data.y_count):
            fig, ax = plt.subplots()
            fig.set_size_inches(10,7)
            ax.set_title(f"{shape} - {x_index+1}x{y_index+1}")
            ax.set_xlabel("Time, us"); ax.set_ylabel("Voltage, mV")
            ax.plot(data.filtered_waveforms_us_mV[f"x{x_index}_y{y_index}"][0],
                    data.filtered_waveforms_us_mV[f"x{x_index}_y{y_index}"][1], color="black", lw=1.5)
            ax.scatter(data.peak_times_us[f"x{x_index}_y{y_index}"],
                       data.peak_voltages_mV[f"x{x_index}_y{y_index}"], 100, color="red")
            data.show_plots()




def main():
    system("cls")
    settings_json = "settings.json"
    shape_list = ( "spherical", "ellipsoid", "conical" )
    
    for shape in shape_list:
        data = Data(settings_json, shape)
        print("\n\n>>>>>>>",data)
        data.process()
        

    
    data.show_plots()


if __name__ == "__main__":
    main()
    
