from os import system
from sys import exit
from data_class import Data
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 



def main():
    system("cls")

    data_folder = "data\\"
    file_namings = ["no_silbione", "with_silbione"]
    number_of_measurements = 16
    sample_thickness_mm = 2.2
    
    data = Data(data_folder, file_namings, number_of_measurements, sample_thickness_mm)

    data.process()
    
    data.show_plots()



if __name__ == "__main__":
    main()
    
