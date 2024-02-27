from os import system
from sys import exit
from data_class import Data
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 



def main():
    system("cls")

    data_folder = "data - changed names\\"
    file_namings = ["no coupling agent", "liquid us gel", "solid us gel", "silbione gel", "overcured silbione gel"]
    number_of_measurements = 10
    
    data = Data(data_folder, file_namings, number_of_measurements)

    data.process()
    
    data.show_plots()



if __name__ == "__main__":
    main()
    
