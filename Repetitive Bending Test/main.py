from os import system
from sys import exit
from data_class import Data



def main():
    system("cls")

    data_folder = "20240119\\" # "bendtestpulseecho\\"
    file_namings = ["normalPE", "afterBending3"]
    number_of_measurements = 25
    
    data = Data(data_folder, file_namings, number_of_measurements)

    data.process_with_error_bars()
    
    data.process_and_compare_averaged()
    
    data.show_plots()



if __name__ == "__main__":
    main()
    
