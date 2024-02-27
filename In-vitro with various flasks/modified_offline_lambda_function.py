
from datetime import datetime
import numpy as np
import ProcessBMDL as bmdl
from sys import exit
from os import system, getcwd

def lambda_handler(data, ants, posts, strats):
    piezonumbers = np.asarray(data[:,3]) // 1000
    timestamps = np.asarray(data[:,5]) * 100 / 32
    
    SUCCESS, volume, r, C, coordinates, refined_timestmaps = bmdl.ProcessData(piezonumbers, timestamps, ants, posts, strats)

    if not SUCCESS:
        return None

    # print(f">>> VOLUME = {volume} mL")
    # print(f">>> Rad = {r} mm, @ {C} mm")
    # print(f">>> {coordinates.shape[0]} Coords")
    # print(f"{coordinates} mm")
    # print(f">>> Selected timestamps =\n{refined_timestmaps.values()} mm \n")

    return volume



if __name__ == "__main__":
    system("cls")

    Shape1 = np.empty(10).reshape(10,1)
    Shape2 = Shape1.copy()
    Shape3 = Shape1.copy()
    Shape4 = Shape1.copy()
    Shape5 = Shape1.copy()
    Shape6 = Shape1.copy()
    Shape7 = Shape1.copy()


    
    ####################################################################### GOOD
   
    # Shape 1 mea 2
    csvfilename = "Shape 1 - armut\\20231221_120554_CO_Sensor_mea2.csv"
    ants = (41, 65); posts = (75, 150)
    strats = ("min", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N
    
    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape1[i] = lambda_handler(partial_data, ants, posts, strats)


    
    ####################################################################### GOOD

    # Shape 2
    csvfilename = "Shape 2 - ellipse\\20231219_61711_PM_CO_Sensor.csv"
    ants = (30, 60); posts = (130, 160)
    strats = ("mean", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N
    
    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape2[i] = lambda_handler(partial_data, ants, posts, strats)



    ####################################################################### GOOD
 
    # Shape 3 mea 2
    csvfilename = "Shape 3 - 250ml ucgenimsi\\20231221_123237_CO_Sensor_mea2.csv"
    ants = (31, 47); posts = (125, 160)
    strats = ("mean", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 5
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape3[i] = lambda_handler(partial_data, ants, posts, strats)

    ####################################################################### GOOD
  
    # Shape 3 mea 3
    csvfilename = "Shape 3 - 250ml ucgenimsi\\20231221_120554_CO_Sensor_mea3.csv"
    ants = (31, 47); posts = (125, 160)
    strats = ("mean", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 5
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape3[i+5] = lambda_handler(partial_data, ants, posts, strats)
    

    
    ####################################################################### 
    # # Shape 3 mea 4 - NOT USED
    # csvfilename = "Shape 3 - 250ml ucgenimsi\\20231221_123237_CO_Sensor_mea4.csv"
        



    ####################################################################### GOOD
  
    # # Shape 4
    csvfilename = "Shape 4 - ucgen\\20231220_125607_PM_CO_Sensor.csv"
    ants = (45, 70); posts = (150, 250)
    strats = ("mean", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape4[i] = lambda_handler(partial_data, ants, posts, strats)


    
    ####################################################################### GOOD
  
    # # Shape 5
    csvfilename = "Shape 5 - 500ml sphere\\20231220_142337_CO_Sensor.csv"
    ants = (50, 100); posts = (150, 250)
    strats = ("max", "mean")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape5[i] = lambda_handler(partial_data, ants, posts, strats)


    
    ####################################################################### GOOD
  
    # # Shape 6
    csvfilename = "Shape 6 - 250ml sphere\\20231220_142337_CO_Sensor.csv"
    ants = (50, 100); posts = (100, 180)
    strats = ("max", "min")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape6[i] = lambda_handler(partial_data, ants, posts, strats)


    ####################################################################### GOOD
  
    # # Shape 7
    csvfilename = "Shape 7 - 100ml sphere\\20231220_142337_CO_Sensor.csv"
    ants = (50, 100); posts = (100, 180)
    strats = ("max", "min")
    data = np.genfromtxt(csvfilename, delimiter=",", skip_header=5)

    N = 10
    length = data.shape[0]//N

    for i in range(N):
        partial_data = data[length*i:length*(i+1) , :]
        # print(f'\n>>> "{csvfilename}" - {i+1}')
        Shape7[i] = lambda_handler(partial_data, ants, posts, strats)


    #######################################################################

    # the order of the flasks shapes will be as follows:
    #       1-7 (100 mL)
    #       2-3-6 (250 mL)
    #       4-5 (500 mL)


    overallData = np.concatenate((Shape1,Shape7,Shape2,Shape3,Shape6,Shape4,Shape5), axis=1)
    
    flask_names = np.array([["100mL Armut", "100mL kure", "250mL armut", "250mL ucgen", "250mL kure",
                            "500mL ucgen", "500mL kure"]]).reshape(7,1)
    
    
    means = np.mean(overallData, axis=0).reshape(7,1)
    stds = np.std(overallData, axis=0).reshape(7,1)
    err = stds / np.sqrt(10)

    print(f"Means (1,7,2,3,6,4,5):\n{means}")
    print(f"Stds (1,7,2,3,6,4,5):\n{stds}")
    mean_std = np.concatenate((means, stds, err), axis=1)


    np.savetxt("n_10 in-vitro results - python export.csv", overallData, delimiter=",", header="Shape1,Shape2,Shape3,Shape4,Shape5,Shape6,Shape7 [all values in mL]")
    np.savetxt("in-vitro mean_std_error - python export.csv", mean_std, delimiter=",", header="Mean,Std,standard error of mean")
