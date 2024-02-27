"""Contains the functions used in the lambda_function file. """


from sys import exit
from os import system
import numpy as np
from scipy.optimize import minimize,differential_evolution
from datetime import datetime



def now() -> str:
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")



def Sphere_func(parameters : np.array, coordinates: np.array) -> float:
    # (Fitted-actual) error function to be minimized 
    x0, y0, z0, r = parameters
    x, y, z = np.transpose(coordinates)
    return np.sum((np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) - r)**2) 



def Calc_Volume_1(coordinates: np.array) -> float:
    # input coordinates: columns are x, y and z. Returns the volume
    x, y, z = np.transpose(coordinates)
    initial_guess = np.array([np.mean(x), np.mean(y), np.mean(z), 50]) # np.std(x) as r
    r_limit = (1e6*(3/4)/np.pi)**(1/3)
    best_fit_sphere = minimize(Sphere_func, initial_guess, coordinates, bounds=[(-20,20), (-20,20), (5,170), (5, r_limit)])
    x0, y0, z0, r = best_fit_sphere.x
    return 4/3 * np.pi * r**3 / 1000, r, (x0,y0,z0)




def noise_free(masked_data: np.array) -> np.array:
    CLOCK_PERIOD = 0.7 # normally it is 0.5 microseconds, set to 0.7 to be safe
    #   Input array is the timestamp array containing the timestamps within the
    #   upper and lower bounds. But, to be counted as a legitimate echo timestamp,
    #   there needs to be at least two timestamps within a time interval, equal to
    #   the clock period. If a timestamp is alone then that timestamp will be disregarded.

    #   mask1: element True, if following element is within CLOCK_PERIOD
    mask1 = np.append(np.abs(masked_data[:-1]-masked_data[1:]) < CLOCK_PERIOD, False)
    #   mask2: element True, if previous element is within CLOCK_PERIOD
    mask2 = np.insert(np.abs(masked_data[1:]-masked_data[:-1]) < CLOCK_PERIOD, 0, False)
    combined_mask = mask1 | mask2

    #combined_mask = np.zeros(masked_data.size)
    #for j in range(masked_data.size):
    #    differences = np.abs(masked_data - masked_data[j])
    #    combined_mask[j] = np.any((differences > 0.3) & (differences < 0.7))
    
    return masked_data[combined_mask]
    #return masked_data




def evaluate_data(data: dict) -> dict:
    # # lower and upper bounds for wall locations [us]
    ant_lim_low = 15
    ant_lim_high = 40
    post_lim_low = 105
    post_lim_high = 120
    


    refined_data = {'1':np.empty(0) , '2':np.empty(0) , '3':np.empty(0) , '4':np.empty(0)}
    
    data_is_enough = False
    counter = 0
    for key in data.keys():
        ant_mask = (data[key]>ant_lim_low) & (data[key]<ant_lim_high)
        post_mask = (data[key]>post_lim_low) & (data[key]<post_lim_high)
        anteriors = data[key] [ant_mask] # within bounds
        posteriors = data[key] [post_mask] # within bounds

        
        # ANTERIOR STRATEGY
        if anteriors.size > 0:
            noise_free_anteriors = noise_free(anteriors) # eliminate noise peaks
            if noise_free_anteriors.size > 0:
                counter += 1
                refined_data[key] = np.append(refined_data[key],np.mean(noise_free_anteriors))
        
        # POSTERIOR STRATEGY
        if posteriors.size > 0:
            noise_free_posteriors = noise_free(posteriors)
            if noise_free_posteriors.size > 0:
                counter += 1
                refined_data[key] = np.append(refined_data[key],np.mean(noise_free_posteriors)) 
  
  
  
        # FUTURE WORK: Here, anterior and posterior timestamps can be seperated!
    print(f'$$$ # of coordinates: {counter}')
    if counter >= 4: data_is_enough = True
    return data_is_enough, refined_data






def ProcessData(csvFileName, csvDATA: np.array):
    #   Each CSV data file will come from a single monitoring action.
    # Therefore, we don't need to compute multiple volumes even if there are
    # multiple complete data packages (multiple 1-2-3-4 full recordings).
    # We can treat them as a single monitoring action and take average which will
    # increase our calculation accuracy.
    print(f'\n$$$ Processing {csvFileName} data (Length:{csvDATA.shape[0]})...')
    __Vwater = 1420 # m/s
    XY = {'1': [0,0] , '2':[13,0], '3':[0,13] , '4':[13,13]}
    
    data = {'1': np.empty((0)) , '2':np.empty((0)), '3':np.empty((0)) , '4':np.empty((0))}

    piezonumbers = csvDATA[:,0]
    timestamps = csvDATA[:,1]

    mask1 = (piezonumbers == 1)
    data['1'] = np.append(data['1'], timestamps[mask1])
    mask2 = (piezonumbers == 2)
    data['2'] = np.append(data['2'], timestamps[mask2])
    mask3 = (piezonumbers == 3)
    data['3'] = np.append(data['3'], timestamps[mask3])
    mask4 = (piezonumbers == 4)
    data['4'] = np.append(data['4'], timestamps[mask4])

    data_is_enough, refined_data = evaluate_data(data)
    if data_is_enough:
        print(f'$$$ {now()} - Successful monitoring session')
        # form x,y,z array from timestamps (columns are x,y,z)
        coordinates = np.empty((0,3))
        for key in refined_data.keys():
            for timestamp in refined_data[key]:
                dist = __Vwater * timestamp * 1E-3 / 2 # in mm
                coordinates = np.append(coordinates, np.array([[XY[key][0],XY[key][1],dist]]),axis=0)           
        # feed it to the function and calculate the volume
        volume, r, C = Calc_Volume_1(coordinates)
        return True, round(volume), r, C, coordinates, refined_data
    else:
        print(f'$$$ {now()} - Not enough data points!')
        volume = 'Empty'
        return False, None, None, None, None, None

