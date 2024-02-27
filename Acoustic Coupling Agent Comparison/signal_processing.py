import numpy as np 
from scipy.signal import butter, filtfilt, hilbert
from scipy.optimize import minimize


def tof_to_distance_mm(time_us):
    return 1e3 * 1480*(time_us*1e-6)/2




def Sphere_func(parameters : np.array, coordinates: np.array) -> float:
    """(Fitted-actual) error function to be minimized."""
    x0, y0, z0, r = parameters
    x, y, z = np.transpose(coordinates)
    return np.sum((np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) - r)**2) 




def calculate_the_volume(coordinates: np.array) -> float:
    """Used to calculate the volume, given a coordinates array in (x,y,z) system."""
    x, y, z = np.transpose(coordinates)
    initial_guess = np.array([np.mean(x), np.mean(y), np.mean(z), 50]) # np.std(x) as r
    r_limit = (1e6*(3/4)/np.pi)**(1/3)
    best_fit_sphere = minimize(Sphere_func, initial_guess, coordinates, bounds=[(-15,15), (-15,15), (10,170), (5, r_limit)])
    x0, y0, z0, r = best_fit_sphere.x
    return 4/3 * np.pi * r**3 / 1000, r, (np.round(x0,2),np.round(y0,2),np.round(z0,2))



def amplify_signal(voltage, gain=5):
    max_amp = np.max(voltage)
    amplified = voltage**gain
    amplified = (amplified - np.min(amplified)) / (np.max(amplified)-np.min(amplified))
    return max_amp*(amplified- np.mean(amplified))




def apply_bandpass_filter(voltage, cutoffs, order, sampling_rate):
	"""Applies a band-pass filter to the inputted array."""
	b,a = butter(order, [cutoffs[0], cutoffs[1]], btype='bandpass', analog=False, output='ba', fs=sampling_rate)
	return filtfilt(b, a, voltage) 




def get_envelope(voltage):
	"""Takes the 1-sided envelope of the inputted signal using Hilbert transform."""
	return np.abs(hilbert(voltage))




def detect_echo_peaks(time:np.array, voltage:np.array, anterior_limits, posterior_limits):
    """Used to detect the anterior and posterior echo locations in a given sectioned data."""
    ant_idx_start = np.argmin( np.abs( time - anterior_limits[0] ) ) # anterior start
    ant_idx_end = np.argmin( np.abs( time - anterior_limits[1] ) ) # anterior stop
    pos_idx_start = np.argmin( np.abs( time - posterior_limits[0] ) ) # posterior start
    pos_idx_end = np.argmin( np.abs( time - posterior_limits[1] ) ) # posterior stop
    # Search for where the maximum voltage occurs in the given intervals
    anterior_peak_index = ant_idx_start+np.argmin( np.abs( voltage[ant_idx_start:ant_idx_end] - np.max(voltage[ant_idx_start:ant_idx_end]) ) )
    posterior_peak_index = pos_idx_start+np.argmin( np.abs( voltage[pos_idx_start:pos_idx_end] - np.max(voltage[pos_idx_start:pos_idx_end]) ) )
    return np.array([time[anterior_peak_index],time[posterior_peak_index]]), np.array([voltage[anterior_peak_index],voltage[posterior_peak_index]]) 




