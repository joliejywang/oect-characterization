# import packages
import numpy as np

def sine(point, offset, amp, points_per_cycle):
    return offset + amp * np.sin(2 * np.pi * (point % points_per_cycle) / points_per_cycle)

# Function for running fits of measured data to get ∆Ig and ∆Id
def fitSineV3(t, amp, phi, I0, f_freq):
    return I0 + amp * np.sin(2 * np.pi * f_freq * t + phi)

# Function for running fits of measured data to get ∆Ig and ∆Id
# def sinFit(t,amp,phi,I0):
    # return I0 + amp*np.sin(2*np.pi*f_freq*t + phi)

