# import packages
import tkinter
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import needed functions from measurement
from measurement_functions import generateLinSweep, runListSweep

# get inst from connect gui
from connect_gui import use_inst
inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
folder = get_folder_path()

# set up gui window
window = tkinter.Tk()

def set_settings():

def run_transfer_curves():

def transfer_measurement():

# create widgets
# create input fields for setting

# run gui
window.mainloop()