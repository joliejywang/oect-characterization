# import packages
import tkinter
from tkinter import *

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