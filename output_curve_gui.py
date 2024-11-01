# import packages
import tkinter
from tkinter import *
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# get inst from connect gui
from connect_gui import use_inst
inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
path = get_folder_path()

# import needed functions from measurement
from measurement_functions import generateLinSweep, runListSweep

# gui window
window = tkinter.Tk()

# create widgets
# create input fields for setting
Label(window, text="Initial Vg: ").grid(row=0, column=0)
vg_init_entry = Entry(window)
vg_init_entry.grid(row=0, column=1)

Label(window, text="Final Vg: ").grid(row=1, column=0)
vg_final_entry = Entry(window)
vg_final_entry.grid(row=1, column=1)

Label(window, text="Vg Step: ").grid(row=2, column=0)
vd_step_entry = Entry(window)
vd_step_entry.grid(row=2, column=1)

Label(window, text="Initial Vd: ").grid(row=3, column=0)
vd_init_entry = Entry(window)
vd_init_entry.grid(row=3, column=1)

Label(window, text="Final Vd: ").grid(row=4, column=0)
vd_final_entry = Entry(window)
vd_final_entry.grid(row=4, column=1)

Label(window, text="Time Per Point: ").grid(row=5, column=0)
time_per_point_entry = Entry(window)
time_per_point_entry.grid(row=5, column=1)

Label(window, text="Points Per Sweep: ").grid(row=6, column=0)
points_per_sweep_entry = Entry(window)
points_per_sweep_entry.grid(row=6, column=1)

reverse_sweep_is_checked = BooleanVar()
Checkbutton(window, text="Reverse Sweep: ", variable=reverse_sweep_is_checked).grid(row=7, columnspan=2)

Label(window, text="Ig Range (in 1e-X format with X between 0 and 7 inclusive): ").grid(row=8, column=0)
ig_range_entry = Entry(window)
ig_range_entry.grid(row=8, column=1)

Label(window, text="Id Range (in 1e-X format with X between 0 and 7 inclusive): ").grid(row=9, column=0)
id_range_entry = Entry(window)
id_range_entry.grid(row=9, column=1)

Label(window, text="Time Hold Before Measurement: ").grid(row=10, column=0)
time_hold_entry = Entry(window)
time_hold_entry.grid(row=10, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=11, column=0)

# run transfer curve button
Button(window, text="Run Transfer Curves: ", command=run_output_curves).grid(row=12, column=0)

# enter filename to save as
Label(window, text="Filename to Save As: ").grid(row=13, column=0)
filename_entry = Entry(window)
filename_entry.grid(row=13, column=1)

# save last transfer measurement
Button(window, text="Save Last Transfer Measurement: ", command=save_output_curve_measurement).grid(row=14, column=0)

# run gui
window.mainloop()