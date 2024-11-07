# import packages
import tkinter
from tkinter import *
from tkinter import messagebox
import pandas as pd

# get inst from connect gui
from connect_gui import use_inst
inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
path = get_folder_path()

# import needed functions
from measurement_functions import generatePulse

# gui window
window = tkinter.Tk()

# get variables
def get_vd():
    Vd = [float(vd_entry.get())]
    return Vd

def get_vg_init_list():
    Vg_init = []
    for vg_init in vg_init_entries:
        Vg_init.append(float(vg_init.get()))
    return Vg_init

def get_vg_final():
    Vg_final = [float(vg_final_entry.get())]
    return Vg_final

def set_settings():
    Vd = get_vd()
    Vg_init = get_vg_init_list()
    Vg_final = get_vg_final()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nVd: {Vd}"
                                  f"\nVg init list: {Vg_init}"
                                  f"\nVg final: {Vg_final}")

def run_rise_time_measurement_series():
    # get settings
    Vd = get_vd()
    Vg_init = get_vg_init_list()
    Vg_final = get_vg_final()

# vd entry
Label(window, text="Vd: ").grid(row=0, column=0)
vd_entry = Entry(window)
vd_entry.grid(row=0, column=1)

# create vg init list entries
Label(window, text="Vg Init List: ").grid(row=1, column=0)
vg_init_entries = []
for i in range(2):
    vg_init_entry = Entry(window)
    vg_init_entry.grid(row=1, column=i+1)
    vg_init_entries.append(vg_init_entry)

# vg final entry
Label(window, text="Vg Final: ").grid(row=2, column=0)
vg_final_entry = Entry(window)
vg_final_entry.grid(row=2, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=3, column=0)


# run gui
window.mainloop()