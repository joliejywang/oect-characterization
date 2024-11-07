# import packages
import tkinter
from tkinter import *
from tkinter import messagebox
import numpy as np
import pickle

# get inst from connect gui
from connect_gui import use_inst
inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
path = get_folder_path()

# import needed functions
from measurement_functions import runListSweep, generateFreqsV2, generateSineSweep

# gui window
window = tkinter.Tk()

# get variable values
def get_gate_voltage_list():
    gate_voltage_list = []
    for gate_voltage in gate_voltage_entries:
        gate_voltage_list.append(float(gate_voltage.get()))
    return gate_voltage_list

def get_drain_voltage_list():
    drain_voltage_list = [float(drain_voltage_entry.get())]
    return drain_voltage_list

def get_min_freq():
    min_freq = int(min_freq_entry.get())
    return min_freq

def get_max_freq():
    max_freq = int(max_freq_entry.get())
    return max_freq

def get_id_range():
    Id_range = float(id_range_entry.get())
    return Id_range

def get_ig_range():
    Ig_range = float(ig_range_entry.get())
    return Ig_range

def get_points_per_decade():
    points_per_decade = int(points_per_decade_entry.get())
    return points_per_decade

def get_cycles():
    cycles = int(cycles_entry.get())
    return cycles

def get_est_capacitance():
    Est_capacitance = float(est_capacitance_entry.get())
    return Est_capacitance

def get_max_meas_time():
    max_meas_time = float(max_meas_time_entry.get())
    return max_meas_time

def get_min_meas_time():
    min_meas_time = float(min_meas_time_entry.get())
    return min_meas_time

def get_vg_amplitude():
    Vg_amplitude = float(vg_amp_entry.get())
    return Vg_amplitude

def get_t_hold():
    t_hold_before_meas = float(t_hold_entry.get())
    return t_hold_before_meas

# get settings
def set_settings():
    gate_voltage_list = get_gate_voltage_list()
    drain_voltage_list = get_drain_voltage_list()
    min_freq = get_min_freq()
    max_freq = get_max_freq()
    Id_range = get_id_range()
    Ig_range = get_ig_range()
    points_per_decade = get_points_per_decade()
    cycles = get_cycles()
    Est_capacitance = get_est_capacitance()
    max_meas_time = get_max_meas_time()
    min_meas_time = get_min_meas_time()
    Vg_amplitude = get_vg_amplitude()
    t_hold_before_meas = get_t_hold()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nGate voltage list: {gate_voltage_list}"
                                  f"\nDrain voltage list: {drain_voltage_list}"
                                  f"\nMin freq: {min_freq}"
                                  f"\nMax freq: {max_freq}"
                                  f"\nId range: {Id_range}"
                                  f"\nIg range: {Ig_range}"
                                  f"\nPoints per decade: {points_per_decade}"
                                  f"\nCycles: {cycles}"
                                  f"\nEstimated capacitance: {Est_capacitance}"
                                  f"\nMax meas time: {max_meas_time}"
                                  f"\nMin meas time: {min_meas_time}"
                                  f"\nVg amplitude: {Vg_amplitude}"
                                  f"\nTime hold before meas: {t_hold_before_meas}")

# create gate voltage list entries
Label(window, text="Gate Voltage List: ").grid(row=0, column=0)
gate_voltage_entries = []
for i in range(2):
    gate_voltage_entry = Entry(window)
    gate_voltage_entry.grid(row=0, column=i+1)
    gate_voltage_entries.append(gate_voltage_entry)

# drain voltage list entry
Label(window, text="Drain Voltage List: ").grid(row=1, column=0)
drain_voltage_entry = Entry(window)
drain_voltage_entry.grid(row=1, column=1)

# minimum frequency entry
Label(window, text="Minimum frequency: ").grid(row=2, column=0)
min_freq_entry = Entry(window)
min_freq_entry.grid(row=2, column=1)

# maximum frequency entry
Label(window, text="Maximum frequency: ").grid(row=3, column=0)
max_freq_entry = Entry(window)
max_freq_entry.grid(row=3, column=1)

# id range entry
Label(window, text="Id Range: ").grid(row=4, column=0)
id_range_entry = Entry(window)
id_range_entry.grid(row=4, column=1)

# ig range entry
Label(window, text="Ig Range: ").grid(row=5, column=0)
ig_range_entry = Entry(window)
ig_range_entry.grid(row=5, column=1)

# points per decade entry
Label(window, text="Points Per Decade: ").grid(row=6, column=0)
points_per_decade_entry = Entry(window)
points_per_decade_entry.grid(row=6, column=1)

# cycles entry
Label(window, text="Cycles: ").grid(row=7, column=0)
cycles_entry = Entry(window)
cycles_entry.grid(row=7, column=1)

# estimated capacitance entry
Label(window, text="Capacitance Estimate (Better to overestimate): ").grid(row=8, column=0)
est_capacitance_entry = Entry(window)
est_capacitance_entry.grid(row=8, column=1)

# max measurement time entry
Label(window, text="Maximum Measurement Time (Should be 1e-3 or lower): ").grid(row=9, column=0)
max_meas_time_entry = Entry(window)
max_meas_time_entry.grid(row=9, column=1)

# min measurement time entry
Label(window, text="Minimum Measurement Time (Minimum is 2e-5): ").grid(row=10, column=0)
min_meas_time_entry = Entry(window)
min_meas_time_entry.grid(row=10, column=1)

# Vg amplitude entry
Label(window, text="Vg Amplitude: ").grid(row=11, column=0)
vg_amp_entry = Entry(window)
vg_amp_entry.grid(row=11, column=1)

# time hold before measurement entry
Label(window, text="Time Hold Before Measurement: ").grid(row=12, column=0)
t_hold_entry = Entry(window)
t_hold_entry.grid(row=12, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=11, column=0)

# run gui
window.mainloop()