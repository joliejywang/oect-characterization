# import packages
import tkinter
from tkinter import *
from tkinter import messagebox

# plot function
from plot import plotPulse

# measurement function
from measurement_functions import generatePulse, runListSweep

# get inst from connect gui
from connect_gui import use_inst
from transfer_curve_gui import vd_init_entry

inst = use_inst()

# get polymer, device, path from input gui
from polymer_device_path_input_gui import get_polymer, get_device, get_folder_path
polymer = get_polymer()
device = get_device()
path = get_folder_path()

# gui window
window = tkinter.Tk()

# get entry values
def get_vd():
    Vd = float(vd_entry.get())
    return Vd

def get_vg_init():
    Vg_init = float(vg_init_entry.get())
    return Vg_init

def get_vg_final():
    Vg_final = float(vg_final_entry.get())
    return Vg_final

def get_id_range():
    Id_range = float(id_range_entry.get())
    return Id_range

def get_ig_range():
    Ig_range = float(ig_range_entry.get())
    return Ig_range

def get_pulse_width():
    pulse_width = float(pulse_width_entry.get())
    return pulse_width

def get_t_hold_before_meas():
    t_hold_before_meas = float(t_hold_before_meas_entry.get())
    return t_hold_before_meas

def get_meas_time():
    meas_time = float(meas_time_entry.get())
    return meas_time

def set_settings():
    Vd = get_vd()
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    Id_range = get_id_range()
    Ig_range = get_ig_range()
    pulse_width = get_pulse_width()
    t_hold_before_meas = get_t_hold_before_meas()
    meas_time = get_meas_time()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nVd: {Vd}"
                                  f"\nVg init: {Vg_init}"
                                  f"\nVg final: {Vg_final}"
                                  f"\nId range: {Id_range}"
                                  f"\nIg range: {Ig_range}"
                                  f"\nPulse width: {pulse_width}"
                                  f"\nTime hold before measurement: {t_hold_before_meas}"
                                  f"\nMeasurement time: {meas_time}")

def run_rise_time_measurement():
    # get variables
    Vd = get_vd()
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    Id_range = get_id_range()
    Ig_range = get_ig_range()
    pulse_width = get_pulse_width()
    t_hold_before_meas = get_t_hold_before_meas()
    meas_time = get_meas_time()

    vg_list = generatePulse(Vg_init, Vg_final, pulse_width, meas_time)
    vd_list = "{:.1E}".format(Vd)

    # parameters format (list of strings, units of seconds, amps, volts):
    # [(0) measurment points,(1) time per point,(2) measurement delay,(3) aquisition time,
    #  (4) gate current range,(5) drain current range,(6) gate points,(7) drain points,
    #  (8) premeasurement voltage hold time]

    parameters = []
    parameters.append("{:.0f}".format(len(vg_list.split(','))))
    parameters.append("{:.1E}".format(meas_time))
    parameters.append("{:.1E}".format(meas_time / 2))
    parameters.append("{:.1E}".format((meas_time / 2) * 0.8))
    parameters.append("{:.0E}".format(Ig_range))
    parameters.append("{:.0E}".format(Id_range))
    parameters.append("{:.0f}".format(len(vg_list.split(','))))
    parameters.append("1")
    parameters.append(t_hold_before_meas)

    pulse = runListSweep(parameters, vg_list, vd_list)

    inst.write(':outp1 off')
    inst.write(':outp2 off')

    plotPulse(pulse)


# vg final entry
Label(window, text="Vd: ").grid(row=0, column=0)
vd_entry = Entry(window)
vd_entry.grid(row=0, column=1)

Label(window, text="Vg init: ").grid(row=1, column=0)
vg_init_entry = Entry(window)
vg_init_entry.grid(row=1, column=1)

Label(window, text="Vg final: ").grid(row=2, column=0)
vg_final_entry = Entry(window)
vg_final_entry.grid(row=2, column=1)

Label(window, text="Id range: ").grid(row=3, column=0)
id_range_entry = Entry(window)
id_range_entry.grid(row=3, column=1)

Label(window, text="Ig range: ").grid(row=4, column=0)
ig_range_entry = Entry(window)
ig_range_entry.grid(row=4, column=1)

Label(window, text="Pulse width: ").grid(row=5, column=0)
pulse_width_entry = Entry(window)
pulse_width_entry.grid(row=5, column=1)

Label(window, text="Time hold before measurement: ").grid(row=6, column=0)
t_hold_before_meas_entry = Entry(window)
t_hold_before_meas_entry.grid(row=6, column=1)

Label(window, text="Measurement time: ").grid(row=7, column=0)
meas_time_entry = Entry(window)
meas_time_entry.grid(row=7, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=8, column=0)

# button to run measurement
Button(window, text="Run Rise Time Measurement: ", command=run_rise_time_measurement).grid(row=9, column=0)

# run gui
window.mainloop()