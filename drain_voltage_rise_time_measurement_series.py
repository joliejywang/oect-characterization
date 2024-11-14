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
from measurement_functions import generatePulse, runListSweep

# gui window
window = tkinter.Tk()

# get variables
def get_vd_final():
    Vd_final = []
    for vd_final in vd_final_entries:
        Vd_final.append(float(vd_final.get()))
    return Vd_final

def get_vg_list():
    Vg = []
    for vg in vg_entries:
        Vg.append(float(vg.get()))
    return Vg

def get_vd_init():
    Vd_init = float(vd_init_entry.get())
    return Vd_init

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
    Vd_final = get_vd_final()
    Vg = get_vg_list()
    Vd_init = get_vd_init()
    Id_range = get_id_range()
    Ig_range = get_ig_range()
    pulse_width = get_pulse_width()
    t_hold_before_meas = get_t_hold_before_meas()
    meas_time = get_meas_time()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nVd final: {Vd_final}"
                                  f"\nVg: {Vg}"
                                  f"\nVd init: {Vd_init}"
                                  f"\nId range: {Id_range}"
                                  f"\nIg range: {Ig_range}"
                                  f"\nPulse width: {pulse_width}"
                                  f"\nTime hold before measurement: {t_hold_before_meas}"
                                  f"\nMeasurement time: {meas_time}")

def run_rise_time_measurement_series():
    # get settings
    Vd_final = get_vd_final()
    Vg = get_vg_list()

    for i in range(0, len(Vg)):
        for j in range(0, len(Vd_final)):
            Vd_init = get_vd_init()
            Id_range = get_id_range()
            Ig_range = get_ig_range()
            pulse_width = get_pulse_width()  # time in seconds
            t_hold_before_meas = get_t_hold_before_meas()

            # ---------------------------------------------------------------------------------------------------

            meas_time = get_meas_time()  # Max sweep steps of 2500

            vd_list = generatePulse(Vd_init, Vd_final[j], pulse_width, meas_time)
            vg_list = "{:.1E}".format(Vg[i])

            # parameters format (list of strings, units of seconds, amps, volts):
            # [(0) measurment points,(1) time per point,(2) measurement delay,(3) aquisition time,
            #  (4) gate current range,(5) drain current range,(6) gate points,(7) drain points,
            #  (8) premeasurement voltage hold time]

            parameters = []
            parameters.append("{:.0f}".format(len(vd_list.split(','))))
            parameters.append("{:.1E}".format(meas_time))
            parameters.append("{:.1E}".format(meas_time / 2))
            parameters.append("{:.1E}".format((meas_time / 2) * 0.8))
            parameters.append("{:.0E}".format(Ig_range))
            parameters.append("{:.0E}".format(Id_range))
            parameters.append("{:.0f}".format(len(vg_list.split(','))))
            parameters.append("{:.0f}".format(len(vd_list.split(','))))
            parameters.append(t_hold_before_meas)

            pulse = runListSweep(parameters, vg_list, vd_list)

            # plotPulse(pulse)

            ###     Save last RISE TIME measurement     ###
            # Select folder and filename where the data will be stored

            # Make sure to change filename, otherwise it will overwrite previously saved meaurement
            filename = "Drain rise time Vg " + str(Vg[i]) + " Vd " + str(Vd_init) + " to " + str(Vd_final[j])

            # Create the path before trying to save
            # path = "C:\\Users\\skeen\\Desktop\\2022-02-17_MM461_1B\\"
            folder = "Device " + str(device) + "\\"

            # ------------------------------------------------------------------------------------------

            outp_data = pd.DataFrame()
            outp_data["Time (s)"] = pulse[:, 1]
            outp_data["Ids (A)"] = pulse[:, 2]
            outp_data["Igs (A)"] = pulse[:, 0]

            outp_data.to_csv(path + folder + filename)

    inst.write(':outp1 off')
    inst.write(':outp2 off')

# vd final entry
Label(window, text="Vd Final List: ").grid(row=0, column=0)
vd_final_entries = []
for i in range (5):
    vd_final_entry = Entry(window)
    vd_final_entry.grid(row=0, column=i+1)
    vd_final_entries.append(vd_final_entry)

# create vg init list entries
Label(window, text="Vg List: ").grid(row=1, column=0)
vg_entries = []
for i in range(13):
    vg_entry = Entry(window)
    vg_entry.grid(row=1, column=i+1)
    vg_entries.append(vg_entry)

# vd init entry
Label(window, text="Vd Init: ").grid(row=2, column=0)
vd_init_entry = Entry(window)
vd_init_entry.grid(row=2, column=1)

# id range entry
Label(window, text="Id range: ").grid(row=3, column=0)
id_range_entry = Entry(window)
id_range_entry.grid(row=3, column=1)

# ig range entry
Label(window, text="Ig range: ").grid(row=4, column=0)
ig_range_entry = Entry(window)
ig_range_entry.grid(row=4, column=1)

# pulse width entry
Label(window, text="Pulse width: ").grid(row=5, column=0)
pulse_width_entry = Entry(window)
pulse_width_entry.grid(row=5, column=1)

# time hold before measurement entry
Label(window, text="Time hold before measurement: ").grid(row=6, column=0)
t_hold_before_meas_entry = Entry(window)
t_hold_before_meas_entry.grid(row=6, column=1)

# measurement time entry
Label(window, text="Measurement time: ").grid(row=7, column=0)
meas_time_entry = Entry(window)
meas_time_entry.grid(row=7, column=1)

# show settings to check
Button(window, text="Show or Check Settings: ", command=set_settings).grid(row=8, column=0)

# button to run measurement
Button(window, text="Run Rise Time Measurement Series: ", command=run_rise_time_measurement_series).grid(row=9, column=0)

# run gui
window.mainloop()