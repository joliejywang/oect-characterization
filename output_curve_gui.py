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

# get individual variables
def get_vg_init():
    Vg_init = float(vg_init_entry.get())
    return Vg_init

def get_vg_final():
    Vg_final = float(vg_final_entry.get())
    return Vg_final

def get_vg_step():
    Vg_step = float(vg_step_entry.get())
    return Vg_step

def get_vd_init():
    Vd_init = float(vd_init_entry.get())
    return Vd_init

def get_vd_final():
    Vd_final = float(vd_final_entry.get())
    return Vd_final

def get_time_per_point():
    time_per_point = float(time_per_point_entry.get())
    return time_per_point

def get_points_per_sweep():
    points_per_sweep = int(points_per_sweep_entry.get())
    return points_per_sweep

def get_reverse_sweep():
    reverse_sweep = reverse_sweep_is_checked.get()
    return reverse_sweep

def get_ig_range():
    Ig_range = float(ig_range_entry.get())
    return Ig_range

def get_id_range():
    Id_range = float(id_range_entry.get())
    return Id_range

def get_t_hold_before_meas():
    t_hold_before_meas = float(time_hold_entry.get())
    return t_hold_before_meas

# check if settings are correct
def set_settings():
    # get entered labels
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    Vg_step = get_vg_step()
    Vd_init = get_vd_init()
    Vd_final = get_vd_final()
    time_per_point = get_time_per_point()
    points_per_sweep = get_points_per_sweep()
    reverse_sweep = get_reverse_sweep()
    Ig_range = get_ig_range()
    Id_range = get_id_range()
    t_hold_before_meas = get_t_hold_before_meas()

    messagebox.showinfo("Inputs", f"Inputs Entered:"
                                  f"\nVg init: {Vg_init}"
                                  f"\nVg final: {Vg_final}"
                                  f"\nVg step: {Vg_step}"
                                  f"\nVd init: {Vd_init}"
                                  f"\nVd final: {Vd_final}"
                                  f"\nTime per point: {time_per_point}"
                                  f"\nPoints per sweep: {points_per_sweep}"
                                  f"\nReverse sweep: {reverse_sweep}"
                                  f"\nIg range: {Ig_range}"
                                  f"\nId range: {Id_range}"
                                  f"\nTime Hold: {t_hold_before_meas}")

def run_output_curves():
    # initialize global variables
    global output, Vd_output, Vg_values_output

    # get entered labels
    Vg_init = get_vg_init()
    Vg_final = get_vg_final()
    Vg_step = get_vg_step()
    Vd_init = get_vd_init()
    Vd_final = get_vd_final()
    time_per_point = get_time_per_point()
    points_per_sweep = get_points_per_sweep()
    reverse_sweep = get_reverse_sweep()
    Ig_range = get_ig_range()
    Id_range = get_id_range()
    t_hold_before_meas = get_t_hold_before_meas()

    # Generate array of gate voltage values based on user input
    Vg_values_output = np.empty((round((Vg_final - Vg_init) / Vg_step) + 1,))
    for i in range(0, len(Vg_values_output)):
        if Vg_init > Vg_final:
            step = -1 * abs(Vg_step)
        else:
            step = abs(Vg_step)
        Vg_values_output[i] = round(Vg_init + step * i, 4)

    # Generate sweep for drain voltage based on user input
    vd_list = generateLinSweep(points_per_sweep, Vd_init, Vd_final, reverse=reverse_sweep)
    vg_list = "{:.1E}".format(Vg_values_output[0])

    # Create measurement parameters
    # parameters format (list of strings, units of seconds, amps, volts):
    # [(0) measurment points,(1) time per point,(2) measurement delay,(3) aquisition time,
    #  (4) gate current range,(5) drain current range,(6) gate points,(7) drain points,
    #  (8) pre-measurement voltage hold time]
    parameters = []
    parameters.append("{:.0f}".format(len(vd_list.split(','))))  # 0
    parameters.append("{:.1E}".format(time_per_point))  # 1
    parameters.append("{:.1E}".format(time_per_point / 2))  # 2
    parameters.append("{:.1E}".format((time_per_point / 2) * 0.8))  # 3
    parameters.append("{:.0E}".format(Ig_range))  # 4
    parameters.append("{:.0E}".format(Id_range))  # 5
    parameters.append("1")  # 6
    parameters.append("{:.0f}".format(len(vd_list.split(','))))  # 7
    parameters.append(t_hold_before_meas)

    # Create empty array to store output data
    output = np.empty((len(Vg_values_output), int(parameters[0]), 4))

    # Run output curve for each gate voltage and plot results
    output[0] = runListSweep(parameters, vg_list, vd_list)
    Vd_output = np.asarray([float(i) for i in vd_list.split(',')])
    plt.plot(Vd_output, output[0, :, 2], label="$V_{GS}$ = " + str(Vg_values_output[0]))

    for i in range(1, len(output)):
        vg_list = "{:.1E}".format(Vg_values_output[i])
        output[i] = runListSweep(parameters, vg_list, vd_list)
        plt.plot(Vd_output, output[i, :, 2], label="$V_{GS}$ = " + str(Vg_values_output[i]))

    plt.gca().invert_yaxis()
    plt.xlabel("$V_{DS}$ (V)")
    plt.ylabel("$I_{DS}$ (A)")
    plt.title("Output curves")
    plt.legend()

    # Turn instrument output off after measurement
    inst.write(":outp1 off")
    inst.write(":outp2 off")

def save_output_curve_measurement():
    # initialize global variables
    global output, Vd_output, Vg_values_output

    # get filename input
    filename = filename_entry.get()

    # Create the path before trying to save
    # path = "C:\\Users\\skeen\\Dropbox\\000_Postdoctoral_Research\\Voltage-dependent mobility with Sophia\\OECT Data\\2023-02-20_pg1T2-g5T2_100mM_NaCl\\"
    folder = polymer + "\\" + device + "\\"

    outp_data = pd.DataFrame()
    outp_data["Time (s)"] = output[0, :, 1]
    outp_data["Vds (V)"] = Vd_output
    for i in range(0, len(output)):
        outp_data["Ids (A) at Vg = " + str(Vg_values_output[i]) + " V"] = output[i, :, 2]
        outp_data["Igs (A) at Vg = " + str(Vg_values_output[i]) + " V"] = output[i, :, 0]

    outp_data.to_csv(path + folder + filename)

# create widgets
# create input fields for setting
Label(window, text="Initial Vg: ").grid(row=0, column=0)
vg_init_entry = Entry(window)
vg_init_entry.grid(row=0, column=1)

Label(window, text="Final Vg: ").grid(row=1, column=0)
vg_final_entry = Entry(window)
vg_final_entry.grid(row=1, column=1)

Label(window, text="Vg Step: ").grid(row=2, column=0)
vg_step_entry = Entry(window)
vg_step_entry.grid(row=2, column=1)

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