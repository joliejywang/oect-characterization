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

    # backend code
    # Create the path before trying to save
    # path = "C:\\Users\\skeen\\Dropbox\\000_Postdoctoral_Research\\Voltage-dependent mobility with Sophia\\OECT Data\\2023-02-20_pg1T2-g5T2_100mM_NaCl\\"
    folder = polymer + "\\" + device + "\\"

    for i in range(0, len(Vd)):
        for j in range(0, len(Vg_init)):
            for k in range(0, len(Vg_final)):
                # Vg_final = -0.65
                Id_range = 1e-2
                Ig_range = 1e-4
                pulse_width = 0.5  # time in seconds
                t_hold_before_meas = 0.5

                # ---------------------------------------------------------------------------------------------------

                meas_time = 2e-4  # Max sweep steps of 2500

                vg_list = generatePulse(Vg_init[j], Vg_final[k], pulse_width, meas_time)
                vd_list = "{:.1E}".format(Vd[i])

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

                filename = "Rise time Vg " + "{:.2f}".format(Vg_init[j]) + " to " + "{:.2f}".format(
                    Vg_final[k]) + " Vd " + str(Vd[i])

                outp_data = pd.DataFrame()
                outp_data["Time (s)"] = pulse[:, 1]
                outp_data["Vds (V)"] = Vd[i]
                outp_data["Ids (A)"] = pulse[:, 2]
                outp_data["Igs (A)"] = pulse[:, 0]

                outp_data.to_csv(path + folder + "Rise time\\" + filename)

    inst.write(':outp1 off')
    inst.write(':outp2 off')




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