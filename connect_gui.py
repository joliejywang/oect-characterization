# import packages
import tkinter
from tkinter import *
from tkinter import messagebox

import pyvisa

window = tkinter.Tk()
rm = pyvisa.ResourceManager()
# inst = None

# view available resources
def available_resources():
    # clear previous lists before usage
    resource_list.delete(0, tkinter.END)

    # detect resources
    resource_tuple = rm.list_resources()

    # # check if resources detected
    # if len(resource_tuple) == 0:
    #     raise Exception("no resources detected")
    # else:
    #     # create list for listbox of resources
    #     for resource in resource_tuple:
    #         resource_list.insert(END, resource)
    for resource in resource_tuple:
        resource_list.insert(END, resource)

def select_resource_from_list():
    chosen = resource_list.get(ACTIVE)
    global inst
    inst = rm.open_resource(chosen)
    inst_id = inst.query("*IDN?")
    messagebox.showinfo("connected", f"connected to {inst_id}")

# for other gui that need to use inst that is connected to import into file
def use_inst():
    global inst
    return inst

def process_timeout_value():
    timeout_value = float(time_val.get())
    global inst
    inst.timeout = timeout_value

# gui parts
# resource list widget
resource_list = Listbox(window, width=25, height=10)
resource_list.pack(pady=5, padx=5)
tkinter.Button(window, text="List Resources", command=available_resources()).pack(pady=5)

# select resource
tkinter.Button(window, text="Connect to Instrument", command=select_resource_from_list).pack(pady=5)

# insert value for time out
tkinter.Label(window, text="Set the instrument timeout value in ms: ").pack()
time_val = Entry(window, width=20)
time_val.pack()
tkinter.Button(window, text="Set timeout value", command=process_timeout_value).pack(pady=5)

# run gui
window.mainloop()