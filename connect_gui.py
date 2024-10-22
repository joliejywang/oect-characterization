# import packages
import tkinter
from tkinter import *
from tkinter import messagebox

import pyvisa

window = tkinter.Tk()
rm = pyvisa.ResourceManager()

# view available resources
def available_resources():
    # clear previous lists before usage
    resource_list.delete(0, tkinter.END)

    # detect resources
    resource_tuple = rm.list_resources()

    # check if resources detected
    if len(resource_tuple) == 0:
        raise Exception("no resources detected")
    else:
        # create list for listbox of resources
        for resource in resource_tuple:
            resource_list.insert(END, resource)

def select_resource_from_list():
    chosen = resource_list.get(ACTIVE)
    inst = rm.open_resource(chosen)
    inst_id = inst.query("*IDN?")
    messagebox.showinfo("connected", f"connected to {inst_id}")

# gui parts
# resource list widget
resource_list = Listbox(window, width=25, height=10)
resource_list.pack(pady=5, padx=5)
tkinter.Button(window, text="List Resources", command=available_resources()).pack(pady=5)

# select resource
tkinter.Button(window, text="Connect to Instrument", command=select_resource_from_list).pack(pady=5)

# run gui
window.mainloop()