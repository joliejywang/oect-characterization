# import needed packages
import tkinter
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox

window = tkinter.Tk()

# declare global variables
# folder_path = ""
# polymer = ""
# device = ""

def fill_path():
    global folder_path
    folder_path = filedialog.askdirectory(initialdir="C:/")

def show_inputs():
    global polymer
    polymer = polymer_entry.get()
    global device
    device = device_entry.get()
    global folder_path

    # throw exception if a value is missing out of three inputs
    if (len(polymer) == 0) or (len(device) == 0) or (len(folder_path) == 0):
        raise Exception("WARNING: missing at least one input")

    messagebox.showinfo("Inputs", f"Inputs Entered:\nPolymer: {polymer}\nDevice: {device}\nPath: {folder_path}")

def get_polymer():
    global polymer
    return polymer

def get_device():
    global device
    return device

def get_folder_path():
    global folder_path
    return folder_path

# create widgets
# polymer entry
Label(window, text="Polymer: ").grid(row=0, column=0)
polymer_entry = Entry(window)
polymer_entry.grid(row=0, column=1)

# device entry
Label(window, text="Device: ").grid(row=1, column=0)
device_entry = Entry(window)
device_entry.grid(row=1, column=1)

# path entry
Label(window, text="Path: ").grid(row=2, column=0)
Button(window, text="select path", command=fill_path).grid(row=2, column=1)

# show
Button(window, text="Show Inputs", command=show_inputs).grid(row=3, column=0)

# run gui
window.mainloop()

