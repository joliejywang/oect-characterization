import numpy as np
import time

# import functions from math
from math_functions import sine

# get inst value from connection
# import connect
# inst = connect.inst

from connect_gui import use_inst
inst = use_inst()

def runListSweep(parameters, ch1_list, ch2_list, wait=True):
    # parameters format (list of strings, units of seconds, amps, volts):
    # [(0) measurment points,(1) time per point,(2) measurement delay,(3) aquisition time,
    #  (4) gate current range,(5) drain current range,(6) gate points,(7) drain points,
    #  (8) premeasurement voltage hold time]

    # Trigger TIM minimum of 2E-5
    # Acquisition time minimum of 8E-6

    # acq delay should be on the order of 5E-5 or greater
    # Measurement time + acq delay < trigger interval by at least 1E-5

    inst.write(":sour1:volt:lev:imm " + ch1_list.split(',')[0])
    inst.write(":sour2:volt:lev:imm " + ch2_list.split(',')[0])

    inst.write(":outp1 on")
    inst.write(":outp2 on")

    if wait:
        time.sleep(parameters[8])

    # Sets the measurement list of voltages for the gate (channel 1)
    inst.write(":sour1:func:mode volt")
    inst.write(":sour1:volt:mode list")
    inst.write(":sour1:list:volt " + ch1_list)
    # Sets the measurement list of voltages for the drain (channel 2)
    inst.write(":sour2:func:mode volt")
    inst.write(":sour2:volt:mode list")
    inst.write(":sour2:list:volt " + ch2_list)

    # Set range and interval for measurement
    inst.write(":sens1:func \"curr\"")
    inst.write(":sens1:curr:rang:auto off")
    inst.write(":sens1:curr:prot " + parameters[4])
    inst.write(":sens1:curr:rang " + parameters[4])
    inst.write(":sens2:func \"curr\"")

    inst.write(":sens2:curr:rang:auto off")
    inst.write(":sens2:curr:rang " + parameters[5])
    inst.write(":sens2:curr:prot " + parameters[5])

    # Source output ranging set to fixed mode
    inst.write(":sour1:volt:rang:auto off")
    inst.write(":sour1:volt:rang 2")
    inst.write(":sour2:volt:rang:auto off")
    inst.write(":sour2:volt:rang 2")
    # Mesurement wait time set to OFF
    inst.write(":sens1:wait off")
    inst.write(":sour1:wait off")
    inst.write(":sens2:wait off")
    inst.write(":sour2:wait off")

    # Set trigger source to the same mode
    inst.write(":trig1:sour tim")
    inst.write("trig1:tim " + parameters[1])
    inst.write("trig1:acq:coun " + parameters[0])
    inst.write(":trig1:acq:del " + parameters[2])
    inst.write("trig1:tran:coun " + parameters[6])
    inst.write(":trig1:tran:del 0")
    inst.write(":trig2:sour tim")
    inst.write("trig2:tim " + parameters[1])
    inst.write("trig2:acq:coun " + parameters[0])
    inst.write(":trig2:acq:del " + parameters[2])
    inst.write("trig2:tran:coun " + parameters[7])
    inst.write(":trig2:tran:del 0")

    # Measurement interval is set to the same value
    inst.write(":sens1:curr:aper " + parameters[3])
    inst.write(":sens2:curr:aper " + parameters[3])

    inst.write(":form:elem:sens curr,time")

    # Runs the measurement
    inst.write(":init (@1,2)")
    # Fetches the measurement data
    t1 = time.time()
    data_out = inst.query(":fetc:arr? (@1,2)")
    t2 = time.time()
    # print(t1-t2)
    # Convert string of data to numpy array
    data = np.asarray([float(i) for i in data_out.split(',')])
    # Return the measurement results

    # inst.write(":syst:beep 800,0.25")

    data = np.reshape(data, (int(parameters[0]), 4))

    return data


# def pulsing(parameters):
    # parameters = [Vd, Vg_low, Vg_high, time per pulse, duty cycle, number of pulses]
    # return null


def findIdRange(Vd, Vg, Ig_range):
    inst.write(":sour1:volt:lev:imm " + "{:.3E}".format(Vg))
    inst.write(":sour2:volt:lev:imm " + "{:.3E}".format(Vd))

    inst.write(":outp1 on")
    inst.write(":outp2 on")

    time.sleep(10)

    inst.write(":sens2:curr:aper:auto on")
    inst.write(":sens2:curr:prot 1")
    inst.write(":sens2:func \"curr\"")

    Id_range = 1e-7
    overload = True

    while overload:
        inst.write(":sens2:curr:rang " + '{:.0E}'.format(Id_range))
        data_out = float(inst.query(":meas:curr? (@2)"))
        if data_out == np.nan or abs(data_out) > 0.8 * Id_range:
            Id_range = Id_range * 10
        else:
            overload = False

    return Id_range


def generateSineSweep(cycles, points_per_cycle, amp, offset):
    sweep = ""
    for i in range(0, int(points_per_cycle * cycles)):
        sweep += "{:.3E}".format(round(sine(i, offset, amp, points_per_cycle), 5)) + ","
    sweep += "{:.3E}".format(round(sine(int(points_per_cycle * cycles), offset, amp, points_per_cycle), 5))
    # sweep += "{:.3E}".format(round(sine(i + 1, offset, amp, points_per_cycle), 5))
    return sweep


def generateLinSweep(points_per_sweep, v_start, v_end, reverse=False):
    sweep = ""
    for i in range(0, points_per_sweep):
        sweep += "{:.3E}".format(round(v_start + (i / points_per_sweep) * (v_end - v_start), 5)) + ","
    # sweep += "{:.3E}".format(round(v_start + ((i + 1) / points_per_sweep) * (v_end - v_start), 5))
    sweep += "{:.3E}".format(round(v_start + (points_per_sweep / points_per_sweep) * (v_end - v_start), 5))
    if reverse:
        sweep += ","
        for j in range(1, points_per_sweep):
            i = points_per_sweep - j
            sweep += "{:.3E}".format(round(v_start + (i / points_per_sweep) * (v_end - v_start), 5)) + ","
        sweep += "{:.3E}".format(round(v_start, 5))
    return sweep


# def generateFreqs(start_freq, end_freq, points_per_dec):
#     decs = np.subtract(np.log10(end_freq), np.log10(start_freq))
#     freqs = np.logspace(np.log10(start_freq), np.log10(end_freq), int(decs * points_per_dec))
#     meas_time = "{:.1E}".format(1 / (12 * freqs[0]))
#     for i in range(1, len(freqs)):
#         meas_time += "," + "{:.1E}".format(1 / (12 * freqs[i]))
#     meas_time = np.asarray([float(i) for i in meas_time.split(',')])
#     freqs = 1 / (12 * meas_time)
#     return freqs, meas_time


def generateFreqsV2(start_freq, end_freq, points_per_dec, max_meas_time=1e-3, min_meas_time=2e-5):
    decs = np.subtract(np.log10(end_freq), np.log10(start_freq))
    freqs = np.logspace(np.log10(start_freq), np.log10(end_freq), int(decs * points_per_dec))
    meas_time = np.empty((len(freqs),))
    points_per_cycle = np.empty((len(freqs),))
    meas_time[0] = max_meas_time
    for i in range(0, len(freqs)):
        if i > 0:
            meas_time[i] = meas_time[i - 1]
        points_per_cycle[i] = np.round(1 / (freqs[i] * meas_time[i]))
        while points_per_cycle[i] < 100 and meas_time[i] > min_meas_time:
            meas_time[i] = float("{:.1E}".format(meas_time[i] / 2))
            points_per_cycle[i] = np.round(1 / (freqs[i] * meas_time[i]))
        if meas_time[i] < min_meas_time:
            meas_time[i] = min_meas_time
            points_per_cycle[i] = np.round(1 / (freqs[i] * meas_time[i]))
        freqs[i] = 1 / (meas_time[i] * points_per_cycle[i])
    return freqs, points_per_cycle, meas_time


def generatePulse(Vi, Vf, width, t):
    v_list = ""
    points = round(width / t)
    for i in range(0, points):
        v_list += "{:.3E}".format(Vi) + ","
    for j in range(0, points):
        v_list += "{:.3E}".format(Vf) + ","
    for k in range(0, points - 1):
        v_list += "{:.3E}".format(Vi) + ","
    v_list += "{:.3E}".format(Vi)
    return v_list