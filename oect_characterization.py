# import needed packages
import tkinter as tk
import pyvisa

# import functions from other files
# connect file
from math_functions import sine, fitSineV3
from measurement_functions import runListSweep, findIdRange, generateSineSweep, generateLinSweep, generateFreqsV2, generatePulse
from plot import plotMobility, plotPulse