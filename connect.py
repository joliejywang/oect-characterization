# import needed packages
import pyvisa

# View the available resources to get instrument address
rm = pyvisa.ResourceManager()
rm.list_resources()

# Establish connection with the SMU

# Copy-paste the instrument address from above list of resources
inst = rm.open_resource('USB0::0x0957::0xCE18::MY51143745::INSTR')
print(inst.query("*IDN?"))

# Set the instrument timeout value in ms
# Must be longer that the longest measurement otherwise you will get an error
inst.timeout = 300*1e3