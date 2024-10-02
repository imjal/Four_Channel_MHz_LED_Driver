from pymeasure.instruments.thorlabs import ThorlabsPM100USB

# will_device_id = 'USB0::0x1313::0x8078::P0015224::INSTR'
# instrum = ThorlabsPM100USB(will_device_id)

# print(instrum.power)
# print(instrum.wavelength)

import pyvisa
from ThorlabsPM100 import ThorlabsPM100
rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0x1313::0x8078::P0015224::INSTR',
                        timeout=1)
power_meter = ThorlabsPM100(inst=inst)
print(power_meter.read)
