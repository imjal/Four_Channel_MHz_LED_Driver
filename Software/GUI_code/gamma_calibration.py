
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import time
from abc import ABC, abstractmethod


import logging

from simple_pid import PID
from pymeasure.instruments.thorlabs import ThorlabsPM100USB

import guiSequence as seq


class CalibrateProjector(ABC):
    def __init__(self, starting_power, log_dir, data_dir, sleep_time=2, wavelength=660, threshold=0.01):
        # PID variables
        self.starting_power = starting_power
        self.sleep_time = sleep_time
        self.threshold = threshold
        
        # configure logging
        self.log_dir = log_dir
        self.data_dir = data_dir
        os.makedirs(self.log_dir, exist_ok=True)
        os.makedirs(self.data_dir, exist_ok=True)
        self.data_filename = self.configureDataFile()

        # configure measurement
        self.measurement_wavelength = wavelength
        self.instrum = self.getInstrument()
        
    
    def configureLogger(self):
        log_filename = os.path.join(self.log_dir, datetime.now().strftime('gamma_calibration_%Y%m%d_%H%M%S.log'))
        logging.basicConfig(filename=log_filename, level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s')
        logging.info(f'Begin Logging at {datetime.now()}')


    def configureDataFile(self):
        data_filename = os.path.join(self.data_dir, datetime.now().strftime('gamma_calibration_data_%Y%m%d_%H%M%S.csv'))
        with open(data_filename, 'w') as file:
            file.write('Timestamp,Level,Power,Control\n')
        
        logging.info(f'PID data will be saved to {data_filename}')
        return data_filename

    def writePIDData(self, level, setpoint, time_elapsed, power, p, i, d, control):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.data_filename, 'a') as file:
            file.write(f'{timestamp},{level},{time_elapsed},{power},{setpoint},{p},{i},{d},{control}\n')

    def getInstrument(self):
        try:
            instrum = ThorlabsPM100USB('USB0::0x1313::0x8078::P0023944::INSTR')
        except:
            logging.error('Could not connect to Thorlabs PM100D')
            return None
        return instrum

    @abstractmethod
    def run_calibration(self):
        pass

    def setUpPlot(self, level, setpoint):
        plt.ion()
        fig = plt.figure(layout='constrained')
        gs = GridSpec(4, 2, figure=fig)
        controlvspower = fig.add_subplot(gs[0:2, 0])
        controlvspower.set_xlabel('Control')
        controlvspower.set_ylabel('Power')
        line, = controlvspower.plot([], [], marker='o', linestyle='')

        timevspower = fig.add_subplot(gs[2:4, 0])
        timevspower.set_xlabel('Time')
        timevspower.set_ylabel('Power (uW)')
        line2, = timevspower.plot([], [], marker='o', linestyle='-', color='r')

        timevcontrol = fig.add_subplot(gs[0, 1])
        timevcontrol.set_xlabel('Time')
        timevcontrol.set_ylabel('Control')
        line3, = timevcontrol.plot([], [], marker='o', linestyle='-', color='orange')

        timevp = fig.add_subplot(gs[1, 1])
        timevp.set_xlabel('Time')
        timevp.set_ylabel('P')
        line4, = timevp.plot([], [], marker='o', linestyle='-', color='g')

        timevi = fig.add_subplot(gs[2, 1])
        timevi.set_xlabel('Time')
        timevi.set_ylabel('I')
        line5, = timevi.plot([], [], marker='o', linestyle='-', color='b')

        timevd = fig.add_subplot(gs[3, 1])
        timevd.set_xlabel('Time')
        timevd.set_ylabel('D')
        line6, = timevd.plot([], [], marker='o', linestyle='-', color='purple')

        fig.suptitle(f'Gamma Calibration for Level {level}, Setpoint: {setpoint}')

        self.fig = fig

    def plotPIDData(self, time_elapsed, power, p, i, d, control):

        data = [[control, power], [time_elapsed, power], [time_elapsed, control], [time_elapsed, p], [time_elapsed, i], [time_elapsed, d]]

        for i, ax in enumerate(self.fig.axes):
            line = ax.lines[0]
            line.set_xdata(np.append(line.get_xdata(), data[i][0]))
            line.set_ydata(np.append(line.get_ydata(), data[i][1]))
            ax.relim()
            ax.autoscale_view()
        plt.draw()
        plt.pause(0.01)

class CalibrateEvenOdd8Bit(CalibrateProjector):

    def __init__(self, starting_power, log_dir, data_dir, sleep_time=2, wavelength=660, debug=False):
        super().__init__(starting_power, log_dir, data_dir, sleep_time, wavelength)

        # PID variables
        self.debug = debug
        levels = [2**i for i in range(8)]
        levels.reverse()
        self.levels = levels
        self.set_points = [self.starting_power * level / 128 for level in levels]

    def run_calibration(self, gui, widget):
        for i, level in enumerate(self.levels):
            self.setUpPlot(level, self.set_points[i])

            # configure PID to the settings in simulation, set the set_point to the intended power level
            div_fact = np.clip(2**(i-1), 1, None)
            pid = PID(0.00139/div_fact, 0.2/div_fact, 0.00000052/div_fact, setpoint=self.set_points[i], sample_time=None)  # works in microwatts
            pid.output_limits = (0, 1)

            start_time = time.time()
            power = 0
            while True:
                control = pid(power, dt=0.01)
                # send control level to the driver
                if self.debug:
                    power = float(200*np.sqrt(control))
                else:
                    # send the sequence to the device TODO: this part is not done yet
                    path = "."
                    seq.loadSequence(gui, widget, path) # load the sequence
                    gui.ser.uploadSyncConfiguration() # upload the sequence to the driver

                    # wait for the sequence to load, and change, and measure
                    time.sleep(self.sleep_time)

                    # measure the power meter
                    power = self.instrum.measure_power(wavelength=self.measurement_wavelength)

                # write the data out to a file
                elapsed_time = time.time() - start_time
                p, intg, d = pid.components
                self.writePIDData(level, self.set_points[i], elapsed_time, power, p, intg, d, control)
                self.plotPIDData(elapsed_time, power, p, intg, d, control)


                # stop the pid loop if the power is within 0.01 of the setpoint
                if abs(pid.setpoint - power) < self.threshold:
                    logging.info(f'Gamma calibration for level {level} complete - Control: {control} Power: {power}')
                    break
            
            # Save the figure in the data directory
            self.fig.savefig(os.path.join(self.data_dir, f'gamma_calibration_{level}.png'))
            plt.close(self.fig)
            

def run_gamma_calibration(gui, widget, debug=False):
    starting_power_at_128 = 130 # TODO: Must set in microwatts
    log_dir = 'calibration_logs'
    data_dir = 'calibration_data'
    calibrator = CalibrateEvenOdd8Bit(starting_power_at_128, log_dir, data_dir, debug=True)
    calibrator.run_calibration(gui, widget)

if __name__ == "__main__":
    run_gamma_calibration(None, None, debug=True)

