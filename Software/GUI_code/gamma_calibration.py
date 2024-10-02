
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import time
from abc import ABC, abstractmethod
import logging

from simple_pid import PID
# from pymeasure.instruments.thorlabs import ThorlabsPM100USB
import pyvisa
from ThorlabsPM100 import ThorlabsPM100

import guiSequence as seq


class CalibrateProjector(ABC):
    def __init__(self, starting_power, calibration_dir, sleep_time=2, wavelength=660, threshold=0.01, debug=False):
        # PID variables
        self.starting_power = starting_power
        self.sleep_time = sleep_time
        self.threshold = threshold
        
        # configure logging
        self.dirname = calibration_dir
        os.makedirs(self.dirname, exist_ok=False)
        self.configureLogger()

        # configure sequence file
        self.seq_filename = os.path.join(self.dirname, 'sequence.csv')
        
        # configure PID data file
        self.pid_data_filename = self.configurePIDDataFile()
        self.plot_dirname = os.path.join(self.dirname, 'plots')
        os.makedirs(self.plot_dirname, exist_ok=True)

        # configure final measurement data file
        self.final_data_filename = os.path.join(self.dirname, 'calibrated_control.csv')
        with open(self.final_data_filename, 'w') as file:
            file.write('Level,Control,Power\n')

        # configure measurement
        self.measurement_wavelength = wavelength
        self.instrum = self.getInstrument() if not debug else None
        if self.instrum is not None:
            self.instrum.sense.correction.wavelength = self.measurement_wavelength
    
    def createSequenceFile(self, control, level, mode='RGB'):
        if mode == 'RGB':
            mapping = [6, 4, 2]
        else: 
            mapping = [5, 3, 1]
        with open(self.seq_filename, 'w') as file:
            file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")

            for j in range(8):
                for i in range(3):
                    if i == 0 and j == level:
                        file.write(f"1, {float(control * 100)}, 70.1, {mapping[i]}\n")
                    else:
                        file.write(f"1, 0, 70.1, {mapping[i]}\n") # set other rows to 0
    
    def configureLogger(self):
        log_filename = os.path.join(self.dirname, datetime.now().strftime('gamma_calibration_%Y%m%d_%H%M%S.log'))
        logging.basicConfig(filename=log_filename, level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s')
        logging.info(f'Begin Logging at {datetime.now()}')


    def configurePIDDataFile(self):
        data_filename = os.path.join(self.dirname, datetime.now().strftime('gamma_calibration_data_%Y%m%d_%H%M%S.csv'))
        with open(data_filename, 'w') as file:
            file.write('Timestamp,Level,TimeElapsed,Power,SetPoint,P,I,D,Control\n')
        
        logging.info(f'PID data will be saved to {data_filename}')
        return data_filename

    def writePIDData(self, level, setpoint, time_elapsed, power, p, i, d, control):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.pid_data_filename, 'a') as file:
            file.write(f'{timestamp},{level},{time_elapsed},{power},{setpoint},{p},{i},{d},{control}\n')

    def writeControlPowerData(self, level, control, power):
        with open(self.final_data_filename, 'a') as file:
            file.write(f'{level},{control},{power}\n')

    def getInstrument(self):
        try:
            # theresa_device_id ='USB0::4883::32888::P0015224::0::INSTR'
            # will_device_id = 'USB0::0x1313::0x8078::P0015224::INSTR'
            # instrum = ThorlabsPM100USB(will_device_id)

            rm = pyvisa.ResourceManager()
            inst = rm.open_resource('USB0::0x1313::0x8078::P0015224::INSTR',
                                    timeout=1)
            instrum = ThorlabsPM100(inst=inst)
            instrum.sense.average.count = 100 # take an average measurement over 100 samples
        except:
            logging.error('Could not connect to Thorlabs PM100D')
            raise Exception('Could not connect to Thorlabs PM100D')
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

    def __init__(self, starting_power, calibration_dir, sleep_time=2, wavelength=660, threshold=0.01, debug=False):
        super().__init__(starting_power, calibration_dir, sleep_time=sleep_time, wavelength=wavelength, threshold=threshold, debug=debug)

        # PID variables
        self.debug = debug
        levels = [2**i for i in range(8)]
        levels.reverse()
        self.levels = levels
        self.set_points = [self.starting_power * level / 128 for level in levels]

    def run_calibration(self, gui):
        for i, level in enumerate(self.levels):
            self.setUpPlot(level, self.set_points[i])

            # configure PID to the settings in simulation, set the set_point to the intended power level
            div_fact = np.clip(2**(i-1), 1, None)
            # old PID settings, not converging fast enough due to characteristics of the percentages
            # pid = PID(0.00139/div_fact, 0.2/div_fact, 0.00000052/div_fact, setpoint=self.set_points[i], sample_time=None)  # works in microwatts
            # converges in 8 iterations
            pid = PID(0.00139, 0.4 * (2**level), 0.00000052, setpoint=power, sample_time=None)
            pid.output_limits = (0, 1)

            start_time = time.time()
            power = 0
            while True:
                control = pid(power, dt=0.01)
                # send control level to the driver
                if self.debug:
                    power = float(250/(2**level)*np.sqrt(control))
                else:
                    # send the sequence to the device
                    self.createSequenceFile(control, i)
                    seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename) # load the sequence
                    seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename) # load the sequence
                    gui.ser.uploadSyncConfiguration() # upload the sequence to the driver

                    # wait for the sequence to load, and change, and measure
                    time.sleep(self.sleep_time)

                    # measure the power meter
                    power = self.instrum.read * 1000000 # convert to microwatts
                    print(control, power, self.set_points[i])

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
            self.writeControlPowerData(level, control, power)
            self.fig.savefig(os.path.join(self.plot_dirname, f'gamma_calibration_{level}.png'))
            plt.close(self.fig)
            

def run_gamma_calibration(gui, widget, debug=False):
    starting_power_at_128 = 213 # TODO: Must set in microwatts
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'calibration_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(starting_power_at_128, calibration_dir, debug=debug)
    calibrator.run_calibration(gui)

if __name__ == "__main__":
    run_gamma_calibration(None, None, debug=True)