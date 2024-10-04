
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import time
from abc import ABC, abstractmethod
import logging
import tkinter as tk

from simple_pid import PID
# from pymeasure.instruments.thorlabs import ThorlabsPM100USB
import pyvisa
from ThorlabsPM100 import ThorlabsPM100

import guiSequence as seq


def _from_rgb(rgb):
    """translates an rgb tuple of int to a tkinter friendly color code
    """
    return "#%02x%02x%02x" % rgb

class CalibrateProjector(ABC):
    def __init__(self, calibration_dir, sleep_time=3, wavelength=660, threshold=0.01, debug=False):
        # PID variables
        self.sleep_time = sleep_time
        self.threshold = threshold

        self.leds = ['R', 'G', 'B', 'O', 'C', 'V']
        self.peak_wavelengths = [660, 550, 450, 590, 510, 410]
        
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
            file.write('LED,Level,Control,Power\n')

        # configure measurement
        self.measurement_wavelength = wavelength
        self.instrum = self.getInstrument() if not debug else None

    
    def createSequenceFile(self, led, control, level, current=100, mode='RGB'):
        if mode == 'RGB':
            mapping = [6, 4, 2]
        else: 
            mapping = [5, 3, 1]
        with open(self.seq_filename, 'w') as file:
            file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")

            for j in range(8):
                for i in range(3):
                    if i == led and j == level:
                        file.write(f"1, {float(control * 100)}, {current}, {mapping[i]}\n")
                    else:
                        file.write(f"1, 0, {current}, {mapping[i]}\n") # set other rows to 0
    
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

    def writeControlPowerData(self, led, level, control, power):
        with open(self.final_data_filename, 'a') as file:
            file.write(f'{led},{level},{control},{power}\n')

    def getInstrument(self):
        try:
            # theresa_device_id ='USB0::4883::32888::P0015224::0::INSTR'
            # will_device_id = 'USB0::0x1313::0x8078::P0015224::INSTR'
            # instrum = ThorlabsPM100USB(will_device_id)
            rm = pyvisa.ResourceManager()
            inst = rm.open_resource('USB0::0x1313::0x8078::P0015224::INSTR',
                                    timeout=1)
            instrum = ThorlabsPM100(inst=inst)
            instrum.sense.average.count = 500 # take an average measurement over 100 samples
        except:
            logging.error('Could not connect to Thorlabs PM100D')
            raise Exception('Could not connect to Thorlabs PM100D')
        return instrum

    @abstractmethod
    def run_calibration(self, gui):
        pass

    # def run_gamma_check(self, gui, channel):
    #
    #     self.sendSequenceTable()


    def setUpPlot(self, led, level, setpoint, isCurrent=False):
        plt.ion()
        fig = plt.figure(layout='constrained')
        gs = GridSpec(4, 2, figure=fig)
        controlvspower = fig.add_subplot(gs[0:2, 0])
        name = 'Current' if isCurrent else 'PWM'
        controlvspower.set_xlabel(name)
        controlvspower.set_ylabel('Power')
        line, = controlvspower.plot([], [], marker='o', linestyle='')

        timevspower = fig.add_subplot(gs[2:4, 0])
        timevspower.set_xlabel('Time')
        timevspower.set_ylabel('Power (uW)')
        line2, = timevspower.plot([], [], marker='o', linestyle='-', color='r')

        timevcontrol = fig.add_subplot(gs[0, 1])
        timevcontrol.set_xlabel('Time')
        timevcontrol.set_ylabel(name)
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

        fig.suptitle(f'Gamma Calibration for LED {led}  Level {level}, Setpoint: {setpoint}')

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

    def __init__(self, gui, calibration_dir, sleep_time=2, wavelength=660, threshold=0.01, debug=False):
        super().__init__(calibration_dir, sleep_time=sleep_time, wavelength=wavelength, threshold=threshold, debug=debug)

        # PID variables
        self.debug = debug
        levels = [2**i for i in range(8)]
        levels.reverse()
        self.levels = levels


    def run_calibration(self, gui, start_led=0, start_level=0):
        self.max_powers = None if self.debug else self.grab_all_max_powers(gui) # what do we want the max led intensity to be?
        print(self.max_powers)

        for led in [0, 3, 4, 5]:
        # for led in [0]:
            if self.instrum is not None:
                self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]
            set_points = [self.max_powers[led] * level/128 for level in self.levels]
            last_control = 0
            for i, level in enumerate(self.levels[start_level:]):
                self.setUpPlot(led, level, set_points[i])

                # configure PID to the settings in simulation, set the set_point to the intended power level
                div_fact = np.clip(2**(i-1), 1, None)
                # old PID settings, not converging fast enough due to characteristics of the percentages
                # pid = PID(0.00139/div_fact, 0.2/div_fact, 0.00000052/div_fact, setpoint=self.set_points[i], sample_time=None)  # works in microwatts
                # converges in 8 iterations
                pid = PID(0.00139, 0.2 * (2**i), 0.00000052, setpoint=set_points[i], sample_time=None)
                pid.output_limits = (0, 1)

                start_time = time.time()
                power = 0.0
                itr = 0
                while True:
                    control = pid(power, dt=0.01)
                    # send control level to the driver
                    if self.debug:
                        power = float(250/(2**i)*np.sqrt(control))
                    else:
                        # send the sequence to the device
                        if led > 2:
                            self.createSequenceFile(led % 3, control, i, mode='OCV')
                        else:
                            self.createSequenceFile(led % 3, control, i)
                        seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename) # load the sequence
                        seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename) # load the sequence
                        gui.ser.uploadSyncConfiguration() # upload the sequence to the driver

                        # wait for the sequence to load, and change, and measure
                        time.sleep(self.sleep_time)

                        # measure the power meter
                        power = self.instrum.read * 1000000.0 # convert to microwatts
                        print(control, power, set_points[i])

                    # write the data out to a file
                    elapsed_time = time.time() - start_time
                    p, intg, d = pid.components
                    self.writePIDData(level, set_points[i], elapsed_time, power, p, intg, d, control)
                    self.plotPIDData(elapsed_time, power, p, intg, d, control)

                    itr = itr + 1
                    # stop the pid loop if the power is within the signficant digits of the settings
                    # self.threshold = 10**(-i-1)
                    if abs(pid.setpoint - power) < self.threshold:
                        logging.info(f'Gamma calibration for led {led} level {level} complete - Control: {control} Power: {power}')
                        break

                    if itr > 100 and abs(control - last_control) < float(1/65535): 
                        logging.info(f'Gamma calibration for led {led} level {level} did not finish - Control: {control} Power: {power}')
                        self.run_finetune_current_calibration(gui, control, led, i)
                        break

                # Save the figure in the data directory
                self.writeControlPowerData(led, level, control, power)
                self.fig.savefig(os.path.join(self.plot_dirname, f'gamma_calibration_{led}_{level}.png'))
                plt.close(self.fig)

    def run_finetune_current_calibration(self, gui, last_pwm_control, led, start_level=0):
        set_points = [self.max_powers[led] * level/128 for level in self.levels]
        for i, level in enumerate(self.levels[start_level:]):
            self.setUpPlot(led, level, set_points[i], isCurrent=True)

            # configure PID to the settings in simulation, set the set_point to the intended power level
            pid = PID(0.00139, 0.2 * (2**i), 0.00000052, setpoint=set_points[i], sample_time=None)
            pid.output_limits = (0, 1)

            start_time = time.time()
            power = 0.0
            itr = 0
            while True:
                control = pid(power, dt=0.01)
                # send control level to the driver
                if self.debug:
                    power = float(250/(2**i)*np.sqrt(control))
                else:
                    # control the CURRENT instead of the PWM in order to get a step down at lower bits
                    if led > 2:
                        self.createSequenceFile(led % 3, last_pwm_control, i, current=control, mode='OCV')
                    else:
                        self.createSequenceFile(led % 3, last_pwm_control, i, current=control)
                    seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename) # load the sequence
                    seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename) # load the sequence
                    gui.ser.uploadSyncConfiguration() # upload the sequence to the driver

                    # wait for the sequence to load, and change, and measure
                    time.sleep(self.sleep_time)

                    # measure the power meter
                    power = self.instrum.read * 1000000.0 # convert to microwatts
                    print(control, power, set_points[i])

                # write the data out to a file
                elapsed_time = time.time() - start_time
                p, intg, d = pid.components
                self.writePIDData(level, set_points[i], elapsed_time, power, p, intg, d, control)
                self.plotPIDData(elapsed_time, power, p, intg, d, control)

                itr = itr + 1
                # stop the pid loop if the power is within the signficant digits of the settings
                # self.threshold = 10**(-i-1)
                if abs(pid.setpoint - power) < self.threshold or itr > 50:
                    logging.info(f'Gamma calibration for led {led} level {level} complete - Control: {control} Power: {power}')
                    break

            # Save the figure in the data directory
            self.writeControlPowerData(led, level, control, power)
            self.fig.savefig(os.path.join(self.plot_dirname, f'gamma_calibration_current_{led}_{level}.png'))
            plt.close(self.fig)

    def grab_all_max_powers(self, gui, percent_of_max=0.8):
        max_powers = []

        for led in range(6):
            if self.instrum is not None:
                self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]

            if led > 2:
                self.createSequenceFile(led % 3, 1,0, mode='OCV')
            else:
                self.createSequenceFile(led, 1, 0)
            seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename)  # load the sequence
            seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename)  # load the sequence
            gui.ser.uploadSyncConfiguration()  # upload the sequence to the driver

            # wait for the sequence to load, and change, and measure
            time.sleep(self.sleep_time)

            # measure the power meter
            max_powers += [self.instrum.read * 1000000.0 * percent_of_max]  # convert to microwatts

        max_power_data_file = os.path.join(self.dirname, 'max_power_data.csv')
        with open(max_power_data_file, 'w') as file:
            file.write('LED,Power\n')
            for i in range(6):
                file.write(f'{i},{max_powers[i]}\n')
        return max_powers
    

    def run_gamma_check(self, step_size=10):

        def record_power(control):
            power = 0 if self.debug else self.instrum.read * 1000000.0
            print(f"{control}, {power}")
            with open(self.gamma_check_power_filename, 'a') as file:
                file.write(f'{control},{power},\n')

        for led in [5]:
            if self.instrum is not None:
                self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]
            
            self.gamma_check_power_filename = os.path.join(self.dirname, f'gamma_check_{led}.csv')
            with open(self.gamma_check_power_filename, 'w') as file:
                file.write('Control,Power\n')

            # Create main window
            root = tk.Tk()
            root.geometry('%dx%d+%d+%d' % (1140, 912, 1920, 0))
            root.configure(background=_from_rgb((0, 0, 0)))
            root.overrideredirect(True)
            root.state("zoomed")
            # root.attributes("-fullscreen", True)
            root.bind("<F11>", lambda event: root.attributes("-fullscreen",
                                                not root.attributes("-fullscreen")))
            root.bind("<Escape>", lambda event: root.attributes("-fullscreen", False))
            # root.geometry("500x500+200+200")
            root.title("Gamma Calibration Screen")
            root.resizable(width = False, height = False)

            def get_colour(index):
                colours = [_from_rgb(tuple([i if j == index else 0 for j in range(3)])) for i in range(0, 256, step_size)]
                values = [tuple([i if j == index else 0 for j in range(3)]) for i in range(0, 256, step_size)]
                for c, v in zip(colours, values):
                    yield c, v
                yield None
            def start():
                out = next(colour_getter)
                if out is None:
                    return
                color, value = out
                root.configure(background=color) # set the colour to the next colour generated
                time.sleep(self.sleep_time)
                record_power(value[led % 3])
                root.after(self.sleep_time * 1000, start) # unit is milliseconds
            
            colour_getter = get_colour(led % 3)

            start()
            root.mainloop()
            # startButton = tk.Button(root,text="START",command=start)
            # startButton.pack()


def run_gamma_calibration(gui, debug=False):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'calibration_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.1)
    calibrator.run_calibration(gui)


def run_gamma_check(gui, debug=True):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'calibration_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.1)
    calibrator.run_gamma_check()

if __name__ == "__main__":
    run_gamma_calibration(None, debug=True)
    # run_gamma_check(None, debug=False)