
from datetime import datetime
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
import numpy as np
import time
from abc import ABC, abstractmethod
import logging
import tkinter as tk
import concurrent.futures

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
            file.write('LED,Level,PWM,Current,Power\n')

        # configure measurement
        self.measurement_wavelength = wavelength
        self.instrum = self.getInstrument() if not debug else None

    
    def createSequenceFile(self, led, control, level, current=1, mode='RGB'):
        if mode == 'RGB':
            mapping = [6, 4, 2]
        else: 
            mapping = [5, 3, 1]
        with open(self.seq_filename, 'w') as file:
            file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")

            for j in range(8):
                for i in range(3):
                    if i == led and j == level:
                        file.write(f"1, {float(control * 100)}, {current * 100}, {mapping[i]}\n")
                    else:
                        file.write(f"1, 0, 0, {mapping[i]}\n") # set other rows to 0
    
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

    def writeControlPowerData(self, led, level, pwm, current, power):
        with open(self.final_data_filename, 'a') as file:
            file.write(f'{led},{level},{pwm},{current},{power}\n')

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

    def measure_power(self): # returns in microwatts
        def record_power():
            power = 0 if self.debug else self.instrum.read * 1000000.0
            return power
        num_tries = 5
        while num_tries > 0: # this function is so faulty so we gotta break off a thread
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(record_power)  # Task that takes 6 seconds
                try:
                    # Attempt to get the result with a 5-second timeout
                    power = future.result(timeout=5)
                    break
                except concurrent.futures.TimeoutError:
                    num_tries -= 1
                    print(f"Measurement Failed. Trying again {num_tries} more times.")
                    self.instrum = self.getInstrument() # untested
        return power

    @abstractmethod
    def run_calibration(self, gui):
        pass

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

    def finetune_specific_mask(self, gui, led, level, calibration_csv_filename=None):
        set_point = 0
        isFinetune = False
        if calibration_csv_filename is not None:
            isFinetune = True
            df = pd.read_csv(calibration_csv_filename)
            row = df[(df['LED'] == led) & (df['Level'] == level*2)]
            set_point = row['Power'].item() /2

        if self.instrum is not None:
            self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]

        # set_points = [self.max_powers[led] * level / 128 for level in self.levels]
        last_control = 0
        last_level_control = 0
        i = self.levels.index(level)
        self.setUpPlot(led, level, set_point)

        # configure PID to the settings in simulation, set the set_point to the intended power level
        div_fact = np.clip(2 ** (i - 1), 1, None)


        starting_out = 0
        if isFinetune:
            row = df[(df['LED'] == led) & (df['Level'] == level)]
            starting_out = row['PWM'].item() / 100
        print(starting_out)
        pid = PID(0.00139, 0.2 * (2 ** i), 0.000000052, setpoint=set_point, sample_time=None,
                  starting_output=starting_out)
        pid.output_limits = (0, 1)

        start_time = time.time()
        self.send_led_bitmask_intensity(gui, led, i, starting_out, 1)
        power = self.measure_power()
        itr = 0
        while True:
            control = pid(power, dt=0.01)
            # send control level to the driver
            if self.debug:
                power = float(250 / (2 ** i) * np.sqrt(control))
            else:
                # send the sequence to the device
                self.send_led_bitmask_intensity(gui, led, i, control, 1)

                # measure the power meter
                power = self.measure_power()
                # power = self.instrum.read * 1000000.0 # convert to microwatts
                print(control, power, set_point)

            # write the data out to a file
            elapsed_time = time.time() - start_time
            p, intg, d = pid.components
            self.writePIDData(level, set_point, elapsed_time, power, p, intg, d, control)
            self.plotPIDData(elapsed_time, power, p, intg, d, control)

            itr = itr + 1
            # stop the pid loop if the power is within the signficant digits of the settings
            # self.threshold = 10**(-i-1)
            # if level < 8:
            #     self.threshold = self.threshold/10
            if abs(pid.setpoint - power) < self.threshold:
                logging.info(
                    f'Gamma calibration for led {led} level {level} complete - Control: {control} Power: {power}')

                # Save the figure in the data directory, needs to be here cuz the other one calls another function
                self.writeControlPowerData(led, level, control * 100, 1 * 100, power)
                self.fig.savefig(os.path.join(self.plot_dirname, f'gamma_calibration_{led}_{level}.png'))
                plt.close(self.fig)

                last_level_control = control
                break

            if abs(control - last_control) <= float(1 / 65535) and itr > 3:
                logging.info(
                    f'Gamma calibration for led {led} level {level} did not finish - Control: {control}, Power: {power}')
                self.run_finetune_current_calibration(gui, last_level_control, led, i)
                break


    def run_calibration(self, gui, start_led=0, start_level=0, calibration_csv_filename=None):
        self.max_powers = None if self.debug else self.grab_all_max_powers(gui) # what do we want the max led intensity to be?
        print(self.max_powers)

        isFinetune = False
        if calibration_csv_filename is not None:
            isFinetune=True
            df = pd.read_csv(calibration_csv_filename)

        for led in [3, 4, 5]:

            if self.instrum is not None:
                self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]

            set_points = [self.max_powers[led] * level / 128 for level in self.levels]
            last_control = 0
            last_level_control = 0
            for i, level in enumerate(self.levels):
                if i < start_level:
                    continue
                self.setUpPlot(led, level, set_points[i])

                # configure PID to the settings in simulation, set the set_point to the intended power level
                div_fact = np.clip(2**(i-1), 1, None)
                # old PID settings, not converging fast enough due to characteristics of the percentages
                # pid = PID(0.00139/div_fact, 0.2/div_fact, 0.00000052/div_fact, setpoint=self.set_points[i], sample_time=None)  # works in microwatts
                # converges in 8 iterations

                starting_out = 0
                if isFinetune:
                    row = df[(df['LED'] == led) & (df['Level'] == level)]
                    starting_out = row['PWM'].item()/100
                print(starting_out)
                pid = PID(0.00139, 0.2 * (2**i), 0.00000052, setpoint=set_points[i], sample_time=None, starting_output=starting_out)
                pid.output_limits = (0, 1)

                start_time = time.time()
                self.send_led_bitmask_intensity(gui, led, i, starting_out, 1)
                power = self.measure_power()
                itr = 0
                while True:
                    control = pid(power, dt=0.01)
                    # send control level to the driver
                    if self.debug:
                        power = float(250/(2**i)*np.sqrt(control))
                    else:
                        # send the sequence to the device
                        self.send_led_bitmask_intensity(gui, led, i, control, 1)

                        # measure the power meter
                        power = self.measure_power()
                        # power = self.instrum.read * 1000000.0 # convert to microwatts
                        print(control, power, set_points[i])

                    # write the data out to a file
                    elapsed_time = time.time() - start_time
                    p, intg, d = pid.components
                    self.writePIDData(level, set_points[i], elapsed_time, power, p, intg, d, control)
                    self.plotPIDData(elapsed_time, power, p, intg, d, control)

                    itr = itr + 1
                    # stop the pid loop if the power is within the signficant digits of the settings
                    # self.threshold = 10**(-i-1)
                    # if level < 8:
                    #     self.threshold = self.threshold/10
                    if abs(pid.setpoint - power) < self.threshold:
                        logging.info(f'Gamma calibration for led {led} level {level} complete - Control: {control} Power: {power}')

                        # Save the figure in the data directory, needs to be here cuz the other one calls another function
                        self.writeControlPowerData(led, level, control * 100, 1 * 100, power)
                        self.fig.savefig(os.path.join(self.plot_dirname, f'gamma_calibration_{led}_{level}.png'))
                        plt.close(self.fig)

                        last_level_control = control
                        break

                    if abs(control - last_control) <= float(1/65535) and itr > 3:
                        logging.info(f'Gamma calibration for led {led} level {level} did not finish - Control: {control}, Power: {power}')
                        self.run_finetune_current_calibration(gui, last_level_control, led, i)
                        break

    def send_led_bitmask_intensity(self, gui, led, level, pwm, current):
        if led > 2:
            self.createSequenceFile(led % 3, pwm, level, current=current, mode='OCV')
        else:
            self.createSequenceFile(led % 3, pwm, level, current=current)
        seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename)  # load the sequence
        seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename)  # load the sequence
        gui.ser.uploadSyncConfiguration()  # upload the sequence to the driver

        # give some time for the hardware to catchup
        time.sleep(self.sleep_time)

    def run_finetune_current_calibration(self, gui, last_pwm_control, led, start_level=0):
        set_points = [self.max_powers[led] * level/128 for level in self.levels]
        for i, level in enumerate(self.levels):
            if i < start_level:
                continue
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
                    logging.info(f'Gamma calibration for led {led} level {level} complete - Current{control}, PWM: {last_pwm_control} Power: {power}')
                    break

            # Save the figure in the data directory
            self.writeControlPowerData(led, level, last_pwm_control * 100, control * 100, power)
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

    def measure_all_bit_masks(self, gui, filename):
        df = pd.read_csv(filename)
        max_powers = []

        for led in [0, 1, 2, 3, 4, 5]:
            for j in range(8):
                level =  2 ** (8 - j - 1)
                row = df[(df['LED'] == led) & (df['Level'] ==level)]
                if self.instrum is not None:
                    self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]

                pwm = row['PWM'].item()
                current = row['Current'].item()
                if led > 2:
                    self.createSequenceFile(led % 3, pwm/100, j, current=current/100, mode='OCV')
                else:
                    self.createSequenceFile(led, pwm/100.0, j, current=current/100)
                seq.loadSequence(gui, gui.sync_digital_low_sequence_table, self.seq_filename)  # load the sequence
                seq.loadSequence(gui, gui.sync_digital_high_sequence_table, self.seq_filename)  # load the sequence
                gui.ser.uploadSyncConfiguration()  # upload the sequence to the driver

                # wait for the sequence to load, and change, and measure
                time.sleep(self.sleep_time)

                # measure the power meter
                power = self.instrum.read * 1000000.0
                self.writeControlPowerData(led, level, pwm, current, power)

    # def run_spectral_measurement(self):
    #     def record_power(control):
    #         power = 0 if self.debug else self.instrum.read * 1000000.0
    #         print(f"{control}, {power}")
    #         with open(self.gamma_check_power_filename, 'a') as file:
    #             file.write(f'{control},{power},\n')
    #
    #     for led in [0, 1, 2, 3]:
    #         if self.instrum is not None:
    #             self.instrum.sense.correction.wavelength = self.peak_wavelengths[led]
    #
    #         self.gamma_check_power_filename = os.path.join(self.dirname, f'gamma_check_{led}.csv')
    #         with open(self.gamma_check_power_filename, 'w') as file:
    #             file.write('Control,Power\n')
    #
    #         # Create main window
    #         root = tk.Tk()
    #         root.geometry('%dx%d+%d+%d' % (1140, 912, 1920, 0))
    #         root.configure(background=_from_rgb((0, 0, 0)))
    #         root.overrideredirect(True)
    #         root.state("zoomed")
    #         # root.attributes("-fullscreen", True)
    #         root.bind("<F11>", lambda event: root.attributes("-fullscreen",
    #                                                          not root.attributes("-fullscreen")))
    #         root.bind("<Escape>", lambda event: root.attributes("-fullscreen", False))
    #         # root.geometry("500x500+200+200")
    #         root.title("Gamma Calibration Screen")
    #         root.resizable(width=False, height=False)
    #
    #         def get_colour(index):
    #             colours = [_from_rgb(tuple([i if j == index else 0 for j in range(3)])) for i in
    #                        [255]]
    #             values = [tuple([i if j == index else 0 for j in range(3)]) for i in  [255]]
    #             for c, v in zip(colours, values):
    #                 yield c, v
    #             yield None
    #
    #         def start():
    #             out = next(colour_getter)
    #             if out is None:
    #                 root.destroy()
    #                 return
    #             color, value = out
    #             root.configure(background=color)  # set the colour to the next colour generated
    #             time.sleep(self.sleep_time)
    #             with concurrent.futures.ThreadPoolExecutor() as executor:
    #                 future = executor.submit(record_power, value[led % 3])  # Task that takes 6 seconds
    #                 try:
    #                     # Attempt to get the result with a 5-second timeout
    #                     result = future.result(timeout=5)
    #                 except concurrent.futures.TimeoutError:
    #                     print(f"Measurement Failed at {value}. Moving on.")
    #             # record_power(value[led % 3])
    #             root.after(self.sleep_time * 1000, start)  # unit is milliseconds
    #
    #         colour_getter = get_colour(led % 3)
    #
    #         start()
    #         root.mainloop()
    #         startButton = tk.Button(root,text="START",command=start)
    #         startButton.pack()

    def run_gamma_check(self, step_size=1):

        def record_power(control):
            power = 0 if self.debug else self.instrum.read * 1000000.0
            print(f"{control}, {power}")
            with open(self.gamma_check_power_filename, 'a') as file:
                file.write(f'{control},{power},\n')

        for led in [0, 1, 2]:
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
                    root.destroy()
                    return
                color, value = out
                root.configure(background=color) # set the colour to the next colour generated
                time.sleep(self.sleep_time)
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(record_power, value[led % 3])  # Task that takes 6 seconds
                    try:
                        # Attempt to get the result with a 5-second timeout
                        result = future.result(timeout=5)
                    except concurrent.futures.TimeoutError:
                        print(f"Measurement Failed at {value}. Moving on.")
                # record_power(value[led % 3])
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


def run_gamma_calibration_finetune(gui, debug=False):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'calibration_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.02)
    calibrator.run_calibration(gui, calibration_csv_filename="calibration_20241023_164950\calibrated_control.csv")

def run_finetune_specific_mask(gui, debug=False):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'calibration_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.02)
    calibrator.finetune_specific_mask(gui, 2, 1, calibration_csv_filename="calibration_20241025_114932\calibrated_control.csv")

def measure_bitmasks(gui, debug=False):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'measure_bitmasks_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.1, sleep_time=3)
    calibrator.measure_all_bit_masks(gui, "calibration_20241025_114932\calibrated_control.csv")



def run_gamma_check(gui, debug=True):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'measurement_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.1, sleep_time=2)
    calibrator.run_gamma_check()

def run_spectral_measurement(gui, debug=True):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    calibration_dir = f'spectral_measurement_{timestamp}'

    calibrator = CalibrateEvenOdd8Bit(gui, calibration_dir, debug=debug, threshold=0.1, sleep_time=2)
    calibrator.run_spectral_measurement()


def create_sequence_file_rocv(dirname, calibration_csv_filename):
    os.makedirs(dirname, exist_ok=True)
    df = pd.read_csv(calibration_csv_filename)

    led_rows = []
    for i in [0, 1, 2, 3]:
        led_rows += [df['LED'] == i]

    mapping = [6, 4, 2]
    even_filename = os.path.join(dirname, "rgb.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3):
                if i == 0:
                    row = df[(df['LED'] == i) & (df['Level'] == 2** (8-j -1))]
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i]}\n")

                else:
                    file.write(f"1, {0}, {0}, {mapping[i]}\n")

    mapping = [5, 3, 1]
    even_filename = os.path.join(dirname, "ocv.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3, 6):
                row = df[(df['LED'] == i) & (df['Level'] == 2 ** (8-j -1))]
                try:
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i-3]}\n")
                except:
                    import pdb; pdb.set_trace()

def create_sequence_file_rgbo(dirname, calibration_csv_filename):
    os.makedirs(dirname, exist_ok=True)
    df = pd.read_csv(calibration_csv_filename)

    led_rows = []
    for i in [0, 1, 2, 3]:
        led_rows += [df['LED'] == i]

    mapping = [6, 4, 2]
    even_filename = os.path.join(dirname, "rgb.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3):
                row = df[(df['LED'] == i) & (df['Level'] == 2 ** (8 - j - 1))]
                try:
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i]}\n")
                except:
                    import pdb; pdb.set_trace()

    mapping = [5, 3, 1]
    even_filename = os.path.join(dirname, "ocv.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3, 6):
                if i == 0:
                    row = df[(df['LED'] == i) & (df['Level'] == 2 ** (8 - j - 1))]
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i-3]}\n")

                else:
                    file.write(f"1, {0}, {0}, {mapping[i-3]}\n")

def create_sequence_file_rgbocv(dirname, calibration_csv_filename):
    os.makedirs(dirname, exist_ok=True)
    df = pd.read_csv(calibration_csv_filename)

    led_rows = []
    for i in [0, 1, 2, 3]:
        led_rows += [df['LED'] == i]

    mapping = [6, 4, 2]
    even_filename = os.path.join(dirname, "rgb.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3):
                row = df[(df['LED'] == i) & (df['Level'] == 2 ** (8 - j - 1))]
                try:
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i]}\n")
                except:
                    import pdb; pdb.set_trace()

    mapping = [5, 3, 1]
    even_filename = os.path.join(dirname, "ocv.csv")
    with open(even_filename, 'w') as file:
        file.write("LED #,LED PWM (%),LED current (%),Duration (s)\n")
        for j in range(8):
            for i in range(3, 6):
                row = df[(df['LED'] == i) & (df['Level'] == 2 ** (8 - j - 1))]
                try:
                    file.write(f"1, {float(row['PWM'].item())}, {row['Current'].item()}, {mapping[i-3]}\n")
                except:
                    import pdb; pdb.set_trace()

if __name__ == "__main__":
    # run_gamma_calibration(None, debug=True)
    # run_gamma_check(None, debug=False)
    # create_sequence_file_rgbo("1025-rgbocv", "calibration_20241025_114932\calibrated_control.csv")
    # create_sequence_file_rocv("1025-rocv", "calibration_20241025_114932\calibrated_control.csv")
    # # run_spectral_measurement(None, debug=False)
    create_sequence_file_rgbocv("1025-rgbocv", "calibration_20241025_114932\calibrated_control.csv")

