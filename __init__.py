"""
OpenSignals Reader
------------------

This package provides an OpenSignal class for easy import of OpenSignals (r)evolution file
containing signals acquired with BITalino (r)evolution with automatic data conversion
using the official BITalino transfer functions.

Visit http://biosignalsplux.com/ and http://bitalino.com for more information about OpenSignals

..

:copyright: (c) 2018 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.


Notes
-----
..  This module is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..  Transfer functions are compatible for conversion data of BITalino (r)evolution
	devices only; older versions of BITalino are not supported

Author
------
..  Pedro Gomes, Master Student, University of Applied Sciences Hamburg

Thesis Supervisors
------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes, PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
29-08-2018

"""
import json
import warnings
import matplotlib.pyplot as plt
import numpy as np
import bitalino as bit


class OpenSignalsReader():
	"""This class reads OpenSignals files (metadata & sensor data) and converts sensor data into original units
	(if required) and includes easy functions to plot sensor signals and to print acquisition data to the console.

	Methods
	-------
	read_file():
		Reads OpenSignals (r)evolution file.
	info():
		Prints acquisition's metadata to the console.
	plot():
		Plots sensor signals (customizable: plots all signals, plots all signals in preferred order or only individual
		signals.

	"""

	def __init__(self, os_file=None, raw=False, show=False):
		# Member variables
		self.signals = {}
		self.samples = {}
		self.sampling_resolution = None
		self.file = None
		self.file_name = None
		self.os_file = None
		self.sensors = None
		self.name = None
		self.column = None
		self.sync_interval = None
		self.time = None
		self.comments = None
		self.mac = None
		self.channels = None
		self.date = None
		self.mode = None
		self.digital = None
		self.device = None
		self.position = None
		self.sampling_rate = None
		self.sampling_resolution = None
		self.labels = None
		self.special = None
		self.raw = False
		self.file_name = os_file
		self.ranges = None
		self.units = None
		self.transfer_functions = None

		# Read files
		if os_file is not None:
			self.read_file(os_file, raw)
			if show:
				self.plot()

	def read_file(self, os_file, raw=False):
		"""Reads OpenSignals file (metadata & sensor data) and converts sensor data into original units (if required).

		Notes
		-----
		.. 	Signals of type 'RAW', 'CUSTOM or 'UNKNOWN' will not be converted; if 'raw' is true )
		.. 	If 'raw' is True, no signals will be converted.

		Parameters
		----------
		os_file : str, file object
			OpenSignals file.
		raw : bool, optional
			If True, sensor data will be converted to original units

		"""
		# Clear data
		self.__init__()

		# Set conversion
		self.raw = raw

		# Check input
		if type(os_file) is str:
			os_file = open(os_file, 'r')
		elif type(os_file) is not file:
			raise TypeError("Incompatible input. Please specify file path or file object.")

		self.file = os_file
		self.file_name = os_file.name

		# Check if OpenSignal file is being used
		if 'OpenSignals' in str(self.file.readline()):
			self.os_file = True
			self._read_metadata()
			self._read_sensordata()
		else:
			self.os_file = False
			warnings.warn("Provided file does not seem to be an OpenSignals (r)evolution file.")

		# Close file
		self.file.close()

	def _read_metadata(self):
		"""Reads metadata from OpenSignals file stored in a string with JSON dictionary.

		"""
		# Load metadata
		data = json.loads(self.file.readline().replace('#', ''))
		data = data[str(data.keys()[0])]

		# Save metadata
		if self.raw:
			self.sensors = ['RAW' for x in data['sensor']]
		else:
			self.sensors = [str(x) for x in data['sensor']]
		self.name = data['device name']
		self.column = [str(x) for x in data['column']]
		self.sync_interval = data['sync interval']
		self.time = data['time']
		self.comments = data['comments'] if data['comments'] != '' else 'n/a'
		self.mac = str(data['device connection'])
		self.channels = data['channels']
		self.date = data['date']
		self.mode = data['mode']
		self.digital = data['digital IO']
		self.firmware = data['firmware version']
		self.device = str(data['device'])
		self.position = data['position']
		self.sampling_rate = data['sampling rate']
		self.sampling_resolution = float(data['resolution'][-1])
		self.labels = [str(x) for x in data['label']]
		self.special = data['special']

		if self.device == 'bitalino':
			self.ranges = bit.ranges
			self.units = bit.units
			self.transfer_functions = bit.transfer_functions

		# TODO add biosignalsplux

	def _read_sensordata(self):
		"""Reads sensor data and calls required functions to convert into original units it required and possible.

		"""
		# Load data
		data = np.loadtxt(self.file)

		n = np.size(self.sensors)
		sensors = []

		for sensor in self.sensors:
			inc = 1

			if sensor in self.signals.keys():
				while (sensor + str(inc)) in self.signals.keys():
					inc += 1
				sensor += str(inc)

			self.signals.update({sensor: self._convert_data(sensor, data[:, -n])})

			n -= 1
			sensors.append(sensor)
		# Update sensor names
		self.sensors = sensors

		# Prepare time vector
		self._time_vector()

		# Update channel labels, dictionary with sensors as keys
		label = self.labels
		self.labels = dict()
		for i, val in enumerate(self.sensors):
			self.labels.update({val: label[i]})

	def _convert_data(self, sensor, samples):
		"""Calls sensor specific transfer functions to convert raw sensor signals to their original units.

		Note
		----
		Signals of type 'RAW', 'CUSTOM or 'UNKNOWN' will returned without conversion.

		Parameters
		----------
		sensor : str
			OpenSignals sensor key (e.g. 'ECG' for Electrocardiography sensor).
		samples : array
			Raw sensor data.

		Returns
		-------
		samples : array
			Raw or converted sensor data.

		"""
		if sensor not in self.transfer_functions.keys():
			return samples
		elif ('RAW' or 'CUSTOM') not in sensor or sensor not in self.transfer_functions.keys():
			sensor = ''.join([x for x in sensor if not x.isdigit()])
			return self.transfer_functions[sensor](samples, self.sampling_resolution)
		else:
			return samples

	def _time_vector(self):
		"""Computes time vector.

		"""
		size = np.size(self.signals.values()[0])
		self.t = np.linspace(0, float(size) / self.sampling_rate, size, 1. / self.sampling_rate)

	def plot(self, sensors=None, show=True, figsize=None, start=None, end=None):
		"""Plots sensor signals.

		Notes
		-----
		..	Signals are plotted in provided order in array (you can define the order yourself).
		.. 	If no input is provided, all available sensor signals will be plotted.

		Parameters
		----------
		sensors : str, array, optional
			Sensor keys of the signals to be plotted. If None, all available signals will be plotted.
		show : bool, optional
			If True, shows the figure of the plotted signals; default: True
		figsize : 2-element array, optional
			Figsize as used in pyplot figures ([width, height]); default: None (will be adjusted dynamically)
		start : int, float, optional
			Start time of visualized interval; default: 0
		end : int, float, optional
			End time of visualized interval; default: max duration
		"""
		# Check plot start time
		if start is None or start < 0:
			start = 0
		elif start >= self.t[-1]:
			start = 0
			warnings.warn('Start time > total signal duration. Start time set to 0.')
		elif start >= end:
			start = 0
			end = self.t.size
			warnings.warn('Start time > end time. Start time set to 0, end time set to max duration.')

		# Check plot end time
		if end is None or end > self.t[-1]:
			end = self.t[-1]

		if sensors is None:
			sensors = self.sensors if np.size(self.sensors) > 1 else self.sensors[0]

		# Create plots
		if np.size(sensors) > 1:
			if sensors is None:
				sensors = self.sensors

			# Plot multiple sensor signals
			if figsize is None:
				figsize = (12, 6)

			rows = np.size(sensors) if np.size(sensors) <= 4 else 4
			columns = 1 if np.size(sensors) <= 4 else 2
			fig, axs = plt.subplots(nrows=rows, ncols=columns, sharex=True, figsize=figsize)

			k = 0
			for i in range(rows):
				if columns == 2:
					for j in range(columns):
						sens = sensors[k]
						key = self._get_key(sens)
						if sensors not in self.transfer_functions.keys():
							key = 'RAW'
						axs[i][j].plot(self.t, self.signals[sens])
						axs[i][j].axis([start, end, self.ranges[key][0], self.ranges[key][1]])
						axs[i][j].set_ylabel('%s (%s)' % (sens, self.units[key]))
						axs[i][j].grid()
						k += 1
					axs[-1][0].set_xlabel('Time(s)')
					axs[-1][1].set_xlabel('Time(s)')
				else:
					sens = sensors[k]
					key = self._get_key(sens)
					axs[i].plot(self.t, self.signals[sens])
					axs[i].axis([start, end, self.ranges[key][0], self.ranges[key][1]])
					axs[i].set_ylabel('%s (%s)' % (sens, self.units[key]))
					axs[i].grid()
					k += 1
					axs[-1].set_xlabel('Time(s)')
		else:
			# Plot single sensor signal
			if figsize is None:
				figsize = (12, 3)
			fig = plt.figure(figsize=figsize)
			ax = fig.add_subplot(111)
			key = self._get_key(sensors)
			ax.plot(self.t, self.signals[sensors])
			ax.axis([start, end, self.ranges[key][0], self.ranges[key][1]])
			ax.set_xlabel('Time (s)')
			ax.set_ylabel('%s (%s)' % (sensors, self.units[key]))
			ax.grid()

		if show:
			plt.show()

	def _get_key(self, sensor):
		key = ''.join([x for x in sensor if not x.isdigit()])
		if key not in self.transfer_functions.keys():
			return'RAW'
		else:
			return key

	def info(self):
		"""
		Prints metadata to the console.
		"""
		if self.file_name is not None:
			print("ACQUISITION INFORMATION")
			print("Device Name: %s" % self.name)
			print("MAC Address: %s" % self.mac)
			print("Device: %s" % self.device)
			print("Firmware: %s" % self.firmware)
			print("File: %s" % self.file_name)
			print("Date: %s" % self.date)
			print("Time: %s" % self.time)
			print("Comments: %s" % self.comments)
			print("Sampling Rate: %s" % self.sampling_rate)
			print("Sampling Resolution: %s" % self.sampling_resolution)
			print("Sensors: %s" % self.sensors)
			print("Channels: %s" % self.channels)
			print("Labels: %s" % self.labels)
			print("Columns: %s" % self.column)
			print("Sync Interval: %s" % self.sync_interval)
			print("Mode: %s" % self.mode)
			print("Position: %s" % self.position)
			print("Special: %s" % self.special)
		else:
			print('No file has been read yet.')


if __name__ == '__main__':
	# Example 1: Create OpenSignals object with automatic signal conversion and plot results
	acq = OpenSignalsReader('SampleECG.txt')

	# Example 2: Read file using the .read_file() method
	acq.read_file('SampleECG.txt')
	acq.plot()

	# Example 3: Create OpenSignals object use only raw sensor signals
	acq = OpenSignalsReader('SampleECG.txt', raw=True)

	# Example 4: Create OpenSignals object and plot automatically
	acq = OpenSignalsReader('SampleECG.txt', show=True)

	# Example 5: Get list of sensors and use the sensor's key to plot and use signal data
	# Output: ['ECG'] <- Key!
	acq.sensors

	# Accessing ECG signal
	acq.signals['ECG']

	# Plots ECG signal (use list of keys when using multiple signals (e.g. acq.plot(['ECG', 'ECG1'])
	acq.plot('ECG')

	# Example 6: Access metadata (for entire list, check the class above)
	acq.sampling_resolution			# Sampling resolution
	acq.sampling_rate				# Sampling rate
	acq.labels						# Labels of each channel
	acq.firmware					# Device firmware
	acq.mac							# MAC address

	# Example 7: Print all metadata to the console
	acq.info()
