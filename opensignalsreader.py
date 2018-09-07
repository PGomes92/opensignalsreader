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
07-09-2018

"""
import json
import warnings
import matplotlib.pyplot as plt
import numpy as np
import bitalino_transfer_functions as bit


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

	def __init__(self, os_file=None, show=False, raw=False):
		# Member variables
		self.file = None
		self.file_name = None
		self.device = None
		self.name = None
		self.mac = None
		self.sensors = None
		self.sampling_resolution = dict()
		self.sampling_rate = None
		self.transfer_functions = None
		self.channels = None
		self.labels = dict()
		self.digital = dict()
		self.info = dict()
		self._converted_signals = dict()
		self._raw_signals = dict()
		self._sensor_channels = dict()
		self._channels_sensors = dict()
		self._resolutions = dict()
		self._labels = []

		# Read files
		if os_file is not None:
			self._read_file(os_file)
			if show:
				self.plot(raw=raw)

	def _read_file(self, os_file):
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
		# Check input
		if type(os_file) is str:
			os_file = open(os_file, 'r')
		elif type(os_file) is not file:
			raise TypeError("Incompatible input. Please specify file path or file object.")

		self.file = os_file
		self.file_name = os_file.name

		# Check if OpenSignal file is being used
		if 'OpenSignals' in str(self.file.readline()):
			self._read_metadata()
			self._read_sensordata()
		else:
			raise TypeError("Provided file does not seem to be an OpenSignals (r)evolution file.")

		# Close file
		self.file.close()

	def _read_metadata(self):
		"""Reads metadata from OpenSignals file stored in a string with JSON dictionary.

		"""
		# Load metadata
		data = json.loads(self.file.readline().replace('#', ''))
		data = data[str(data.keys()[0])]

		# Save metadata
		self.sensors = [str(x) for x in data['sensor']]
		for i, sens in enumerate(self.sensors):
			if 'BITREV' in sens:
				self.sensors[i] = sens.replace('BITREV', '')
			elif 'BIT' in sens:
				self.sensors[i] = sens.replace('BIT', '')

		# Save metadata in info dictionary
		self.info = dict()
		self.info.update({'device name': data['device name']})
		self.info.update({'column': [str(x) for x in data['column']]})
		self.info.update({'time': data['time']})
		self.info.update({'date': data['date']})
		self.info.update({'comments': data['comments']})
		self.info.update({'mac': data['device connection']})
		self.info.update({'channels': data['channels']})
		self.info.update({'firmware': data['firmware version']})
		self.info.update({'device': data['device']})
		self.info.update({'sampling rate': data['sampling rate']})
		self.info.update({'resolution': data['resolution']})
		self.info.update({'label': data['label']})

		# Save critical data
		self.device = self.info['device']
		self.name = self.info['device name']
		self.mac = self.info['mac']
		self.sampling_rate = self.info['sampling rate']
		self.channels = self.info['channels']
		self._labels = [str(x) for x in self.info['label']]

		if self.device == 'bitalino' or 'bitalino_rev':
			self.ranges = bit.ranges
			self.units = bit.units
			self.transfer_functions = bit.transfer_functions

		# TODO add biosignalsplux
		if self.device == 'biosignalsplux':
			pass

	def _read_sensordata(self):
		"""Reads sensor data and calls required functions to convert into original units it required and possible.

		"""
		# Load data
		data = np.loadtxt(self.file)

		# Prepare data
		n = np.size(self.sensors)
		sensors = []
		self._resolutions = [int(x) for x in self.info['resolution'] if int(x) in [6, 8, 10, 12, 16]]

		for sensor in self.sensors:
			inc = 1

			if sensor in self._converted_signals.keys():
				while (sensor + str(inc)) in self._converted_signals.keys():
					inc += 1
				sensor += str(inc)
			self.sampling_resolution.update({sensor: self._resolutions[-n]})
			self._raw_signals.update({sensor: data[:, -n]})
			self._converted_signals.update({sensor: self._convert_data(sensor, data[:, -n])})

			n -= 1
			sensors.append(sensor)

		# Update sensor names
		self.sensors = sensors

		# Prepare time vector
		self._time_vector()

		# Update channel labels, dictionary with sensors as keys
		for i, sensor in enumerate(self.sensors):
			self.labels.update({sensor: self._labels[i]})
			self._sensor_channels.update({self.channels[i]: sensor})
			self._channels_sensors.update({sensor: self.channels[i]})

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
		sensor = ''.join([x for x in sensor if not x.isdigit()])
		if 'RAW' in sensor:
			return samples
		elif 'CUSTOM' in sensor:
			return samples
		elif sensor in self.transfer_functions.keys():
			output = self.transfer_functions[sensor](samples, self.sampling_resolution[sensor])
			return output
		else:
			return samples

	def _time_vector(self):
		"""Computes time vector.

		"""
		size = np.size(self._converted_signals.values()[0])
		self.t = np.linspace(0, float(size) / self.sampling_rate, size, 1. / self.sampling_rate)

	def signal(self, sensors=None):
		"""Returns converted signal data from selected sensor(s) using the respective label(s) or channel(s).

		Parameters
		----------
		sensors : str, int, list, array
			Sensor label(s) or channel(s).

		Return
		------
		signal : array, dict
			Converted sensor data of the selected label(s) or channel(s). Returns array for single sensor signal
			or dictionary when returning multiple signals.

		Raises
		------
		ValueError
			If channel numbers do not exist (array input).
		ValueError
			If sensor labels do not exist (array input).
		TypeError
			If sensors contain a mix of sensors labels and channel numbers (array input).
		ValueError
			If channel number does not exist (int or str input).
		ValueError
			If sensor label does not exist (int or str input).
		"""
		# Return all signals
		if sensors is None:
			output = self._raw_signals

		# Return all signals of a list in a wrapped in a dictionary
		elif type(sensors) is list and len(sensors) > 1:
			# Get signals for channel number list
			if all(isinstance(x, int) for x in sensors):
				# Check if channels exists
				_na_sensors = [x for x in sensors if x not in self._sensor_channels.keys()]
				if _na_sensors:
						raise ValueError("Could not find channel(s) %s in available channels." % str(_na_sensors))
				keys = [self._sensor_channels[key] for key in sensors]
				output = dict()

				# Prepare output
				for key in keys:
					output.update({key: self._raw_signals[key]})

			# Get signals for sensor label ist
			elif all(isinstance(x, str) for x in sensors):
				# Check if labels exists
				_na_sensors = [x for x in sensors if x not in self.sensors]
				if _na_sensors:
					raise ValueError("Could not find '%s' in available sensor data." % str(_na_sensors))

				output = dict()
				# Prepare output
				for key in sensors:
					output.update({key: self._raw_signals[key]})
			# Raise error if list contains mixed channel numbers and sensor labels
			else:
				raise TypeError('Please provide only channel numbers (int) or channel labels (str).')

		# Return single sensor data wrapped in numpy array
		elif (type(sensors) is list and len(sensors) == 1) or (type(sensors) is int) or (type(sensors) is str):
			if type(sensors) is list:
				sensors = sensors[0]

			# Channel number
			if type(sensors) is int:
				if sensors in self._sensor_channels.keys():
					output = self._raw_signals[self._sensor_channels[sensors]]
				else:
					raise ValueError("Could not find channel %i in available channels." % sensors)

			# Sensor label
			if type(sensors) is str:
				if sensors in self.sensors:
					output = self._raw_signals[sensors]
				else:
					raise ValueError("Could not find '%s' in available sensor data." % sensors)

		return output

	def raw(self, sensors=None):
		"""Returns raw digital signal data from selected sensor(s) using the respective label(s) or channel(s).

		Parameters
		----------
		sensors : str, int, list, array
			Sensor label(s) or channel(s).

		Return
		------
		signal : array, dict
			Raw sensor data of the selected label(s) or channel(s). Returns array for single sensor signal
			or dictionary when returning multiple signals.

		Raises
		------
		ValueError
			If channel numbers do not exist (array input).
		ValueError
			If sensor labels do not exist (array input).
		TypeError
			If sensors contain a mix of sensors labels and channel numbers (array input).
		ValueError
			If channel number does not exist (int or str input).
		ValueError
			If sensor label does not exist (int or str input).
		"""
		# Return all signals
		if sensors is None:
			output = self._raw_signals

		# Return all signals of a list in a wrapped in a dictionary
		elif type(sensors) is list and len(sensors) > 1:
			# Get signals for channel number list
			if all(isinstance(x, int) for x in sensors):
				# Check if channels exists
				_na_sensors = [x for x in sensors if x not in self._sensor_channels.keys()]
				if _na_sensors:
						raise ValueError("Could not find channel(s) %s in available channels." % str(_na_sensors))
				keys = [self._sensor_channels[key] for key in sensors]
				output = dict()

				# Prepare output
				for key in keys:
					output.update({key: self._raw_signals[key]})

			# Get signals for sensor label ist
			elif all(isinstance(x, str) for x in sensors):
				# Check if labels exists
				_na_sensors = [x for x in sensors if x not in self.sensors]
				if _na_sensors:
					raise ValueError("Could not find '%s' in available sensor data." % str(_na_sensors))

				output = dict()
				# Prepare output
				for key in sensors:
					output.update({key: self._raw_signals[key]})
			# Raise error if list contains mixed channel numbers and sensor labels
			else:
				raise TypeError('Please provide only channel numbers (int) or channel labels (str).')

		# Return single sensor data wrapped in numpy array
		elif (type(sensors) is list and len(sensors) == 1) or (type(sensors) is int) or (type(sensors) is str):
			if type(sensors) is list:
				sensors = sensors[0]

			# Channel number
			if type(sensors) is int:
				if sensors in self._sensor_channels.keys():
					output = self._raw_signals[self._sensor_channels[sensors]]
				else:
					raise ValueError("Could not find channel %i in available channels." % sensors)

			# Sensor label
			if type(sensors) is str:
				if sensors in self.sensors:
					output = self._raw_signals[sensors]
				else:
					raise ValueError("Could not find '%s' in available sensor data." % sensors)

		return output

	def plot(self, sensors=None, raw=False, interval=None, show=True, figsize=None):
		"""Plots sensor signals.

		Notes
		-----
		..	Signals are plotted in provided order in array (you can define the order yourself).
		.. 	If no input is provided, all available sensor signals will be plotted.

		Parameters
		----------
		sensors : str, array, optional
			Sensor keys of the signals to be plotted. If None, all available signals will be plotted.
		raw : bool, optional
			If true, plot raw sensor data, otherwise plot converted sensor signals.
		interval : list, array (2-elements)
			Visualization interval [x_min, x_max]
		show : bool, optional
			If True, shows the figure of the plotted signals; default: True
		figsize : 2-element array, optional
			Figsize as used in pyplot figures ([width, height]); default: None (will be adjusted dynamically)

		Raises
		------
		ValueError
			If channel number does not exist.
		ValueError
			If sensor label does not exist.
		TypeError
			If sensors contain a mix of sensors labels and channel numbers.
		"""
		# Check plot start time
		if sensors is None:
			sensors = self.sensors if np.size(self.sensors) > 1 else self.sensors[0]
		elif sensors is not None:
			if type(sensors) is list and len(sensors) > 1:
				if all(isinstance(x, int) for x in sensors):
					_na_sensors = [x for x in sensors if x not in self._sensor_channels.keys()]
					if _na_sensors:
							raise ValueError("Could not find channel(s) %s in available channels." % str(_na_sensors))
					sensors = [self._sensor_channels[key] for key in sensors]
				elif all(isinstance(x, str) for x in sensors):
					_na_sensors = [x for x in sensors if x not in self.sensors]
					if _na_sensors:
						raise ValueError("Could not find '%s' in available sensor data." % str(_na_sensors))
					sensors = [x for x in sensors if x in self.sensors]
				else:
					raise TypeError('Please provide only channel numbers (int) or channel labels (str).')
			elif (type(sensors) is list and len(sensors) == 1) or (type(sensors) is int) or (type(sensors) is str):
				if type(sensors) is list:
					sensors = sensors[0] if type(sensors[0]) is str else self._sensor_channels[sensors[0]]
				elif type(sensors) is int:
					sensors = self._sensor_channels[sensors]

		if interval is None:
			interval = [0, self.t[-1]]
		elif interval[0] >= interval[1] or interval[0] < 0 or interval[1] < 0:
			interval = [0, self.t[-1]]
			warnings.warn('Invalid interval. Interval set to [0, max_signal_duration].')

		# Create plots
		if np.size(sensors) > 1:
			if sensors is None:
				sensors = self.sensors

			# Plot multiple sensor signals
			if figsize is None:
				figsize = (12, 6)

			n_plots = len(sensors)

			# Prepare figure and plots
			if len(sensors) in [1,2,3]:
				rows, columns = len(sensors), 1
			elif len(sensors) == 4:
				rows, columns = 2, 2
			elif len(sensors) in [5, 6]:
				rows, columns = 3, 2

			fig, axs = plt.subplots(nrows=rows, ncols=columns, sharex=True, figsize=figsize)
			fig.suptitle('Sensor Signals\n(OpenSignals (r)evolution file: %s)' % self.file_name)

			# Plot signals
			k = 0
			column = 0
			row = 0

			if len(sensors) == 5:
				axs[2][1].set_visible(False)
			for n, sens in enumerate(sensors):
				key = self._get_key(sens)
				channel = self._channels_sensors[sens]

				if sens not in self.transfer_functions.keys():
					key = 'RAW'

				if raw:
					signal = self._raw_signals[sens]
					units = 'RAW'
					ranges = [0, 2**self.sampling_resolution[sens]]
				else:
					signal = self._converted_signals[sens]
					units = self.units[key]
					ranges = self.ranges[key]
				if columns == 2:
					axs[row][column].plot(self.t, signal)
					axs[row][column].axis([interval[0], interval[1], ranges[0], ranges[1]])
					axs[row][column].set_ylabel('CH%i - %s (%s)' % (channel, sens, units))
					axs[row][column].grid()
					if row == rows - 1:
						axs[row][column].set_xlabel('Time (s)')
				else:
					axs[row].plot(self.t, signal)
					axs[row].axis([interval[0], interval[1], ranges[0], ranges[1]])
					axs[row].set_ylabel('CH%i - %s (%s)' % (channel, sens, units))
					axs[row].grid()
					if row == rows - 1:
						axs[row].set_xlabel('Time (s)')
				# Manage rows and columns
				if row < rows - 1:
					row += 1
				else:
					row = 0
					column += 1
		else:
			# Plot single sensor signal
			channel = self._channels_sensors[sensors]
			if figsize is None:
				figsize = (12, 3)
			if raw:
				signal = self._raw_signals[sensors]
				units = 'RAW'
				ranges = [0, 2 ** self.sampling_resolution[sensors]]
			else:
				signal = self._converted_signals[sensors]
				units = self.units[sensors]
				ranges = self.ranges[sensors]
			fig = plt.figure(figsize=figsize)
			fig.suptitle('Sensor Signals - OpenSignals (r)evolution file: %s' % self.file_name)
			ax = fig.add_subplot(111)
			key = self._get_key(sensors)
			ax.plot(self.t, signal)
			ax.axis([interval[0], interval[1], ranges[0], ranges[1]])
			ax.set_xlabel('Time (s)')
			ax.set_ylabel('CH%i - %s (%s)' % (channel, sensors, units))
			ax.grid()

		if show:
			plt.show()

	def _get_key(self, sensor):
		key = ''.join([x for x in sensor if not x.isdigit()])
		if key not in self.transfer_functions.keys():
			return'RAW'
		else:
			return key


if __name__ == '__main__':
	# Create OpenSignals object and read file with automatic signal conversion and plot results
	acq = OpenSignalsReader('SampleECG.txt')

	# Plot all signals
	acq.plot()

	# Plot all raw signals
	acq.plot(raw=True)

	# Plot ECG signal using channel number
	acq.plot(2)

	# Plot ECG signal using channel label
	acq.plot('ECG')

	# Access signal using channel number
	acq.signal(1)

	# Access signal using channel label
	acq.signal('ECG')

	# Access raw signal using channel number
	acq.raw(1)

	# Access raw signal using channel label
	acq.raw('ECG')
