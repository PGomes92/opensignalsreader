# -*- coding: utf-8 -*-
"""
BITalino (r)evolution Transfer Functions
----------------------------------------
Visit http://bitalino.com/

Provides functions to convert raw BITalino (r)evolution sensor samples (single samples
or entire datsets) into original physical units and supports the following sensors:

	* ECG		* EEG		* EMG 		* EDA 		* ACC	* GMR
	* TEMP		* NTC		* LUX 		* OSL		* BPR	* EOG
	* EGG

Visit the following page for more detailed information about the transfer functions:
http://bitalino.com/en/learn/documentation

Visit the following page for BITalino (r)evolution sample signals:
https://github.com/BITalinoWorld/revolution-sample-data

..

:copyright: (c) 2018 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.


Notes
-----
..  Transfer functions are compatible for conversion data of BITalino (r)evolution
	devices only; older versions of BITalino are not supported

Author
------
..  Pedro Gomes, pgomes92@gmail.com

Notes Up To Version 0.2.1
-------------------------
..  Versions 0.2.1 and earlier of this package were part of the master thesis
    "Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	Thesis Author: Pedro Gomes, Master Student, University of Applied Sciences Hamburg
..  1st Supervisor: Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg
..  2nd Supervisor: Hugo Silva, PhD, Instituto de Telecomunicacoes, PLUX wireless biosignals S.A.

Last Update
-----------
06.06.2018

"""
import warnings
import numpy as np


def ecg(samples=None, resolution=10):
	"""Converts raw ECG values into original physical unit (mV).

	BITalino Sensor: Electrocardiography (ECG)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_ECG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	ECG(mv) = ((((ADC / (2 ** n)) - 0.5) * VCC) / G_ECG) * 1000

		- ADC   Value sampled from the channel
		- n     Number of bits of the acquiring channel
		- VCC   Operating voltage (3.3V)
		- G_ECG Sensor gain (1000)

	[RANGE]
	-1.5mV to 1.5mV

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	ecg_samples : array of floats
		Converted ECG samples (mV).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check output
	ecg_samples = np.asarray([((((float(s) / (2 ** resolution)) - 0.5) * 3.3) / 1100) * 1000 for s in samples])
	_check_ranges(ecg_samples, 'ECG')
	return ecg_samples


def eeg(samples, resolution=10):
	"""Converts raw EEG values into original physical unit (µV).

	BITalino Sensor: Electroencephalography (EEG)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_EEG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	EEG(µV) = ((((ADC/2^n) - 0.5) * VCC) / GEEG) * (10^6)

		EEG(µV)		Sample value in µV.
		ADC			Value sampled from the channel (raw value)
		n 			Number of bits of the channel
		VCC			Operating voltage; 3.3V
		GEEG		Sensor gain (40 000)

	[RANGE]
	-41.25µV to 41.25µV

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	eeg_samples : array of floats
		Converted EEG sample.

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check output
	eeg_samples = np.asarray([((((float(s) / 2 ** resolution) - 0.5) * 3.3) / 40000) * (10 ** 6) for s in samples])
	_check_ranges(eeg_samples, 'EEG')
	return eeg_samples


def emg(samples=None, resolution=10):
	"""Converts raw EMG values into original physical unit (mV).

	BITalino Sensor: Electromyography (EMG)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_EMG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	EMG(mV) = ((((ADC/2^n) - 0.5) * VCC) / GEMG) * 1000

		EMG(mV)		Sample value in mV.
		ADC			Value sampled from the channel (raw value)
		n 			Number of bits of the channel
		VCC			Operating voltage; 3.3V
		GEMG		Sensor gain (1009)

	[RANGE]
	-1.65mV to 1.65mV

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	emg_samples : array of floats
		Converted EMG samples (mV).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check output
	emg_samples = np.asarray([((((float(s) / 2 ** resolution) - 0.5) * 3.3) / 1009) * 1000 for s in samples])
	_check_ranges(emg_samples, 'EMG')
	return emg_samples


def eda(samples=None, resolution=10):
	"""Converts raw EDA values into original physical unit (µS).

	BITalino Sensor: Electrodermal Activity (EDA)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_EDA_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	EDA(µS) = ((((ADC/2^n) * VCC) - 0.574) / 0.132)

		EMG(mV)		Sample value in mV.
		ADC			Value sampled from the channel (raw value)
		n 			Number of bits of the channel
		VCC			Operating voltage; 3.3V
		GEMG		Sensor gain (1009)

	[RANGE]
	-4.4µS to 21µS

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit; default: 10-bit).

	Returns
	-------
	eda_samples : array of floats
		Converted EDA samples (µS).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check input
	eda_samples = np.asarray([((float(s) / (2 ** resolution)) * 3.3) / 0.132 for s in samples])
	_check_ranges(eda_samples, 'EDA')
	return _check_ranges


def acc(samples=None, resolution=10, c_min=400, c_max=600):
	"""Converts raw ACC values into original units.

	BITalino Sensor: Accelerometer (ACC)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_ACC_Sensor_Datasheet.pdf

	Transfer Function
	-----------------
	ACC(g) = ((ADC - Cmin)/(Cmax - Cmin)) * 2 - 1

		- ACC(g)	Sample value in g.
		- ADC		Value sampled from the channel (raw value)
		- Cmin		Minimum calibration value when performing slow 360° rotation
		- Cmax 		Maximum calibration value when performing slow 360° rotation

	[RANGE]
	-3g to 3g

	Notes
	-----
	Cmin and Cmax should be calibrated for each sensors as indicated in the datasheet.

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	c_min : int
		Minimum calibration value.
	c_max : int
		Maximum calibration value.

	Returns
	-------
	acc_values : array of floats
		Converted ACC samples (g).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Conversion
	acc_samples = np.asarray([2 * ((float(s) - c_min) / (c_max - c_min) - .5) for s in samples])
	_check_ranges(acc_samples, 'ACC')
	return acc_samples


def temp(samples=None, resolution=10, unit='c'):
	"""Converts raw TMP values into original physical units (°C, °F or °K).

	BITalino Sensor: Temperature (TMP)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_TMP_Sensor_Datasheet.pdf

	[Transfer Function]
	TMP(°C) = (ADC/2^n - 0.5) * 100
	TMP(°F) = TMP(°C) * 9/5 +32
	TMP(°K)  = TMP(°C) - 273.15

		- TMP(°C)   TMP value in degrees Celsius (°C)
		- TMP(°F)   TMP value in degrees Fahrenheit (°F)
		- TMP(°K)    TMP value in Kelvin
		- ADC       Value sampled from the channel
		- n         Number of bits of the acquiring channel

	[RANGE]
	-40°C to 125°C

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution (6-bit or 10-bit; default: 10-bit).
	unit : char
		Unit of the output (c for degrees Celcius; f for degrees Fahrenheit; k for Kelvin)

	Returns
	-------
	temp_samples : array of floats
		Converted temperature values (°C, °F or °K).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute in degrees Celcius
	temp = np.asarray([(float(s) / 2 ** resolution - 0.5) * 100 for s in samples])

	if unit is 'f':
		# Convert to Fahrenheit
		temp = temp * 9. / 5 + 32
	elif unit is 'k':
		# Convert to Kelvin
		temp = temp - 273.15

	# Check samples
	_check_ranges(temp, 'temp')
	return temp


def ntc(samples=None, resolution=10, unit='c'):
	"""Converts raw NTC values into original physical units (°C, °F or K).

	BITalino Sensor: High-Definition Temperature Sensor (NTC)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/NTC_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	NTC(V) = ADC * VCC / 2^n
	NTC(ohm) = (1 * 10^4 * NTC(V)) / (VCC - NTC(V))
	TMP(°K) = 1 / (a0 + a1 * log(NTC(ohm)) + a2 * log(NTC(ohm))^3)
	TMP(°C) = NTC(°K) - 273.15
	TMP(°F) = NTC(°C) * 9/5 + 32

		- NTC(V)    NTC output in Volt (V)
		- NTC(ohm)  NTC output in Ohm (ohm)
		- TMP(°K)   Temperature value in Kelvin (°K)
		- TMP(°C)   Temperature value in Celcius
		- TMP(°F)   Temperature value in Fahrenheit
		- ADC       Value sampled from the channel
		- n         Number of bits of the channel
		- VCC       Operating voltage (3.3V)

		- a0 = 1.12764514 * 10^3
		- a1 = 2.34282709 * 10^4
		- a2 = 8.77303013 * 10^8

	[RANGE]
	0°C to 50°C

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).
	unit : char
		Unit of the output (c for degrees Celcius; f for degrees Fahrenheit; k for Kelvin).

	Returns
	-------
	temp_samples : array of floats
		Converted temperature samples (°C, °F or °K).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	ntc = samples * 3.3 / 2 ** resolution
	ntc = 10 ** 4 * ntc / (3.3 - ntc)

	a0 = 1.12764514 * 10 ** 3
	a1 = 1.34282709 * 10 ** 4
	a2 = 8.77303013 * 10 ** 8

	# Compute temperature in Kelvin
	temp = np.asarray(1 / (a0 + a1 * np.log(ntc) + a2 * np.log(ntc) ** 3))

	if unit is 'c':
		# Convert to Celcius
		temp = temp - 273.15
	elif unit is 'f':
		# Convert to Fahrenheit
		temp = (temp - 273.15) * 9. / 5 + 32

	_check_ranges(temp, 'ntc')
	return temp


def lux(samples=None, resolution=10):
	"""Converts raw LUX value into original unit (%).

	BITalino Sensor: Light (LUX)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_ECG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	LUX(%) = ADC / 2**n * 100

		- ADC   Value sampled from the channel
		- n     Number of bits of the acquiring channel

	[RANGE]
	0% to 100%

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	lux_samples : array of floats
		Converted LUX samples (mV).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check samples
	lux_samples = np.asarray([s * 100. / 2 ** resolution for s in samples])
	_check_ranges(lux_samples, 'LUX')
	return lux_samples


def osl(samples=None, resolution=10):
	"""Converts single, raw SpO2 Reader (OSL) value into original physical units
	(oxygen saturation in % or heart rate in bpm).

	BITalino Sensor: SpO2 Reader (OSL)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/telemedicine/CMS-50D_Plus.pdf

	[TRANSFER FUNCTION]
	SpO2(%) = 0.25 * 2**(10-n) * ADC - 0.8
	HR(bpm) = 0.25 * 2**(10-n) * ADC - 0.8

	[RANGE]
	SpO2: 0% to 100%
	HR: 30bpm to 250bpm

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	osl_samples : array of floats
		Converted SpO2 samples (%) or heart rate sample (bpm).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check input
	osl_samples = np.asarray([0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples])
	_check_ranges(osl_samples, 'OSL')
	return osl_samples


def bpr(samples=None, resolution=10):
	"""Converts single, raw Blood Pressure Reader (BPR) value into original physical units
	(blood pressure in mmHg).

	BITalino Sensor: Blood Pressure Reader (BPR)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/telemedicine/Goodwill_Studio_Ver1.pdf

	[TRANSFER FUNCTION]
	BPR(mmHg) = 0.25 * 2**(10-n) * ADC - 0.8


	[RANGE]
	0mmHg to 300mmHg

	--

	Parameters
	----------
	samples: int
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	bpr_samples : float
		Converted blood pressure samples (mmHg).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check samples
	bpr_samples = np.asarray([0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples])
	_check_ranges(bpr_samples, 'BPR')
	return bpr_samples


def gmr(samples=None, resolution=10):
	"""Converts raw Glucose Meter Reader (GMR) samples into original physical units
	(Glucose level in mg/dL).

	BITalino Sensor: Glucose Meter Reader (GMR)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/telemedicine/caresens_ii.pdf

	[TRANSFER FUNCTION]
	GMR(mg/dL) = 0.25 * 2**(10-n) * ADC - 0.8

	[RANGE]
	20mg/dL to 600mg/dL

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	gmr_samples : array of floats
		Glucose level samples (mg/dL).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check samples
	gmr_samples = np.asarray([0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples])
	_check_ranges(gmr_samples, 'GMR')
	return gmr_samples


def eog(samples=None, resolution=10):
	"""Converts raw EOG values into original physical unit (mV).

	BITalino Sensor: Electrooculography (EOG)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_EOG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	ECG(mV) = ((((ADC / (2 ** n)) - 0.5) * VCC) / G_EOG) * 1000

		- ADC   Value sampled from the channel
		- n     Number of bits of the acquiring channel
		- VCC   Operating voltage (3.3V)
		- G_EOG Sensor gain (2040)

	[RANGE]
	-0.81mV to 0.81mV

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor values.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	eog_samples : array of floats
		Converted EOG samples (mV).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check samples
	eog_samples = np.asarray([((((float(s) / (2 ** resolution)) - 0.5) * 3.3) / 2040) * 1000 for s in samples])
	_check_ranges(eog_samples, 'EOG')
	return eog_samples


def egg(samples=None, resolution=10):
	"""Converts raw EGG values into original physical unit (mV).

	BITalino Sensor: Electrogastrography (EGG)

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_EGG_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	ECG(mV) = ((((ADC / (2 ** n)) - 0.5) * VCC) / G_EGG) * 1000

		- ADC   Value sampled from the channel
		- n     Number of bits of the acquiring channel
		- VCC   Operating voltage (3.3V)
		- G_EGG Sensor gain (6110)

	[RANGE]
	-0.27mV to 0.27mV

	--

	Parameters
	----------
	samples : array of ints
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	egg_sample : array of floats
		Converted EGG samples (mV).

	Raises (via _check_input())
	---------------------------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute & check samples
	egg_samples = np.asarray([((((float(s) / (2 ** resolution)) - 0.5) * 3.3) / 6110) * 1000 for s in samples])
	_check_ranges(egg_samples, 'EGG')
	return egg_samples


def _check_input(samples, resolution=None):
	"""Checks if provided input data is valid.

	Parameters
	----------
	samples : array
		Raw sensor values.
	resolution : sampling resolution
		Sampling resolution.
	Returns
	-------
	samples : np.ndarray
		Raw sensors values in a NumPy array.

	Raises
	-----
	TypeError
		If no samples provided in input data.
	TypeError
		If unsupported data format provided (no list or np.ndarray)
	ValueError
		Unsupported sampling resolution (not 6-bit or 10-bit).
	"""
	# Check samples
	if samples is None:
		raise TypeError("No sensor samples provided.")
	elif type(samples) is list:
		samples = np.asarray(samples)
	elif type(samples) is not np.ndarray:
		raise TypeError("Unsupported data format")

	# Check resolution
	if resolution is not 6 and not 10:
		raise ValueError("Unsupported sampling resolution for BITalino devices (6-bit or 10-bit supported)")

	return np.asarray(samples)


def _check_ranges(signal, sensor):
	"""Checks if the computed sensor values lie within the sensor's range and raises warnings if otherwise.

	Parameters
	----------
	samples : array
		Raw sensor values.
	resolution : sampling resolution
		Sampling resolution.

	"""
	if np.min(signal) < ranges[sensor][0]:
		warnings.warn("Minimum value of the %s sensor data is smaller than the min range of this sensor (%f < %f)."
			% (sensor, np.min(signal), units[sensor], ranges[sensor][1], units[sensor]), stacklevel=2)
	if np.max(signal) > ranges[sensor][1]:
		warnings.warn("Maximum value of the %s sensor data is greater than the max range of this sensor (%f%s > %f%s)."
			 % (sensor, np.max(signal), units[sensor], ranges[sensor][1], units[sensor]), stacklevel=2)


# DICTIONARIES
units = {
	'RAW': '-',
	'ECG': 'mV',
	'EEG': u'µV',
	'EMG': 'mV',
	'EDA': u'µS',
	'ACC': 'g',
	'GLUC': 'mg/dL',
	'TEMP': u'°C',
	'NTC': u'°C',
	'LUX': '%',
	'HR': 'bpm',
	'OXI': '%',
	'SYS': 'mmHg',
	'DIA': 'mmHg',
	'EOG': 'mV',
	'EGG': 'mV'
}

ranges = {
	'RAW': [0, 2**10],
	'ECG': [-1.5, 1.5],
	'EEG': [-41.24, 41.25],
	'EMG': [-1.65, 1.65],
	'EDA': [-4.4, 21],
	'ACC': [-3, 3],
	'GLUC': [20, 600],
	'TEMP': [-40, 125],
	'NTC': [0, 50],
	'LUX': [0, 100],
	'HR': [30, 250],
	'OXI': [0, 100],
	'SYS': [0, 300],
	'DIA': [0, 300],
	'EOG': [-0.81, 0.81],
	'EGG': [-0.27, 0.27]
}

transfer_functions = {
	'ECG': ecg,
	'EEG': eeg,
	'EMG': emg,
	'EDA': eda,
	'ACC': acc,
	'GLUC': gmr,
	'TEMP': temp,
	'NTC': ntc,
	'LUX': lux,
	'HR': osl,
	'OXI': osl,
	'SYS': gmr,
	'DIA': gmr,
	'EOG': eog,
	'EGG': egg
}


if __name__ == '__main__':
	"""Example Script

	For newest sample data visit:
		https://github.com/BITalinoWorld/revolution-sample-data

	For more information & documentation visit:
		http://bitalino.com/en/learn/documentation
	"""
	# Load data using NumPy
	signal = np.loadtxt('SampleECG.txt')[:, -1]

	# Convert raw ECG data to mV using the transfer function at 10-bit resolution
	ecg_signal = ecg(signal, 10)

	# Alternatively (10-bit set as default resolution)
	ecg_signal = ecg(signal)

	# Convert raw ECG data to mV using the transfer function at 6-bit resolution
	ecg_signal = ecg(signal, 6)
