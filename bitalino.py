#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BITalino (r)evolution Transfer Functions
----------------------------------------
Visit http://bitalino.com/

Provides functions to convert raw BITalino (r)evolution sensor samples (single samples
or entire datsets) into original physical units and supports the following sensors:

	* ECG		* EEG		* EMG 		* EDA 		* ACC	* GMR
	* TEMP		* NTC		* LUX 		* OSL		* BPR

Visit the following page for more detailed information about the transfer functions:
http://bitalino.com/en/learn/documentation

Visit the following page for BITalino (r)evolution sample signals:
https://github.com/BITalinoWorld/revolution-sample-data

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
import warnings
import numpy as np


def ecg(samples=None, resolution=10):
	"""Converts raw ECG values into original physical unit (mV).

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
	samples : array
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	ecg_sample : float
		Converted ECG sample (mV).

	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute output
	return [((((float(s) / (2 ** resolution)) - 0.5) * 3.3) / 1100) * 1000 for s in samples]


def eeg(samples, resolution=10):
	"""Converts single, raw EEG value into original physical unit (µV).

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/EEG_Sensor_Datasheet.pdf

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
	samples : array
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [((((float(s) / 2 ** resolution) - 0.5) * 3.3) / 40000) * (10 ** 6) for s in samples]


def emg(samples=None, resolution=10):
	"""Converts single raw EMG value into original physical unit (mV).

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
	samples : int
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	emg_sample : float
		Converted EMG sample.

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [((((float(s) / 2 ** resolution) - 0.5) * 3.3) / 1009) * 1000 for s in samples]


def eda(samples=None, resolution=10):
	"""Converts raw EDA values into original physical unit (µS).

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
	samples : int
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit; default: 10-bit).

	Returns
	-------
	eda_sample : float
		Converted EDA sample.

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [(((float(s) / (2 ** resolution)) * 3.3) - 0.574) / 0.132 for s in samples]


def acc(samples=None, resolution=10, c_min=28000, c_max=38000):
	"""Converts raw acc value(s) into original units.

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/REVOLUTION_ACC_Sensor_Datasheet.pdf

	[TRANSFER FUNCTION]
	ACC(g) = ((ADC - Cmin)/(Cmax - Cmin)) * 2 - 1

		- ACC(g)	Sample value in g.
		- ADC		Value sampled from the channel (raw value)
		- Cmin		Minimum calibration value when performing slow 360° rotation
		- Cmax 		Maximum calibration value when performing slow 360° rotation

	[RANGE]
	-3g to 3g

	--

	Parameters
	----------
	samples : int
		Raw digital sensor value.
	c_min : int
		Minimum calibration value.
	c_max : int
		Maximum calibration value.

	Returns
	-------
	acc_value : float
		Converted ACC sample (g).

	Raises
	------
	TypeError
		If no input data is provided.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	if c_min is 28000 or c_max is 38000:
		warnings.warn("Using standard calibration values. \
		Sensor calibration is recommended.")

	# Conversion
	return [((float(s) - c_min) / (c_max - c_min)) * 2 - 1 for s in samples]


def temp(samples=None, resolution=10, unit='c'):
	"""Converts single, raw TMP value into original physical units (°C, °F or °K).

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
	samples : int
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution (6-bit or 10-bit; default: 10-bit).
	unit : char
		Unit of the output (c for degrees Celcius; f for degrees Fahrenheit; k for Kelvin)

	Returns
	-------
	temp : float
		Converted temperature value (°C, °F or °K).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	# Compute in degrees Celcius
	temp = [(float(s) / 2 ** resolution - 0.5) * 100 for s in samples]

	if unit is 'f':
		# Convert to Fahrenheit
		temp = temp * 9. / 5 + 32
	elif unit is 'k':
		# Convert to Kelvin
		temp = temp - 273.15

	return temp


def ntc(samples=None, resolution=10, unit='c'):
	"""Converts single, raw NTC value into original physical units (°C, °F or K).

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
	samples : int
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).
	unit : char
		Unit of the output (c for degrees Celcius; f for degrees Fahrenheit; k for Kelvin).

	Returns
	-------
	temp : float
		Converted temperature sample (°C, °F or °K).

	Raises
	------
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
	temp = 1 / (a0 + a1 * np.log(ntc) + a2 * np.log(ntc) ** 3)

	if unit is 'c':
		# Convert to Celcius
		temp = temp - 273.15
	elif unit is 'f':
		# Convert to Fahrenheit
		temp = (temp - 273.15) * 9. / 5 + 32

	return temp


def lux(samples=None, resolution=10):
	"""Converts single, raw LUX value into original unit (%).

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
	samples : array
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	lux_sample : float
		Converted LUX sample (mV).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [s * 100. / 2 ** resolution for s in samples]


def osl(samples=None, resolution=10):
	"""Converts single, raw SpO2 Reader (OSL) value into original physical units
	(oxygen saturation in % or heart rate in bpm).

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
	samples : array
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	converted_sample : float
		Converted SpO2 sample (%) or heart rate sample (bpm).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples]


def bpr(samples=None, resolution=10):
	"""Converts single, raw Blood Pressure Reader (BPR) value into original physical units
	(blood pressure in mmHg).

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
	bpr : float
		Converted blood pressure samples (mmHg).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples]


def gmr(samples=None, resolution=10):
	"""Converts raw Glucose Meter Reader (GMR) samples into original physical units
	(Glucose level in mg/dL).

	See sensor datasheet for more information:
	http://bitalino.com/datasheets/telemedicine/caresens_ii.pdf

	[TRANSFER FUNCTION]
	GMR(mg/dL) = 0.25 * 2**(10-n) * ADC - 0.8

	[RANGE]
	20mg/dL to 600mg/dL

	--

	Parameters
	----------
	samples : int
		Raw digital sensor value.
	resolution : int, optional
		Sampling resolution used during acquisition (6-bit or 10-bit).

	Returns
	-------
	gmr : float
		Glucose level (mg/dL).

	Raises
	------
	TypeError
		If no input data is provided.
	Warnings
		If sampling resolution other than 6-bit or 10-bit is being used.
	"""
	# Check input
	samples = _check_input(samples, resolution)

	return [0.25 * 2 ** (10 - resolution) * float(s) - 0.8 for s in samples]


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


# DICTIONARIES
units = {
	'RAW': '-',
	'ECG': 'mV',
	'EEG': u'µV',
	'EMG': 'mV',
	'EDA': u'µmS',
	'ACC': 'g',
	'GLUC': 'mg/dL',
	'TEMP': u'°C',
	'NTC': u'°C',
	'LUX': '%',
	'HR': 'bpm',
	'OXI': '%',
	'SYS': 'mmHg',
	'DIA': 'mmHg'
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
	'DIA': [0, 300]
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
	'DIA': gmr
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
