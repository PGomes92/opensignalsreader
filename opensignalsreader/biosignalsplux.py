

# DICTIONARIES
units = {
	'RAW': '-',
	'ECG': 'mV',
	'EEG': u'µV',
	'EMG': 'mV',
	'EDA': u'µmS',
	'ACC': 'g',
	'PZT': '%',
	'FSR': 'lb',
	'RIP': '%',
	'GON': u'°',
	'LOAD': '-',
	'BVP': u'µA',
	'PLATFORM': 'kgf',
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
	'PZT': [0, 100],
	'FSR': [0, 150],
	'RIP': [0, 100],
	'GON': [-180, 180],
	'LOAD': [0, 110],
	'BVP': [0, 20],
	'PLATFORM': [0, 200],
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
	'PZT': pzt,
	'FSR': fsr,
	'RIP': rip,
	'GON': gon,
	'LOAD': load,
	'BVP': bvp,
	'PLATFORM': platform,
	'GLUC': gmr,
	'TEMP': temp,
	'NTC': ntc,
	'LUX': lux,
	'HR': osl,
	'OXI': osl,
	'SYS': gmr,
	'DIA': gmr
}
