# opensignalsreader
Python package to read OpenSignals (r)evolution files and automatic sensor data conversion for BITalino (r)evolution sensor data.

### Examples & How to Use this Package
Reading OpenSignals (r)evolution files.
```python
from opensignalsreader import OpenSignalsReader
acq = OpenSignalsReader('SampleECG.txt')

# alternatively
acq = OpenSignalsReader()
# then later in code
acq.read_file('SampleECG.txt')
```
If required, no data conversion will be conducted, instead, the original raw digital values are stored.
```python
acq = OpenSignalsReader('SampleECG.txt', raw=True)
```
The OpenSignalsReader class comes with plotting features which figure size does automatically adjust to the number of plotted signals at a time.
```python
acq = OpenSignalsReader('SampleECG.txt', show=())

# alternatively
acq = OpenSignalsReader('SampleECG.txt)
acq.plot()
```
If you want to plot specific signals, use the _.sensors_ variable to get the list of available sensor signal and pass the required sensor keys (str) to the _.plot()_ method.
```python
# Read Data
acq = OpenSignalsReader('opensignalsfile.txt')

# Get sensor keys
sensors = acq.sensors   # Example Return: ['ECG', 'EMG', 'EDA']

# Plot all available signals
acq.plot()

# Plot only ECG signal
acq.plot('ECG')

# Plot ECG & EDA signal
acq.plot(['ECG', 'EDA])
```

### BITalino (r)evolution Transfer Functions
This package includes the _bitalino_ module which contains all available transfer functions of the current BITalino (r)evolution sensors.
It is used by the _OpenSignalsReader_ class to convert raw signal samples imported from OpenSignals files.

This package can also be useful if you want to convert sensor signals within your own software when not importing signals from the OpenSignals files.

BITalino sample series can be converted into their original units using the sensor's transfer function. See below how to use the functions of the _bitalino_ module on the example of the ECG sensor.
```python
import numpy as np
from opensignalsreader.bitalino import ecg
signal = np.loadtxt('SampleECG.txt', 'r')[:, -1]

# Convert signal
ecg_signal = ecg(signal)

# Convert signal acquired with 6-bit sampling resolution
ecg_signal = ecg(signal, 6)
```

### Useful Links
Detailed documentation about BITalin (r)evolution sensors can be found here:
http://bitalino.com/en/learn/documentation

BITalino (r)evolution sample files can be found here:
https://github.com/BITalinoWorld/revolution-sample-data 

### Dependencies
- matplotlib
- numpy

### Context of this Work
This package is part of the master thesis	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)" at the University of Applied Sciences Hamburg, Germany.
