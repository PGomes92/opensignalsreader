#### Update v.0.2.2 - 06.06.2019
- improved support for Python 3 (fixes #2)
- fixed pip install issues under Python 3
- added setup.py
- restructured package
	- transfer functions are now found in the submodule 'transfer_functions'
	- 'bitalino_tf.py' renamed to 'bit.py'
	- 'biosignalsplux.py' renamed to 'bio.py'

#### Update v.0.2.1 - 13.10.2018
- added changelog!
- added installation using `pip
- fixed unsuitable fig size for plot figures with only one plot (previous size did not visualize the x-label correctly)
- added support for new BITalino sensors
  - Electrooculography (EOG)
  - Electrogastrography (EGG)
- added datasheet link for EEG sensor
- added `_check_ranges()` function that checks if the computed converted sensor samples lie within the technical ranges of the sensor being used (warnings are triggered if this is not the case)
- changed License type (see LICENSE.txt)