#### Update 13.10.2018
- added changelog!
- added installation using `pip
- fixed unsuitable fig size for plot figures with only one plot (previous size did not visualize the x-label correctly)
- added support for new BITalino sensors
  - Electrooculography (EOG)
  - Electrogastrography (EGG)
- added datasheet link for EEG sensor
- added `_check_ranges()` function that checks if the computed converted sensor samples lie within the technical ranges of the sensor being used (warnings are triggered if this is not the case)
- changed License type (see LICENSE.txt)