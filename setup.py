# OPENSIGNALS SETUP SCRIPT

import setuptools
from opensignalsreader import __author__, __version__, __email__, name, description
with open("README.md", "r") as fh:
	long_description = fh.read()

# Create setup
setuptools.setup(
	name=name,
	version=__version__,
	author=__author__,
	author_email=__email__,
	description=description,
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/PGomes92/opensignalsreader",
	keywords=['opensignals', 'opensignalsreader', 'physiological signals', 'BITalino', 'biosignalsplux'],
	setup_requires=[
		'numpy',
		'matplotlib'
	],
	install_requires=[
		'numpy',
		'matplotlib'
	],
	packages=setuptools.find_packages(),
	package_data={
		"opensignalsreader": [
			"files/*",
			"README.md"]},
	include_package_data=True,
	classifiers=[
		'Intended Audience :: Developers',
		'Intended Audience :: Education',
		'Intended Audience :: Science/Research',
		'Natural Language :: English',
		'License :: OSI Approved :: BSD License',
		'Programming Language :: Python',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Operating System :: OS Independent',
	],
)
