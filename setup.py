from setuptools import setup, find_packages

def requirements_from_file(file_name):
	return open(file_name).read().splitlines()

setup(
	names="ThetaRecov",
	version="0.0.1",
	packages=find_packages(),
	install_requires=requirements_from_file('requirements.txt'),
	entry_points={
		'console_scripts':[
			'tajimaD_windows=ThetaRecov.cli_windows:main',
			'tajimaD_overall=ThetaRecov.cli_overall:main',
		],
	},
	description="A package for correct computation of theta and Tajima's D under missing data",
	author = 'Akira S. Hirao',
	author_email='akihirao@gmail.com',
	url='https://github.com/akihirao/ThetaRecov'
)
