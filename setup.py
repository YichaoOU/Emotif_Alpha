from setuptools import setup, find_packages

setup(
    name='Emotif_alpha',
    version='1.0',
    author='Yichao Li',
    author_email='yl079811@ohio.edu',
    packages=['Emotif'],
	scripts=['Emotif_alpha'],
    url='https://github.com/YichaoOU/Emotif_Alpha',
    license='LICENSE',
	package_data={'': ["*","algo/*","_dataset/*","_templates/*","utils/*"]},
	include_package_data=True,
    description='Ensemble motif discovery',
)