from setuptools import setup, find_packages

setup(
    name='marco_stuff',
    version='1.0.0',
    author='Marco Larranaga',
    description='Marco stuff',
    packages=find_packages(),    
    install_requires=['numpy>=1.18.1',
        'netCDF4>=1.5.3'],
)
