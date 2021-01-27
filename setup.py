from distutils.core import setup
import setuptools

setup(
    name='pyloads-wind-turbine',  # How you named your package folder (MyLib)
    packages=setuptools.find_packages(),  # Chose the same as "name"
    version='1.0.1',  # Start with a small number and increase it with every change you make
    license='MIT',  # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    description='Loads and deflection for wind turbine blade.',  # Give a short description about your library
    author='Gonzalo Mazzini',  # Type in your name
    author_email='gmazzini@itba.edu.ar',  # Type in your E-Mail
    url='https://github.com/GonMazzini/pyloads/tree/v1.0.1',  # Provide either the link to your github or to your website
    download_url='https://github.com/GonMazzini/pyloads/archive/v1.0.1.tar.gz',  # I explain this later on
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)