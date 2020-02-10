import setuptools

setuptools.setup(
    name="emsplot-martin",  # local package
    version="1.0.0",
    author="Martin",
    description="A small visualization package for NetCDF files from the EMS",
    packages=setuptools.find_packages(),
    license='WTFPL',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=["cartopy",  # numpy is a dependency of cartopy
                      "matplotlib",
                      "xarray"],  # actually we could manage with just NetCDF4
    python_requires='>=3.5',
)
