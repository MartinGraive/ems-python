import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="emsplot-martin",  # local package
    version="1.1.1",
    author="Martin",
    description="A small visualization package for NetCDF files from the EMS",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    license='WTFPL',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=["cartopy",  # numpy is a dependency of cartopy
                      "matplotlib",
                      "netCDF4",
                      "xarray"],  # actually we could manage with just NetCDF4
    python_requires='>=3',
)
