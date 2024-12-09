from setuptools import setup, find_packages

REQUIREMENTS = ["cadquery", "numpy==1.24.4", "nlopt==2.7.1", "spyder", "spyder-kernels==2.5.*", "h5py"]

setup(
      name="mesh_from_CAD",
      version="0.1.0",
      author="Soner Bilyaz",
      author_email="un24029@bristol.ac.uk",
      description="Python Code which generates structured meshes with the cadquery library",
      long_description=open("README.md").read(),                                # Describe it in detail with a README file
      packages=find_packages(),                                                 # Automatically find and include all packages and subpackages
      install_requires=REQUIREMENTS,                                            # List of dependencies
      python_requires=">=3.10",                                                 # Specify supported Python versions
      
      ### Extra specifications ###
      url="",                                       # Repository url
      
      classifiers=[                                     # IT IS JUST REQUIRED FOR UPLOADING IT TO PyPI. Has a certain syntax format, but when you type wrong PyPI doesnt give an error, just ignores it.
        "Programming Language :: Python :: 3.10",       # Compatible Python version
        "License :: OSI Approved :: MIT License",       # License type
        "Operating System :: OS Independent",           # OS compatibility
    ],
      package_data={
    "package_tutorial": [],        # Data that are required to run the package
    },
      entry_points={
        "console_scripts": [
            "mesh=Cadquery_Get_coords_from_CAD",      # Define a console command to mesh the step file
        ],
      }
      
      )

