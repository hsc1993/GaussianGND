# Gaussian GND [![DOI](https://zenodo.org/badge/470862080.svg)](https://zenodo.org/badge/latestdoi/470862080)

This repository superimposes atomistic data for dislocations from MD simulation onto spatial tessellations, calculating mesh size dependent GND (geometrically necessary dislocation) signal. 

The project consists of three parts: 
1. Analyzing **atomistic data** using OVITO and export TXT files containing dislocation information.
2. Generating **spatial tessellation** using open-source library.
3. Using **cell-edge detection process** to overlay dislocation information with spatial tessellation to generate GND signal.


## Table of Contents

- [Install](#install)
- [Usage](#usage)
- [Examples](#example)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)


## Install

This project consists of code developed in [Python](https://www.python.org/) and [MATLAB](https://www.mathworks.com/products/matlab.html). Commercial Software [OVITO](https://www.ovito.org/) is used for data processing. 


## Usage

### Dislocation extraction
Starting from the [sample MD result]() of Fe grain boundary with a prismatic loop, we use OVITO to export XYZ file that contains coordinates of all particles. The set up for OVITO and the dislocation extraction image are as shown below:

<img src="ovito_setup.png" width="300" height="700">     <img src="dislocation.png" width="300" height="700"> 

The output file should look like the [dislocation.XYZ](). Then the python script [dxa_analysis.py](dxa_analysis.py), which is built upon [OVITO's Python interface](https://docs.ovito.org/python/) extracts dislocation information from the XYZ file and outputs dislocation data for the types '100', '110', '111' and other into txt files. 

### Spatial tessellation
There are two types of meshes available: [Voronoi Tessellation](tessellation_voronoi.m) and [regular Hexahedron Tesselation](tessellation_cubic.m). To construct spatial tessellation, run the corresponding MATLAB code and the coordinates of all volume mesh vertices will be saved into text files in designated directories.

### Gaussian GND calculation
The [main.py](main.py) calculates the GND density by combining the dislocation information and meshed volumes. The output from this script includes: 
1. Simulation parameters like the mesh size and the dimension of the simulation volume; 
2. Volume of each meshed volume element; 
3. GND signal intensity of each meshed volume element; 
4. Vertices of dislocation segments with the truncating effect of meshed volume elements; 
5. The correspondence between each dislocation segment and the meshed volume element that fully contains it.

### Analyze GND calculation

## Example



## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by research group at UCLA, Johns Hopkins University, Hongkong City University and Pennsylvania University.


## License

[MIT](LICENSE) © Sicong He
