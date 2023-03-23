# Gaussian GND

This repository superimposes atomistic data from MD simulation onto spatial tessellations, calculating mesh size dependent GND (geometrically necessary dislocation) signal. 

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
Once atomistic results are obtained from MD, use OVITO to export XYZ file that contains coordinates of all particles. Then the python script [dxa_analysis.py](dxa_analysis.py), which is built upon, [OVITO's Python interface](https://ovito.org/manual/python/introduction/running.html) extracts dislocation information from the XYZ file. 

### Spatial tessellation
There are two types of meshes available: [Voronoi Tessellation](full_voronoi_random) and [regular Hexahedron Tesselation](full_voronoi_edgevariate). To construct spatial tessellation, run the MATLAB code and the coordinates of all volume mesh vertices will be saved into text files.

### Gaussian GND calculation
The [main.py](main.py) calculates the GND density by combining the dislocation information and meshed volumes. The output from this script includes: 
1. Simulation parameters like the mesh size and the dimension of the simulation volume; 
2. Volume of each meshed volume element; 
3. GND signal intensity of each meshed volume element; 
4. Vertices of dislocation segments with the truncating effect of meshed volume elements; 
5. The correspondence between each dislocation segment and the meshed volume element that fully contains it.


## Example
https://user-images.githubusercontent.com/56003395/227333746-c0249328-7a25-4444-899a-64329008e719.mp4



## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by research group at UCLA, Johns Hopkins University, Hongkong City University and Pennsylvania University.


## License

[MIT](LICENSE) © Sicong He
