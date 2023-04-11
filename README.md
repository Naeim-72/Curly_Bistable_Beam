PyAnsys Finite Element Bistability Analysis for Curly Beams
This repository contains a finite element script based on the PyAnsys package for performing bistability analysis of curly beams. Bistability refers to the phenomenon where a structure can have two stable equilibrium states, and curly beams are known for exhibiting bistability behavior due to their unique geometric shape.

The script in this repository provides a computational approach for analyzing the bistability of curly beams using finite element method (FEM) simulations. PyAnsys is a Python library that interfaces with the Ansys Mechanical software, a popular tool for finite element analysis (FEA). The script takes advantage of PyAnsys' capabilities to generate and solve finite element models of curly beams, and it provides a convenient way to analyze the structural behavior of these beams under different loading conditions.

Dependencies
To use the script, you need to have the following dependencies installed:

Python 3.8
PyAnsys package (version 0.59.4 or higher)
Ansys Mechanical software (version 2019 or higher)
Installation
Clone the repository to your local machine using git clone.
Make sure you have Python 3.x installed. If not, you can download it from the official Python website (https://www.python.org/).
Install PyAnsys package by running pip install pyansys in your terminal or command prompt.
Make sure you have Ansys Mechanical software installed on your machine. If not, you can download it from the Ansys website (https://www.ansys.com/).
Open the script in your preferred Python environment (e.g., Jupyter Notebook, Spyder, etc.).
Usage
Update the curly beam geometry and material properties in the script as needed.
Define the boundary conditions and loads applied to the curly beam.
Run the script to generate the finite element model and solve for the displacements and reactions of the curly beam.
Analyze the results to determine the bistability behavior of the curly beam, such as the critical load for bistability and the deformation of the beam in different equilibrium states.
Modify the parameters in the script to perform sensitivity analyses or parametric studies.
Examples
The repository includes example scripts that demonstrate how to use the PyAnsys finite element script for bistability analysis of curly beams. These examples provide step-by-step instructions on how to modify the script for different beam geometries, material properties, and loading conditions, and how to analyze the results to obtain meaningful insights into the bistability behavior of the curly beams.

Contributing
If you would like to contribute to this repository, you are welcome to submit pull requests with bug fixes, feature enhancements, or other improvements. Please provide clear documentation and explanation of your changes to facilitate the review process.

License
This repository is released under the MIT License, which means you are free to use, modify, and distribute the code for both commercial and non-commercial purposes. However, please note that the PyAnsys package and Ansys Mechanical software may have their own licensing requirements, and you should comply with their terms and conditions when using this script.

Acknowledgements
This work was inspired by the research on bistability of curly beams, and it would not have been possible without the PyAnsys package and Ansys Mechanical software. We acknowledge the developers of these tools for providing valuable resources that enable researchers and engineers to study and understand the behavior of complex structures like curly beams.

References
Reference 1: Ghavidelnia, Naeim and Yin, Kaiyang and Cao, Bo and Eberl, Christoph, Curly Beam with Programmable Bistability. Available at SSRN: https://ssrn.com/abstract=4381506 or http://dx.doi.org/10.2139/ssrn.4381506
Reference 2: [Paper title or website link]
Please cite the relevant references when using or referencing this script in your
