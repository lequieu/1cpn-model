1CPN Model
====================

Welcome to the LAMMPS implementation of the 1CPN model of chromatin. 

The first thing you'll want to do is familiarize youself with the 1CPN manuscript:

Lequieu, Cordoba, Moller, de Pablo "1CPN: A coarse-grained multi-scale model of chromatin" (2019) `J. Chem. Phys. 150, 215102 <https://doi.org/10.1063/1.5092976>`_


Getting Started
--------------------
The best way to get started with 1CPN is go to the `documentation <https://1cpn-model.readthedocs.io/>`_, and follow the "Quick Start" tutorial. This section will explain everything you need to know to use 1CPN model including:

* Compiling LAMMPS with 1CPN
* Generating an initial 1CPN configuration and input files
* Running a Simple Simulation 
* Vizualizing and Analyzing Results



What's Included
--------------------
* **src/**
    * **lammps/** - Core of 1CPN model, as implemented into LAMMPS
    * **viz/** - Tools for vizualizing 1CPN in VMD
    * **analysis/** - Tools for analyzing a 1CPN trajectory
    * **include/** - Includes for .cpp code 
* **init/** - Scripts to generate initial configuration for 1CPN
* **inputs/** - LAMMPS input scripts for 1CPN model
* **doc/** - Documentation written in .rst and Sphinx (for building, see below)
* **test/** - Testing Code
* **utils/** - Utilities for 1CPN model
* **bin/** - Directory where built binary files are written


Issues? 
--------------------

If you run into any issues don't hesitate to reach out. We'd love to help you get up and running.
