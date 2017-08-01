1CPN Model
====================

This is the 1CPN model for Chromatin. Yay!

If you haven't yet read the 1CPN manuscript (CITE). Go read it now.


Getting Started
--------------------
The best way to get started with 1CPN is to build the documentation (see below), and follow the "Quick Start" tutorial. This section will explain everything you need to know to use 1CPN model including:

* Compiling LAMMPS with 1CPN
* Generating an initial 1CPN configuration and input files
* Running a Simple Simulation 
* Vizualizing and Analyzing Results

If you're impatient, and don't want to build the documentation, you can just open up ``./docs/quickstart.rst`` in your favorite text editor and follow the "Quick Start" turotial from there.


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



Build Documentation
--------------------

First you'll need to install `Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`. Which should be included in your Linux distributions repository, and can be installed as simply as ::

  # apt-get install python-sphinx

If this doesn't work for you, visit the Sphinx website for information about setting Sphinx up on your machine. 

Next, making the documentation is a breeze. Simply execute ::

  $ cd doc
  $ make html
  $ firefox _build/html/index.html


