
.. _label-quickstart:

Quick Start
==================

This a quick start guide to getting going quickly with the 1CPN model.


Getting 1CPN and compiling with LAMMPS 
---------------------------------------

1. Install git

::

    sudo apt-get install git

2. Navigate to a directory where you would like to install the 1CPN Model.

::

    mkdir 1cpn
    cd 1cpn
    echo "export D_1CPN=`pwd`" >> ~/.bashrc
    source ~/.bashrc


3. Download a fresh copy of LAMMPS and checkout the `stable_16Mar2018` tag. Even if you already have a copy of LAMMPS, its recommended that you download a new copy specifically for 1CPN.  
::

  git clone -b stable https://github.com/lammps/lammps.git lammps-1cpn
  git checkout stable_16Mar2018


4. Clone the 1CPN model

:: 

  git clone https://github.com/lequieu/1cpn-model.git

5. Link 1CPN with LAMMPS

::

  cd $D_1CPN/1cpn-model/src/lammps
  make link

The Makefile assumes that the LAMMPS src code is located at ${D_1CPN}/lammps-1cpn. If you have LAMMPS located at another location, you can specify by redefining the ``LAMMPS_SRC`` variable.

::

  make link LAMMPS_SRC=<path to LAMMPS src>


6. Build LAMMPS with the `ASPHERE` and `MOLECULE` package.

::

  cd $D_1CPN/lammps-1cpn/src
  make yes-ASPHERE
  make yes-MOLECULE
  make serial



Running your first simulation
-------------------------------

Now that LAMMPS has been linked and compuled with the 1CPN-model, we're ready to setup and analyze our first 1CPN simulation. First lets make a new directory called `example` where the input and output files from this simulation will be stored. After the directory is made, navigate into it.

::

    cd ${D_1CPN}/example
    mkdir -p ${D_1CPN}/example


Generating 1CPN Input Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. The next step is to generate the initial configuration of the 1CPN-model. This is achieved using the `${D_1CPN}/init/init_1cpn.py` script. To start out, we'll generate a section of chromatin consisiting of 20 nucleosomes, each with a nucleosome repeat length (NRL) of 187 base pairs. From the `${D_1CPN}/example` directory, issue the following command:

::

    ${D_1CPN}/1cpn-model/init/init_1cpn.py -n 20 -nrl 187

The full list of arguements that can be passed to ``init_1cpn.py`` see :ref:`label-initialization` or issue ``init_1cpn.py`` with the ``-h`` flag:

::

    ${D_1CPN}/1cpn-model/init/init_1cpn.py -h

Running a Simulation
^^^^^^^^^^^^^^^^^^^^^^

2. Next we need to copy the LAMMPS input files necessary for setting up and running a simulation. These files are located at `${D_1CPN}/inputs`. Copy them to our `example` directory.

::

    cp ${D_1CPN}/1cpn-model/inputs/in.* .

3. Now we're ready to run a simple simulation. This is as simple as executing

::

    ${D_1CPN}/lammps-1cpn/src/lmp_serial -i in.1cpn

Users unfamilair with LAMMPS are referred to the `LAMMPS docmentation <https://lammps.sandia.gov/doc/Manual.html>`_ for descriptions on LAMMPS output, and how to run different simulations.


Congratulations! You've just performed (hopefully) your first simulation with 1CPN!


Vizualizing the Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that you've run your simulation. You'll probably want to visualize it. 1CPN comes packaged with several different vizualization options. To learn more, check out :ref:`label-viz`.

Analyzing the Simulation
---------------------------

Visualizing a simulation is fun, but you'll probably want to perform some sort of analysis on it. 
For example, you might want to compute the the end-to-end distance if a single chromatin fiber, or the distance (or angle) between two nucleosomes.

To perform these sorts of analysis, 1CPN comes packaged with a variety of analysis scripts. See more at :ref:`label-analysis`.

