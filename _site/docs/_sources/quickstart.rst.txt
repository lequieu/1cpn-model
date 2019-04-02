
Quick Start
==================

This a quick start guide to getting going quickly with the 1CPN model.

Getting 1CPN 
---------------------------
1. Navigate to a directory where you would like to install the 1CPN Model. Optionally, you can create add a new bash variable `D_1CPN` so that 1CPN codes can be easily accessed. The path `D_1CPN` will be used throughout this tutorial.

::

    mkdir 1cpn
    cd 1cpn
    echo "export D_1CPN=`pwd`" >> ~/.bashrc
    source ~/.bashrc

2. Install `git` (if you don't have it) and clone the 1CPN-model.

::

    sudo apt-get install git
    git clone https://lequieu@bitbucket.org/lequieu/1cpn-model.git
    

Compiling LAMMPS with 1CPN
---------------------------

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
  git checkout -b stable_16Mar2018


4. Clone the 1CPN model

:: 

  git clone https://lequieu@bitbucket.org/lequieu/1cpn-model.git
  git clone https://github.com/lequieu/1cpn-model.git

5. Link 1CPN with LAMMPS

::

  cd $D_1CPN/1cpn-model/src/lammps
  make link

6. Build LAMMPS

::

  cd $D_1CPN/lammps-1cpn/src
  make serial


::

  git clone -b stable https://github.com/lammps/lammps.git lammps-1cpn
  cd lammps-1cpn
  git checkout -b stable_16Mar2018


4. Clone the 1CPN model.

:: 

  git clone https://lequieu@bitbucket.org/lequieu/1cpn-model.git
  git clone https://github.com/lequieu/1cpn-model.git

5. Link 1CPN with LAMMPS.

::

  cd $D_1CPN/1cpn-model/src/lammps
  make link

The Makefile assumes that the LAMMPS src code is located at ${D_1CPN}/lammps-1cpn. If you have LAMMPS located at another location, you can specify by redefining the `LAMMPS_SRC` variable.

::

  make link LAMMPS_SRC=<path to LAMMPS src>


6. Build LAMMPS with the `ASPHERE` and `MOLECULE` package.

::

  cd $D_1CPN/lammps-1cpn/src
  make yes-ASPHERE
  make yes-MOLECULE
  make serial



Running your first simulation
----------------------------

Now that LAMMPS has been linked and compuled with the 1CPN-model, we're ready to setup and analyze our first 1CPN simulation. First lets make a new directory called `example` where the input and output files from this simulation will be stored. After the directory is made, navigate into it.

::

    cd ${D_1CPN}/example
    mkdir -p ${D_1CPN}/example


Generating 1CPN Input Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. The next step is to generate the initial configuration of the 1CPN-model. This is achieved using the `${D_1CPN}/init/init_1cpn.py` script. To start out, we'll generate a section of chromatin consisiting of 20 nucleosomes, each with a nucleosome repeat length (NRL) of 187 base pairs. From the `${D_1CPN}/example` directory, issue the following command:

::

    ${D_1CPN}/1cpn-model/init/init_1cpn.py -n 20 -nrl 187

The full list of arguements that can be passed to `init_1cpn.py` are described LINK TO ME.

Running a Simulation
^^^^^^^^^^^^^^^^^^^^^^

2. Next we need to copy the LAMMPS input files necessary for setting up and running a simulation. These files are located at `${D_1CPN}/inputs`. Copy them to our `example` directory.

::
    cp ${D_1CPN}/1cpn-model/inputs/in.* .

3. Now we're ready to run a simple simulation. This is as simple as executing

::

    ${D_1CPN}/lammps-1cpn/src/lmp_serial -i in.1cpn

Users unfamilair with LAMMPS are referred to the LAMMPS docmentation LiNK ME, for descriptions on LAMMPS output, and how to run different simulations.


Vizualizing the Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :ref:`label-viz`.

Analyzing the Simulation
---------------------------

