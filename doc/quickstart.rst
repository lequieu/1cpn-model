
Quick Start
==================

This a quick start guide to getting going quickly with the 1CPN model.


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


3. Download a fresh copy of LAMMPS and checkout the `stable_31Mar2017` tag. Even if you already have a copy of LAMMPS, its recommended that you download a new copy specifically for 1CPN.

::

  git clone -b stable https://github.com/lammps/lammps.git lammps-1cpn
  git checkout -b stable_31Mar2017


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


Generating 1CPN Input Files
---------------------------


Running a Simple Simulation
----------------------------


Vizualizing the Simulation
---------------------------
