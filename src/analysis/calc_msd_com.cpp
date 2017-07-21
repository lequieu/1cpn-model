#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"

#define TWID 500

struct AngleInfo{
    std::string vectname;
    int sitea;
    int siteb;
    double angle;
};

int main(int argc, char**argv){
    if (argc<3){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <com output file name>"<<std::endl;
        exit(1);
    }

    int itraj, ntraj;
    ntraj = std::atoi(argv[1]);

        //MASSES are hardcoded!
    std::vector<double> masses(3+1);
    masses[1] = 196666.0000;
    masses[2] = 1950.000000;  
    masses[3] = 19500.00000;  

    std::string filename;
     


    int nframes, natom, dumpfreq;
    int dumpfreqprev;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
      
    //get the min number of frames        
    int minframes = 1e9;
    for (int itraj = 0; itraj<ntraj; itraj++){
      filename = "./" + std::to_string(itraj) + "/traj.dump";
      parser.load_dump(filename.c_str()); 
      nframes = parser.get_numFrames();
      if (nframes < minframes){
         minframes = nframes;
      }
      parser.close();
    }

    std::vector<double> sumsqdisp(minframes,0);
    std::vector<double> meansqdisp(minframes,0);


    for (int itraj = 0; itraj<ntraj; itraj++){
      filename = "./" + std::to_string(itraj) + "/traj.dump";
      //Load the trajectory into the parser
      parser.load_dump(filename.c_str()); 
      //Get the vector for the types of atoms
      atom_types = parser.get_type(); 
      dumpfreq = parser.get_dumpfreq(); 
      if (itraj==0){
        dumpfreqprev = dumpfreq;
      }
      else{
        if (dumpfreq != dumpfreqprev){
            std::cerr << "Error! dumpfreq is inconsistent between trajectories!" << std::endl;
            exit(1);
        }
      }
      //Get the number of atoms
      natom = parser.get_numAtoms();
      box_dim = parser.get_boxDim();

      double com0[3];

      bool firstframe = true;
      //Loop through the dump file using the parser
      for(size_t i=0; i<minframes; i++) {

          //Write new dump file
          int type,prevtype;

          //The actual functions from the parser
          atoms = parser.get_coord();
          quats = parser.get_quat();
          vects_f = parser.get_vect(quats,'f');
          vects_u = parser.get_vect(quats,'u');
          vects_v = parser.get_vect(quats,'u');
          
          //calculate COM
          double com[3], totalmass;
          com[0] =0.0;
          com[1] =0.0;
          com[2] =0.0;
          totalmass = 0;
          for (size_t iatom=0; iatom < natom; iatom++){
            type = atom_types[iatom];
            com[0] += atoms[iatom][0]*masses[type];
            com[1] += atoms[iatom][1]*masses[type];
            com[2] += atoms[iatom][2]*masses[type];
            totalmass += masses[type];
          }
          com[0] /= totalmass;
          com[1] /= totalmass;
          com[2] /= totalmass;

          if (firstframe){
            com0[0] = com[0];
            com0[1] = com[1];
            com0[2] = com[2];
          }
          double sqdisp,dx,dy,dz;
          dx = com[0]- com0[0];
          dy = com[1]- com0[1];
          dz = com[2]- com0[2];
          sumsqdisp[i] += (dx*dx + dy*dy + dz*dz);

          parser.next_frame();
          if (firstframe) firstframe = false;
      }  

      parser.close();
    }

    for (size_t i=0;i<minframes+1;i++){
        meansqdisp[i] = sumsqdisp[i]/(minframes+1);
    }

    //setup file to write to
    std::ofstream ofile;
    std::string outfile;
    outfile = argv[2];
    ofile.open(outfile);
    for (size_t i=0;i<minframes+1;i++){
        ofile << i*dumpfreq << " " << meansqdisp[i] << std::endl;
    }
    ofile.close();

}
