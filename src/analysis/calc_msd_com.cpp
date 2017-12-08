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


    std::string filename;
     
    int nframes, natom, dumpfreq;
    int dumpfreqprev;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    std::vector<float> box_dim;
      
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

      std::vector<double> com0(3);
      std::vector<double> com(3,0);

      bool firstframe = true;
      //Loop through the dump file using the parser
      for(size_t i=0; i<minframes; i++) {

          //Write new dump file
          int type,prevtype;

          //The actual functions from the parser
          parser.next_frame();

          //calculate COM
          com = parser.get_com();

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
