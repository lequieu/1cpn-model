#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"



int main(int argc, char**argv){

    if (argc != 3){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <output file> "<<std::endl;
        exit(1);
    }

   
    // Load dump
    long long ntimestep,t;
    int natoms;

    std::string dumpfilename;
    std::string outfilename;
    dumpfilename = argv[1];
    outfilename = argv[2];

    std::ofstream ofile;
    ofile.open(outfilename.c_str());
    ofile << "# <timestep> <S20,w>" << std::endl;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<float> box_dim;
    natoms = parser.get_numAtoms();
    ntimestep = parser.get_numFrames();

    int nnucl; 
    std::vector<int> nucl_ids;
    double Rnucl = 54.6; //Arya2006 
    double S1 = 11.1; //Arya2006 eq 32
    double S20w;
       
    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<ntimestep; i++) {
       
        //Get the frame timestep and move to next frame
        parser.next_frame();
        t = parser.get_current_timestep();

        if (firstframe){
            nnucl = 0;
            std::vector<int> types = parser.get_types();
            for (size_t j=0;j<natoms;j++){
              if (types[j] == 1){ //is nucleosome
                nucl_ids.push_back(j+1);
                nnucl++;
              }
            }
        }
  
        //unwrap the coordinates
        parser.unwrap_coords();

        //compute sedimentation coeff 
        double sum=0;
        double dx,dy,dz,dr;
        int nuclj, nuclk;
        for(size_t j=0;j<nnucl-1;j++){
          nuclj = nucl_ids[j];
          for(size_t k=j+1;k<nnucl;k++){
            nuclk = nucl_ids[k];
            dr = parser.get_dist(nuclk,nuclj);
            sum += 1.0/dr;
          }
        }
        S20w = S1 * (1 + 2.* Rnucl / nnucl * sum);
        ofile << t<< "\t" <<S20w << "\t" << sum << std::endl;
    
        if (firstframe) firstframe = false;
    }  

    ofile.close();
   
}
