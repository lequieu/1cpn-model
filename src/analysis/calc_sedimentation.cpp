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

    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
    natoms = parser.get_numAtoms();
    atom_types = parser.get_type(); 
    ntimestep = parser.get_numFrames();

    int nnucl; 
    std::vector<int> nucl_ids;
    double Rnucl = 50; //this is just a guess
    double S1 = 11.1; //Arya2006 eq 32
    double S20w;
       
    bool firstframe = true;
    //Loop through the dump file using the parser
    //for(size_t i=0; i<timestep; i++) {
    for(size_t i=0; i<ntimestep; i++) {
       
        parser.next_frame();

        vects_f = parser.get_vect('f');
        //vects_v = parser.get_vect(quats,'v');
        vects_u = parser.get_vect('u');
        t = parser.get_current_timestep();

        if (firstframe){
            nnucl = 0;
            for (size_t j=0;j<natoms;j++){
              if (atom_types[j] == 1){ //is nucleosome
                nucl_ids.push_back(j);
                nnucl++;
              }
            }

        }
  
        //compute sedimentation coeff 
        double sum=0;
        double dx,dy,dz,dr;
        int nuclj, nuclk;
        for(size_t j=0;j<nnucl-1;j++){
          nuclj = nucl_ids[j];
          for(size_t k=j+1;k<nnucl;k++){
            nuclk = nucl_ids[k];
            dx = parser.coords_[nuclk][0] - parser.coords_[nuclj][0];
            dy = parser.coords_[nuclk][1] - parser.coords_[nuclj][1];
            dz = parser.coords_[nuclk][2] - parser.coords_[nuclj][2];
            dr = sqrt(dx*dx+dy*dy+dz*dz);
            sum += 1.0/dr;
          }
        }
        S20w = S1 * (1 + 2.* Rnucl / nnucl * sum);
        ofile << t<< "\t" <<S20w << "\t" << sum << std::endl;
    
        if (firstframe) firstframe = false;
    }  

    ofile.close();
   
}
