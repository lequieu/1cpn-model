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
    ofile << "# <timestep> <total nucl average angle>  <theta(s)> ....." << std::endl;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<float> box_dim;
    natoms = parser.get_numAtoms();
    ntimestep = parser.get_numFrames();

    int nnucl; 
    std::vector<int> nucl_ids;
    std::vector<double> angles;
      
    bool firstframe = true;
    //Loop through the dump file using the parser
    //for(size_t i=0; i<timestep; i++) {
    for(size_t i=0; i<ntimestep; i++) {
       
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
            //resize the angles array
            angles.resize(nucl_ids.size()-1);
        }
  
        //compute sedimentation coeff 
        double sum=0, angle = 0;
        int nuclj, nuclk, nucll;
        for(size_t j=1;j<nnucl-1;j++){
          nuclj = nucl_ids[j-1];
          nuclk = nucl_ids[j];
          nucll = nucl_ids[j+1];
          angle = parser.get_angleSites(nuclj,nuclk,nucll);
          angles[j] = angle;
          sum += angle;
        }

        //calculate the average angle of nucleosomes
        sum *= 1.0/(nnucl-2.0);

        //write out the angles to the output file
        ofile << t << "\t" << sum << "\t";
        for(auto& angle : angles)
            ofile << angle << "\t";
        ofile << std::endl;
    
        if (firstframe) firstframe = false;
    }  

    ofile.close();
   
}
