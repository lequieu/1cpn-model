#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"

//Simple scalar product function for ease of use
double scalar_product(std::vector<double> a, std::vector<double> b)
{
    if( a.size() != b.size() ) // error check
    {
        puts( "Error a's size not equal to b's size" ) ;
        return -1 ;  // not defined
    }

    // compute
    double product = 0;
    for (int i = 0; i <= a.size()-1; i++)
       product += (a[i])*(b[i]); // += means add to product
    return product;
}

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
    ofile << "# <timestep> <rouse modes> ....." << std::endl;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<float> box_dim;
    natoms = parser.get_numAtoms();
    ntimestep = parser.get_numFrames();

    int nnucl; 
    std::vector<int> nucl_ids;
    std::vector<double> rmodes;
      
    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<ntimestep; i++) {
      
        //Change frame and calculate time step
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
            //resize the rouse modes array
            rmodes.resize(nucl_ids.size()-1);
        }
 
        //Initialize some starting variables
        double invN = 1.0/float(nnucl); //inverse of number of nucleosomes
        double prefactor = sqrt(2.0*invN); //prefactor for rouse mode calculation

        //compute each rouse mode 
        for(size_t p=1;p<nnucl;p++){ //iterate through possible rouse modes
          std::vector<double> rmodeVect = {0,0,0};
          for(size_t j=0;j<nnucl;j++){ //iterate through all nucleosomes
            for(size_t k=0;k<3;k++) { //iterate through 3 dimensions
              rmodeVect[k] += prefactor*parser.coords_[nucl_ids[j]][k]*cos(float(p)*M_PI*invN*(float(j)+0.5)); //Normally its j-0.5, but for this we start at 0 instead of 1
            }
          }
          rmodes[p-1] = sqrt(scalar_product(rmodeVect,rmodeVect));
        }

        //write out the rouse modes to the output file
        ofile << t << "\t\t";
        for(auto& rmode : rmodes)
            ofile << rmode << "\t\t";
        ofile << std::endl;
    
        if (firstframe) firstframe = false;
    }  
    ofile.close();   
}
