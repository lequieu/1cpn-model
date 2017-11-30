#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"

#define TWID 500

struct DistInfo{
    int sitea;
    int siteb;
    double dist;
};

int main(int argc, char**argv){
    
    bool badinputflag = false;     
    if (argc < 5) { //must be at least 4 args
        badinputflag = true;
    }
    else if (((argc-3)%2)!=0){ //must have multiple of 2 args (excluding dump file)
        badinputflag = true;
    }

    if (badinputflag){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <output filename> <site 1a> <site 1b> <site 2a> <site 2b> etc..."<<std::endl;
        std::cout<<"Note indexing starts at 1" << std::endl;
        exit(1);
    }
    
    // Load dump
    int timestep, natom;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
       
    //Load the trajectory into the parser
    parser.load_dump(argv[1]); 
    //Get the number of frames
    timestep = parser.get_numFrames();
    std::cout<<timestep<<std::endl;
    //Get the vector for the types of atoms
    atom_types = parser.get_type(); 
    //Get the number of atoms
    natom = parser.get_numAtoms();
    box_dim = parser.get_boxDim();

    //setup file to write to
    std::ofstream ofile;
    ofile.open(argv[2]);

    // figure out what I need to calculate
    int ndist = (argc-3)/2;
    std::vector<DistInfo> distinfo(ndist);
   
    int pos = 3;
    for (size_t i=0; i< ndist; i++){
        distinfo[i].sitea = std::atoi(argv[pos]);
        distinfo[i].siteb = std::atoi(argv[pos+1]);
        
        //check inputs
        if ((distinfo[i].sitea > natom)||(distinfo[i].sitea <= 0)){
            std::cout<<"Invalid cluster" <<  i << "sitea: " << distinfo[i].sitea << std::endl;
            exit(1);
        }
        if ((distinfo[i].siteb > natom)||(distinfo[i].siteb <= 0)){
            std::cout<<"Invalid cluster" <<  i << "siteb: " << distinfo[i].siteb << std::endl;
            exit(1);
        }
        pos +=2;
    }

    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<timestep; i++) {

        //Write new dump file
        int type,prevtype;

        //The actual functions from the parser
        parser.next_frame();
        //quats = parser.get_quat();

        //Make sure to move to the next frame
        //Since everything after this is processing the data we can put the next frame option here

        ofile << i << "\t";
    
        double dx,dy,dz,dist,angle;
        for (size_t j = 0; j<ndist;j++){
                        
            std::vector<double> vectA(3);
            std::vector<double> vectB(3);
            vectA = parser.coords_[distinfo[j].sitea-1]; 
            vectB = parser.coords_[distinfo[j].siteb-1]; 
            
            dx = vectB[0] - vectA[0];
            dy = vectB[1] - vectA[1];
            dz = vectB[2] - vectA[2];
            dist = sqrt(dx*dx + dy*dy + dz*dz);
            distinfo[j].dist = dist;

            ofile << distinfo[j].dist << "\t";
        } 
        ofile << std::endl;

        if (firstframe) firstframe = false;
    }  
    ofile.close();
}
