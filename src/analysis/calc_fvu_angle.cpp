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
    /*
        Script that parses a dump file and computes the angle between the f,v or u vectors of some number of sites
        The script takes several arguments so the angle between many sites can be calculated in one go
    
    */
    
    bool badinputflag = false;     
    if (argc < 5) { //must be at least 4 args
        badinputflag = true;
    }
    else if (((argc-3)%3)!=0){ //must have multiple of 3 args (excluding dump file)
        badinputflag = true;
    }

    if (badinputflag){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <output filename> <f,v,u vector> <site 1a> <site 1b> <fvu> <site 2a> <site 2b> etc..."<<std::endl;
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
    int nangle = (argc-1)/3;
    std::vector<AngleInfo> angleinfo(nangle);
   
    int pos = 3;
    for (size_t i=0; i< nangle; i++){
        angleinfo[i].vectname = argv[pos];
        angleinfo[i].sitea = std::atoi(argv[pos+1]);
        angleinfo[i].siteb = std::atoi(argv[pos+2]);
        
        //check inputs
        if ((angleinfo[i].vectname.compare("f") != 0) &&
            (angleinfo[i].vectname.compare("v") != 0) &&
            (angleinfo[i].vectname.compare("u") != 0)){
            std::cout << "Invalid fvu character: "<<argv[pos] << std::endl;
            exit(1);
        }
        if ((angleinfo[i].sitea > natom)||(angleinfo[i].sitea <= 0)){
            std::cout<<"Invalid cluster" <<  i << "sitea: " << angleinfo[i].sitea << std::endl;
            exit(1);
        }
        if ((angleinfo[i].siteb > natom)||(angleinfo[i].siteb <= 0)){
            std::cout<<"Invalid cluster" <<  i << "siteb: " << angleinfo[i].siteb << std::endl;
            exit(1);
        }
    
        pos +=3;
    }

    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<timestep; i++) {

        //Write new dump file
        int type,prevtype;

        //The actual functions from the parser
        parser.next_frame();
        vects_f = parser.get_vect('f');
        vects_v = parser.get_vect('v');
        vects_u = parser.get_vect('u');
        //Make sure to move to the next frame
        //Since everything after this is processing the data we can put the next frame option here

        ofile << i << "\t";
    
        double dot,angle;
        for (size_t j = 0; j<nangle;j++){
            std::vector<double> vectA(3);
            std::vector<double> vectB(3);
            if (angleinfo[j].vectname.compare("f") == 0){
                vectA = vects_f[angleinfo[j].sitea-1]; 
                vectB = vects_f[angleinfo[j].siteb-1]; 
            }
            else if (angleinfo[j].vectname.compare("v") == 0){
                vectA = vects_v[angleinfo[j].sitea-1]; 
                vectB = vects_v[angleinfo[j].siteb-1]; 
            }
            else if (angleinfo[j].vectname.compare("u") == 0){
                vectA = vects_u[angleinfo[j].sitea-1]; 
                vectB = vects_u[angleinfo[j].siteb-1]; 
            }
            
            dot = vectA[0]*vectB[0] +  vectA[1]*vectB[1] +  vectA[2]*vectB[2];
            if (dot > 1) dot =1;
            if (dot <-1) dot =-1;
            angleinfo[j].angle = acos(dot)*180./M_PI;
            ofile << angleinfo[j].angle << "\t";
        } 
        ofile << std::endl;

        if (firstframe) firstframe = false;
    }  
    ofile.close();
}
