#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include "math_vector.h"

//For quaternion use only (probably try to do it better later on)
using namespace LAMMPS_NS;

//This is a standard c++ class that we plan to use for 1cpn analysis
//The TrajectoryIterator contains a load->get format which will be run from a separate c++ script
class TrajectoryIterator {

    private:
        std::ifstream dumpFile_;
        std::string filename_;
        int numAtoms_;                  //initial numAtoms
        int numAtomsPrev_;              //prev numAtoms
        std::vector<float> boxDim_;     //initial boxDim
        std::vector<float> boxDimPrev_; //prev boxDim
        std::streampos pos_;
        std::streampos posPrev_;
        long long timestep_; 
        long long timestepPrev_; 
        
        void split(const std::string&, char, std::vector<std::string>&);
    public:
        //Initialize variables of the class here
        TrajectoryIterator(){};
        TrajectoryIterator(int numAtoms):numAtoms_(numAtoms){};

        //Declare both the number of atoms here
        void clear_file_errors();
        void close();
        void reset();
        void load_dump(const char *);
        std::vector<std::vector<double>> get_coord(void);
        std::vector<std::vector<double>> get_quat(void);
        std::vector<std::vector<double>> get_vect(std::vector<std::vector<double>>, char);
        std::vector<int> get_type(void);
        std::vector<float> get_boxDim(void);
        long long get_current_timestep(void); 
        int get_current_natoms(void); 
        int get_numAtoms(void);
        int get_numFrames(void);
        int get_dumpfreq(void);
        int next_frame(void); 
        void previous_frame(void);        
        void append_current_frame_to_file(std::string);
};

//resets file to the first frame
void TrajectoryIterator::reset() {
    close();
    std::string myfilename = filename_;
    load_dump(myfilename.c_str());
}

void TrajectoryIterator::close() {
    dumpFile_.close();
}
//Function that loads the specified dump file
void TrajectoryIterator::load_dump(const char *fName) {
    //Open the file here and make sure it exists
    std::string str(fName);
    filename_ = str;
    dumpFile_.open(str, std::ifstream::in);
    if(!dumpFile_.is_open()) {std::cerr<<"Error! "<<fName<<" does not exist!"<<std::endl; exit(1);}
  
    //Size the boxdimensions to three dimensions
    boxDim_.resize(6);
    boxDimPrev_.resize(6);

    //Only store the number of atoms at the beginning.
    //This makes all useable loops of the form get_fxns->next_frame 
    next_frame();

    timestepPrev_ = timestep_;
    numAtomsPrev_ = numAtoms_;
    for (int i=0;i<6;i++){ 
      boxDimPrev_[i] = boxDim_[i];
    }
    
};

int TrajectoryIterator::get_numFrames() {
    if(!dumpFile_.is_open()) {std::cerr<<"Error! Function get_frames has been called before load_dump"<<std::endl; exit(1);}

    std::string line;
    int numFrames = 1; //start at one since loading the file parses first frame
    while (std::getline(dumpFile_,line)) {
        if(line.find("TIMESTEP") != std::string::npos) {
            numFrames++;
        }
    }
    dumpFile_.clear();
    dumpFile_.seekg(pos_, dumpFile_.beg);
    return numFrames;
}

int TrajectoryIterator::get_dumpfreq() {
    if(!dumpFile_.is_open()) {std::cerr<<"Error! Function get_dumpfreq has been called before load_dump"<<std::endl; exit(1);}

    std::string line;
    int numFrames = 0;
    int tprev,delta,dumpfreq;
    long long t;
    bool firstframe= true,secondframe=false;;
    while (std::getline(dumpFile_,line)) {
        if(line.find("TIMESTEP") != std::string::npos) {
            std::getline(dumpFile_,line);

            t = std::stoll(line);
            if (firstframe){
                firstframe=false;
                secondframe=true;
            }
            else{
                delta = t - tprev;
                if (!secondframe){
                  if (dumpfreq != delta){
                      std::cerr << "Error! Uneven dumpfrequencey!" <<std::endl;
                      exit(1);
                  }
                }
                dumpfreq = delta;
                secondframe=false;
            }
            tprev = t;
        }
    }
    dumpFile_.clear();
    dumpFile_.seekg(pos_, dumpFile_.beg);
    return dumpfreq;
}


//Change to next snapshot frame
//FIXME: check that timestep, natoms, etc are populated correclty
//FIXME: also check if correct return 1 if at end of file
int TrajectoryIterator::next_frame(void) {
    if(!dumpFile_.is_open()) {std::cerr<<"Error! trying to read next_frame() but file not open!"<<std::endl; exit(1);}

    //copy to previous
    timestepPrev_ = timestep_;
    numAtomsPrev_ = numAtoms_;
    for (int i=0;i<6;i++){ 
      boxDimPrev_[i] = boxDim_[i];
    }

    //Define a subline variable that maybe holds dump info
    std::string line;
    std::vector<std::string> l;
    //Find the next TIMESTEP of the function
    while (!std::getline(dumpFile_,line).eof()) {
        if (line.find("ITEM: ATOMS") != std::string::npos) {
            posPrev_ = pos_;
            pos_ = dumpFile_.tellg();

            return 0;
        }
        //populate values pertaining to ATOMS
        else if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            std::getline(dumpFile_,line);
            timestep_ = std::stoll(line);
        }
        else if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos) {
            std::getline(dumpFile_,line);
            numAtoms_ = std::stoi(line);
        }
        else if (line.find("ITEM: BOX BOUNDS") != std::string::npos) {
            for(int i=0; i<3; i++){
                std::getline(dumpFile_,line);
                split(line,' ',l); 
                boxDim_[2*i] = std::stof(l[0]);
                boxDim_[2*i+1] = std::stof(l[1]);
            }
        }
    }
    //if you get here, must be at end of file
    //pos_ = dumpFile_.tellg();
    return 1; 
};

//Change to previous snapshot frame
void TrajectoryIterator::previous_frame(void) {
    //Set the point for the input file
    pos_ = posPrev_;

    timestep_ = timestepPrev_;
    numAtoms_ = numAtomsPrev_;
    for (int i=0;i<6;i++){ 
      boxDim_[i] = boxDimPrev_[i];
    }
    
};

//Return the coordinates as a vector of vectors
//FIXME I dont like the double parsing of the file when coords and quat are gotten, when next frame is called, I think internal vectors should be populated for type, coord, quat
//FIXME then the get_coord will only copy that data over?, or send a pointer to the TrajectoryIterator object?
std::vector<std::vector<double>> TrajectoryIterator::get_coord(void) {
    std::vector<std::vector<double>> atom_pos;

    //Initialize and size the atom vectors
    atom_pos.resize(numAtoms_); 
    for(size_t i=0; i<numAtoms_; i++) {atom_pos[i].resize(3);}
    
    //Set the point for the input file
    dumpFile_.seekg(pos_);

    int index = 0;
    int type;
    double x[3];
    std::string line;
    for(size_t i=0; i<numAtoms_; i++) {
        std::getline(dumpFile_,line);
        if (!line.compare("")){
          std::cerr << "Error! Line is empty while reading coords, and it shouldn't be!" << std::endl;
          exit(1);
        }
        std::stringstream sin(line);
        sin >> index >> type;
        index -= 1;
        for(size_t j=0; j<3; j++) {
            sin >> x[j];
            atom_pos[index][j] = x[j];
        }
    }
    return atom_pos;
};

//Return the quaterions as a vector of vectors
std::vector<std::vector<double>> TrajectoryIterator::get_quat() {
    std::vector<std::vector<double>> atom_quat;

    //Initialize and size the quat vectors to num atoms
    atom_quat.resize(numAtoms_);
    for(size_t i=0; i<numAtoms_; i++) {atom_quat[i].resize(4);}

    //Set the point for the input file
    dumpFile_.seekg(pos_);

    //Loop through the snapshot to get the quaternion info
    std::string line;
    int index = 0, type;
    double x[3], q[4];
    for(size_t i=0; i<numAtoms_; i++) {
        std::getline(dumpFile_,line);
        if (!line.compare("")){
          std::cerr << "Error! Line is empty while reading quats, and it shouldn't be!" << std::endl;
          exit(1);
        }
        std::stringstream sin(line);
        sin >> index >> type;
        index -= 1;
        for(size_t j=0; j<3; j++) {sin >> x[j];}
        for(size_t j=0; j<4; j++) {
            sin >> q[j];
            atom_quat[index][j] = q[j];
        }
    }
    return atom_quat;
};

//Convert the quaternions to f,v,u vectors
std::vector<std::vector<double>> TrajectoryIterator::get_vect(std::vector<std::vector<double>> quats, char type) {
    //Initialize the orient vector
    std::vector<std::vector<double>> vects;
    vects.resize(numAtoms_);
    for (size_t i=0; i<numAtoms_; i++) {vects[i].resize(3);}

    //Check to see if the type is correct
    double fvu0[3];
    double fvu[3];
    if(type == 'f') { fvu0[0] = 1; fvu0[1] = 0; fvu0[2] = 0; }
    else if(type == 'v') { fvu0[0] = 0; fvu0[1] = 1; fvu0[2] = 0; }
    else if(type == 'u') { fvu0[0] = 0; fvu0[1] = 0; fvu0[2] = 1; }
    else { std::cout<<"Warning: "<<type<<" does not name a vector type"<<std::endl; }

    //Construct the LAMMPS formalism for quaternion calculations
    quaternion rquat;
    double norm;
    for(size_t i=0; i<quats.size(); i++) {
        norm = sqrt(quats[i][0]*quats[i][0] +  quats[i][1]*quats[i][1] + quats[i][2]*quats[i][2] + quats[i][3]*quats[i][3]);
        if (fabs(norm - 1.0) > 1e-3){
            std::cerr<<"Norm of quat is too big ("<<norm<<")! You've screwed something up with quat math!"<<std::endl;
            exit(1);
        }
        rquat[0] = quats[i][0];
        rquat[1] = quats[i][1];
        rquat[2] = quats[i][2];
        rquat[3] = quats[i][3];

        //Using the lammps equations (for now!) calculate the vector of choice
        quat_vec_rot(fvu,fvu0,rquat);

        //Repack into the vector structure
        vects[i][0] = fvu[0];
        vects[i][1] = fvu[1];
        vects[i][2] = fvu[2];
    }
    return vects;
};



//Get the type of each of the atoms
//This should only really be called once
std::vector<int> TrajectoryIterator::get_type() {
    std::vector<int> atom_type;

    //Initialize and size the vector
    atom_type.resize(numAtoms_);

    //Set the point for the input file
    dumpFile_.seekg(pos_);

    //Loop and only find the typing of the atom
    int index, type;
    std::string line;
    for(size_t i=0; i<numAtoms_; i++) {
        std::getline(dumpFile_,line);
        std::stringstream sin(line);
        sin >> index >> type;
        index -= 1;
        atom_type[index] = type;
    }
    return atom_type;
};

int TrajectoryIterator::get_numAtoms() {
    return numAtoms_;
};

std::vector<float> TrajectoryIterator::get_boxDim() {
    return boxDim_;
}
void TrajectoryIterator::append_current_frame_to_file(std::string filename){
    //Set the point for the input file
    //dumpFile_.seekg(pos_);

    std::ofstream file(filename, std::ofstream::app); //FIXME what flags?
    file << "ITEM: TIMESTEP" << std::endl;
    file << timestep_ << std::endl;
    file << "ITEM: NUMBER OF ATOMS" << std::endl;
    file << numAtoms_ << std::endl;
    file << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
    file << boxDim_[0] << " " << boxDim_[1] << std::endl;
    file << boxDim_[2] << " " << boxDim_[3] << std::endl;
    file << boxDim_[4] << " " << boxDim_[5] << std::endl;
    file << "ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4]" << std::endl;
    std::vector<int> types;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    atoms = get_coord();
    quats = get_quat();
    types = get_type();
    for (int i = 0; i < numAtoms_; i++){
        file << i+1 << " " << types[i] << " ";
        file << atoms[i][0] << " " << atoms[i][1] << " " << atoms[i][2] << " ";
        file << quats[i][0] << " " << quats[i][1] << " " << quats[i][2] << " " << quats[i][3];
        file << std::endl;
    }
    file.close();

}

/*
    Private method to split string by given delimiter, returns into elems
*/
void TrajectoryIterator::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    elems.clear(); //empty elems
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) //ignore empty entries
          elems.push_back(item);
    }
    //return elems;
}
long long TrajectoryIterator::get_current_timestep(){
    return timestep_;
}   
int TrajectoryIterator::get_current_natoms(){
    return numAtoms_;
}   
void TrajectoryIterator::clear_file_errors(){
    dumpFile_.clear();
} 
