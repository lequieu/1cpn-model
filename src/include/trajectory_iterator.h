#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>


// A few useful functions copied from LAMMPS
typedef double vector[3];                // 0:x  1:y  2:z
typedef double quaternion[4];                // quaternion
typedef double form[6];                        // 0:xx 1:yy 2:zz 3:zy 4:zx 5:yx
inline void quat_vec_rot(vector &dest, vector &src, quaternion &q) {
  quaternion aa={q[0]*q[0], q[1]*q[1], q[2]*q[2], q[3]*q[3]};
  form ab={q[0]*q[1], q[0]*q[2], q[0]*q[3], q[1]*q[2], q[1]*q[3], q[2]*q[3]};
  dest[0] = (aa[0]+aa[1]-aa[2]-aa[3])*src[0]+
            ((ab[3]-ab[2])*src[1]+(ab[1]+ab[4])*src[2])*2.0;
  dest[1] = (aa[0]-aa[1]+aa[2]-aa[3])*src[1]+
            ((ab[2]+ab[3])*src[0]+(ab[5]-ab[0])*src[2])*2.0;
  dest[2] = (aa[0]-aa[1]-aa[2]+aa[3])*src[2]+
            ((ab[4]-ab[1])*src[0]+(ab[0]+ab[5])*src[1])*2.0;
}

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
        std::vector<float> halfBox_;
        std::vector<int> types_;
        std::streampos pos_;
        std::streampos posPrev_;
		bool crash_ = false;   //Checks for a crash
        bool firstFrame_ = true;
        long long timestep_; 
        long long timestepPrev_; 

		bool get_crash(void);
		bool check_crash(std::string);
        double check_pbc(double,int);
        void get_type(void); 
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
        void get_info(void);
        std::vector<std::vector<double>> coords_;
        std::vector<std::vector<double>> quats_;
        std::vector<std::vector<double>> get_vect(char);
        std::vector<int> get_types(void);
        std::vector<float> get_boxDim(void);
        std::vector<double> get_com(void);
        std::vector<double> get_distVect(int,int);
        double get_dist(int,int);
        double get_angleSites(int,int,int);
        long long get_current_timestep(void); 
        int get_current_natoms(void); 
        int get_numAtoms(void);
        int get_numFrames(void);
        int get_dumpfreq(void);
        int next_frame(void); 	
		bool isFloat(std::string);
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

//Reads string before it is cast as a float to make sure it is a float
bool TrajectoryIterator::isFloat( std::string myString ) {
	std::istringstream iss(myString);
	float f;
	iss >> std::noskipws >> f; //noskipws considers leading whitespace invalid
	return iss.eof() && !iss.fail();
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
    halfBox_.resize(3);

    //Only store the number of atoms at the beginning.
    //This makes all useable loops of the form get_fxns->next_frame 
    next_frame();

    timestepPrev_ = timestep_;
    numAtomsPrev_ = numAtoms_;

    //Initialize and size the quat vectors to num atoms
    quats_.resize(numAtoms_);
    coords_.resize(numAtoms_);
    for(size_t i=0; i<numAtoms_; i++) {
        quats_[i].resize(4);
        coords_[i].resize(3);
    }

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
            //Populate the internal vectors here if the system hasn't crashed
            if(firstFrame_) {
                get_type();
                firstFrame_ = false;
                return 0;
            }
            get_info();
            if(!crash_) { return 0; }
            else {crash_ = false;}
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
                halfBox_[i] = 0.5*(boxDim_[2*i+1]-boxDim_[2*i]);
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

//This checks to see if the trajectory has become in any way corrupted
bool TrajectoryIterator::check_crash(std::string myString) {
    std::vector<std::string> l;
	split(myString,' ',l); 
	for(size_t j=1; j<l.size(); j++) {
		if(!isFloat(l[j])) {
			std::cerr << "Warning! Read error occurred at timestep: " <<timestep_<<std::endl;
			crash_ = true;
			return true;
		}
	}
	crash_ = false;
	return false;
}

//Get the information of the atom vectors and quaternions
//FIXME I dont like the double parsing of the file when coords and quat are gotten, when next frame is called, I think internal vectors should be populated for type, coord, quat
//FIXME then the get_coord will only copy that data over?, or send a pointer to the TrajectoryIterator object?
void TrajectoryIterator::get_info() {
    //std::vector<std::vector<double>> atom_quat;

    if (numAtoms_ <= 0){
        std::cout << "Error! Trying to get_quat() but numAtoms <= 0. Could the traj file be empty?" << std::endl;
        exit(1);
    }

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
		if (!check_crash(line)) {
            std::stringstream sin(line);
            sin >> index >> type;
            index -= 1;
            for(size_t j=0; j<3; j++) {
                sin >> x[j];
				coords_[index][j] = x[j]; 
            }
            for(size_t j=0; j<4; j++) {
                sin >> q[j];
                quats_[index][j] = q[j];
            }
        }
        else {
            return;
        }
    }
    return;
};

//Convert the quaternions to f,v,u vectors
std::vector<std::vector<double>> TrajectoryIterator::get_vect(char type) {
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
    for(size_t i=0; i<quats_.size(); i++) {
        norm = sqrt(quats_[i][0]*quats_[i][0] +  quats_[i][1]*quats_[i][1] + quats_[i][2]*quats_[i][2] + quats_[i][3]*quats_[i][3]);
        if (fabs(norm - 1.0) > 1e-3){
            std::cerr<<"Norm of quat is too big ("<<norm<<")! You've screwed something up with quat math!"<<std::endl;
            exit(1);
        }
        rquat[0] = quats_[i][0];
        rquat[1] = quats_[i][1];
        rquat[2] = quats_[i][2];
        rquat[3] = quats_[i][3];

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
void TrajectoryIterator::get_type() {
    if (numAtoms_ <= 0){
        std::cout << "Error! Trying to get_type() but numAtoms <= 0. Could the traj file be empty?" << std::endl;
        exit(1);
    }
    //Initialize and size the vector
    types_.resize(numAtoms_);

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
        types_[index] = type;
    }
    return;
};

std::vector<double> TrajectoryIterator::get_com() {
    std::vector<double> com(3);
    std::vector<double> masses(3);
    
    //The masses are hard-coded here for now
    masses[0] = 196666.0000;
    masses[1] = 1950.000000;  
    masses[2] = 19500.00000;  

    int type;
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;
    double totalmass = 0;
    for (size_t iatom=0; iatom < numAtoms_; iatom++){
        type = types_[iatom]-1;
        com[0] += coords_[iatom][0]*masses[type];
        com[1] += coords_[iatom][1]*masses[type];
        com[2] += coords_[iatom][2]*masses[type];
        totalmass += masses[type];
    }
    com[0] /= totalmass;
    com[1] /= totalmass;
    com[2] /= totalmass;

    return com;
};

double TrajectoryIterator::check_pbc(double dist, int dim) {
    if (dist > halfBox_[dim]) {return dist - halfBox_[dim];}
    else if (dist < -halfBox_[dim]) {return dist + halfBox_[dim];}
    else {return dist;}
};

// Returns the Euclidean distance between siteA and siteB
// Note: the sites that are input are the same as that of the in.lammps file!
double TrajectoryIterator::get_dist(int siteA, int siteB) {
    double dist = 0.0, dx = 0.0;
    for (size_t i = 0; i < 3; i++) {
        dx = coords_[siteA-1][i]-coords_[siteB-1][i];
        dx = check_pbc(dx,i);
        dist += dx*dx;
    }
    return sqrt(dist);
};

// Returns the angle site ABC (B is the middle site)
// Note: the sites that are input are the same as that of the in.lammps file!
double TrajectoryIterator::get_angleSites(int siteA, int siteB, int siteC) {
    double angle = 0.0, magA = 0, magB = 0; 
    std::vector<double> vectA;
    std::vector<double> vectB;

    vectA = get_distVect(siteA,siteB);
    vectB = get_distVect(siteB,siteC);
    magA = get_dist(siteA,siteB);
    magB = get_dist(siteB,siteC);

    for (size_t i = 0; i < 3; i++) {
        vectA[i] /= magA;
        vectB[i] /= magB;
        angle += vectA[i]*vectB[i];
    }

    if (angle > 1) angle = 1;
    if (angle <-1) angle =-1;

    return acosf(angle)*180./M_PI;
};

// Returns the vector between siteA and siteB
// Note: the tail of the vector begins at siteB and the head is at siteA!
// Note: the sites that are input are the same as that of the in.lammps file!
std::vector<double> TrajectoryIterator::get_distVect(int siteA, int siteB) {
    std::vector<double> distVect(3);
    double dx = 0;
    for (size_t i = 0; i < 3; i++) {
        dx = coords_[siteA-1][i]-coords_[siteB-1][i];
        dx = check_pbc(dx,i);
        distVect[i] = dx;
    }
    return distVect;
};

int TrajectoryIterator::get_numAtoms() {
    return numAtoms_;
};

std::vector<float> TrajectoryIterator::get_boxDim() {
    return boxDim_;
};
void TrajectoryIterator::append_current_frame_to_file(std::string filename){
    //Set the point for the input file
    //dumpFile_.seekg(pos_);

    if(!crash_) {
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
        get_info();
        for (int i = 0; i < numAtoms_; i++){
            file << i+1 << " " << types[i] << " ";
            file << coords_[i][0] << " " << coords_[i][1] << " " << coords_[i][2] << " ";
            file << quats_[i][0] << " " << quats_[i][1] << " " << quats_[i][2] << " " << quats_[i][3];
            file << std::endl;
        }
        file.close();
    }
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
bool TrajectoryIterator::get_crash(){
	return crash_;
}
std::vector<int> TrajectoryIterator::get_types(){
    get_type();
    return types_;
}
