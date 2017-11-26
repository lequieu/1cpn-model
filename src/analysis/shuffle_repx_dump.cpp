#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"


//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    elems.clear(); //empty elems
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) //ignore empty entries
          elems.push_back(item);
    }
    //return elems;
}


int main(int argc, char**argv){

    if (argc != 4){
        std::cout<<"Usage: "<<argv[0]<<" <temper log file> <temper dump file prefix (e.g.) traj.temper> <new dump prefix> "<<std::endl;
        exit(1);
    }
    
    // Load dump
    int timestep, natom;

    std::string logfilename;
    std::string dumpfilenameprefix;
    std::string newdumpfilenameprefix;
    logfilename = argv[1];
    dumpfilenameprefix = argv[2];
    newdumpfilenameprefix = argv[3];

    //parse log file
    std::ifstream logfile(logfilename);
    if(!logfile.is_open()) {std::cerr<<"Error! "<<logfilename<<" does not exist!"<<std::endl; exit(1);}
   
    //get number of replicas
    std::string line;
    std::vector<std::string> l; 
    int nreplicas,flag=0;
    do { 
      flag = std::getline(logfile,line).eof();
      //std::cout << line << std::endl;
    }
    while(!(line.find("Step T0") != std::string::npos) && !flag);
    if (flag){
        std::cerr << "Error! logfile doesn't look like a temper logfile" <<std::endl;
        exit(1);
    }
    split(line,' ',l);
    nreplicas = l.size() - 1;
    std::cout << "nreplicas: " << nreplicas << std::endl;
    

    //Set up all vectors needed for the trajectory parser class
    std::vector<TrajectoryIterator> parsers(nreplicas);
    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
    
    //max frames
    bool uneven_length_flag = false;
    int ntimesteps = 0;
    int maxtimesteps = 0;

    //Load the trajectories into the parsers
    parsers.resize(nreplicas);
    std::cout << "nframes: ";
    for (int i=0;i < nreplicas; i++ ){
        std::string dumpfilename;
        dumpfilename = dumpfilenameprefix + "." +  std::to_string(i) + ".dump";
        parsers[i].load_dump(dumpfilename.c_str());

        //check that all files are same size
        int nframes, nframes_prev;
        nframes = parsers[i].get_numFrames();
        std::cout << nframes << " ";
        if ((i) && (nframes != nframes_prev)){
            uneven_length_flag = true;
            //calc min of prev and current
            maxtimesteps = ((nframes_prev < nframes)? nframes_prev:nframes); 
        }
        nframes_prev = nframes; 
        
        //clear newdumpfiles (in case they exist)
        std::string newdumpfilename;
        newdumpfilename = newdumpfilenameprefix + "." + std::to_string(i) + ".dump";
        std::ofstream file(newdumpfilename, std::ofstream::out); 
    }
    std::cout << std::endl;

    if (uneven_length_flag){
      std::cerr << "Warning! dump files do not have equal number of frames. Only reading first " << maxtimesteps << " frames" << std::endl;
    }

    
    //parse log file
    long long log_timestep, dump_timestep;
    while (!std::getline(logfile,line).eof()){
      split(line,' ',l);
      //each iteration of while loop is different line of log.lammps
      
      //make sure line is correct length, and l[0] isn't "Step" (the header of a temper simulation)
      if (l.size() == (nreplicas+1) && (l[0].compare("Step"))){
        log_timestep = std::stoll(l[0]);
        
        //if dump file's are uneven, only use the first n frames
        if ((uneven_length_flag) &&  (ntimesteps > maxtimesteps)) continue;
        ntimesteps++;
        
        //for each temperature (for each timestep)
        for(int itemp = 0; itemp<nreplicas; itemp++){
           int idump; // the dumpfile that contains that temperature
           for(idump = 0; idump < nreplicas; idump++){
                if (!l[idump+1].compare(std::to_string(itemp))) break;
           } 
           if (idump == nreplicas+1){
             std::cerr << "Error! itemp " << itemp << " not found in line " << line << std::endl;
             exit(1);
           }

           // find  the spot in that dumpfile where log_timestep == dump_timestep;
           int flag;
           do{
             dump_timestep = parsers[idump].get_current_timestep();
             flag = parsers[idump].next_frame(); //return 1 at eof
           }while(!flag && (log_timestep != dump_timestep));

           if (flag){
             if (log_timestep != dump_timestep){
              //on last frame, flag can == 0, but if timesteps are equal its okay
              std::cerr << "Error! log_timestep " << log_timestep << " is not found in dumpfile " << idump << std::endl;
              exit(1);
             }
             else{ //flip eof bit
                parsers[idump].clear_file_errors();
             }
           }
           else{ 
             //dump_timestep corresponds to prev frame, but only if not eof
             parsers[idump].previous_frame(); 
           }
            
           // write that frame into the new dump file
           std::string newdumpfilename;
           newdumpfilename = newdumpfilenameprefix + "." + std::to_string(itemp) + ".dump";
           parsers[idump].append_current_frame_to_file(newdumpfilename); 
        }
      }
    }
    logfile.close(); 
    
}
