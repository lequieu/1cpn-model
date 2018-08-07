#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"



int main(int argc, char**argv){
    /*
        Script to analyze a two nucleosome configuration
    */
    if (argc != 3){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <timeseries output file>"<<std::endl;
        exit(1);
    }

    // Load dump
    long long ntimestep,t;
    int natoms;
    int nuclA, nuclB;
    std::ofstream ofile;

    std::string dumpfilename;
    std::string outfilename;
    dumpfilename = argv[1];
    outfilename = argv[2];

    ofile.open(outfilename.c_str()); 

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<float> box_dim;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
    natoms = parser.get_numAtoms();
    ntimestep = parser.get_numFrames();

    double halfbox[3];

    // define vectors
    std::vector<double> r(3),rhat(3),fA(3),fB(3);
    double thetaA,thetaB,phi;
    double rnorm,fAdotfB, rdotfA, rdotfB;

    // define histograms
    std::vector<long long> count_in_angle_region(4,0);
    double thetaA_saddle=90., phi_saddle = 90.; // saddle points to define regions
    std::vector<long long> count_in_r_region(2,0);
    double rnorm_saddle = 140; //only bin values into histogram less than this distance

    bool firstframe = true;
    for(size_t i=0; i<ntimestep; i++) {
        if (firstframe){
          int nnucl = 0;
          std::vector<int> types = parser.get_types();
          for (size_t j=0; j<natoms; j++){
            if (types[j] == 1){
              if (nnucl == 0) nuclA = j;
              else if (nnucl == 1) nuclB = j;
              else{
                std::cout << "Warning! More than two nucleosomes in dump file! Are you sure you're doing what you expect?" << std::endl;
              }
              nnucl ++;
            }
          }
          //std::cout << "Nucl ID are: " << nuclA << " " << nuclB << std::endl;

        }
        //The actual functions from the parser
        parser.next_frame();
        t = parser.get_current_timestep();
        vects_f = parser.get_vect('f');
        box_dim = parser.get_boxDim();

        
        fA = vects_f[nuclA];
        fB = vects_f[nuclB];
       
        r = parser.get_distVect(nuclB,nuclA);
        ////r[0] = parser.coords_[nuclB][0] - parser.coords_[nuclA][0];
        ////r[1] = parser.coords_[nuclB][1] - parser.coords_[nuclA][1];
        ////r[2] = parser.coords_[nuclB][2] - parser.coords_[nuclA][2];

        ////apply pbc
        //halfbox[0] = 0.5* (box_dim[1] - box_dim[0]);
        //halfbox[1] = 0.5* (box_dim[3] - box_dim[2]);
        //halfbox[2] = 0.5* (box_dim[5] - box_dim[4]);

        //if (r[0] >  halfbox[0]) r[0] -= halfbox[0]; 
        //if (r[0] <- halfbox[0]) r[0] += halfbox[0]; 
        //if (r[1] >  halfbox[1]) r[1] -= halfbox[1]; 
        //if (r[1] <- halfbox[1]) r[1] += halfbox[1]; 
        //if (r[2] >  halfbox[2]) r[2] -= halfbox[2]; 
        //if (r[2] <- halfbox[2]) r[2] += halfbox[2]; 

        rnorm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        rhat[0] = r[0]/rnorm;
        rhat[1] = r[1]/rnorm;
        rhat[2] = r[2]/rnorm;

        rdotfA = rhat[0]*fA[0] +  rhat[1]*fA[1] +  rhat[2]*fA[2];
        if (rdotfA > 1) rdotfA =1;
        if (rdotfA <-1) rdotfA =-1;
        thetaA= acos(rdotfA) * 180. / M_PI;

        rdotfB = rhat[0]*fB[0] +  rhat[1]*fB[1] +  rhat[2]*fB[2];
        if (rdotfB > 1) rdotfB =1;
        if (rdotfB <-1) rdotfB =-1;
        thetaB= acos(rdotfB) * 180. / M_PI;

        fAdotfB = fA[0]*fB[0] +  fA[1]*fB[1] +  fA[2]*fB[2];
        if (fAdotfB > 1) fAdotfB =1;
        if (fAdotfB <-1) fAdotfB =-1;
        phi = acos(fAdotfB) * 180. / M_PI;
        

        ofile << t << " " <<rnorm << " " << phi << " " << thetaA << " " << thetaB << std::endl;
        //threshold rnorm
        if (rnorm < rnorm_saddle){ //bound r region
          //ofile << t << " " <<rnorm << " " << fAdotfB << " " << rdotfA << " " << rdotfB << std::endl;
          count_in_r_region[0]++;

          //now determine region
          if ((phi <  phi_saddle) && (thetaA >= thetaA_saddle)) count_in_angle_region[0]++;
          if ((phi <  phi_saddle) && (thetaA <  thetaA_saddle)) count_in_angle_region[1]++;
          if ((phi >= phi_saddle) && (thetaA >= thetaA_saddle)) count_in_angle_region[2]++;
          if ((phi >= phi_saddle) && (thetaA <  thetaA_saddle)) count_in_angle_region[3]++;
        }
        else{ //unbound r region
          count_in_r_region[1]++;
        }

        if (firstframe) firstframe = false;
    }  
    ofile.close();

    long long sum=0; 
    for (size_t i=0;i<2;i++) sum += count_in_r_region[i];
    for (size_t i=0;i<2;i++){ 
      double prob = (double)count_in_r_region[i]/sum;
      std::cout << prob << " ";
    }
    sum=0; 
    std::vector<double> prob_in_region(4);
    for (size_t i=0;i<4;i++) sum += count_in_angle_region[i];
    for (size_t i=0;i<4;i++){ 
      prob_in_region[i] = (double)count_in_angle_region[i]/sum;
    }
    std::cout << prob_in_region[0] << " " <<  prob_in_region[1] << " " <<  prob_in_region[2]+prob_in_region[3] <<std::endl;
    //std::cout << std::endl;

}
