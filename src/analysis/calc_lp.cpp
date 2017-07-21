#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"
#include "math_vector.h"


std::vector<double> rotate(std::vector<double> in, std::vector<double> u, const double theta){
    LAMMPS_NS::vector myin, myout;
    std::vector<double> out(3);
    LAMMPS_NS::quaternion q;

    for (size_t i=0;i<3;i++){
      myin[i] = in[i];
    }
    double norm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    

    q[0] = cos(theta/2.0);
    q[1] = u[0]/norm  * sin(theta/2.0);
    q[2] = u[1]/norm  * sin(theta/2.0);
    q[3] = u[2]/norm  * sin(theta/2.0);

    LAMMPS_NS::quat_vec_rot(myout, myin, q);
    out[0] = myout[0]; 
    out[1] = myout[1]; 
    out[2] = myout[2]; 
    return out;

}

int main(int argc, char**argv){

    if (argc != 6){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <lp output> <start frame> <end frame> <frame usage freq>"<<std::endl;
        exit(1);
    }

    ////check rotate
    //std::vector<double> a(3);
    //std::vector<double> b(3); 
    //std::vector<double> c(3); 
    //a[0] = 1; a[1] = 0; a[2] = 0;
    //b[0] = 0; b[1] = 1; b[2] = 0;
    //double theta = 90. * M_PI / 180.;
    //c = rotate(a,b,theta);
    //std::cout << c[0] << " " << c[1] << " " <<c[2] << std::endl;
    //exit(1);
    
    // Load dump
    long long timestep;
    int natoms;

    std::string dumpfilename;
    std::string outfilename;
    dumpfilename = argv[1];
    outfilename = argv[2];
    long long startframe, endframe, freq;
    startframe = std::atoll(argv[3]);
    endframe = std::atoll(argv[4]);
    freq = std::atoll(argv[5]);


    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
    natoms = parser.get_numAtoms();
    atom_types = parser.get_type(); 
    timestep = parser.get_numFrames();

    if ((endframe == -1) || (endframe > timestep)){
        endframe = timestep;
    }
    if (freq < 1) freq = 1;
    std::cout << "Reading frames (start: "<<startframe<<" end: "<<endframe<< " freq: "<<freq<<")"<<std::endl;
    
    int nbins, nignore;
    nignore = 10; //beads from both ends to ignore, possible input
    nbins = natoms - 2*nignore-1; //could change nbins to user defined input if desired

    // for bend angle Lp
    std::vector<double> lpdotsum;
    std::vector<long long> lpdotcount;
    lpdotsum.resize(nbins,0);
    lpdotcount.resize(nbins,0);

    // for end-to-end dist Lp
    std::vector<double> lpr2sum;
    std::vector<long long> lpr2count;
    lpr2sum.resize(nbins,0);
    lpr2count.resize(nbins,0);


    // for twist angle Lt
    std::vector<double> ltdotsum;
    std::vector<long long> ltdotcount;
    ltdotsum.resize(nbins,0);
    ltdotcount.resize(nbins,0);

    double l0_sum=0, l0_avg; //9.9
    long long l0_count=0;

    double w0 = 108. * M_PI / 180.;

    //make sure all atoms are dna type
    for (size_t j = 0; j<natoms;j++){
        if (atom_types[j] != 1){
            std::cerr << "Error! All atom_types must be 1 to compute lp" <<std::endl;
            exit(1);
        }
    }

       
    bool firstframe = true;
    //Loop through the dump file using the parser
    //for(size_t i=0; i<timestep; i++) {
    for(size_t i=startframe; i<endframe; i++) {
        if ((i % freq) != 0){
          parser.next_frame();
          continue;
        }


        //The actual functions from the parser
        atoms = parser.get_coord();
        quats = parser.get_quat();
        vects_f = parser.get_vect(quats,'f');
        //vects_v = parser.get_vect(quats,'v');
        vects_u = parser.get_vect(quats,'u');
    
        double lpdot,ltdot, bin, normA, normB;
        double dx,dy,dz,dr2;
        std::vector<double> vectA(3);
        std::vector<double> vectB(3);
        std::vector<double> vectC(3);
        std::vector<double> vectD(3);
        for (size_t j = nignore; j<natoms-nignore-1; j++){
          //vect A bond vector from j to j+1
          vectA[0] = atoms[j+1][0] - atoms[j][0];
          vectA[1] = atoms[j+1][1] - atoms[j][1];
          vectA[2] = atoms[j+1][2] - atoms[j][2];
          normA = sqrt(vectA[0]*vectA[0] + vectA[1]*vectA[1] + vectA[2]*vectA[2]);
          vectA[0] /= normA;
          vectA[1] /= normA;
          vectA[2] /= normA;

          l0_sum += normA;
          l0_count ++;

          //vect C f vector from j 
          vectC[0] = vects_f[j][0];
          vectC[1] = vects_f[j][1];
          vectC[2] = vects_f[j][2];


          for (size_t k = j; k<natoms-nignore-1; k++){
            //vect B bond vector from k to k+1
            vectB[0] = atoms[k+1][0] - atoms[k][0]; 
            vectB[1] = atoms[k+1][1] - atoms[k][1]; 
            vectB[2] = atoms[k+1][2] - atoms[k][2]; 
            normB = sqrt(vectB[0]*vectB[0] + vectB[1]*vectB[1] + vectB[2]*vectB[2]);
            vectB[0] /= normB;
            vectB[1] /= normB;
            vectB[2] /= normB;

            lpdot = vectA[0]*vectB[0] +  vectA[1]*vectB[1] +  vectA[2]*vectB[2];
            if (lpdot > 1) lpdot =1;
            if (lpdot <-1) lpdot =-1;

            //r2 for lp
            dx = atoms[k][0] - atoms[j][0];
            dy = atoms[k][1] - atoms[j][1];
            dz = atoms[k][2] - atoms[j][2];
            dr2 = dx*dx + dy*dy + dz*dz;

            //vect D f vector from j 
            vectD[0] = vects_f[k][0];
            vectD[1] = vects_f[k][1];
            vectD[2] = vects_f[k][2];

            //correct ltdot by rotating vectD around vector B by w0?
            //vectD  = rotate(vects_f[k], vects_u[k],-w0);
            //vectD  = rotate(vects_f[k], vectB,-w0);

            ltdot = vectC[0]*vectD[0] +  vectC[1]*vectD[1] +  vectC[2]*vectD[2];
            if (ltdot > 1) ltdot =1;
            if (ltdot <-1) ltdot =-1;


            //angle = acos(dot)*180./M_PI;
            //ofile << angleinfo[j].angle << "\t";
            bin = k-j;
            if ((bin >= 0) && (bin < nbins)){
              lpdotsum[bin] += lpdot;
              lpdotcount[bin] ++;

              lpr2sum[bin] += dr2;
              lpr2count[bin] ++;

              ltdotsum[bin] += ltdot;
              ltdotcount[bin] ++;
            }
            else{
                std::cout <<" FREAK OUT" << std::endl;
                exit(1);
            }
          }
        } 

        parser.next_frame();
        if (firstframe) firstframe = false;
    }  

    //setup file to write to
    std::ofstream ofile;
    ofile.open(outfilename.c_str());
    ofile << "# <contour length> <bond correlation> <twist correlation>" << std::endl;
    double ltavg, lpavg,lpr2avg;
    l0_avg = (double) l0_sum / l0_count;
    for (size_t i = 0; i < nbins; i++){
      lpavg = (double) lpdotsum[i] / lpdotcount[i];
      lpr2avg = (double) lpr2sum[i] / lpr2count[i];
      ltavg = (double) ltdotsum[i] / ltdotcount[i];
      ofile << i*l0_avg << "\t";
      ofile << lpavg << "\t";
      ofile << lpr2avg << "\t";
      ofile << ltavg << std::endl;
    }
    ofile.close();
   
}
