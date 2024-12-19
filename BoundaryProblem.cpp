// Dongjin Lee, GIST, 2024 fall
// last updated: 2024-12-14
// Scientific programming course
// Probabilistic potential theory with brownian motion simulation
// Configuration: Starting at a point apart from horizontal plane with distance d=1.
// Implementation: Direction is determined by exponential function with complex power.
// Assummed symmetric case, so that binning only considers length from origin to reached position.

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#define SIMPLE_SPRNG

#include <mpi.h>
#include "/home/sprng5/include/sprng_cpp.h"

using namespace std;

// Parameters
double d = 1.0;
double epsilon = 1e-10;
int initial = 5;
int log10_N = 7;
double bin_d_ratio = 0.1;

// Function to simulate a Brownian motion step
double random_walk(double d, double epsilon) {
    double z = d;
    complex<double> xy(0.0, 0.0);

    while (z > epsilon) {
        complex<double> tmp_direction = exp(M_PI * sprng() * complex<double>(0, 1));
        if (imag(tmp_direction) > sprng()) //Metric for surface area on sphere contains sine term, which should be considered by this way.
        {
        xy += z * imag(tmp_direction) * exp(2 * M_PI * sprng() * complex<double>(0, 1));
        z += z * real(tmp_direction);
        }
    }
    return std::abs(xy);
}



// Function for simulating and binning.
vector<long long> binning_running(double d, double epsilon, int bin_num, int num) {
    double bin_width = bin_d_ratio * d;
    vector<long long> bins(bin_num);
    for (int i = 0; i < num; i++) {
        int bin_index = min(static_cast<int>(abs(random_walk(d, epsilon)) / bin_width), bin_num-1);
        bins[bin_index]++;
    }
    return bins;
}

int main(int argc, char *argv[]) {
    int myid, nprocs;
    int gtype, bin_num;
    MPI::Init(argc, argv);       /* Initialize MPI                          */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);	/* find process id                 */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* find number of processes      */
    
    if(myid == 0)
    {
      printf("Type in a generator type (integers: 0,1,2,3,4,5):  ");
      scanf("%d", &gtype);
      printf("Type in the number of bins (e.g., 1000 is recommended):  ");
      scanf("%d", &bin_num);
    }
    MPI::COMM_WORLD.Bcast(&gtype,1,MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&bin_num,1,MPI::INT,0);
    
    // Verifying distinct seeds are set.
    sleep(myid);
    init_sprng(make_sprng_seed(),SPRNG_DEFAULT, gtype);
    printf("\n\nProcess %d, print information about stream:\n", myid);
    print_sprng();
    
    string file_name;
    for (int n=initial;n<=log10_N;n++){
      // Binning data
      vector<long long> binned_data = binning_running(d, epsilon, bin_num, pow(10, n));
      file_name = "./results/"+to_string(myid)+"_"+to_string(n)+".txt";
      // Saving data.
      ofstream fout(file_name);
      for (int j=0;j<bin_num;j++){
        fout << binned_data[j] << "\n";
      }
      fout.close();
    }
    // Saving information on parameters.
    if (myid==0){
      ofstream fout("./results/info.txt");
      fout << d << "\n";
      fout << epsilon << "\n";
      fout << initial << "\n";
      fout << log10_N << "\n";
      fout << bin_num << "\n";
      fout << nprocs << "\n";
      fout << bin_d_ratio << "\n";
      fout.close();
    }
    
    MPI::Finalize();		/* Terminate MPI                           */

    return 0;
}
