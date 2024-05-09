/* ------------
 * The code is adapted from the XSEDE online course Applications of Parallel Computing. 
 * The copyright belongs to all the XSEDE and the University of California Berkeley staff
 * that moderate this online course as well as the University of Toronto CSC367 staff.
 * This code is provided solely for the use of students taking the CSC367 course at 
 * the University of Toronto.
 * Copying for purposes other than this use is expressly prohibited. 
 * All forms of distribution of this code, whether as given or with 
 * any changes, are expressly prohibited. 
 * -------------
*/
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <vector>
#include <cmath>

#define density 0.0005
#define cutoff  0.01

int main( int argc, char **argv ) {    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //  process command line parameters
    if( find_option( argc, argv, "-h" ) >= 0 ) {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    //  set up MPI
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //  allocate generic resources
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    

    double size = sqrt(density * n);
    int num_bins = ceil(size / cutoff);
    int total_bins = num_bins * num_bins;

    int bins_per_proc = (total_bins + n_proc - 1) / n_proc;

    int *bin_offsets = (int *)malloc((n_proc) * sizeof(int));
    for (int i = 0; i < n_proc; ++i) {
        bin_offsets[i] = min(i * bins_per_proc, total_bins);
    }

    int *bin_counts = (int *)malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; ++i) {
        bin_counts[i] = min(bins_per_proc, total_bins - bin_offsets[i]);
    }

    // create vector for each bin
    std::vector<particle_t> *bins[total_bins];
    for (int i = 0; i < total_bins; ++i) {
        bins[i] = new std::vector<particle_t>();
        bins[i]->reserve(n); // reserve space for particles
    }

    std::vector<particle_t> all_particles;
    all_particles.reserve(n);
    int sizes[n_proc];

    set_size( n );
    if( rank == 0 ) {
        init_particles( n, particles );

        for (int i = 0; i < n; ++i) {
            int x = int(particles[i].x / cutoff);
            int y = int(particles[i].y / cutoff);
            bins[y * num_bins + x]->push_back(particles[i]);
        }
    }
    
    // broadcast all bins
    for (int i = 0; i < total_bins; ++i) {
        int size = bins[i]->size();
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bins[i]->resize(size);
        MPI_Bcast(bins[i]->data(), size, PARTICLE, 0, MPI_COMM_WORLD);
    }

    std::vector<particle_t> local;
    local.reserve(n);

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for( int step = 0; step < NSTEPS; ++step ) {
        navg = 0, dmin = 1.0, davg = 0.0;
        
        //  save current step if necessary
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );

        double time = read_timer();
        // compute all forces
        for (int i = 0; i < bin_counts[rank]; ++i) {
            int offset = bin_offsets[rank] + i;
            for (particle_t &p : *bins[offset]) {
                p.ax = p.ay = 0;
                int x = int(p.x / cutoff);
                int y = int(p.y / cutoff);
                for (int jx = x - 1; jx <= x + 1; ++jx) {
                    for (int jy = y - 1; jy <= y + 1; ++jy) {
                        if (jx >= 0 && jx < num_bins && jy >= 0 && jy < num_bins) {
                            for (particle_t &j : *bins[jy * num_bins + jx]) {
                                apply_force(p, j, &dmin, &davg, &navg);
                            }
                        }
                    }
                }
                local.push_back(p);
            }
        }

        if( find_option( argc, argv, "-no" ) == -1 ) {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

          if (rank == 0) {
            // Computing statistical data
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        // clear bins
        for (int i = 0; i < total_bins; ++i) {
            bins[i]->resize(0);
        }

        //  move particles
        for (particle_t &p : local) {
            move(p);
        }

        all_particles.resize(n);

        int size = local.size();
        MPI_Allgather(&size, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

        int offsets[n_proc];
        offsets[0] = 0;
        for (int i = 1; i < n_proc; ++i) {
            offsets[i] = offsets[i - 1] + sizes[i - 1];
        }

        MPI_Allgatherv(local.data(), local.size(), PARTICLE, all_particles.data(), sizes, offsets, PARTICLE, MPI_COMM_WORLD);
        local.resize(0);

        for (particle_t &p : all_particles) {
            int x = int(p.x / cutoff);
            int y = int(p.y / cutoff);
            bins[y * num_bins + x]->push_back(p);
        }
    }
    simulation_time = read_timer() - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

      if( find_option( argc, argv, "-no" ) == -1 )
      {
        if (nabsavg) absavg /= nabsavg;
      //  -the minimum distance absmin between 2 particles during the run of the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
      //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
      if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");     
        
      // Printing summary data
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //  release resources
    if ( fsum )
        fclose( fsum );
    free( bin_offsets );
    free( bin_counts );
    // free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
