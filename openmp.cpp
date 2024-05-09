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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define dt      0.0005
#define min_r   (cutoff / 100)

struct Node {
    particle_t* particle;
    Node* next;
};

void push(Node** head, particle_t* particle) {
    Node* new_node = new Node;
    new_node->particle = particle;
    new_node->next = *head;
    *head = new_node;
}   


//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    // Setup 

    double size = sqrt(density * n);
    double bin_size = cutoff;
    int bin_count = ceil(size / bin_size);

    Node** bins = new Node*[bin_count * bin_count];
    for (int i = 0; i < bin_count * bin_count; i++) {
        bins[i] = nullptr;
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin) 
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < 1000; step++ )
    {
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for (int i = 0; i < bin_count * bin_count; i++) {
            Node* current = bins[i];
            while (current != nullptr) {
                current->particle->ax = current->particle->ay = 0;
                int x = int(current->particle->x / bin_size);
                int y = int(current->particle->y / bin_size);

                int jx_start = max(0, x - 1);
                int jx_end = min(bin_count - 1, x + 1);
                int jy_start = max(0, y - 1);
                int jy_end = min(bin_count - 1, y + 1);

                for (int jx = jx_start; jx <= jx_end; jx++) {
                    for (int jy = jy_start; jy <= jy_end; jy++) {
                        Node* current2 = bins[jy * bin_count + jx];
                        while (current2 != nullptr) {
                            apply_force(*(current->particle), *(current2->particle), &dmin, &davg, &navg);
                            current2 = current2->next;
                        }
                    }
                }
                current = current->next;
            }
        }


        //
        // Delete Node objects in each bin
        //
        #pragma omp for
        for (int i = 0; i < bin_count * bin_count; i++) {
            Node* current = bins[i];
            while (current != nullptr) {
                Node* next = current->next;
                current->particle = nullptr;
                current->next = nullptr;
                current = next;
            }
            bins[i] = nullptr;
        }
        //
        //  move particles
        //
        #pragma omp for
        for (int i = 0; i < n; i++) {
            move(particles[i]);
            push(&bins[int(particles[i].y / bin_size) * bin_count + int(particles[i].x / bin_size)], &particles[i]);
        }
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
