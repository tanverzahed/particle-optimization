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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "common.h"
#include <stack>

#define density 0.0005 
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

struct Node {
    particle_t* particle;
    Node* next;
};
std::stack<Node*> availableNodes;

inline void push(Node** head, particle_t* particle) {
    Node* new_node = availableNodes.top();
    availableNodes.pop();
    new_node->particle = particle;
    new_node->next = *head;
    *head = new_node;
} 

int main(int argc, char** argv) {
    int navg, nabsavg = 0;
    double davg, dmin, absmin = 1.0, absavg = 0.0;

    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
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

    double size = sqrt(density * n);
    double bin_size = cutoff;
    Node* nodePool = new Node[n];
    for (int i = 0; i < n; i++) {
        availableNodes.push(&nodePool[i]);
    }
    int bin_count = ceil(size / bin_size);

    Node** bins = new Node*[bin_count * bin_count];
    for (int i = 0; i < bin_count * bin_count; i++) {
        bins[i] = nullptr;
    }

    int step, i, jx_start, jx_end, jy_start, jy_end, jx, jy, x, y;
    Node* current;

    //  simulate a number of time steps
    double simulation_time = read_timer( );
	
    for(step = 0; step < NSTEPS; ++step) {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        // compute forces

        for (i = 0; i < bin_count * bin_count; ++i) {
            current = bins[i];
            while (current != nullptr) {
                current->particle->ax = current->particle->ay = 0;
                x = current->particle->x / bin_size;
                y = current->particle->y / bin_size;

                jx_start = max(0, x - 1);
                jx_end = min(bin_count - 1, x + 1);
                jy_start = max(0, y - 1);
                jy_end = min(bin_count - 1, y + 1);

                for (jx = jx_start; jx <= jx_end; ++jx) {
                    for (jy = jy_start; jy <= jy_end; ++jy) {
                        Node *current2 = bins[jy * bin_count + jx];
                        while (current2 != nullptr) {
                            apply_force(*(current->particle), *(current2->particle), &dmin, &davg, &navg);
                            current2 = current2->next;
                        }
                    }
                }
                current = current->next;
            }
        }

        // Delete Node objects in each bin
        for (i = 0; i < bin_count * bin_count; ++i) {
            current = bins[i];
            while (current != nullptr) {
                Node* next = current->next;
                availableNodes.push(current);
                current = next;
            }
            bins[i] = nullptr;
        }

        // Move particles
        for (i = 0; i < n; ++i) {
            move(particles[i]);
            push(&bins[int(particles[i].y / bin_size) * bin_count + int(particles[i].x / bin_size)], &particles[i]);
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          // Computing statistical data
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //  save if necessary
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
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
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    // Clearing space
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
