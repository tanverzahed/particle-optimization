Serial implementation:
Our implementation of Serial has a runtime of O(n). We accomplish this by applying the Binning
strategy, which involves us taking the array of particles and splitting them based on its location.
We set up each bin to contain a neighborhood of particles. In our implementation, we decided to
use a linked list instead of an array. This was to reduce the time it takes to insert a particle into a
bin, since for an array we would need to “guess” the maximum number of particles in any bin. In
addition, if our guess was wrong then we would need to reallocate the array, which would be
more expensive. We also tried using vectors but no matter how we implemented it we were
unable to get a better time than the linked list. We believe this is due to the overhead for
inserting.
Once we set up the bins we then compute forces for each particle by looping through each bin
then applying force. We also loop through the surrounding 8 neighboring bins to apply force.
Afterwards, we loop through all the particles and move them all. In our initial implementation,
we would delete all nodes and then recreate them in order to rebin for the next simulation.
However, we realized that we could get much faster time by having a stack to reuse our bins.

With these optimization we were able to achieve the following times:
![image](https://github.com/tanverzahed/particle-optimization/assets/113176044/2bbc303c-6121-4a20-b558-a276fe3b115b)
