Serial implementation:
<br/>

Our implementation of Serial has a runtime of O(n). We accomplish this by applying the Binning
strategy, which involves us taking the array of particles and splitting them based on its location.
We set up each bin to contain a neighborhood of particles. In our implementation, we decided to
use a linked list instead of an array. This was to reduce the time it takes to insert a particle into a
bin, since for an array we would need to “guess” the maximum number of particles in any bin. In
addition, if our guess was wrong then we would need to reallocate the array, which would be
more expensive. We also tried using vectors but no matter how we implemented it we were
unable to get a better time than the linked list. We believe this is due to the overhead for
inserting.
<br/>

Once we set up the bins we then compute forces for each particle by looping through each bin
then applying force. We also loop through the surrounding 8 neighboring bins to apply force.
Afterwards, we loop through all the particles and move them all. In our initial implementation,
we would delete all nodes and then recreate them in order to rebin for the next simulation.
However, we realized that we could get much faster time by having a stack to reuse our bins.
<br/>

With these optimization we were able to achieve the following times:
<br/>
![image](https://github.com/tanverzahed/particle-optimization/assets/113176044/2bbc303c-6121-4a20-b558-a276fe3b115b)
<br/>

Evidently our slope is ~1 meaning that our runtime is O(n). Also below is a log-log plot that
shows our linear runtime.
![image](https://github.com/tanverzahed/particle-optimization/assets/113176044/b300e1ce-154a-41dd-a890-94ceb1de7eb8)
<br/>


Drawbacks:
One major flaw in our serial implementation is that we do not utilize spatial locality fully. While,
our bin is an array of linked lists. Each linked list does not store contiguous memory meaning
that when we loop over the particles in a bin we are not using spatial locality at all. One solution
would be to do an array but as we discussed above that does not turn out to be any faster.

---
OpenMP implementation: <br/>

![image](https://github.com/tanverzahed/particle-optimization/assets/113176044/850a89e7-9d7b-4f08-87d5-bf35263d1ec7)
<br/>

Our OpenMP implementation closely mirrors the serial implementation, with notable differences
introduced to accommodate parallel execution. Chief among these alterations are the inclusion of
OpenMP pragmas and the abandonment of stack-based node reuse due to concurrency issues
arising from multiple thread access.
<br/>


Synchronization Strategy:
<br/>

Ensuring data consistency during parallel updates is paramount for correctness. However, this
necessitates synchronization mechanisms that can potentially introduce performance bottlenecks.
To address this concern, our strategy involves minimizing the number of critical sections and
leveraging reduction clauses in OpenMP loops. This approach aims to strike a balance between
correctness and performance, optimizing parallel execution without sacrificing data integrity.
After some profiling we saw that the average number of particles per bin was within the range of
0 to 6. And were pretty evenly distributed, as such our implementation does not spend too much
time trying to distribute the work as it is mostly balanced as it is.
<br/>


Challenges and Solutions:<br/>

The efficiency plots indicate that the code scales well up to 2 threads, after which the efficiency
decreases. This is likely due to the overhead of thread management and synchronization., which
can significantly degrade performance. Initially we thought this issue to be false sharing, so to
mitigate this issue, we propose padding the particle struct to align it precisely with cache line
boundaries, typically 64 bytes. By doing so, we aimed to eliminate unnecessary cache line
contention, thereby enhancing overall performance.
<br/>

Padding can be achieve with an implementation shown in the figure below <br/>

![image](https://github.com/tanverzahed/particle-optimization/assets/113176044/fc280c06-0780-4c27-a23a-c2268b8ed79f)
<br/>

Implementing the proposed padding solution, as discussed previously, did not yield a noticeable
improvement. Initially, the unmodified struct occupies 48 bytes, allowing approximately 1.5
particle_t instances to fit within a single cache line. Consequently, if two particles, sharing the
same half of the cache line, are subject to a write operation, false sharing may occur. However,
the likelihood of this scenario transpiring is diminished due to the random distribution of
particles and the uncorrelated assignment of bins with their initial memory storage. Hence, the
probability of false sharing happening for this reason remains low.
