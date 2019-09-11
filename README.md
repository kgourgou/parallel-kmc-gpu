# Experimental implementation of parallel KMC with CUDA-C

Part of an initial research project was to explore 
how parallel Kinetic Monte Carlo can be mapped to the architecture 
of an NVidia GPU. 

This code contains my initial experiments to accomplish this. It supports
multiple GPUs and, up to my benchmarks at least, it offers a modest speed-up
even though it is unoptimized. It is **not** production-level code; at this
point, it can only simulate an 1D Ising system (I also have a version for 2D). It's been a while since I wrote it, 
so I can't guarantee it will run on your system. I only share it
as a proof of concept. I have since moved on to implementing code
in [SPPARKS](http://spparks.sandia.gov); see also [2]. 

Because it was written with the "fractional-step" point-of-view in mind (see [1,3]), i.e.,
that an asynchronous parallel kMC algorithm is really defined by the lattice
decomposition + computation schedule, I put a lot of focus into making schedules
easy to change. This makes the implementation of the Strang splitting a
relatively simple addition. I've been writing code to replicate this functionality
in SPPARKS.   




## References
1. Arampatzis, G., Katsoulakis, M.A., Plecháč, P., Taufer, M. and Xu, L., 2012.
   Hierarchical fractional-step approximations and parallel kinetic Monte Carlo
   algorithms. Journal of Computational Physics, 231(23), pp.7795-7814.
2. Plimpton, S., Battaile, C., Chandross, M., Holm, L., Thompson, A., Tikare,
   V., Wagner, G., Webb, E., Zhou, X., Cardona, C.G. and Slepoy, A., 2009.
   Crossing the mesoscale no-man’s land via parallel kinetic Monte Carlo. Sandia
   Report SAND2009-6226.
3. Arampatzis, G., Katsoulakis, M.A. and Plechác, P., 2014. Parallelization,
   processor communication and error analysis in lattice kinetic Monte Carlo.
   SIAM Journal on Numerical Analysis, 52(3), pp.1156-1182. 
