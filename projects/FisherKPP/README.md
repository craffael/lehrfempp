This folder contains the code developed in the context of the [bachelor thesis](http://www.sam.math.ethz.ch/%7Ehiptmair/StudentProjects/Loher.Amelie/BScThesis_LoherAmelie_FisherKPP.pdf) "Solving the Fisher/KPP equation on the globe" by Amélie Loher under the advisorship of Prof. Dr. Ralf Hiptmair. The content encloses the main simulation, the simulation for the model problems on simplified domains, as well as the implementation for the convergence analysis.

**strangsplitting.h**: Header file for the implementation of the Strang splitting method applied to the FisherKPP equation. 

**strangsplitting.cc**: Implementation of the Strang splitting method applied to the FisherKPP equation. The discretisation is based on the method of lines approach consisting of the spatial Galerkin discretisation combined with suitable timestepping methods. 

**strangsplitting_main.cc**: Main file for the simulation on the earth's surface. It assembles the carrying capacity maps, and computes the population density. 

**modelproblem_circle_main.cc**: Model problem posed on a two-dimensional ball. Computes the solution based on the numerical scheme implemented in **strangsplitting.cc**. The initial data is either originating from one source only, or from two sources. 

**modelproblem_threecircles_main.cc**: Model problem posed on three balls at positive distance from eachother. Computes the solution based on the numerical scheme implemented in **strangsplitting.cc**. The boundary conditions are either homogeneous Neumann boundary conditions, or else non-local flux conditions based on decaying function kernels. There are three different function kernels to choose from.

**modelproblem_island_main.cc**: Model problem posed on two arbitrarily shaped islands at positive distance from eachother. Computes the solution based on the numerical scheme implemented in **strangsplitting.cc**. The diffusion coefficient can either be set to normal or very low. 

**convergence_meshwidth_main.cc**: File for convergence studies. The solution is computed on the two arbitrarily shaped islands, as in **modelproblem_island_main.cc**. Non-local flux boundary conditions are chosen. Here, the mesh width is refined, and the time step size is kept constant. Solutions on coarses meshes are interpolated onto the finest mesh. The L2-norm of the error is computed. 

**convergence_tau_main.cc**: File for convergence studies. The solution is computed on the two arbitrarily shaped islands, as in **modelproblem_island_main.cc**. Non-local flux boundary conditions are chosen. Here, the time step size is varied, and the mesh width is kept constant. The L2-norm of the error is computed. 

**convergence_linked_main.cc**:  File for convergence studies. The solution is computed on the two arbitrarily shaped islands, as in **modelproblem_island_main.cc**. Non-local flux boundary conditions are chosen. Here, the mesh width is refined, and the time step size is varied accordingly. Solutions on coarses meshes are interpolated onto the finest mesh. The L2-norm of the error is computed. 

**convergence_scripts**: This folder contains Matlab scripts to plot the error norms computed with the files **convergence_meshwidth_main.cc**, **convergence_tau_main.cc**, and **convergence_linked_main.cc**. Each file has its corresponing Matlab script.

**meshes**: Contains all the underlying meshes used in the files above.
