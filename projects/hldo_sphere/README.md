# Hodge Laplacian and Dirac Operators on the surface of the 3-Sphere

This folder contains the code developed in the context of the bachelor 
thesis *"Hodge Laplacian and Dirac Operators on the surface of 
the 3-Sphere"* by Nico Graf under the supervision of 
Prof. Dr. Ralf Hiptmair.

>The thesis discusses the h-convergence of the discrete solutions 
>for the Dirac operator and the Hodge Laplacians on the surface of 
>the 3-sphere. At first, it describes the discretization with finite 
>element methods based on the 
>Lehrfem++. The h-convergence experiments 
>then show algebraic convergence with rate one for all source problems.
>Moreover, the Dirac operator turned out to be solvable numerically with 
>the Laplace Operator and a suitable modification of the load function.

The full thesis can be found in **THESIS**

---

To reproduce the results in the thesis, use the main.sh script.
The script requires an existing `build` folder in the projects
root folder containing
the compiled code. The script by default runs all the experiments in
the *experiment* section of the **THESIS**, this might take a day or
two. The results are then available in the corresponding experiment folders
`build/projects/hldo_sphere/experiments/<experiment>/results`.

Folders called `k_x_yz` contain the `vkt` files for all refinement 
level in the corresponding
`<experiment>` executed with the value `k = x.yz`. The `vtk` files can
be viewed with e.g. `paraview` and show various measurements on 
the meshes.

Files called `result_*_x_yz-u_vw.csv` contain the `<experiment>`
for a list of refinement levels and values $k \in [x.yz, u.vw]$,
tabulated in `csv` format and they can be displayed using the
scripts `build/projects/hldo_sphere/experiments/<experiment>/plot*.py`
Examples of such calls are done in the script main.sh.

To generate only partial results (run experiments only on coarse meshes e.g. up
to refinement level 4 instead of 7)
pass the option `-s | --short` to the script.
With this option the code will run
significantly faster.\
To reproduce all the experiment including the ones in the *debugging* section
of the **THESIS** use the option `-a | --all`. The two option can be
combined to generate only short versions of all experiments.\
The option `-p | --plot` allows for running the plotting scripts for the hodge laplace
experiments. This option only works together with `-s`.\
For more plots plotting scripts are provided in the experiment build folders.

