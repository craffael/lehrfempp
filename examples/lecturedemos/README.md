# LehrFEM++ demo codes for Chapter 2 of the Course "Numerical Methods for PDEs"

The files in this folder demonstrate the use of various facilities offered by the
LehrFEM++ library. These codes are discussed in detail in Section 2.7 of the 
"Numerical Methods for Partial Differential Equation" lecture notes. 

Remark: Parts of the codes are included in this document as snippets. Thus some of the
source codes feature lines like `/* SAM_LISTING_BEGIN_1 */`. These control the inclusion
of the C++ code into LaTeX files. 

- `lecturedemomesh`: shows how to use the container capabilities of `lf::mesh::Mesh`
- `lecturedemodof`: teaches how to access d.o.f. information via
  `lf::assemble::DofHandler`
- `lecturedemoassemble`: demonstrates the use of local provider objects for assembly of
  Galerkin matrices and right-hand side vectors. 
- `lecturedemoquad`: introduces the use of the local quadrature rules in LehrFEM++
- `lecturedemorefine`: shows how to use (local) mesh refinement
- `lecturedemotwonorm`: illustrates different ways of computing the L2 norm of a finite
  element function 
  
  


