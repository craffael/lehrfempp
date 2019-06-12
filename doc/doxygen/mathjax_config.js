// This file defines a number of macros which can be used in 
// formulas of the doxygen documentation.

MathJax.Hub.Config({
  TeX: {
    Macros: {
      vec: ["\\boldsymbol{#1}", 1],
      grad: "\\vec{\\operatorname{grad}}",
      norm: ["\\left\\lVert#1\\right\\rVert", 1]
    }
  }
});