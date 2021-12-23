# The HotFill Algorithm

This repository has the source code of Hotfill. Hotfill is an algorithm, based on dynamic programming paradigm, that finds a tool-path for solid raster infill that keeps the cooling time intervals between depositions below user-specified limits. The development of Hotfill comes from the observation that one of the factors responsible for the mechanical strength of objects manufactured by thermoplastic extrusion 3D printers is the adhesion between adjacent material beads in the same layer. This, in turn, depends on the time lapsed between their deposition. We obversed that Hotfill is fast and produces tool-paths that have low fabrication time while respecting the cooling time constraints. We also included a set of representative slices of STL models for several objects where low mechanical strength could be an issue, such as human prostheses, mechanical parts, and coat hangers. 

For illustration purposes, we show below four examples of outputs for HotFill with different cooling time limits, where the green horizontal lines are the raster lines, the dashed black lines represents the machine jumps, the purple and green small lines represents the raster links, and the red dashed lines the coolest contact for that path.

![](https://github.com/ecassiana/hotfill/blob/main/hotfillfig.jpg)

# Video

# Citation

If you want more information about the algorithm or find this work useful for your research, please read and cite our paper:

```
@article{Elis2022,
  title = {HotFill: A Cooling Time Constrained Raster-Fill Planning Algorithm for Extrusion Additive Manufacturing}
  author = {Elis C. Nakonetchnei, Jorge Stolfi, Neri Volpato, Ricardo D. da Silva, and Rodrigo Minetto},
  journal = {?},
  volume = {?},
  pages = {1-?},
  year = {2022},
  doi = {?},
}
