# The HotFill Algorithm

This repository has the source code of Hotfill. Hotfill is an algorithm, based on dynamic programming paradigm, that finds a tool-path for solid raster infill that keeps the cooling time intervals between depositions below user-specified limits. The development of Hotfill comes from the observation that one of the factors responsible for the mechanical strength of objects manufactured by thermoplastic extrusion 3D printers is the adhesion between adjacent material beads in the same layer. This, in turn, depends on the time lapsed between their deposition. We obversed that Hotfill is fast and produces tool-paths that have low fabrication time while respecting the cooling time constraints. We also included a set of representative slices of STL models for several objects where low mechanical strength could be an issue, such as human prostheses, mechanical parts, and coat hangers. 

For illustration purposes, we show below four examples of outputs for HotFill with different cooling time limits, where the green horizontal lines are the raster lines, the dashed black lines represents the machine jumps, the purple and green small lines represents the raster links, and the red dashed lines the coolest contact for that path.

![](https://github.com/ecassiana/hotfill/blob/main/hotfillfig.jpg)

# Video

# Installation (Linux)

01) chmod u+x *.sh
02) conda create -n hotfillenv
03) source activate hotfillenv
04) conda install pillow
05) conda install pandas
06) conda install -c bioconda wkhtmltopdf
07) pip install shapely
08) pip install pdfkit
09) pip install PyX==0.15
10) sudo apt-get install texlive-fonts-extra
11) sudo apt-get install -y texlive-latex-extra
12) sudo apt install ghostscript

# Installation (Windows 10)
You can use the Linux environment on Windows 10 to run the Hotfill tests. 
01) Install [Ubuntu on Windows 10](https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)
02) Install [Anaconda on Windows Ubuntu Terminal](https://gist.github.com/kauffmanes/5e74916617f9993bc3479f401dfec7da)
03) Follow [Installation (Linux)](https://github.com/ecassiana/hotfill#installation-linux) steps.

# Running
01) ./run.sh
This script will run a batch of tests using distinct cooling tines and bandpath sizes.
For each test it will save in a folder tests/out the imagens and information about
the fabrication time, air time, largest cooling contact, etc.

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
