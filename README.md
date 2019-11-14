## ariaDNE: a robustly implemented algorithm for Dirichlet Energy of the Normal
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1465949.svg)](https://doi.org/10.5281/zenodo.1465949)

![teaser](./images/teaser.jpg?raw=true)

**DNE background:**
Dirichlet Normal Energy (DNE) quantifies local curvature in 3D surfaces; a discretized version for 3D meshes is effective for biological studies of morphological surfaces. For a continuous 2D manifold embedded in 3D, DNE is defined by,

![continuous](https://latex.codecogs.com/gif.latex?\int_{\Omega}&space;\mathscr{K}_1^2&space;&plus;&space;K_2^2&space;~&space;dA)

where *K*<sub>1</sub>, *K*<sub>2</sub> are the principal curvatures. The discrete version is defined by suming local curvature of each vertex across a surface mesh, weighted by the area distributed to that vertex,

![discrete](https://latex.codecogs.com/gif.latex?\sum_{V_j}&space;(K_{j,&space;1}^2&space;&plus;&space;K_{j,2}^2&space;)&space;\cdot&space;\mbox{Area}(V_j))

Unlike other shape metrics, DNE is landmark-free and independent of the surface’s initial position, orientation and scale. However, recent studies found that the currently published implementation of DNE is sensitive to various surface preparation procedures, raising concerns regarding comparability and objectivity of utilizing DNE in morphological research. This package provides a robustly implemented algorithm for DNE (ariaDNE).

**ariaDNE algorithm:** We leverage the observation that local curvature can be stably estimated via a weighted PCA approach. The procedure is outlined as follows:

1. Apply a weighted PCA localized around the query point by means of the Gaussian kernel fucntion.
2. Find the principal component that is cloest to the vertex normal, and estimate the local curvature by its principal score
![curvature](https://latex.codecogs.com/gif.latex?\frac{\sigma_{\mbox{chosen}}}{\sigma_1&space;&plus;&space;\sigma_2&space;&plus;&space;\sigma_3})
3. Summing the local curvature across the mesh and obtain the glocal DNE.

For the exact procedure, see this document ([pdf](./images/algorithm.pdf)).

**Code:**
The function 'ariaDNE.m' computes ariaDNE of a mesh surface.

**Demos:**
This package includes a demo for computing ariaDNE on a 3D mesh representing a second mandibular tooth from *Callicebus*. Download the package and run example.m in Matlab.

[2019/10] We provide code for plotting the curvature field of ariaDNE as in the paper. Run plotexample.m in Matlab. 

[2019/11] We provide a demo for computing ariaDNE for a folder of 3D meshes saved in '.ply' files. Run exampleDNEfolder.m in Matlab. 


**Compatibility and dependencies:**
- Developed and tested with Matlab 2017b on Mac and Linux. 

**Disclaimer:** The code is for acadmic use only. Please contact the author [Shan Shan](https://sshanshans.github.io) for any questions or bugs.

**License:**
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.



