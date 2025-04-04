# SysidMatrix

## Transfer Matrix Measurements using a recursive least-squares algorithm

Detailed discussion is in the CEBAF tech-note [JLAB-TN-25-14](https://jlabdoc.jlab.org/docushare/dsweb/Get/Document-296114/25-014.pdf).

The Matlab file tmm_2D.m reconstructs 2x2 transfer matrices and tmm_4D does the same 
for 4x4 transfer matrices which can be arbitrarily coupled. 

I extended the method to handle nonlinear maps which is described in  CEBAF tech-note [JLAB-TN-25-16](https://jlabdoc.jlab.org/docushare/dsweb/Get/Document-296130/25-016.pdf).
The Matlab script with the example is called tmm_2d_nonlin.m.
