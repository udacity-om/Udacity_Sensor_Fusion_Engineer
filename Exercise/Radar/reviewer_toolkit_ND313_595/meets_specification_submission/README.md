# Radar Target Generation and Detection

## 2D CFAR Steps

To implement the 2D CFAR proecss, based on experimentation, I selected 16 training cells in the range and 12 training cells in the doppler dimension. In each dimension, I used 4 guard cells, and a threshold offset of 6 db, also mostly based on experimentation.

For the non-thresholded cells at the edges, I set the signal values to 0 to suppress effects from cells which were outside of the CUT cell matrix.