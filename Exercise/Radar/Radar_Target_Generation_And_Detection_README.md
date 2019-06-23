# Radar Target Generation and Detection

## Selection of Training, Guard cells and offset
Below are the image of Range-Amplitude and Doppler-Ampltiude plot. 
![alt text](Images/Range_Ampltude_Plot.jpg)
![alt text](Images/Doppler_Ampltude_Plot.jpg)
As we can see, the signal spread is more in Doppler axis and very less in Range axis. So to not include the signal tail in noise estimation, the number of Guard Cells should be more in Doppler axis and few sufficient Guard Cells in Range axis. Hence the numbers 3 and 6 for Guard Cells in Range and Doppler axis respectively.

## Implementation of 2D CFAR
## Steps taken to suppress the non-thresholded cells at the edges
