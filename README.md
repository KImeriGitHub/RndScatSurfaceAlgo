# RndScatSurfaceAlgo
Algorithm to paper 'Realizing Machine Learning through Scattering Surfaces with Rapid Wave Propagation'

Sorry, those files are not well updated and commented.

These are all matlab files.
For chapter 3 you might want to run MNISTTest.m . For the plot use RndMatGen.m

For chapter 4 you might want to run testoptDiskLoc . 

---------------------------------------------------
Further explanation for some files:
---------------------------------------------------
MNISTTest:
- Loads digit images for the training and testing
- Lower length of loopvariable ii to not do it too often. I.e. use 'ii=1:1'.
- gets a random matrix W by placing disks on the plane and sending waves through the disks. Done in getRndScatterMat
---------------------------------------------------
RndMatGen:
- Places disks randomly on the plane. Several configurations are given to place sources and receivers on the plane. Comment out that one you desire.
- The plotting works by splitting up the process on the different cores. Lower xPlot, yPlot length to achieve faster computaiton.
- Gives different Statistical evaluation on the end.
---------------------------------------------------
PlotGen:
- Several configurations are given to place disks on the plane. Comment out that one you desire.
- The plotting works by splitting up the process on the different cores. Lower xPlot, yPlot length to achieve faster computaiton
---------------------------------------------------
testoptDiskLoc:
- loads digitimages data for testing and training
- reduces data for speed
- lower kRange and NRange to make the program faster
---------------------------------------------------
optimalDiskLoc:
- minmizes U N - Y to find best disk location to get the ML algorithm. The disks are infinitesimally small.
---------------------------------------------------
testgetOptCircles
- Similar to optimalDiskLoc, but uses not infinitesimally small disks instead tries different sized disks.
- Results are not very satisfying.
---------------------------------------------------
