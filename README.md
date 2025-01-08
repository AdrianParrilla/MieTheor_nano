## Mie_theoryEV

A Matlab script to derive nanoparticle size by flow cytometry applying Mie theory.

It implements the theoretical model presented by [De Rond et al. (2018)](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpcy.43) to describe nanoparticle light scattering. Using Mie theory, it correlates scattering cross section (σ<sub>s</sub>) to the size and refractive index (RI) of the particles. The model enables the calibration of arbitraty scattering units from the flow cytometer to size units, facilitating data interpretation and comparison between different instruments. 

The main funcitons of the script are based on the work done by [Jan Schäfer](https://de.mathworks.com/matlabcentral/fileexchange/36831-matscat) on the MatScat package. The code was adapted to consider only the collection angles of the detector and a specific numerical aperture (NA). 

The test data folder contains the theoretical scattering intensities for polystyrene and silica beads as well as for extracellular vesicles (EVs) considering 3 different internal RI. It also contains a FCS file of EVs from Jurkat cells that were stain with a membrane dye to improve detection, and the output obtained after interpolating the data into the model with R.
