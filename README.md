# M-ISIC
A proposed blind extracting algorithm of hidden data embedded by the multi-carrier spread-spectrum scheme. See Details in  "Lattice-Aided Extraction of Spread-Spectrum Hidden Data".

# Main function
Compare different extracting algorithms including:
1. MMSE(Minimum mean-square error)
2. MF(Matched Filter)
3. JADE(Jean-Franccois Cardoso,et al."Blind beamforming for non Gaussian signals","High-order contrasts for independent component analysis")
4. M_IGLS(M. Li, et al. “Extracting spread-spectrum hidden data from digital media”)
5. M_ISIC(Proposed method)

This program simulates a data hiding scenario. Firstly,a random watermark <img src="http://latex.codecogs.com/gif.latex?B" /> will be embeded in the host.After that <img src="http://latex.codecogs.com/gif.latex?B'" /> will be extracted using different algorithm.Then calculate BER of <img src="http://latex.codecogs.com/gif.latex?B'" />.
The differences between the algorithms mentioned above are mainly shown by BER in different image distortion.  