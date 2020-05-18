# blue-whale-migration

This repository contains code and data necessary to generate the analyses and figures associated with the manuscript **"Animal-borne metrics enable acoustic detection of a dispersed migration"** (currently in submission). Manuscript authors: William K. Oestreich, James A. Fahlbush, David E. Cade, John Calambokidis, Tetyana Margolina, John Joseph, Ari S. Friedlaender, Megan F. McKenna, Alison K. Stimpert, Brandon L. Southall, Jeremy A. Goldbogen, and John P. Ryan.

The repository is organized as follows:
* **Main directory:** all scripts and miscellaneous files (e.g., .shp file for central coast of CA).
* **Sub-directory "tag_data":** all tag-associated data 
  * tag-on-whale image 
  * data from individual tag deployments
  * analyzed/aggregated tag data.
* **Sub-directory "acoustic_data":** all acoustics data 
  * daily power spectral density from MBARI's MARS hydrophone 
  * Aug - Dec (peak blue whale song seasong) 50th and 90th percentile daily mean spectrum levels
  * modeled received level results for the MARS hydrophone (modeled as in Collins, M. D. A split-step Pade solution for the parabolic equation method. J. Acoust. Soc. Am. 93, 1736-1742 (1993). 
  * .mat file containing ~30 minutes of MARS data with a prominent example blue whale song
  * .mat file containing geographic location of the MARS hydrophone
  
*Full acoustic dataset from the MARS hydrophone (Aug 2015 - present) is in the process of being organized, staged, and uploaded to a public repository through a collaboration between the Monterey Bay Aquarium Research Institute (MBARI) and Amazon Web Services (AWS). Acoustic data necessary for the full analysis of blue whale song in the present paper is provided here.*

* **Sub-directory "m_map":** m_map mapping scripts (see citation below)

All main text figures (1-4) and extended data figure 2 are generated using the scripts found in the main directory. Extended data figure 1 can be generated from the full-resolution tag data for deployment Bm181021-TDR11. This full-resolution dataset is far too large for storage via github, but is available upon request to the authors.

A variety of publicly-available software packages are used here and cited below:

brewermap.m: Brewer, Cynthia A., 2002. http://www.ColorBrewer.org, accessed May 17, 2020.

SolarAzEl.m: Koblick, D. Vectorized Solar Azimuth and Elevation Estimation (https://www.mathworks.com/matlabcentral/fileexchange/23051-vectorized-solar-azimuth-and-elevation-estimation). MATLAB Cent. File Exch. (2020).
  
m_map: Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

cmocean.m: Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M., & DiMarco, S. F. (2016). True colors of oceanography. Oceanography, 29(3), 10.
