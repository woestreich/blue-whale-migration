# blue-whale-migration

This repository contains code and data necessary to generate the analyses and figures associated with the manuscript **"Animal-borne metrics enable acoustic detection of a dispersed migration"** (currently in submission).

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



  
