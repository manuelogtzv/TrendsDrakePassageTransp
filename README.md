# Trends in Drake Passage Transport and properties
Repository with data sets and code for recreating results of Nature Communications paper titled "Compensating transport trends in the Drake Passage frontal regions yield no acceleration in net transport". Directories Codes and Datasets contain the auxiliary codes and ancillary datasets to run the main routines. Main routines calculate all variables of interest and plot the results. Main routines are:

1) plotFig1_map.m:  Code for plotting Figure 1 - Area of study, plus mean temperature, salinity and velocity transects. 
2) plotFig2_trendDPtransp.m: Code for plotting Figure 2 - Time series and trends in net total, geostrophic and reference Drake Passage transport.
3) plotFig3_trendsVel.m: Code for plotting Figure 3 - Trends in the cross-transect velocity as a function of depth and distance across the Passage.
4) plotFig4_tranps_timeseries.m: Code for plotting Figure 4 - Time series, means and trends of transport per distance bin for the geostrophic, total and reference components.
5) plotFig5_streamwiseSv_mean.m: Code for plotting Figure 5 - Time series, mean and trends of total, geostrophic and reference transport per pair of SSH streamlines.
6) plotFig6_trend_spiceheave.m: Code for plotting Figure 6 - Trends in the temperature and salinity as a function of gamma (neutral density) and distance across the passage.
7) plotFig7_streamwise_eke_uv.m: Code for plotting Figure 7 - Time series, means and trends of eddy kinetic energy and eddy momentum flux as a function of SSH streamlines.
8) plotFigS1_trendDPtransp.m: Code for plotting supplementary Figure S1 - Trends for total, geostrophic and reference transport using different number of transects to test sensitivity.
9) plotFigS2_trendcurl.m: Code for plotting supplementary Figure S2 - Trends in wind stress curl.


* Note: Figure 7 was created using Adobe Illustrator.
* Note: Figure S2 requires access to a large dataset that is not suitable for GitHub. The dataset is available upon request.

**IMPORTANT**
The code shared in this repository needs the additional MATLAB libraries:

1) m_map mapping toolbox. Code available at: https://www.eoas.ubc.ca/~rich/map.html
2) PreTEOS-10 Neutral density. Code available at: http://www.teos-10.org/preteos10_software/neutral_density.html
3) TEOS-10 Seawater. Code available at: https://www.teos-10.org/software.htm#1
