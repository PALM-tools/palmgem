# Leaf Area Density (LAD)
Utilization of open source dataset with Leaf Area Index (LAI) and Tree Canopy Height into PALM gridded LAD. 
## LAI
LAI is obtained from source **Add missing source** in raster format. Input is uploaded via import shell script into input schema. 

## Canopy Height
Raster data with canopy height can be downloaded from **Add missing source**. Input is uploaded via import shell script into input postgresql schema. 

## Processing
Due to spatial resolution 10m, we decided to use the simplest solution based on homogenous LAD distribution. The LAD is calculated as \
LAD(z,y,x) = LAI(y,x) / Canopy_height(y,x) * max(min((grid_bottom_height(z,y,x) - Canopy_height) / dz), 1, 0) \
More advanced approach based on Markkanen et al. (2003) will be added in future releases.