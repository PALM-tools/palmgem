# Surface Fraction
The implementation of surface fraction is done using grid fraction approach. 
In current version, surface fraction is not supported with lod2. Only surface fraction related to palm types are available.
Surface fraction is calculated as AREA(type) / AREA(grid cell). Limit of 5% is set to as minimum fraction. 
The rest is nullified. 
Building grid cell are cell with building area majority. 