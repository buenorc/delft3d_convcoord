## Delft3d ConvCoord 

This algorithm converts input data used in Delftd (.ldb, .grd, and .xyz) from Decimal Degrees to Standard UTM 
The algorithm the geodetic datum WGS84 (but can be easily modified to another geodetic datum)

Example: 

1) grid80.grd                     - Delft3D Grid
2) sau_reservoir_samples_mod.xyz  - Sample points defining water depth 
3) sau-reservoir.ldb              - Landboudnary 

Output: Same files but all of them converted to UTM
