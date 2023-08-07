# AR_Detection_Algorithm
A rewrite/reimplementation of an AR detection algorithm by Kennett (2021) (see https://github.com/daemonkennett/ar_detection). This version uses xarray instead of iris, and updates the AR axis algorithm to follow the maximum IVT subject to direction constraints. 

The 'header_aralgorithm', 'functions_aralgorithm', and 'main_aralgorithm' need to be put in the same folder.




INPUTS:
- netCDF files containing monthly zonal and meridional IVT components. Files must follow the filenaming convention 'yyyymm.nc' (e.g. '201601.nc' for January 2016), and they must contain the following netCDF variables:
    - 'lon': variable containing the longitudes of the grid
    - 'lat': variable containing the latitudes of the grid
    - 'IVTu_yyyymmddhh': variable(s) containing the zonal component of IVT in the region of interest, at the year (yyyy), month (mm), day (dd), and UTC hour (hh) of interest (e.g. 'IVTu_2016010100' for January 1st, 2016, at 00 UTC)
    - 'IVTv_yyyymmddhh': variable(s) containing the meridional component of IVT in the region of interest, at the year (yyyy), month (mm), day (dd), and UTC hour (hh) of interest (e.g. 'IVTv_2016010100' for January 1st, 2016, at 00 UTC)
- A .nc file containing a mask for land (nonzero values) and sea (zero value). An example file for the ERA5 model can be found in this repository.


 OUTPUTS:
 - netCDF files containing the object shape and axis of all valid atmospheric rivers in the time and geographic period of interest. One netCDF file per object is generated, with the filenaming convention 'yyyy_mm_[time]_[object number]' where 'time' is the time step of the input data and 'object_number' is the number of the object as labeled in the code by ndimage.label (so for example, '2016_01_12_43.nc' would be the data for the object labeled 43 at the 12th time step in January 2016). See an example file in this repository.
     - If the 'write_valid_objects_stats_to_netcdf' flag is set to True, then the file will also contain some stats about the object:
         - The length of the calculated axis
         - The end-to-end straight-line distance of the calculated object
         - The ratio of the previous 2 quantities, which gives a measure of how "curvy" the AR is
         - The poleward mean IVT of the object in kg * m**-1 * s**-1
         - The mean direction of the object's IVT in degrees
         - The percentage of the object's cells which deviated more than 45 degrees from this mean direction
         - The deviation between the object's mean IVT direction and its orientation (defined as the straight line between the axis endpoints) in degrees
    - If the 'regrid_final_results' flag is set to True, then the objects' shapes and axes will be regridded to the uder's defined grid using xESMF
