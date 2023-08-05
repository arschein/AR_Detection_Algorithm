from header_aralgorithm import *
from functions_aralgorithm import *


#%% All scalar vars
restrict_landfalling = 0 #set to 1 if only care about landfalling ARs
restrict_geography = 1 #set to 1 if only care about objects in a geographic region defined by the user
[llcr_lat, llcr_lon, urcr_lat, urcr_lon] = [18,-179.75,54,-116] #coordinates of the lower left corner and upper right corner of the region to care about. MUST BE VALID COORDINATES IN THE MODEL USED (e.g. must be multiples of 0.25 deg if using ERA5 grid). MUST use 180 deg convention, west negative!!!
min_ivt_threshold = 150 
min_poleward_ivt = 50 
min_obj_size = 60
max_obj_size = 200000 #maximum number of grid cells an object can occupy
min_length = 1000 # Minimum length of AR axis in km
min_span = 1000 # Minimum required straight-line distance between AR axis start/end points in km
min_aspect_ratio = 2 # Minimum length/width ratio that an AR must surpass
min_deviation = 45 #units = degrees. ARs must have at least half their cells' IVT directions aligned within this much of the object's mean IVT direction

save_list_of_valid_objects = 1 #if true, will write a .json list of the indices of valid objects which can be used for other purposes

write_valid_objects_to_netcdf = 1 #if true, then all objects that make it through the filters will be written to netcdf 
write_valid_objects_stats_to_netcdf = 1 #if true, then object's stats will be written to netcdf along with object axis and mask

regrid_final_results = 1 #if true, the final results will be regridded to a grid of your choice (see relevant section at the end to adjust the parameters)

#%% All relevant directory paths - must be input as raw strings, at least on Windows

filepath_uv = 'filepath needed'
filepath_landseamask = 'filepath needed' #must also contain the landseamask netCDF filename
filepath_IVT_magnitude = 'filepath needed'
filepath_IVT_direction = 'filepath needed'
filepath_IVT_threshold = 'filepath needed'

filepath_object_axes = 'filepath needed'
filepath_object_axes_validregiononly = 'filepath needed'

filepath_object_mean_ivt_directions = 'filepath needed'
filepath_object_mean_ivt_directions_validregiononly = 'filepath needed'

filepath_save_final_results = 'filepath needed'
filepath_save_final_results_validregiononly = 'filepath needed'

filepath_save_regridded_results = 'filepath needed'
filepath_save_regridded_results_validregiononly = 'filepath needed'

#%% Load preliminary data

allfiles = os.listdir(filepath_uv)

ivt_u_list = [xr.open_dataarray(filepath_uv+'\\'+file, cache=False) for file in allfiles if 'IVTu' in file]
ivt_v_list = [xr.open_dataarray(filepath_uv+'\\'+file, cache=False) for file in allfiles if 'IVTv' in file]

lats = ivt_u_list[0]['lat'].data
lons = ivt_v_list[0]['lon'].data

landseamask = xr.open_dataset(filepath_landseamask, cache=False)
landseamask = landseamask['LSM'][0]
landseamask.data.round() #very crude - rounds cells with >50% land to be land and the rest to be ocean. Good enough for our purposes
landseamask.data = landseamask.data.astype(int)

if not restrict_landfalling:
    landseamask.data = ((landseamask.data * 0) + 1)

# Create a mask of the same shape as the input data, to indicate Northern and Southern Hemispheres (1 for NH and -1 for SH to properly align poleward IVT calculation):
hemisphere = landseamask.copy() 
hemisphere.data[0:(int(len(lats)/2))] = 1
hemisphere.data[(int(len(lats)/2))+1:(int(len(lats)))] = -1 

#make geographic mask
geographic_mask = np.zeros(np.shape(landseamask.data))

#create valid region - make all cells there 1
valid_lat_indexes = np.where((lats < urcr_lat) & (lats > llcr_lat))[0]
if llcr_lon <= 0 and urcr_lon <= 0: #region is entirely in western hemisphere
    valid_lon_indexes = np.where((lons < 360+urcr_lon) & (lons > 360+llcr_lon))[0]
elif llcr_lon >= 0 and urcr_lon >= 0: #region is entirely in eastern hemisphere
    valid_lon_indexes = np.where((lons < urcr_lon) & (lons > llcr_lon))[0]
else: #region straddles eastern and western hemispheres - need to convert to 360 degree system and get values from there - but not relevant to what I need this to do right now so implement this later
    sys.exit('Valid region crosses hemispheres - need to implement a fix!')

if restrict_geography:    
    geographic_mask[np.ix_(valid_lat_indexes, valid_lon_indexes)] = 1
else:
    geographic_mask = np.ones(np.shape(landseamask.data))

#%% Check for, or calculate, magnitude data
master_datelist = ['201511','201512','201601','201602','201603','201604','201611','201612','201701','201702','201703','201704','201711','201712','201801','201802','201803','201804']
threshold_datelist_years = ['2016','2017','2018']
threshold_datelist_months = ['01','02']

try: #Check to see if all IVT magnitude files exist on disk
    for date in master_datelist:
        xr.open_dataarray(filepath_IVT_magnitude+'\\'+date+'.nc', cache=False)
    print('IVT magnitudes for all months exist on disk')

except: #compute IVT magnitudes if not all of them exist on disk. 
    print('Computing IVT magnitudes')
    for i in range(len(ivt_u_list)):
        start_time = timeit.default_timer() 
        calcIVTMagnitudeAndWriteToNetCDF(ivt_u_list[i], ivt_v_list[i], filepath_IVT_magnitude)
        elapsed = timeit.default_timer() - start_time
        print('Time '+str(i+1)+' of '+str(len(ivt_u_list))+' done. Time elapsed: '+str(elapsed)[0:5])

    print('IVT magnitudes for all months calculated and written to disk')
    
#%% Check for, or calculate, direction data

try: #Check to see if all IVT direction files exist on disk. Otherwise, calculate and write them to file for future use
    for date in master_datelist:
        xr.open_dataarray(filepath_IVT_direction+'\\'+date+'.nc', cache=False)
    print('IVT directions for all months exist on disk')

except:
    print('Computing IVT directions')
    for i in range(len(ivt_u_list)):
        start_time = timeit.default_timer() 
        calcIVTDirectionAndWriteToNetCDF(ivt_u_list[i], ivt_v_list[i], filepath_IVT_direction)
        elapsed = timeit.default_timer() - start_time
        print('Time '+str(i+1)+' of '+str(len(ivt_u_list))+' done. Time elapsed: '+str(elapsed)[0:5])
    
    print('IVT directions for all months calculated and written to disk')
    
#%% Load or calculate threshold data
# NOTE: thresholds only exist for predefined months with 2 months of data before and after the month in question!

try:
    for year in threshold_datelist_years:
        for month in threshold_datelist_months:
            xr.open_dataarray(filepath_IVT_threshold+'\\'+year+month+'.nc', cache=False)
    print('IVT thresholds for all relevant months exist on disk')
            
except: 
    ivt_magnitude_list = []
    for date in master_datelist:
        ivt_magnitude_list.append(xr.open_dataarray(filepath_IVT_magnitude+'\\'+date+'.nc', cache=False))
    print('Magnitudes read from file')
    
    print('Computing IVT thresholds')
    for year in threshold_datelist_years:
        for month in threshold_datelist_months:
            start_time = timeit.default_timer()
            idx_center = [array.long_name[-6:] for array in ivt_u_list].index(year+month) #magnitude doesn't have long_name to center the month upon
            ivt_threshold_5months = xr.concat(ivt_magnitude_list[idx_center-2 : idx_center+3], dim='time') #note -2 but +3 as slicing clips on the right but not the left. Arrays must be in sequential time order for this to work
            ivt_threshold = ivt_threshold_5months.quantile(0.85, dim='time')
            ivt_threshold = np.maximum(ivt_threshold, min_ivt_threshold)
            ivt_threshold.to_netcdf(filepath_IVT_threshold+'\\'+year+month+'.nc')
            elapsed = timeit.default_timer() - start_time
            print(year+month+' done. Time elapsed: '+str(elapsed)[0:5])
    print('IVT threshold for all months done')
    

    
#%% Start main loop 
# Compute objects and landfall + size + geographic filters

master_valid_list = []
master_datelist_years = ['2016','2017','2018']
master_datelist_months = ['01','02']

for yyyy in master_datelist_years:
    for mm in master_datelist_months:
        
        start_time_master = timeit.default_timer()

        ivt_magnitude_currentmonth = xr.open_dataarray(filepath_IVT_magnitude+'\\'+yyyy+mm+'.nc', cache=False)
        ivt_direction_currentmonth = xr.open_dataarray(filepath_IVT_direction+'\\'+yyyy+mm+'.nc', cache=False)
        ivt_threshold_currentmonth = xr.open_dataarray(filepath_IVT_threshold+'\\'+yyyy+mm+'.nc', cache=False)
        
        object_mask = np.greater(ivt_magnitude_currentmonth, ivt_threshold_currentmonth) 
        object_mask.data = object_mask.data.astype(int)
        
        num_objs_list = []
        objs_list = []
        object_mask.data, num_objs = ndimage.label(object_mask.data) #I have no idea why ndimage.label won't function on a binary array time-by-time - need to have this first 
        for i in range(object_mask.shape[0]):
            objs, num_objs = ndimage.label(object_mask[i].data)
            objs_list.append(objs) #necessary so that the label indexes reset at each time instead of being summed up over time and running into the thousands, making it a headache to identify object sizes using the below code
            num_objs_list.append(num_objs) 
        print('Potential AR object identification done ('+yyyy+mm+')')
        
        landfall_filter = []
        if restrict_landfalling:
            for i in range(len(num_objs_list)):
                x = []
                y = landseamask.data * objs_list[i]
                for j in range(num_objs_list[i]):
                    if (j+1) not in y: 
                        x.append(j+1)
                landfall_filter.append(x)
                del y #memory intensive
        else: #needed for later code
            for i in range(len(num_objs_list)):
                landfall_filter.append([])
        print('Landfall filter done ('+yyyy+mm+')')
        
        size_filter = []
        for i in range(len(num_objs_list)):
            x = []
            object_sizes = ndimage.sum((objs_list[i] != 0).astype(int), objs_list[i], list(range(1, num_objs_list[i]+1)))
            for j in range(num_objs_list[i]):
                if (j+1) not in landfall_filter[i]:
                    if object_sizes[j] <= min_obj_size or object_sizes[j] >= max_obj_size:
                        x.append(j+1)
            size_filter.append(x)
        print('Size filter done ('+yyyy+mm+')')
        
        geographic_filter = []
        for i in range(len(num_objs_list)):
            x = []
            y = objs_list[i] * geographic_mask #potential objects at a given timestep are all disjoint, so can do one multiplication at the start and then check intersections in the loop
            for j in range(num_objs_list[i]):
                if (j+1) not in landfall_filter[i]+size_filter[i]:
                    if (j+1) not in y: #does the object in question NOT intersect the valid region?
                        x.append(j+1)
            geographic_filter.append(x)
        print('Geographic filter done ('+yyyy+mm+')')
        
          
        #%% Get all AR axes for the current month from disk, or calculate them
        
        try: 
            print('Reading axes from disk ('+yyyy+mm+')')
            object_axis_list = []
            axis_coords_list = []
            for i, num_objs in enumerate(num_objs_list):
                object_axis_list_currenttime = []
                axis_coords_list_currenttime = []
                for j in range(num_objs):
                    if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                        if restrict_geography == 1:
                            filepath_and_filename = filepath_object_axes_validregiononly+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        else:
                            filepath_and_filename = filepath_object_axes+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        object_axis_currenttime_currentobject = xr.open_dataarray(filepath_and_filename, cache=False) #too big to cache in memory
                        object_axis_list_currenttime.append(object_axis_currenttime_currentobject)
                        
                        axis_coords_list_currenttime_currentobject = np.argwhere(object_axis_currenttime_currentobject.data != 0)
                        axis_coords_list_currenttime.append(axis_coords_list_currenttime_currentobject)
                    else:
                        object_axis_list_currenttime.append([])
                        axis_coords_list_currenttime.append([])
                object_axis_list.append(object_axis_list_currenttime)
                axis_coords_list.append(axis_coords_list_currenttime)
                
                printPercentageComplete(i, num_objs_list)
            print('All axes read from disk ('+yyyy+mm+')') 
            
        except: #axis data not on disk, write it  
            print('Read failed. Computing axes ('+yyyy+mm+')') 
            zero_array = 0*object_mask[0].copy()
            for i in range(len(num_objs_list)):
                object_axis_list_currenttime, axis_coords_list_currenttime = calculateAllAxesAtGivenTime(num_objs_list[i], objs_list[i], landfall_filter[i], size_filter[i], geographic_filter[i], ivt_magnitude_currentmonth[i], ivt_direction_currentmonth[i], zero_array, lats, lons)
                for j in range(len(object_axis_list_currenttime)): #This only writes the axis of valid objects! When reading these back in from file and appending them to a list, make sure to append [] for objects that don't exist!
                    if type(object_axis_list_currenttime[j]) is not list: #gross but fast
                        if restrict_geography:
                            filepath_and_filename = filepath_object_axes_validregiononly+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        else:
                            filepath_and_filename = filepath_object_axes+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        
                        #doing this manually instead of with 'to_netcdf()' saves about 4 MB per axis... well worth it for disk space savings considering we have tens of thousands of objects
                        with Dataset(filepath_and_filename, 'w') as f:
                            f.createDimension('lon',len(lons))
                            lonvar = f.createVariable('lon', 'float32',('lon'))
                            lonvar[:] = lons
                            
                            f.createDimension('lat', len(lats))
                            latvar = f.createVariable('lat', 'float32', ('lat'))
                            latvar[:] = lats
                            
                            object_axis_var = f.createVariable('object_axis', 'i4', ('lat','lon'), compression = 'zlib', complevel=9)
                            object_axis_var[:] = object_axis_list_currenttime[j]
                            
                print('All axes for time '+str(i+1)+' of '+str(len(num_objs_list))+' computed ('+yyyy+mm+')')
            print('All axes calculated ('+yyyy+mm+')')
            
            print('Reading axes from disk ('+yyyy+mm+')')
            object_axis_list = []
            axis_coords_list = []
            for i, num_objs in enumerate(num_objs_list):
                object_axis_list_currenttime = []
                axis_coords_list_currenttime = []
                for j in range(num_objs):
                    if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                        if restrict_geography:
                            filepath_and_filename = filepath_object_axes_validregiononly+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        else:
                            filepath_and_filename = filepath_object_axes+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        object_axis_currenttime_currentobject = xr.open_dataarray(filepath_and_filename, cache=False) #too big to cache in memory
                        object_axis_list_currenttime.append(object_axis_currenttime_currentobject)
                        
                        axis_coords_list_currenttime_currentobject = np.argwhere(object_axis_currenttime_currentobject.data != 0)
                        axis_coords_list_currenttime.append(axis_coords_list_currenttime_currentobject)
                    else:
                        object_axis_list_currenttime.append([])
                        axis_coords_list_currenttime.append([])
                object_axis_list.append(object_axis_list_currenttime)
                axis_coords_list.append(axis_coords_list_currenttime)
                printPercentageComplete(i, num_objs_list)
            print('All axes read from disk ('+yyyy+mm+')') 

        #%% Calculate axis lengths
        
        print('Calculating axis lengths ('+yyyy+mm+')')
        
        def calcAxisLength(axis_coords):
            length = 0 
            if np.ndim(axis_coords) > 1: #account for edge case where axis is only 1 cell (happens sometimes)
                for k in range(len(axis_coords)-1):
                    lat1 = lats[axis_coords[k][0]]
                    lat2 = lats[axis_coords[k+1][0]]
                    lon1 = lons[axis_coords[k][1]]
                    lon2 = lons[axis_coords[k+1][1]]
                    length += geopy.distance.distance((lat1,lon1), (lat2,lon2)).km
            return length
        
        axis_length_list = []
        for i, num_obj in enumerate(num_objs_list):
            step_length_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                    sorted_axis = calcProperAxis(axis_coords_list[i][j]) #needed to put the axis in proper order for length calculation (order matters!)
                    step_length_list.append(calcAxisLength(sorted_axis))
                else:
                    step_length_list.append(0)
            axis_length_list.append(step_length_list)
            printPercentageComplete(i, num_objs_list)
        print('All axis lengths calculated ('+yyyy+mm+')') 

        #%% Calculate surface area
        # Grid resolution is usually given in degrees - but this is not the same distance everywhere on the globe! 
        # Calculating the area of each cell is extremely expensive, though. So for each object, get the min lon/lat, max lon/lat, calculate the average grid cell distance, then multiply that by the number of grid points in the object. Should be good enough - it's not super important to get this exact
       
        print('Calculating object surface areas ('+yyyy+mm+')')
        object_area_list = []
        for i, num_obj in enumerate(num_objs_list):
            area_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                    grid_points_list = np.argwhere(objs_list[i] == (j+1)) 
                    lat_idx_min = min([x[0] for x in grid_points_list])
                    lon_idx_min = min([x[1] for x in grid_points_list])
                    lat_idx_max = max([x[0] for x in grid_points_list])
                    lon_idx_max = max([x[1] for x in grid_points_list])
                    lat1 = lats[lat_idx_min]
                    lon1 = lons[lon_idx_min]
                    lat2 = lats[lat_idx_max]
                    lon2 = lons[lon_idx_max]
                    grid_lat_resolution = geopy.distance.distance((lat1, lon1), (lat2, lon1)).km / (len(lats[lat_idx_min : lat_idx_max+1]) - 1)
                    grid_lon_resolution = geopy.distance.distance((lat1, lon1), (lat1, lon2)).km / (len(lons[lon_idx_min : lon_idx_max+1]) - 1)
                    area_list.append(grid_points_list.shape[0] * grid_lat_resolution * grid_lon_resolution) #num points in object * average cell size in km**2
                else:
                    area_list.append(0)
            object_area_list.append(area_list)
            printPercentageComplete(i, num_objs_list)
        print('All object surface areas calculated ('+yyyy+mm+')')

        #%% Calculate object width
        #Width := surface area / axis length
       
        print('Calculating object widths ('+yyyy+mm+')')
        object_width_list = []
        for i, num_obj in enumerate(num_objs_list):
            width_list = []
            for j in range(num_obj):
                if (axis_length_list[i][j] > 0):
                    width_list.append(object_area_list[i][j] / axis_length_list[i][j])
                else:
                    width_list.append(0)
            object_width_list.append(width_list)
            printPercentageComplete(i, num_objs_list)
        print('Object widths calculated ('+yyyy+mm+')')

        #%% Calculate object orientations and distances between axis start and end (note the latter is a straight-line distance, while axis_length_list sums up the whole axis length, i.e. it will be >= distance from axis start to end)
       
        print('Calculating object orientations and axis distances ('+yyyy+mm+')')
        object_orientation_list = []
        axis_distance_list = []
        for i, num_obj in enumerate(num_objs_list):
            dir_list = []
            ax_dist_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                    if axis_length_list[i][j]>0:
                        sorted_axis = calcProperAxis(axis_coords_list[i][j])
                        coords1 = sorted_axis[0]
                        coords2 = sorted_axis[-1]
                        lat1 = lats[coords1[0]]
                        lon1 = lons[coords1[1]]
                        lat2 = lats[coords2[0]]
                        lon2 = lons[coords2[1]]
                    else:
                        coords1 = [0,0]
                        coords2 = [0,0]
                        lat1 = 0
                        lon1 = 0
                        lat2 = 0
                        lon2 = 0
                    dir_list.append((np.arctan2((lons[coords1[1]] - lons[coords2[1]]), (lats[coords1[0]] - lats[coords2[0]])) * 180/np.pi +180) %360)
                    ax_dist_list.append(geopy.distance.distance( (lat1, lon1), (lat2, lon2) ).km)
                else:
                    dir_list.append([])   
                    ax_dist_list.append([])
            object_orientation_list.append(dir_list)
            axis_distance_list.append(ax_dist_list)
            printPercentageComplete(i, num_objs_list)
        print('Object orientations calculated ('+yyyy+mm+')')
        print('Distances between axis starts and ends calculated ('+yyyy+mm+')')

        #%% AR criteria check 1: length check
        
        print('Computing AR length filter ('+yyyy+mm+')')
        length_filter = []
        for i, num_obj in enumerate(num_objs_list):
            filter_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]:
                    if axis_length_list[i][j] <= min_length:
                        filter_list.append(j+1)
            length_filter.append(filter_list)
            printPercentageComplete(i, num_objs_list)
        print('AR length filter done ('+yyyy+mm+')')

        #%% AR criteria check 2: narrowness check
        
        print('Computing AR narrowness filter ('+yyyy+mm+')')
        narrowness_filter = []
        for i, num_obj in enumerate(num_objs_list):
            filter_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]:
                    if (axis_length_list[i][j] / object_width_list[i][j]) <= min_aspect_ratio:
                        filter_list.append(j+1)
            narrowness_filter.append(filter_list)
            printPercentageComplete(i, num_objs_list)
        print('AR narrowness filter done ('+yyyy+mm+')')

        #%% AR criteria check 3: mean IVT poleward component > min threshold
        #much more efficient to preprocess the poleward data for the month and then just do a mean check in the loop
        
        print('Computing AR poleward IVT filter ('+yyyy+mm+')')
        poleward_ivt_filter = []
        valid_objects_poleward_ivt_list = []
        valid_objects_poleward_ivt_list_indexes = []
        poleward_ivt = ivt_v_list[[array.long_name[-6:] for array in ivt_v_list].index(yyyy+mm)] * hemisphere.data
        for i, num_obj in enumerate(num_objs_list):
            filter_list = []
            valid_list = []
            index_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]:
                    poleward_ivt_mean = np.mean(poleward_ivt[i].data[np.where(objs_list[i] == j+1)])
                    if poleward_ivt_mean <= min_poleward_ivt:
                        filter_list.append(j+1)
                    else:
                        valid_list.append(poleward_ivt_mean)
                        index_list.append(j+1)
            
            poleward_ivt_filter.append(filter_list)
            valid_objects_poleward_ivt_list.append(valid_list)
            valid_objects_poleward_ivt_list_indexes.append(index_list)
            printPercentageComplete(i, num_objs_list)
        print('AR poleward IVT filter done ('+yyyy+mm+')')

        #%% Calculate mean object IVT direction
        # need to compute the direction for the mean IVT components, NOT the mean of the IVT directions - they are not quite the same! 
        # read from disk; otherwise, create and save a txt file with the mean directions
        #txt file format: [yyyy] [mm] [t] [(j+1)] [mean dir]
        
        try: #if file exists and is complete, read mean directions from it
            print('Reading AR mean directions from disk ('+yyyy+mm+')')
            if restrict_geography:
                filepath_ivt_mean_dir = filepath_object_mean_ivt_directions_validregiononly
            else:
                filepath_ivt_mean_dir = filepath_object_mean_ivt_directions
            with open(filepath_ivt_mean_dir+'\\'+yyyy+mm+'.txt', 'r') as file:
                ivt_direction_mean_list = []
                for i, num_obj in enumerate(num_objs_list):
                    mean_list = []
                    for j in range(num_obj):
                        if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]+poleward_ivt_filter[i]:
                            x = file.readline()[:-1] #gets rid of \n
                            y = x.split('_')
                            mean_list.append(float(y[-1]))
                            # print(str(i)+' | '+str(j+1)+' | '+y[-2])
                        else:
                            mean_list.append([])
                    ivt_direction_mean_list.append(mean_list)
                print('All IVT mean directions read from file ('+yyyy+mm+')')
        

        except: #mean direction file doesn't exist, calculate them. 
        # Uses a specific setup to greatly speed up mean calculation compared to previous AR algorithm
            print('Read from disk failed. Computing mean directions ('+yyyy+mm+')')
            if restrict_geography:
                filepath_ivt_mean_dir = filepath_object_mean_ivt_directions_validregiononly
            else:
                filepath_ivt_mean_dir = filepath_object_mean_ivt_directions
            with open(filepath_ivt_mean_dir+'\\'+yyyy+mm+'.txt', 'w') as file:
                ivt_u_currentmonth = xr.open_dataarray(filepath_uv+'\\'+'era5_IVTu_'+yyyy+mm+'.nc', cache=False)
                ivt_v_currentmonth = xr.open_dataarray(filepath_uv+'\\'+'era5_IVTv_'+yyyy+mm+'.nc', cache=False)
                
                #loading data into memory is expensive but well worth it - speeds up calculation by almost three orders of magnitude!
                print('Loading u, v data into memory ('+yyyy+mm+')')
                ivt_u_currentmonth = ivt_u_currentmonth.data
                ivt_v_currentmonth = ivt_v_currentmonth.data
                
                print('Calculating mean directions ('+yyyy+mm+')')
                ivt_direction_mean_list = []
                for i, num_obj in enumerate(num_objs_list):
                    mean_list = []
                    for j in range(num_obj):
                        if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]+poleward_ivt_filter[i]:
                            mean_u_ivt = np.mean(ivt_u_currentmonth[i][np.where(objs_list[i]==j+1)])
                            mean_v_ivt = np.mean(ivt_v_currentmonth[i][np.where(objs_list[i]==j+1)])
                            mean_direction = (np.arctan2(mean_u_ivt, mean_v_ivt) * 180/np.pi + 180) % 360
                            mean_list.append(mean_direction)
                            writestr = yyyy +'_'+ mm +'_'+ str(i) +'_'+ str(j+1) +'_'+ str(mean_direction)+'\n'
                            file.write(writestr)
                        else:
                            mean_list.append([])
                    ivt_direction_mean_list.append(mean_list)
                    printPercentageComplete(i, num_objs_list)
                print('All IVT mean directions calculated and written to disk ('+yyyy+mm+')')
                
                del ivt_u_currentmonth, ivt_v_currentmonth #maybe release memory, who knows
                
        #%% AR criteria check 4: coherence check
        # If > 50% of the object's cells have IVT deviating >= 45 deg from the object's mean IVT direction, the object is filtered
        
        print('Computing IVT direction coherence filter ('+yyyy+mm+')')
        coherence_filter = []
        valid_objects_ivt_percentage_deviant_list = []
        valid_objects_ivt_percentage_deviant_list_indexes = []
        for i, num_obj in enumerate(num_objs_list):
            filter_list = []
            valid_list = []
            index_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]+poleward_ivt_filter[i]:
                    deviation = np.abs(ivt_direction_mean_list[i][j] - ivt_direction_currentmonth[i].data[np.where(objs_list[i] == j+1)])
                    if np.mean(deviation > 45) > 0.5: 
                        filter_list.append(j+1)
                    else: 
                        valid_list.append(np.mean(deviation > 45)*100)
                        index_list.append(j+1)
            
            coherence_filter.append(filter_list)
            valid_objects_ivt_percentage_deviant_list.append(valid_list)
            valid_objects_ivt_percentage_deviant_list_indexes.append(index_list)
            printPercentageComplete(i, num_objs_list)
        print('IVT direction coherence filter done ('+yyyy+mm+')')
        
        #%% AR criteria check 5: consistency between axis direction and mean IVT
        # if the object's orientation (as defined in object_orientation_list previously) deviates from the mean IVT direction by >45 degrees, the object is filtered
        # the straight-line distance between first and last axis gridpoints (as defined in axis_distance_list previously) is < min_span
        
        print('Computing IVT direction consistency filter ('+yyyy+mm+')')
        consistency_filter = []
        valid_objects_ivt_direction_deviation_list = []
        valid_objects_ivt_direction_deviation_list_indexes = []
        for i, num_obj in enumerate(num_objs_list):
            filter_list = []
            valid_list = []
            index_list = []
            for j in range(num_obj):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]+poleward_ivt_filter[i]+coherence_filter[i]:
                    # object orientation is calculated as a directed quantity, but really doesn't have a direction - it needs to be tested against ivt mean direction (which really is directed) in both possible orientations. 
                    # both plus and minus 180 degrees need to be considered because we don't know if rotating clockwise or counterclockwise is better, so take both
                    deviation = np.min([np.abs(ivt_direction_mean_list[i][j] - object_orientation_list[i][j]), 
                                        np.abs(ivt_direction_mean_list[i][j] - (object_orientation_list[i][j]-180)), 
                                        np.abs(ivt_direction_mean_list[i][j] - (object_orientation_list[i][j]+180))])
                    if (deviation > 45) or (axis_distance_list[i][j] < min_span):
                        filter_list.append(j+1)
                    else:
                        valid_list.append(deviation)
                        index_list.append(j+1)
                        
                        
            consistency_filter.append(filter_list)
            valid_objects_ivt_direction_deviation_list.append(valid_list)
            valid_objects_ivt_direction_deviation_list_indexes.append(index_list)
            printPercentageComplete(i, num_objs_list)
        print('IVT direction consistency filter done ('+yyyy+mm+')')
        
        #%% Get all valid AR objects
        
        valid_objects = []
        for i, num_obj in enumerate(num_objs_list):
            x = []
            for j in range(num_objs_list[i]):
                if (j+1) not in landfall_filter[i]+size_filter[i]+geographic_filter[i]+length_filter[i]+narrowness_filter[i]+poleward_ivt_filter[i]+coherence_filter[i]+consistency_filter[i]:
                    x.append(j+1)
            valid_objects.append(x)

        if save_list_of_valid_objects:
            master_valid_list.append(valid_objects)
        
        
        #%% Write data for valid objects to netCDF, with all relevant quantities if selected
        
        if write_valid_objects_to_netcdf:
            print("Writing valid objects' stats to netCDF")
            if restrict_geography:
                filepath = filepath_save_final_results_validregiononly
            else:
                filepath = filepath_save_final_results
            
            for i, num_obj in enumerate(num_objs_list):
                for j in range(num_obj):
                    if (j+1) in valid_objects[i]:
                        
                        currentobject_mask = np.zeros(np.shape(objs_list[i]))
                        currentobject_mask[np.where(objs_list[i] == j+1)] = 1
                        
                        filepath_and_filename = filepath+'\\'+yyyy+'_'+mm+'_'+str(i)+'_'+str(j+1)+'.nc'
                        with Dataset(filepath_and_filename, 'w') as f:
                            f.createDimension('lon',len(lons))
                            lonvar = f.createVariable('lon', 'float32',('lon'))
                            lonvar[:] = lons
                            
                            f.createDimension('lat', len(lats))
                            latvar = f.createVariable('lat', 'float32', ('lat'))
                            latvar[:] = lats
                            
                            object_axis_var = f.createVariable('object_axis', 'i4', ('lat','lon'), compression = 'zlib', complevel=9)
                            object_axis_var[:] = object_axis_list[i][j]
                            
                            object_mask_var = f.createVariable('object_mask', 'i4', ('lat','lon'), compression = 'zlib', complevel=9)
                            object_mask_var[:] = currentobject_mask
                            
                            if write_valid_objects_stats_to_netcdf:
                                axis_length_var = f.createVariable('axis_length', 'float32', compression = 'zlib', complevel=9)
                                axis_length_var[:] = axis_length_list[i][j]
                                
                                axis_distance_var = f.createVariable('axis_distance', 'float32', compression = 'zlib', complevel=9)
                                axis_distance_var[:] = axis_distance_list[i][j]
                                
                                axis_ratio_var = f.createVariable('axis_ratio', 'float32', compression = 'zlib', complevel=9)
                                axis_ratio_var[:] = axis_length_list[i][j] / axis_distance_list[i][j]
                                
                                poleward_ivt_mean_var = f.createVariable('ivt_poleward_mean', 'float32', compression = 'zlib', complevel=9)
                                poleward_ivt_mean_var[:] = valid_objects_poleward_ivt_list[i][valid_objects_poleward_ivt_list_indexes[i].index(j+1)]
                                
                                ivt_direction_mean_var = f.createVariable('ivt_direction_mean', 'float32', compression = 'zlib', complevel=9)
                                ivt_direction_mean_var[:] = ivt_direction_mean_list[i][j]
                                
                                ivt_deviant_percentage_var = f.createVariable('ivt_deviant_percentage', 'float32', compression = 'zlib', complevel=9)
                                ivt_deviant_percentage_var[:] = valid_objects_ivt_percentage_deviant_list[i][valid_objects_ivt_percentage_deviant_list_indexes[i].index(j+1)]
                                
                                ivt_direction_deviation_var = f.createVariable('ivt_direction_deviation', 'float32', compression = 'zlib', complevel=9)
                                ivt_direction_deviation_var[:] = valid_objects_ivt_direction_deviation_list[i][valid_objects_ivt_direction_deviation_list_indexes[i].index(j+1)]
                            
                printPercentageComplete(i, num_objs_list)
                
        elapsed_master = timeit.default_timer() - start_time_master        
        print('!!!!!!!!!!!!! '+yyyy+' '+mm+' done!!!!!!!!!!!!! Time elapsed: '+str(elapsed_master)[0:6])

#%% Save list of master objects for easy access for testing. Saved into the same directory as this main .py file

if save_list_of_valid_objects:
    with open('list_of_valid_objects.json','w') as f:
        json.dump(master_valid_list, f, indent=2)
        
#%% Regrid final results (if selected)

if regrid_final_results:
    print('Regridding final AR masks and axes')
    
    if restrict_geography:
        filepath = filepath_save_final_results_validregiononly
    else:
        filepath = filepath_save_final_results
        
    allfiles = os.listdir(filepath)
    
    lats = getVarFromNetCDF(filepath+'\\'+allfiles[0], 'lat')
    lons = getVarFromNetCDF(filepath+'\\'+allfiles[0], 'lon')
    
    x = xr.open_dataset(filepath+'\\'+allfiles[0], cache=False)
	#Adjust lats, lons, and grid spacing here
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], np.arange(-90, 90, 0.08), {"units": "degrees_north"}),
            "lon": (["lon"], np.arange(0, 360, 0.08), {"units": "degrees_east"}),
        }
    )
    
    start_time = timeit.default_timer()
    print('Making regridder')
    regridder = xesmf.Regridder(x, ds_out, 'nearest_s2d') #nearest_s2d is the one to use for binary data!
    elapsed = timeit.default_timer() - start_time
    print('Regridder done (elapsed: '+str(elapsed)[0:6]+')')
    
    if restrict_geography:
        filepath_regridded = filepath_save_regridded_results_validregiononly
    else:
        filepath_regridded = filepath_save_regridded_results

    for i, filename in enumerate(allfiles):
        x = xr.open_dataset(filepath+'\\'+filename)
        dr_out = regridder(x, keep_attrs=True)
        xr.Dataset.to_netcdf(dr_out, 'temp.nc') #save to temp netcdf file, then re-open it and write custom file to save almost 4 mb of disk space per file 
        
        regridded_lats = np.array(getVarFromNetCDF('temp.nc', 'lat'))
        regridded_lons = np.array(getVarFromNetCDF('temp.nc', 'lon'))
        regridded_axis = np.array(getVarFromNetCDF('temp.nc', 'object_axis'))
        regridded_mask = np.array(getVarFromNetCDF('temp.nc', 'object_mask'))

        #saves ~72 mb of disk space per object... well worth it!
        with Dataset(filepath_regridded+'\\'+filename, 'w') as f:
            f.createDimension('lon',len(regridded_lons))
            lonvar = f.createVariable('lon', 'float32',('lon'))
            lonvar[:] = np.array(getVarFromNetCDF('temp.nc', 'lon'))
            
            f.createDimension('lat', len(regridded_lats))
            latvar = f.createVariable('lat', 'float32', ('lat'))
            latvar[:] = np.array(getVarFromNetCDF('temp.nc', 'lat'))
            
            object_axis_var = f.createVariable('object_axis', 'i4', ('lat','lon'), compression = 'zlib', complevel=9)
            object_axis_var[:] = np.array(getVarFromNetCDF('temp.nc', 'object_axis'))
            
            object_mask_var = f.createVariable('object_mask', 'i4', ('lat','lon'), compression = 'zlib', complevel=9)
            object_mask_var[:] = np.array(getVarFromNetCDF('temp.nc', 'object_mask'))
            
        printPercentageComplete(i, allfiles)