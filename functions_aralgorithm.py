from header_aralgorithm import *

#%%

def getVarFromNetCDF(filepath, varname):
    #input: filepath (string), varname (string) (name of variable contained in file)
    #will return a foo string if the variable is not in the file, so can easily condition on this being returned
    with Dataset(filepath, 'r', format='NETCDF4') as f:
        allvars = list(f.variables)
        if varname in allvars:
            if varname == 'time':
                return f.variables[varname]
            else:
                return f.variables[varname][:]
        else:
            return 'NOT_IN_FILE'

        
#------------------------------------------------

def calcIVTMagnitudeAndWriteToNetCDF(u_ivt, v_ivt, filepath_to_save_to):
    #input: u and v IVT xarray DataArrays for one month, with dimensions (time x lat x lon), with long_name ending in 'yyyymm'
    #input: filepath to save to 
    magnitude = np.power(np.power(u_ivt, 2) + np.power(v_ivt, 2), 1/2)
    magnitude.to_netcdf(filepath_to_save_to+'\\'+u_ivt.long_name[-6:]+'.nc', mode='w')
    return 


def calcIVTDirectionAndWriteToNetCDF(u_ivt, v_ivt, filepath_to_save_to):
    #input: u and v IVT xarray DataArrays for one month, with dimensions (time x lat x lon), with long_name ending in 'yyyymm'
    #input: filepath to save to 
    direction = (np.arctan2(u_ivt, v_ivt) * 180 / np.pi + 180) % 360
    direction.to_netcdf(filepath_to_save_to+'\\'+u_ivt.long_name[-6:]+'.nc', mode='w')
    return

def calculateAllAxesAtGivenTime(num_objs, objs_list_currenttime, landfall_filter_currenttime, size_filter_currenttime, geographic_filter_currenttime, ivt_magnitude_currenttime, ivt_direction_currenttime, zero_array, lats, lons):
    num_consecutive_bad_steps_allowed = 3
    object_axes_list_currenttime = []
    axis_coords_list_currenttime = []
    for j in range(num_objs):
        if (j+1) not in landfall_filter_currenttime+size_filter_currenttime+geographic_filter_currenttime:
            
            axis_coords_list_currenttime_currentobject = []
            
            x = np.equal(objs_list_currenttime, (zero_array+(j+1))) #extract a binary mask of the current object
            ivt_magnitude_currenttime_currentobject = np.multiply(x, ivt_magnitude_currenttime)
            ivt_direction_currenttime_currentobject = np.multiply(x, ivt_direction_currenttime)
            ivt_magnitude_max_coords = np.where(ivt_magnitude_currenttime_currentobject.data == np.max(ivt_magnitude_currenttime_currentobject.data)) #== may be dangerous here - revisit if there's problems
            ivt_magnitude_center_coords = ivt_magnitude_max_coords #used when resetting the algorithm to the center point again to calculate the second direction
            max_lat_coord = ivt_magnitude_max_coords[0][0]
            max_lon_coord = ivt_magnitude_max_coords[1][0] 
            
            axis_array_currenttime_currentobject = zero_array.copy()
            
            axis_array_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 1 #sets the axis center on the cell with the max IVT in the object - now to run the western and eastern half of the algorithm
           
           #use the direction of IVT at this center point to decide which direction to go in first. 
           #note the direction up to here is the direction the wind is coming FROM!!! Need to rotate by 180 degrees to get the direction the wind is blowing TO!! Also note 0 = due north before rotation, i.e. 0 = due south when rotated, AND direction increases CLOCKWISE!
           
           #     directions:
           #         N
           #    315  0   45
           #      \  |  /
           #       \ | / 
           #        \|/
           #270 ----------- 90
           #        /|\
           #       / | \
           #      /  |  \
           #    225 180 135 
           
           # again, non-rotated = where the wind is coming FROM, rotated = where the wind is blowing TO
    
            ivt_direction_rotated_currenttime_currentobject = ivt_direction_currenttime_currentobject.copy()
            ivt_direction_rotated_currenttime_currentobject.data = (ivt_direction_rotated_currenttime_currentobject.data - 180) % 360 #gives the direction the wind is blowing towards
           
            
            if 0 < ivt_direction_rotated_currenttime_currentobject.data[max_lat_coord, max_lon_coord] <= 180: #direction = to the east
            #Run east then west   
                cell_record = [i*0 for i in range(num_consecutive_bad_steps_allowed)]
                #Run algorithm east 
                while ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] > 0:
                    if (0 < max_lat_coord < len(lats)-1) and (0 < max_lon_coord < len(lons)-1): #culls points on the boundary of the grid, as ndimage has issues with them
                        ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0 #ensures the algorithm will move to a new point
                       
                        #Define search directions. Note here we only consider N/NE/E/SE/S!
                        north = (np.array(max_lat_coord-1), np.array(max_lon_coord))
                        north_east = (np.array(max_lat_coord-1), np.array(max_lon_coord+1))
                        east = (np.array(max_lat_coord), np.array(max_lon_coord+1))
                        south_east = (np.array(max_lat_coord+1), np.array(max_lon_coord+1))
                        south = (np.array(max_lat_coord+1), np.array(max_lon_coord))
                                  
                        downstream_direction = ivt_direction_rotated_currenttime_currentobject.data[max_lat_coord, max_lon_coord] % 360
                       
                        if downstream_direction > 337.5 or downstream_direction <= 22.5:
                            candidate_cells = [north, north, north_east] #note only two unique candidate cells here, as we have no western direction, but make a list of 3 so the code to check direction below works
                        elif 22.5 < downstream_direction <= 67.5:
                            candidate_cells = [north, north_east, east]
                        elif 67.5 < downstream_direction <= 112.5:
                            candidate_cells = [north_east, east, south_east]
                        elif 112.5 < downstream_direction <= 157.5 :
                            candidate_cells = [east, south_east, south]
                        elif 157.5 < downstream_direction <= 202.5:
                            candidate_cells = [south_east, south_east, south] #again, only 2 unique candidates
                        else: #terminate axis here
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                            
                        #If the first candidate (i.e. adjacent cell with the greatest IVT) has a direction opposite that of the direction of the first center cell, then go to the adjacent cell with the second greatest magnitude and check; if that fails, then go to the third cell; if that fails, then give up and go to the cell with largest magnitude, regardless of direction
                        candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                        candidate_cells_direction = [ivt_direction_rotated_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells] #note these are directions the wind is blowing towards
                        index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                       
                        if not any(candidate_cells_ivt): #reached the edge of the object - all adjacent cells in the direction we care about are 0 and we can move on
                            foo = 1 #dummy condition to move algorithm along
                        elif 0 < candidate_cells_direction[index] < 180: #cell with max ivt is good, set as new center and proceed
                            ivt_magnitude_max_coords = candidate_cells[index]
                            max_lat_coord = ivt_magnitude_max_coords[0]
                            max_lon_coord = ivt_magnitude_max_coords[1]
                            cell_record.append(0)
                        else: #move to next greatest cell
                            candidate_cells_ivt[index] = 0
                            index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                            if 0 < candidate_cells_direction[index] < 180: #cell with second-largest ivt is good, set as new center and proceed
                                ivt_magnitude_max_coords = candidate_cells[index]
                                max_lat_coord = ivt_magnitude_max_coords[0]
                                max_lon_coord = ivt_magnitude_max_coords[1]
                                cell_record.append(0)
                            else: #move to third-greatest cell
                                candidate_cells_ivt[index] = 0
                                index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                if 0 < candidate_cells_direction[index] < 180: #cell with third-largest ivt is good, set as new center and proceed
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(0)
                                else: #terminate, as something is wrong  
                                    candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                                    index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(1)
                       
                        if np.sum(cell_record[-num_consecutive_bad_steps_allowed:]) == num_consecutive_bad_steps_allowed: #If the algorithm had to take X bad steps in a row, terminate this direction of axis calculation
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                            
                        else:
                            axis_array_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 1
                            axis_coords_list_currenttime_currentobject.append(ivt_magnitude_max_coords)
                        
                    else: #point on the boundary, cull it
                        ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 0
                       
                       
                #return start of search to the center
                max_lat_coord = ivt_magnitude_center_coords[0][0]
                max_lon_coord = ivt_magnitude_center_coords[1][0]
                ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 1 #to get the algorithm going again
            
                cell_record = [i*0 for i in range(num_consecutive_bad_steps_allowed)]
                #Now run algorithm west
                while ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] > 0:
                    if (0 < max_lat_coord < len(lats)-1) and (0 < max_lon_coord < len(lons)-1): #culls points on the boundary of the grid, as ndimage has issues with them
                        ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0 #ensures the algorithm will move to a new point
                        
                        #Define search directions. Note here we only consider N/NW/W/SW/S!
                        north = (np.array(max_lat_coord-1), np.array(max_lon_coord))
                        north_west = (np.array(max_lat_coord-1), np.array(max_lon_coord-1))
                        west = (np.array(max_lat_coord), np.array(max_lon_coord-1))
                        south_west = (np.array(max_lat_coord+1), np.array(max_lon_coord-1))
                        south = (np.array(max_lat_coord+1), np.array(max_lon_coord))
                                  
                        upstream_direction = ivt_direction_currenttime_currentobject.data[max_lat_coord, max_lon_coord] % 360 #note difference to previous direction - now we're looking upstream to the west
                        
                        if upstream_direction > 337.5 or upstream_direction <= 22.5:
                            candidate_cells = [north, north, north_west] #note only two unique candidate cells here, as we have no eastern direction, but make a list of 3 so the code to check direction below works
                        elif 292.5 < upstream_direction <= 337.5:
                            candidate_cells = [north, north_west, west]
                        elif 247.5 < upstream_direction <= 292.5:
                            candidate_cells = [north_west, west, south_west]
                        elif 202.5 < upstream_direction <= 247.5 :
                            candidate_cells = [west, south_west, south]
                        elif 157.5 < upstream_direction <= 202.5:
                            candidate_cells =  [south_west, south, south] #again, only 2 unique candidates
                        else: 
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                           
                        #If the first candidate (i.e. adjacent cell with the greatest IVT) has a direction opposite that of the direction of the first center cell, then go to the adjacent cell with the second greatest magnitude and check; if that fails, then go to the third cell; if that fails, then give up and go to the cell with largest magnitude, regardless of direction
                        candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                        candidate_cells_direction = [ivt_direction_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells] #note these are directions the wind is blowing FROM (different than the previous iteration)
                        index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                       
                        if not any(candidate_cells_ivt): #reached the edge of the object - all adjacent cells in the direction we care about are 0 and we can move on
                            foo = 1 #dummy condition to move algorithm along
                        elif 180 < candidate_cells_direction[index] < 359.99: #cell with max ivt is good, set as new center and proceed
                            ivt_magnitude_max_coords = candidate_cells[index]
                            max_lat_coord = ivt_magnitude_max_coords[0]
                            max_lon_coord = ivt_magnitude_max_coords[1]
                            cell_record.append(0)
                        else: #move to next greatest cell
                            candidate_cells_ivt[index] = 0
                            index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                            if 180 < candidate_cells_direction[index] < 359.99: #cell with second-largest ivt is good, set as new center and proceed
                                ivt_magnitude_max_coords = candidate_cells[index]
                                max_lat_coord = ivt_magnitude_max_coords[0]
                                max_lon_coord = ivt_magnitude_max_coords[1]
                                cell_record.append(0)
                            else: #move to third-greatest cell
                                candidate_cells_ivt[index] = 0
                                index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                if 180 < candidate_cells_direction[index] < 359.99: #cell with third-largest ivt is good, set as new center and proceed
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(0)
                                else: #terminate, as something is wrong
                                    candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                                    index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(1)
               
                        if np.sum(cell_record[-num_consecutive_bad_steps_allowed:]) == num_consecutive_bad_steps_allowed: #If the algorithm had to take X bad steps in a row, terminate this direction of axis calculation
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                            
                        else:
                            axis_array_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 1
                            axis_coords_list_currenttime_currentobject.append(ivt_magnitude_max_coords)
                           
                    else: #point on the boundary, cull it
                         ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 0      
           
            else: #object direction = to the west
            #Run west then east
                cell_record = [i*0 for i in range(num_consecutive_bad_steps_allowed)]
                #Run algorithm west
                while ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] > 0:
                    if (0 < max_lat_coord < len(lats)-1) and (0 < max_lon_coord < len(lons)-1): #culls points on the boundary of the grid, as ndimage has issues with them
                        ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0 #ensures the algorithm will move to a new point
                        
                        #Define search directions. Note here we only consider N/NW/W/SW/S!
                        north = (np.array(max_lat_coord-1), np.array(max_lon_coord))
                        north_west = (np.array(max_lat_coord-1), np.array(max_lon_coord-1))
                        west = (np.array(max_lat_coord), np.array(max_lon_coord-1))
                        south_west = (np.array(max_lat_coord+1), np.array(max_lon_coord-1))
                        south = (np.array(max_lat_coord+1), np.array(max_lon_coord))
                                  
                        downstream_direction = ivt_direction_rotated_currenttime_currentobject.data[max_lat_coord, max_lon_coord] % 360
                        
                        if downstream_direction > 337.5 or downstream_direction <= 22.5:
                             candidate_cells = [north, north, north_west] #note only two unique candidate cells here, as we have no eastern direction, but make a list of 3 so the code to check direction below works
                        elif 292.5 < downstream_direction <= 337.5:
                            candidate_cells = [north, north_west, west]
                        elif 247.5 < downstream_direction <= 292.5:
                            candidate_cells = [north_west, west, south_west]
                        elif 202.5 < downstream_direction <= 247.5 :
                            candidate_cells = [west, south_west, south]
                        elif 157.5 < downstream_direction <= 202.5:
                            candidate_cells =  [south_west, south, south] #again, only 2 unique candidates
                        else: #terminate axis here
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                           
                        #If the first candidate (i.e. adjacent cell with the greatest IVT) has a direction opposite that of the direction of the first center cell, then go to the adjacent cell with the second greatest magnitude and check; if that fails, then go to the third cell; if that fails, then give up and go to the cell with largest magnitude, regardless of direction
                        candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                        candidate_cells_direction = [ivt_direction_rotated_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells] #note these are directions the wind is blowing towards
                        index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                        
                        if not any(candidate_cells_ivt): #reached the edge of the object - all adjacent cells in the direction we care about are 0 and we can move on
                            foo = 1 #dummy condition to move algorithm along
                        elif 0 < candidate_cells_direction[index] < 180: #cell with max ivt is good, set as new center and proceed
                            ivt_magnitude_max_coords = candidate_cells[index]
                            max_lat_coord = ivt_magnitude_max_coords[0]
                            max_lon_coord = ivt_magnitude_max_coords[1]
                            cell_record.append(0)
                        else: #move to next greatest cell
                            candidate_cells_ivt[index] = 0
                            index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                            if 0 < candidate_cells_direction[index] < 180: #cell with second-largest ivt is good, set as new center and proceed
                                ivt_magnitude_max_coords = candidate_cells[index]
                                max_lat_coord = ivt_magnitude_max_coords[0]
                                max_lon_coord = ivt_magnitude_max_coords[1]
                                cell_record.append(0)
                            else: #move to third-greatest cell
                                candidate_cells_ivt[index] = 0
                                index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                if 0 < candidate_cells_direction[index] < 180: #cell with third-largest ivt is good, set as new center and proceed
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(0)
                                else: #terminate, as something is wrong
                                    candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                                    index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(1)
            
                        if np.sum(cell_record[-num_consecutive_bad_steps_allowed:]) == num_consecutive_bad_steps_allowed: #If the algorithm had to take X bad steps in a row, terminate this direction of axis calculation
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                            
                        else:
                            axis_array_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 1
                            axis_coords_list_currenttime_currentobject.append(ivt_magnitude_max_coords)
                            
                    else: #point on the boundary, cull it
                        ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 0
                        
                #return start of search to the center
                max_lat_coord = ivt_magnitude_center_coords[0][0]
                max_lon_coord = ivt_magnitude_center_coords[1][0]
                ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 1 #to get the algorithm going again
            
                cell_record = [i*0 for i in range(num_consecutive_bad_steps_allowed)]
                #Run algorithm east
                while ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] > 0:
                    if (0 < max_lat_coord < len(lats)-1) and (0 < max_lon_coord < len(lons)-1): #culls points on the boundary of the grid, as ndimage has issues with them
                        ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0 #ensures the algorithm will move to a new point
                        
                        #Define search directions. Note here we only consider N/NE/E/SE/S!
                        north = (np.array(max_lat_coord-1), np.array(max_lon_coord))
                        north_east = (np.array(max_lat_coord-1), np.array(max_lon_coord+1))
                        east = (np.array(max_lat_coord), np.array(max_lon_coord+1))
                        south_east = (np.array(max_lat_coord+1), np.array(max_lon_coord+1))
                        south = (np.array(max_lat_coord+1), np.array(max_lon_coord))
                                  
                        upstream_direction = ivt_direction_currenttime_currentobject.data[max_lat_coord, max_lon_coord] % 360 #note difference to previous direction - now we're looking upstream to the east
                        
                        if upstream_direction > 337.5 or upstream_direction <= 22.5:
                            candidate_cells = [north, north, north_east] #note only two unique candidate cells here, as we have no western direction, but make a list of 3 so the code to check direction below works
                        elif 22.5 < upstream_direction <= 67.5:
                            candidate_cells = [north, north_east, east]
                        elif 67.5 < upstream_direction <= 112.5:
                            candidate_cells = [north_east, east, south_east]
                        elif 112.5 < upstream_direction <= 157.5 :
                            candidate_cells = [east, south_east, south]
                        elif 157.5 < upstream_direction <= 202.5:
                            candidate_cells = [south_east, south_east, south] #again, only 2 unique candidates
                        else: #terminate axis here
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                           
                        #If the first candidate (i.e. adjacent cell with the greatest IVT) has a direction opposite that of the direction of the first center cell, then go to the adjacent cell with the second greatest magnitude and check; if that fails, then go to the third cell; if that fails, then give up and go to the cell with largest magnitude, regardless of direction
                        candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                        candidate_cells_direction = [ivt_direction_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells] #note these are directions the wind is blowing FROM (different than the previous iteration)
                        index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                        
                        if not any(candidate_cells_ivt): #reached the edge of the object - all adjacent cells in the direction we care about are 0 and we can move on
                            foo = 1 #dummy condition to move algorithm along
                        elif 0 < candidate_cells_direction[index] < 180: #cell with max ivt is good, set as new center and proceed
                            ivt_magnitude_max_coords = candidate_cells[index]
                            max_lat_coord = ivt_magnitude_max_coords[0]
                            max_lon_coord = ivt_magnitude_max_coords[1]
                            cell_record.append(0)
                        else: #move to next greatest cell
                            candidate_cells_ivt[index] = 0
                            index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                            if 0 < candidate_cells_direction[index] < 180: #cell with second-largest ivt is good, set as new center and proceed
                                ivt_magnitude_max_coords = candidate_cells[index]
                                max_lat_coord = ivt_magnitude_max_coords[0]
                                max_lon_coord = ivt_magnitude_max_coords[1]
                                cell_record.append(0)
                            else: #move to third-greatest cell
                                candidate_cells_ivt[index] = 0
                                index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                if 0 < candidate_cells_direction[index] < 180: #cell with third-largest ivt is good, set as new center and proceed
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(0)
                                else: #terminate, as something is wrong
                                    candidate_cells_ivt = [ivt_magnitude_currenttime_currentobject.data[x[0], x[1]] for x in candidate_cells]
                                    index, value = max(enumerate(candidate_cells_ivt), key=operator.itemgetter(1))
                                    ivt_magnitude_max_coords = candidate_cells[index]
                                    max_lat_coord = ivt_magnitude_max_coords[0]
                                    max_lon_coord = ivt_magnitude_max_coords[1]
                                    cell_record.append(1)
                
                        if np.sum(cell_record[-num_consecutive_bad_steps_allowed:]) == num_consecutive_bad_steps_allowed: #If the algorithm had to take X bad steps in a row, terminate this direction of axis calculation
                            ivt_magnitude_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 0
                            
                        else:
                            axis_array_currenttime_currentobject.data[max_lat_coord, max_lon_coord] = 1
                            axis_coords_list_currenttime_currentobject.append(ivt_magnitude_max_coords)
                            
                    else: #point on the boundary, cull it
                         ivt_magnitude_currenttime_currentobject[max_lat_coord, max_lon_coord] = 0     
        
            object_axes_list_currenttime.append(axis_array_currenttime_currentobject)
            axis_coords_list_currenttime.append(axis_coords_list_currenttime_currentobject)                
        
        else: #bad object, give back null data for axis locations so that the list this function output appends to has the correct number of objects for each time (each time should have an entry for each object at the time, even for the bad objects which are represented by [])
            object_axes_list_currenttime.append([])    
            axis_coords_list_currenttime.append([]) 
           
    return object_axes_list_currenttime, axis_coords_list_currenttime
                       


def calcProperAxis(axis_coords_currentobject):
    
    sorted_axis_coords = []
    
    axis_coords_currentobject = axis_coords_currentobject[axis_coords_currentobject[:,1].argsort()] #sorts from west to east, which is the direction we'll be going in
    axis_coords_currentobject_listform = [list(coordinate) for coordinate in axis_coords_currentobject] #operations below only work in list form
    
    #find the origin cell, i.e. the westernmost cell with only one neighbor. Once it's found, append it to the sorted list, delete it from the original list, and go to algorithm below
    # can't just rerun the neighbor search algorithm because cardinal direction + diagonal is valid but results in a neighbor count of 2
    #however, this fails in edge cases where the start has 2 cells adjacent! 
    try:
        neighbor_count = 999
        current_index = 0
        while neighbor_count > 1:
            neighbor_count=0
            
            current_cell = axis_coords_currentobject_listform[current_index]
            
            if [current_cell[0]-1, current_cell[1]-1] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0], current_cell[1]-1] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0]-1, current_cell[1]] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0]-1, current_cell[1]+1] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0]+1, current_cell[1]-1] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0]+1, current_cell[1]] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0]+1, current_cell[1]+1] in axis_coords_currentobject_listform:
                neighbor_count += 1
            if [current_cell[0], current_cell[1]+1] in axis_coords_currentobject_listform:
                neighbor_count += 1

            current_index += 1
                
        origin_coord = axis_coords_currentobject_listform.pop(current_index-1)
        
    except: #select the westernmost cell. This has the possibility to be flawed by a cell or two, but in the context of length calculation this is practically never relevant - still a design flaw that should be fixed! 
    #idea for fixing: if above fails, go back to first cell and look only in cardinal directions. If there's a neighbor there, it's the next cell; if none, then the diagonal is the cell (though I don't see how this could happen - the cardinal direction search should suffice...)
        origin_coord = axis_coords_currentobject_listform.pop(0) 
    
    sorted_axis_coords.append(origin_coord)
    axis_coords_currentobject = np.array([np.array(coord) for coord in axis_coords_currentobject_listform])
    
    #now compute list of L1 norms from origin cell, select the minimum, make that the origin, repeat until list is 1 entry long
    while len(axis_coords_currentobject_listform) > 1:
    
        x = np.abs(axis_coords_currentobject - np.array(origin_coord))
        l1_norms = [np.sum(y) for y in x]
        origin_coord = axis_coords_currentobject_listform.pop(np.argmin(l1_norms))
        sorted_axis_coords.append(origin_coord)
        axis_coords_currentobject = np.array([np.array(coord) for coord in axis_coords_currentobject_listform])
        
    return np.array([np.array(coord) for coord in sorted_axis_coords])                



def printPercentageComplete(i, num_objs_list):
    if i == int(len(num_objs_list)/4):
        print('25% done')
    elif i == int(len(num_objs_list)/2):
        print('50% done')
    elif i == int(3*len(num_objs_list)/4):
        print('75% done')


    