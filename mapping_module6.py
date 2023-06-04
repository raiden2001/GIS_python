# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:35:39 2019

@author: SALVAIJ
"""

shp_2018Q3_path = 'G:/Jason Ma/2022Q3 HERE shapefiles/Corridors/Old corridors/Gatineau1_N.shp'
#shp_2018Q3_AADT_path = 'N:/HERE_GIS_CANADA/2018Q3/2018Q3 corridors/ALL_AADT.shp'
#global wu_oh
attributionstring = 'Data from HERE Technologies. Basemap by OpenStreetMaps.'

zoom_lims = [10000.,
             25000.,
             50000.,
             90000.,
             ]

fraction_shift = 0.01
zoom_shift = 0


def calculate_map_extent(gdf_city):
    import shapely.ops
    
    bounds = shapely.ops.cascaded_union(gdf_city.geometry).bounds
    dx = bounds[2]-bounds[0]
    dy = bounds[3]-bounds[1]
    #map_extent = max([dx,dy,zoom_lims[0]])   
    map_extent = max([dx,dy])   

    return map_extent

def calculate_latitude(gdf_city):
    import shapely.ops
    
    bounds = shapely.ops.cascaded_union(gdf_city.to_crs(epsg=4326).geometry).bounds
    
    latitude = (bounds[3]+bounds[1])/2
    
    return latitude


def calculate_zoom(map_extent):
    if map_extent <= zoom_lims[0]:
        zoom=14
    elif zoom_lims[0] < map_extent <= zoom_lims[1]:
        zoom=13
    elif zoom_lims[1] < map_extent <= zoom_lims[2]:
        zoom=12
    elif zoom_lims[2] < map_extent <= zoom_lims[3]:
        zoom=11
    elif zoom_lims[3] < map_extent:
        zoom=10
    
    return zoom+zoom_shift


    #zoom_lon = np.ceil(np.log2(360 * 2.0 / lon_length))
    #zoom_lat = np.ceil(np.log2(360 * 2.0 / lat_length))
    #zoom = np.max([zoom_lon, zoom_lat])




def part_roads_corridor(corridor,df_corr,map_extent):
    import shapely.geometry, shapely.affinity

    
    offset_scale=map_extent*fraction_shift
    
    start_end = df_corr.sort_values(by='POS_ALONG').iloc[[0,-1],:].reset_index(drop=True)

    start_point = start_end.loc[0,'geometry'].centroid
    end_point   = start_end.loc[1,'geometry'].centroid

    d = end_point.distance(start_point)

    unit_vec = shapely.geometry.LineString( [(0,0),
        ((end_point.x-start_point.x)/d,
         (end_point.y-start_point.y)/d)])

    theta = -90
        
    offset_vec = shapely.affinity.rotate(unit_vec, theta, origin=(0,0))
    offset_x,offset_y = list(offset_vec.coords[1])
    offset_x*=offset_scale
    offset_y*=offset_scale
    
    def translate_row(row):
        geometry = row.geometry
        translated = shapely.affinity.translate(geometry,xoff=offset_x, yoff=offset_y, zoff=0.0)
        return translated
    
    shifted_geometry = df_corr.apply(translate_row,axis=1)
    return shifted_geometry





def add_corr_num(ax,gdf_corr,map_extent,corr_num,legend_fontsize):
    import shapely.geometry
    import matplotlib.patches as mpatches
    from numpy import sqrt
    offset_scale=map_extent*fraction_shift*4
    N = len(gdf_corr)
    midway = int(N/2)
    
    mid_start = gdf_corr.sort_values(by='POS_ALONG').iloc[  midway].geometry.centroid
    mid_end   = gdf_corr.sort_values(by='POS_ALONG').iloc[midway+1].geometry.centroid

    mid_start = shapely.geometry.Point(mid_start)
    mid_end = shapely.geometry.Point(mid_end)

    d = mid_end.distance(mid_start)
    unit_vec_x,unit_vec_y = (mid_end.x-mid_start.x)/d, (mid_end.y-mid_start.y)/d
    unit_vec = shapely.geometry.LineString( [(0,0),(unit_vec_x,unit_vec_y)])

    theta = -90
    offset_vec = shapely.affinity.rotate(unit_vec, theta, origin=(0,0))
    offset_x,offset_y = list(offset_vec.coords[1])
    offset_x*=offset_scale
    offset_y*=offset_scale

    arrow_base_x,arrow_base_y = mid_start.x+offset_x,mid_start.y+offset_y
    arrow_dx,arrow_dy = offset_scale*unit_vec_x,offset_scale*unit_vec_y

    r = sqrt(arrow_dx**2 + arrow_dy**2)/2
    circle = mpatches.Circle((arrow_base_x,arrow_base_y),
                             r,
                             facecolor='white',
                             edgecolor='black',
                             linewidth=0.5, 
                            )
    ax.add_artist(circle)

    ax.text(arrow_base_x,arrow_base_y,
            corr_num,
            verticalalignment='center',
            horizontalalignment='center',
            fontsize=legend_fontsize-1
            )


def add_direction_arrow(ax,gdf_corr,map_extent,corridor):
    import shapely.geometry, shapely.affinity
    
    problem_corridors = ['Montreal3_E',
                         'Montreal3_W',
                         'Trois-Rivieres2_E',
                         'Trois-Rivieres2_W',
                         'Trois-Rivieres2_N',
                         'Trois-Rivieres2_S',
                             ]
    
    offset_scale=map_extent*fraction_shift*4
    N = len(gdf_corr)
    midway = int(N/2)
    
    mid_start = gdf_corr.sort_values(by='POS_ALONG').iloc[  midway].geometry.centroid
    mid_end   = gdf_corr.sort_values(by='POS_ALONG').iloc[midway+1].geometry.centroid
    
    #mid_start,mid_end = midrows.iloc[0].centroid,midrows.iloc[1].centroid

    mid_start = shapely.geometry.Point(mid_start)
    mid_end = shapely.geometry.Point(mid_end)

    d = mid_end.distance(mid_start)
    unit_vec_x,unit_vec_y = (mid_end.x-mid_start.x)/d, (mid_end.y-mid_start.y)/d
    unit_vec = shapely.geometry.LineString( [(0,0),(unit_vec_x,unit_vec_y)])

    if any([corridor==problem_corridor for problem_corridor in problem_corridors]):
        theta = +90
    else:
        theta = -90

    offset_vec = shapely.affinity.rotate(unit_vec, theta, origin=(0,0))
    offset_x,offset_y = list(offset_vec.coords[1])
    offset_x*=offset_scale
    offset_y*=offset_scale

    arrow_base_x,arrow_base_y = mid_start.x+offset_x,mid_start.y+offset_y
    arrow_dx,arrow_dy = offset_scale*unit_vec_x,offset_scale*unit_vec_y

    ax.arrow(arrow_base_x,arrow_base_y, arrow_dx,arrow_dy,
             width=offset_scale*.1,
             fc='k', ec='k')




def add_basemap(ax_or_axes,
                map_extent,
                latitude,
                zoom='auto',
                url='http://a.tile.openstreetmap.org/{z}/{x}/{y}.png',
                # another good one is 'http://tile.stamen.com/toner-lite/{z}/{x}/{y}.png'
                alpha=1,
                aspect='auto',
                scalebar=True,
               ):
    from contextily import bounds2img
    
    if type(ax_or_axes)==list:
        ax=ax_or_axes[0]
        axes = ax_or_axes
    else:
        ax=ax_or_axes
    
    xmin, xmax, ymin, ymax = ax.axis()
    
    centre_loc_x,centre_loc_y=(xmax+xmin)/2.,(ymax+ymin)/2.
    
    floored_map_extent = map_extent = max([map_extent,zoom_lims[0]])  
    
    if aspect == 'square':
        xmin,xmax = centre_loc_x - floored_map_extent/2.,centre_loc_x + floored_map_extent/2.
        ymin,ymax = centre_loc_y - floored_map_extent/2.,centre_loc_y + floored_map_extent/2.
    elif aspect == 'auto':
        pass
    else:
        print('aspect should be "square" or "auto". Using "auto".')
    
    if zoom=='auto':
        zoom = calculate_zoom(floored_map_extent)
    else:
        pass
    
    basemap, extent = bounds2img(xmin, ymin, xmax, ymax, zoom=zoom, url=url)

    if type(ax_or_axes)==list:
        for ax in axes:
            ax.imshow(basemap, extent=extent, interpolation='bilinear',alpha=alpha)
            
            # restore original x/y limits
            ax.axis((xmin, xmax, ymin, ymax))
    else:
        ax.imshow(basemap, extent=extent, interpolation='bilinear',alpha=alpha)
        
        # restore original x/y limits
        ax.axis((xmin, xmax, ymin, ymax))


def add_scalebar(ax,latitude):
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from numpy import cos,pi
    
    xmin, xmax, ymin, ymax = ax.axis()
    
    dx= xmax - xmin
    
    ground_meters = dx*cos(latitude*pi/180.)
    
    ground_kms = ground_meters/1000
    
    bar_kms = max([1, ground_kms//7])
    
    bar_length = bar_kms*1000/cos(latitude*pi/180.)
    
    #print(ground_meters,ground_kms,bar_length)
    
    bar_height = 100
    
    bar = AnchoredSizeBar(ax.transData, #transform
                          bar_length, #Horizontal length of the size bar, given in coordinates of transform.
                          '%d km'%(bar_kms), #Label to display.
                          4,#location
                          pad=0.5,#Padding around the label and size bar, in fraction of the font size. Defaults to 0.1.
                          sep=5, #Seperation between the label and the size bar, in points. Defaults to 2.
                          borderpad=0.5,#Border padding, in fraction of the font size. Defaults to 0.1.
                          frameon=True,#If True, draw a box around the horizontal bar and label. Defaults to True.
                          size_vertical=bar_length/13,#Vertical length of the size bar, given in coordinates of transform. Defaults to 0.
                          #color='black',#Defaults to black
                          #fontproperties=fontprops, #use FontManager
                         )
    ax.add_artist(bar)
    ''' Locations:
    'upper right'  : 1,
    'upper left'   : 2,
    'lower left'   : 3,
    'lower right'  : 4,
    'right'        : 5,
    'center left'  : 6,
    'center right' : 7,
    'lower center' : 8,
    'upper center' : 9,
    'center'       : 10
    '''




def add_cbar(vmin,vmax,cmap,fig,ax):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from pylab import cm, Normalize
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    sm = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cbar = fig.colorbar(sm, cax=cax)
    return cbar,cax



def other_node_type(node_type):
    if node_type == 'REF_IN_ID':
        return 'NREF_IN_ID'
    elif node_type == 'NREF_IN_ID':
        return 'REF_IN_ID'
    
    else:
        print('%s isnt a valid node type'%(node_type))
        raise ValueError('argument node_type must be either "REF_IN_ID" or "NREF_IN_ID". Instead receieved %s'%(node_type))
        
def pop_row(df,index):
    #MAKES THE CHANGE TO df IN PLACE AND RETURNS THE REMOVED ROW
    out_row = df.loc[index,:]
    df_row_removed = df.loc[df.index.difference([index]),:]
    return df_row_removed,out_row

def Order_the_Links(corridor_name,corr_df):
    import pandas as pd
    corridor_direction = corridor_name[-1]
    #print("corridors_direction:",corridor_direction) 
    node_types = ['REF_IN_ID', 'NREF_IN_ID']
   # print(node_types)
    df_corridor_unsorted = corr_df.copy()
    df_corridor_sorted = pd.DataFrame(columns = corr_df.columns.to_list()) #initialize the empty df, to which rows will be added
    #print("df_corridor_unsorted:",df_corridor_unsorted)
    #print("df_corridor_sorted:",df_corridor_sorted)
    #print(df_corridor_unsorted)
    N_rows = len(df_corridor_unsorted)
    #print("N_rows:",N_rows)
    #find the NODE IDs that are unique to both rows
    #in most cases, there will be 2: a start and an end point
    #these may be REF_IN_ID or NREF_IN_ID
    #in toronto because of express and collector, there are more than 2. deal with that later
    all_nodes = pd.concat([df_corridor_unsorted[node_type] for node_type in node_types],
                           keys=node_types
                         ).reset_index(level=0).rename(columns={'level_0':'NODE_COL',0:'NODE_ID'})
    
    start_end_nodes = all_nodes.drop_duplicates(subset='NODE_ID',keep=False)
    #print("start_end_nodes:",start_end_nodes)
    #print(all_nodes)
    #if start_end_nodes.empty:
        #print("start_end_nodes dataframe is empty. Please check the input data.")
        #return df_corridor_sorted
    #else:
    index_coord_A = start_end_nodes.index.values[0]
    index_coord_B = start_end_nodes.index.values[1]
    #print("index_coord_A",index_coord_A)
    #print("index_coord_B",index_coord_B)
        
    coord_A = list(df_corridor_unsorted.loc[index_coord_A,'geometry'].coords)[0]
    coord_B = list(df_corridor_unsorted.loc[index_coord_B,'geometry'].coords)[0]
    #print("coord_A",coord_A)
    #print("coord_B",coord_B)
    
    if corridor_direction == 'N':
        if coord_A[1]-coord_B[1] > 0:
            start_end_iter = 1
        else:
            start_end_iter = 0
    elif corridor_direction == 'E':
        if coord_A[0]-coord_B[0] > 0:
            start_end_iter = 1
        else:
            start_end_iter = 0
    elif corridor_direction == 'S':
        if coord_A[1]-coord_B[1] > 0:
            start_end_iter = 0
        else:
            start_end_iter = 1        
    elif corridor_direction == 'W':
        if coord_A[0]-coord_B[0] > 0:
            start_end_iter = 0
        else:
            start_end_iter = 1
           
    else:
        wuh_oh
         
    #print("coord_A:", coord_A)
    #print("coord_B:", coord_B)
    if len(start_end_nodes) == 2:
        
        unsorted_index = start_end_nodes.index.values[start_end_iter]
        node_in_type,node_in_id  = start_end_nodes.values[start_end_iter]
        '''
        print("start_end_nodes1",start_end_nodes)
        print("node_in_type1:",node_in_type)
        print("node_in_id:",node_in_id)
        print("unsorted_index :",unsorted_index )
        print("start_end_iter",start_end_iter)
        '''
        last_node_type,last_node_id  = start_end_nodes.values[(start_end_iter+1)%2]
        #print("last_node_type",last_node_type)
        #print("last_node_id",last_node_id)
    #elif len(start_end_nodes) >= 2:
       # print("eroor: could not determine the start and end nodes")
        #raise Exception("could not determine the start end nodes")
        
         
        for i in range(0,N_rows):
            #print(i,len(df_corridor_sorted),len(df_corridor_unsorted))
            
            
            current_link_endpoints = [list(df_corridor_unsorted.loc[unsorted_index,'geometry'].coords)[ 0],
                          list(df_corridor_unsorted.loc[unsorted_index,'geometry'].coords)[-1]
                       ]
            #print("current_link_endpoints",current_link_endpoints)
            
            #in one line, remove the row from UNSORTED and add it to SORTED
            df_corridor_unsorted,df_corridor_sorted.loc[i] = pop_row(df_corridor_unsorted,unsorted_index)
            #print("df_corridor_unsorted:",df_corridor_unsorted)
            #print("unsorted_index:",unsorted_index)
            node_out_type = other_node_type(node_in_type)
            #print("node_out_type2:",node_out_type)
            node_out_id = df_corridor_sorted.loc[i,node_out_type]
            #print("node_out_id:",node_out_id)
            
            #print("df_corridor_unsorted:", df_corridor_unsorted)
            #print("df_corridor_sorted",df_corridor_sorted.loc[i])
            #print("node_out_type:",node_out_type)
            #Tries to call the function and pass in the correct argument
            node_type = node_out_type
            other_node = other_node_type(node_type)
            #print("node_type1.1",node_type)
            #print("other_node1",other_node) # outputs the results
            #if node_out_type in ['REF_IN_ID','NREF_IN_ID']:
                
              
            #else:
                #raise ValueError("node_out_type must be either 'REF_IN_ID' or 'NREF_IN_ID, not%s"%node_out_type)
            #print(node_out_type)
            #print(node_out_id)

            
            if i == N_rows-1:
                #end of algorithm
                is_sorted = last_node_id == node_out_id
                print (corridor_name,is_sorted) #they should match
                
                #print("df_corridor_sorted1",df_corridor_sorted)
                if is_sorted:
                    print("sorting is successful")
                    df_corridor_sorted.index+=1
                    df_corridor_sorted = df_corridor_sorted.reset_index().rename(columns={'index':'POS_ALONG_CORR'})
                    return (1,"sorting successful", df_corridor_sorted)
                else:
                    print("It is not sorted correctly, based on the result of shapefiles")
                    print("sorting unsuccessful")
                    return(0,"sorting not successful due to incorrect shaepfile",None)
              
            else:
                #update list of nodes
                all_remaining_nodes = pd.concat([df_corridor_unsorted[node_type] for node_type in node_types],
                                                   keys=node_types
                                                 ).reset_index(level=0
                                                 ).rename(columns={'level_0':'NODE_COL',0:'NODE_ID'})            
                ### update variables for next iteration of loop
                #filtered_nodes = all_remaining_nodes[all_remaining_nodes['NODE_ID'] == node_out_id]
                #if not filtered_nodes.empty:
                  #  filtered_nodes = filtered_nodes.reset_index(drop=True)
                  #  unsorted_index = filtered_nodes.index.values[0]
                  #  node_in_type,node_in_id = filtered_nodes.values[0]
               # else:
               #     print("No matching NODE_ID found for node_out_id:{}, skipping...".format(node_out_id))
               #     continue
                unsorted_index = all_remaining_nodes[all_remaining_nodes['NODE_ID']==node_out_id].index.values[0]
                node_in_type,node_in_id = all_remaining_nodes[all_remaining_nodes['NODE_ID']==node_out_id].values[0]
                #print(df_corridor_sorted)
                '''
                print("all_remaining_nodes",all_remaining_nodes)
                print("unsorted_index", unsorted_index)
                print("node_in_type2", node_in_type) ###### problem with getting wrong value
                print("node_in_id",node_in_id)
                '''
                next_link_endpoints = [list(df_corridor_unsorted.loc[unsorted_index,'geometry'].coords)[ 0],
                                      list(df_corridor_unsorted.loc[unsorted_index,'geometry'].coords)[-1]
                                      ]
                
                if df_corridor_sorted.loc[i,'DIR_TRAVEL']=='B':
                    for j in [0,1]:
                        for k in [0,1]:
                            if current_link_endpoints[j]==next_link_endpoints[k]:
                                diffs = [a-b for (a,b) in zip(current_link_endpoints[j],current_link_endpoints[(j+1)%2])]
                                if diffs[1] > 0:
                                    dir_travel = 'F'
                                elif diffs[1] < 0:
                                    dir_travel = 'T'
                                elif diffs[1] == 0:
                                    if diffs[0] > 0:
                                        dir_travel = 'F'
                                    elif diffs[0] > 0:
                                        dir_travel = 'T'
                
                    df_corridor_sorted.loc[i,'DIR_TRAVEL']=dir_travel
                    df_corridor_sorted.loc[i,'LINK_DIR']='%d%s'%(df_corridor_sorted.loc[i,'LINK_ID'],dir_travel)


            
    else:
        print(corridor_name,'unsorted. Found',len(start_end_nodes),"start/end candidates. Should be only 2")
        df_corridor_unsorted['POS_ALONG_CORR']=pd.Series([0 for i in range(0,len(df_corridor_unsorted))])
        
        cols = df_corridor_unsorted.columns.tolist()
        df_corridor_unsorted = df_corridor_unsorted[[cols[-1],] + cols[0:-1]]
        
        return (0,"sorting not successful due to df_corridor_unsorted",df_corridor_unsorted)
        #print("sorry it not sorted")


def quick_map_show(gdf,column=None):
    import geopandas as gpd
    import pylab as P
    fig=P.figure(figsize=(10,10))
    ax=fig.add_subplot(111)
    
    if column:
        gdf.plot(ax=ax,
                 column=column
                )
    else:
        gdf.plot(ax=ax)
    ax.set_axis_off()
    return fig,ax
    
