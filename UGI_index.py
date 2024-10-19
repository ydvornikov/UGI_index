def m_index(polygons,id_field):
    
    def subset(polygons):
        centroids = polygons.centroid
        polygons_sub = gpd.sjoin(polygons.set_geometry(centroids), town, how='inner')[list(polygons.columns)]
        polygons_sub[polygons.geometry.name] = polygons[polygons.geometry.name]
        
        polygons = polygons_sub.dissolve(by=id_field).reset_index(drop=False)
        polygons['area'] = round(polygons.area, 2)
        
        return polygons
        
    def access(population):
        for i in [2,5,15,20]:
            iso = accessibility[accessibility.time == i].union(accessibility[accessibility.time == i])
            population['acc_'+str(i)] = population[population.geometry.name].map(lambda x: 
                                                                                 1 if iso.contains(x).any() == True else 0)
            population['acc'+str(i)+'_index'] = population['acc_'+str(i)] * population['inhabitant']
        population['acc_index'] = population['acc_2'] + population['acc_5'] + population['acc_15'] + population['acc_20']
        population['pop_index'] = population['inhabitant']*population['acc_index']
        
        return population
    
    def point_stat(polygons, points, point_field):
        polygons_stat = pd.merge(polygons, 
                                 pd.DataFrame(
                                     gpd.sjoin(points,
                                               polygons,
                                               predicate ='within').groupby(id_field)[point_field].sum()),
                                 on = id_field, how ='left')
        polygons_stat[point_field] = polygons_stat[point_field].fillna(0)
        return polygons_stat

    landcover = rasterio.open(ras_path)
    lc = landcover.read(1)
    lc_names = ['nodata','water','sealed','soils','forest','shrubs','lawns']

    def rast_stat_sl_green(polygons):
        for i in [2,4,5,6]:
            polygons[lc_names[i]+"_perc"] = polygons.geometry.apply(lambda x :
                                                                    pd.DataFrame(
                                                                        rasterstats.zonal_stats(x,
                                                                                                np.where(lc == i,1,0),
                                                                                                affine = landcover.transform,
                                                                                                nodata = 0,
                                                                                                stats = "sum"))['sum']*100)
            polygons[lc_names[i]+"_perc"] = round((polygons[lc_names[i]+"_perc"].fillna(0) * 100) / polygons['area'], 2)
        return polygons

    def zone_extract(zones, polygons):
        zn = gpd.GeoDataFrame()
        for x in polygons[id_field].unique():
            sub_pol = polygons[polygons[id_field] == x]
            zone = gpd.clip(zones, sub_pol)
            if len(zone) == 1:
                zone[id_field] = sub_pol[id_field].values[0]
                zn = pd.concat([zn, zone]).reset_index(drop = True)
        zn['zone_area'] = zn.area
        return zn
    
    def green_metrics(zones, polygons, green_stands, green_areas):
        polygons = pd.merge(polygons,
                            zones[[id_field,'zone_area']],
                            on = id_field,
                            how ='left')
        stands = gstands.sjoin(zones[['geometry',id_field]],
                               how="left")
        stands = stands[~stands[id_field].isnull()]
        stands.geometry = stands.geometry.buffer(3.3)

        zones_area = zone_extract(green_areas, zones[['geometry',id_field]])
        g_sp_san = pd.concat([zones_area[[geom_field, id_field]], stands[[geom_field, id_field]]]).dissolve(
            by = id_field).reset_index(drop=False)
        g_sp_san['canopy_area'] = round(g_sp_san.area,2)
        polygons['canopy_area'] = pd.merge(polygons,
                                           g_sp_san[[id_field, 'canopy_area']],
                                           on = id_field,
                                           how ='left')['canopy_area'].fillna(0)
        polygons['dni_full'] = round(polygons['canopy_area'] * 100 / (polygons['zone_area']), 2)
        polygons = polygons.rename(columns={'zone_area': 'sanitary_area'})
        polygons['sanitary_area'] = round(polygons['sanitary_area'], 2).fillna(0)
        
        return polygons
    
    def green_rast(roads, polygons):
        zones = zone_extract(roads, polygons)[['geometry','zone_area',id_field]]
        zones['green_road_area'] = pd.DataFrame(rasterstats.zonal_stats(zones,
                                                                        np.where(lc >= 4, 1, 0),
                                                                        nodata = 0,
                                                                        affine = landcover.transform,
                                                                        stats = 'sum'))['sum']*100
        zones['gri'] = round(((zones['green_road_area']*100) / (zones['zone_area']/2)), 2).fillna(0)
        zones['gri'] = zones.apply(lambda x: 100 if x['gri']>100 else x['gri'], axis=1) #for grid values exceed 100
        zones = zones.rename(columns={'zone_area': 'road_area'})
        polygons = pd.merge(polygons,
                            zones[[id_field, 'road_area', 'gri']],
                            on = id_field,
                            how ='left')
        polygons['road_area'] = round(polygons['road_area'].fillna(0),2)
        polygons['sanitary_area_corrected'] = polygons['sanitary_area'] - (polygons['road_area']/2)
        polygons['dni'] = round(polygons['canopy_area']*100 / polygons['sanitary_area_corrected'],
                                2)
        
        return polygons
    
    def hydrology(polygons, raster1, raster2, raster3):
        polygons['RLT_m'] = round(pd.DataFrame(
            rasterstats.zonal_stats(polygons, raster1.read(1),
                                    nodata = -9999,
                                    affine = raster1.transform,
                                    all_touched = False,
                                    stats = 'sum'))['sum'], 2)
        polygons['RLR_m'] = round(pd.DataFrame(
            rasterstats.zonal_stats(polygons, raster2.read(1),
                                    nodata = -9999,
                                    affine = raster2.transform,
                                    all_touched = False,
                                    stats = 'sum'))['sum'], 2)
        polygons['sri'] = round(100 - ((polygons['RLT_m']/polygons['RLR_m'])*100),2)
        polygons['water_abs'] = round(pd.DataFrame(
            rasterstats.zonal_stats(polygons, raster3.read(1),
                                    nodata = -3.4028234663852886e+38,
                                    affine = raster3.transform,
                                    all_touched = False,
                                    stats = 'sum'))['sum'] * 1000 / polygons['area'] * 10000, 2)
        polygons['delta_runoff'] = polygons['RLR_m'] - polygons['RLT_m']
        return polygons

    def cooling_stat(polygons, cool_raster):
        polygons['ci'] = round(pd.DataFrame(rasterstats.zonal_stats(polygons,
                                                                    np.where(
                                                                        cool_raster.read(1) > 0,
                                                                        cool_raster.read(1),
                                                                        np.nan),
                                                                    nodata = cool_raster.meta['nodata'],
                                                                    affine = cool_raster.transform,
                                                                    all_touched = False,
                                                                    stats = 'mean'))['mean'], 2)
        polygons['ci'] = polygons['ci'].fillna(0.1)
        return polygons

    polygons = point_stat(rast_stat_sl_green(subset(polygons)), access(pop),
                          ['inhabitant','pop_index','acc2_index',
                           'acc5_index','acc15_index','acc20_index'])
    polygons['green_perc'] = polygons['forest_perc']+polygons['shrubs_perc']+polygons['lawns_perc']
    for i in [2,5,15,20]:
        polygons['acc'+str(i)+'_index'] = round(polygons['acc'+str(i)+'_index']/polygons['inhabitant'], 3)
        polygons['gai'] = round(polygons['pop_index']/polygons['inhabitant'], 2)
    
    if id_field in ['name_ru','h3_id']:
        polygons_green = green_metrics(zone_extract(sanitary_zones, polygons), polygons, gstands, g_areas)
        polygons = hydrology(cooling_stat(green_rast(roads, polygons_green), cooling), water_real, water_repl, water_delta)
        polygons['green_perc'] = polygons.apply(
            lambda x: 100 if x['green_perc']>100 else x['green_perc'], axis=1) #for grid values exceed 100
        return polygons
    
    else:
        polygons = polygons[polygons.inhabitant > 0]
        polygons['gpi'] = round((polygons['green_perc']/100*polygons['area'])/polygons['inhabitant'], 2)
        polygons_green = green_metrics(zone_extract(sanitary_zones, polygons), polygons, gstands, g_areas)
        polygons_green = polygons_green.rename(columns={'dni_full': 'dni'})
        polygons = hydrology(cooling_stat(polygons_green, cooling), water_real, water_repl, water_delta)
        return polygons