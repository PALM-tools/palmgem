db=name_of_database
schema=name_of_schema
owner=palm
base_path=path_2_working_dir
log_file=${base_path}/import_process.log
export PGUSER=your_pg_user
export PGPASSWORD=your_pg_password
export PGHOST=localhost
export PGPORT=5432
export PGDATABASE=$db
export QUIET=1

echo "Importing data to database: ${db}, schema: ${schema}, under schema owner: ${owner}" > ${log_file} 2>&1

psql -d $db -c "create schema if not exists \"$schema\"" >> ${log_file} 2>&1
psql -d $db -c "alter schema \"$schema\" owner to ${owner}" >> ${log_file} 2>&1

# vector layers
vshapes=('name_of_UrbanAtlas_shapefile.shp' 'name_of_OpenStreetMap_shapefile.shp')
vlayers=('imported_landcover' 'imported_streetmap')

# raster layers
rshapes=('Name_of_BuildingsHeights_raster.tif' 'name_of_EU_DEM_raster.tif')
rlayers=('imported_buildings' 'imported_dem')

# copy vector layers
for i in ${!vlayers[*]}; do
  echo $i, ${vshapes[i]}, ${vlayers[i]} >> ${log_file} 2>&1

  sql="drop table if exists \"$schema\".\"${vlayers[i]}\" cascade"
  psql -d $db -c "$sql" >> ${log_file} 2>&1

  ogr2ogr -nln ${vlayers[i]} -nlt PROMOTE_TO_MULTI -lco GEOMETRY_NAME=geom -lco SCHEMA=$schema -lco FID=gid -lco PRECISION=NO -overwrite Pg:"dbname=$db host=$PGHOST user=$PGUSER port=$PGPORT" $base_path/${vshapes[i]} >> ${log_file} 2>&1
  psql -d $db -c "alter table \"$schema\".\"${vlayers[i]}\" owner to ${owner}" >> ${log_file} 2>&1

done

# copy raster layers
for i in ${!rlayers[*]}; do
  echo $i, ${rshapes[i]}, ${rlayers[i]} >> ${log_file} 2>&1

  sql="drop table if exists \"$schema\".\"${rlayers[i]}\" cascade"
  psql -d $db -c "$sql" >> ${log_file} 2>&1

  raster2pgsql -I -C -M -t auto $base_path/${rshapes[i]} -q $schema.${rlayers[i]} | psql -q -b -d $db >> ${log_file} 2>&1
  psql -d $db -c "alter table \"$schema\".\"${rlayers[i]}\" owner to ${owner}" >> ${log_file} 2>&1

done