#!/bin/bash

## Documentation

# wget_GMRT_corners queries GMRT servers and downloads bathymetry within a region defined by a NE and SW corner (lat/long)

# corners.xyz input should be of the form: 
#	name, NE corner long, NE corner lat, SW corner long, SW corner lat
#	name ...
# 								and so on

# wget_GMRT_corners will download bathymetry within rectangle defined by corners.xyz, in format $format at resolution $res and 
# save the file in the specified directory $save_direct as filename $name from the corners.xyz file
# GMRT offers masked or unmasked data with GEBCO08 fill between high resolution data

## ***** record *****
# 30 oct 14 | first draft | T.Morrow
#
# ***** ***** *****

## Code

# inputs

# corner file (csv formatted as above)
corners="CORNERS.xy"

# desired resolution (low/default, med, high, max)
res="high"

# output file format (netcdf, esriascii, geotiff)
format="netcdf"

# desired layer (topo, topo-mask)
layer="topo"

#filename construction
PRE=""
POST=".grd"

#read in file to variables
while read -r line; do
    FILENAME=$(echo $line | awk -F, '{print $1;}')
    NElong=$(echo $line | awk -F, '{print $2;}')
    NElat=$(echo $line | awk -F, '{print $3;}')
    SWlong=$(echo $line | awk -F, '{print $4;}')
    SWlat=$(echo $line | awk -F, '{print $5;}')

#build URL
    URL1="https://www.gmrt.org/services/GridServer?minlongitude=" #"http://www.marine-geo.org:80/tools/GridServer?minlongitude="  
    URL2="&maxlongitude="
    URL3="&minlatitude="
    URL4="&maxlatitude="
    URL5="&format="
    URL6="&resolution="
    URL7="&layer="

#assemble URL
WGETURL="$URL1$SWlong$URL2$NElong$URL3$SWlat$URL4$NElat$URL5$format$URL6$res$URL7$layer"

#wget URL
echo $FILENAME

# run wget operation
wget --progress=bar:force --no-check-certificate -O $PRE$FILENAME$POST $WGETURL
    
done < "$corners"
