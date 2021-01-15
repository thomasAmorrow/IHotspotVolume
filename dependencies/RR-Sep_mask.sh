#!/bin/bash

# RR-Sep | regional-residual separation for hotspot tracks
# 	2017, July 25

# Thomas Morrow, University of Idaho, Moscow, Idaho, USA
# contact: morr4998@vandals.uidaho.edu

# This script separates regional and residual bathymetry
#	(typically along hotspot tracks, but maybe useful in other locales)

#	1) Using the DiM filter, area of interest is median filtered, then smoothed.
#	2) Residual (OBS-FILT) volume and area are determined, mean amplitude calculated
#	3) Script selects maximized mean amplitude filter length (ORS)
#	4) Optimized residual and regional bathymetry output as grid files

# Script modified from Kim and Wessel's dim.template.sh in the GMT packages
# 		see Kim, S.-S., and Wessel, P. (2008), "Directional Median Filtering
# 			for Regional-Residual Separation of Bathymetry, Geochem. Geophys. Geosyst.,
# 			9(Q03005), doi:10.1029/2007GC001850.

# script requires GMT functions

# ---- In/out ---- #

echo "Inputs:"
echo "  Bathy file $1"
echo "  Region     $2"
echo
echo "Filter parameters:"
echo "  Min width  $3 km"
echo "  Max width  $4 km"
echo "  Step       $5 km"
echo "  Level      $6 m"

# 0. System defaults (Change if you know what you are doing):
dim_dist=2		# How we compute distances on the grid [Flat Earth approximation]
dim_sectors=8		# Number of sectors to use [8]
dim_filter=m		# Primary filter [m for median]
dim_quantity=l		# Secondary filter [l for lower]
dim_smooth_type=m	# Smoothing filter type [m for median]
dim_smooth_width=50	# Smoothing filter width, in km [50]

#    To prevent edge effects, the input grid domain must be
#    larger than the area of interest.
box=$2 #-R-159/-121/41/59	# Area of interest, a subset of data domain

# Input/output files
bathy=$1 #COBB.masked.grd #neg.submask.grav.grd #HAW_EMP.MASKED.grd # Input bathymetry grid file for the entire data domain
ors=ORStable.txt # Intermediate Optimal Robust Separator analysis results (table)
orsout=ORSout # ORS output work folder
dim=RR_output.grd # Final output DiM-based regional grid
#err=RR_output_err.grd # Final output DiM-based MAD uncertainty grid

gmtset PROJ_ELLIPSOID Sphere

# create a temporary working directory
if [ ! -d ORStmp ]; then

	mkdir ORStmp

fi

# ---- ORS analysis ---- #

if [ ! -f $ors ]; then

	mkdir -p $orsout

fi
	#gmt grdsample $bathy -I2m -G$bathy
	gmt grdsample $bathy $box -GORStmp/trim.t.grd -V # the area of interest

	#trim mask
	gmt grdsample MASK.grd -RORStmp/trim.t.grd -GORStmp/MASK.grd -V # the area of interest


	# A.1. Set filter parameters for an equidistant set of filters:
	minW=$3 # Minimum filter width candidate for ORS  (e.g., 60) in km
	maxW=$4 # Maximum filter width candidate for ORS  (e.g., 600) in km
	intW=$5 # Filter width step (e.g., 20) in km
	#maxcont=500 # maximum contour height for residual volume/area
	level=$6 #300 # step for base contour calculations

			
	#------stop A.1. editing here--------------------------------------

	#STEP=seq $minw $intW $maxW
	#STEP=`gmtmath -T$minW/$maxW/$intW -N1/0 =`

	for width in $(seq "$minW" "$intW" "$maxW")
	#for (( width = $minW; width <= $maxW; width += $intW ))
	#while [ $width -lt $maxW ];
	do
		echo "W = $width km"
		gmt dimfilter $bathy $box -G./ORStmp/${width}.dim.grd -F${dim_filter}${width} -D${dim_dist} -N${dim_quantity}${dim_sectors} # DiM filter
		
		gmt grdfilter ./ORStmp/${width}.dim.grd $box -G$orsout/dim.${width}.grd -F${dim_smooth_type}${dim_smooth_width} -D${dim_dist} # smoothing - TM added $box

		gmt grdmath ORStmp/trim.t.grd $orsout/dim.${width}.grd SUB ORStmp/MASK.grd MUL = ./ORStmp/${width}.resid.grd # residual from DiM
		gmt grdvolume ORStmp/${width}.resid.grd -Sk -C$level -Vl | awk '{print r,$1,$2,$3,$4}' r=${width} >> $ors  # ORS from DiM
		#gmt grdvolume ./tmp/$$.sd.grd -Sk -C0/200/$level -T -Vl #| awk '{print r,$2,$3,$4}' r=${width} >> $ors  # ORS from DiM

			# ORS table output format	(see Wessel 1998 eqn [4,5,6])
			#  filt W | area (A) | volume (V) | max mean height (h=V/A) 

		#width=`expr $width + $intW`


	done

echo "~~DONE!~~"

