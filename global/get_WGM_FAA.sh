#!/bin/bash
​
## Documentation
# retrieves the WGM-2012 FAA data and converts for +/-180 and 0-360 longitude

wget https://cdn.obs-mip.fr/bgi/wgm_grid/data/WGM2012_Freeair_ponc_2min.grd

gmt grdsample WGM2012_Freeair_ponc_2min.grd -Rg -GWGM2012_Freeair_ponc_2min_360.grd
