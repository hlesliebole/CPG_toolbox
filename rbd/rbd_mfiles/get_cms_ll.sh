#!/bin/sh

eval `cat cms.dat`
cmsvar1=Lon
cmsvar2=Lat
#echo Spot: $id
echo ${!cmsvar1} ${!cmsvar2} 
echo ${!cmsvar1} ${!cmsvar2} > spot.dat

