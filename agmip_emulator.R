agmip_emulator <- function(stnlat, stnlon, tasmax = NULL, tasmin = NULL, pr = NULL, year = NULL, chirps=F) {

##############################################################################
# agmip_emulator
# 
# SUMMARY: This script takes three vectors of daily climate input variables 
# for a single pixel and year (maximum temperature, minimum temperature and 
# precipitation) and converts them into the 40 possible input variables for 
# the AgMIP-GGCMI seasonal emulator. The output is a single value of emulated 
# rainfed maize yield anomaly for that pixel and growing season/year. User can
# either input climate data from an external source, or use any year of data 
# from a combined CHIRTS / CHIRPS dataset from 1983-2016.
#
# INPUTS: latitude in degrees N, longitude in degrees E, vector of daily maximum
# temperatures in degrees C for year of interest (optional), vector of daily 
# minimum temperatures in degrees C for year of interest (optional), vector of daily 
# precipitation rates in mm/day for year of interest (optional). If tasmax, tasmin
# and pr are not provided, need to provide a year of interest between 1983 and 2016
# and set chirps=T to use CHIRTS/CHIRPS dataset.
#
# OUTPUTS: emulated rainfed maize yield anomaly for year of interest from 
# 1980-2010 baseline
#
# author:   Meridel Phillips
#           mmp2192@columbia.edu
# created:  August 31, 2021
# 
#
# DEBUG USING INPUT DATA:
#   stnlat = 5.7693
#   stnlon = 34.3989
#   tasmax = sample(-3:45, 365, replace=T)
#   tasmin = sample(-15:30, 365, replace=T)
#   pr = sample(0:15, 365, replace=T)
#   
# DEBUG USING CHIRTS/CHIRPS:
#   stnlat = 5.7693
#   stnlon = 34.3989
#   year = 1983
#   chirps = T
#
#
##############################################################################

#### required libraries

library(R.matlab)
library(RANN)

chirps = as.logical(chirps)

#### make sure either chirps flag = T or climate input is provided

if ((is.null(tasmax)) && (is.null(tasmin)) && (is.null(pr)) && (chirps==F)){
	stop("ERROR: No climate input data. Please provide climate input variables or set chirps=T to use CHIRTS/CHIRPS data.")
}

#### make sure a good year is selected for CHIRTS/CHIRPS

if (chirps==T) {
	if (is.null(year)){
		stop("ERROR: No year provided for CHIRTS/CHIRPS.")
	}
    if ((year<1983) || (year>2016)){
	    stop("ERROR: Year is outside the available 1983-2016 range for CHIRTS/CHIRPS.")
    }
    if ((!is.null(tasmax)) || (!is.null(tasmin)) || (!is.null(pr))) {
    	stop("ERROR: Selected to use CHIRTS/CHIRPS data but one or more climate variables is also provided. Please choose to input three climate variables (tasmax, tasmin and pr) OR set chirps=T, not both.")
    }
}

#### import global emulator variable indices and coefficients

globalcoeffs = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/rfMaize_global_coefficients.mat")
globalcoeffs = globalcoeffs$coefficients
globalvars = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/rfMaize_global_variables.mat")
globalvars = globalvars$variables

#### import latitude and longitude data

lat = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/AgGRIDlat.mat")
lat = lat$lat
lon = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/AgGRIDlon.mat")
lon = lon$lon

### import growing season dates

globalpday = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/rfMaize_pday.mat")
globalpday = globalpday$pday
globalgslength = readMat("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/Emulator_Parameters/rfMaize_gslength.mat")
globalgslength = globalgslength$gslength

#### find indices of stnlat and stnlon

if (stnlon>180){
	stnlon=stnlon-180
}

lats=lat[,1]
lons=lon[1,]
jj = nn2(lats,stnlat)
jj = jj$nn.idx[1]
ii = nn2(lons,stnlon)
ii = ii$nn.idx[1]

#### make sure it's a land pixel

if (is.na(globalpday[jj,ii])){
	stop("ERROR: This location is over the ocean in source data. Please input a land-based location.")
}

#### isolate emulator input variables and coefficients for this pixel

coeffs = globalcoeffs[jj,ii,]
vars = globalvars[jj,ii,]

### isolate growing season dates for this pixel

pday = round(globalpday[jj,ii])
gslength = round(globalgslength[jj,ii])
hday = pday+gslength
if (hday>365){
	hday=hday-365
}

if (hday<pday){
	seasonidx = c(pday:365, 1:hday)
}else{
	seasonidx = c(pday:hday)
}

#### calculate planting window and anthesis time periods 

# planting window: 11-day window around planting date
plantingidx = c((pday-5):(pday+5))
plantingidx[plantingidx>365] = plantingidx[plantingidx>365]-365
# before anthesis: planting date to beginning of anthesis
beforeidx = c(pday:(round(pday+(gslength/2))-5))
beforeidx[beforeidx>365] = beforeidx[beforeidx>365]-365
# during anthesis: 11-day window around middle of growing season
duringidx = c((round(pday+(gslength/2))-5):(round(pday+(gslength/2))+5))
duringidx[duringidx>365] = duringidx[duringidx>365]-365
# after anthesis: end of anthesis to end of growing season
afteridx = c((round(pday+(gslength/2))+5):(round(pday+gslength)))
afteridx[afteridx>365] = afteridx[afteridx>365]-365

#### calculate input variables from tasmax, tasmin, pr

if (chirps==T) {

   globaltasmax = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/CHIRTS_CHIRPS/chirts_05deg_daily_tasmax_", year, ".mat", sep=""))
   globaltasmax = globaltasmax$chirtstmax
   globaltasmin = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/CHIRTS_CHIRPS/chirts_05deg_daily_tasmin_", year, ".mat", sep=""))
   globaltasmin = globaltasmin$chirtstmin
   globalpr = readMat(paste("https://portal.nccs.nasa.gov/datashare/GISS/Impacts/AgMIP/WorldModelers/CHIRTS_CHIRPS/chirps_05deg_daily_pr_", year, ".mat", sep=""))
   globalpr = globalpr$chirpspr

   #### isolate pixel from CHIRTS/CHIRPS data
   tasmax = globaltasmax[jj,ii,]
   tasmin = globaltasmin[jj,ii,]
   pr = globalpr[jj,ii,]
   
}

inputs = array(data=NA, dim=c(1,40))

tas = (tasmax+tasmin)/2

# mean growing season temperature
if (is.element(1,vars)){
	inputs[1,1] = mean(tas[seasonidx])
}
# mean growing season precipitation
if (is.element(2,vars)){
	inputs[1,2] = mean(pr[seasonidx])
}
# number of growing season days where tasmax exceeds 30C
if (is.element(3,vars)){
	inputs[1,3] = length(which((tasmax[seasonidx])>30))
}
# number of growing season days where tasmax exceeds 35C
if (is.element(4,vars)){
	inputs[1,4] = length(which((tasmax[seasonidx])>35))
}
# number of growing season days where tasmin is less than 0C
if (is.element(5,vars)){
	inputs[1,5] = length(which((tasmin[seasonidx])<0))
}
# number of growing season days where tasmin is less than 5C
if (is.element(6,vars)){
	inputs[1,6] = length(which((tasmin[seasonidx])<5))
}
# number of growing season days where pr exceeds 1mm
if (is.element(7,vars)){
	inputs[1,7] = length(which((pr[seasonidx])>1))
}
# maximum consecutive growing season days where pr<0.01mm
if (is.element(8,vars)){
	subpr = pr[seasonidx]
	streak=numeric(0)
	if (subpr[1]<0.01){
		streak[1] = 1
	} else {
		streak[1] = 0
	}
	for (dd in c(2:length(subpr))){
		if (subpr[dd]<0.01){
			streak[dd] = streak[dd-1] + 1
			streak[dd-1] = 0
		} else {
			streak[dd] = 0
		}
	}
	inputs[1,8] = max(streak)
}
# mean temperature during planting window
if (is.element(9,vars)){
	inputs[1,9] = mean(tas[plantingidx])
}
# mean temperature before anthesis
if (is.element(10,vars)){
	inputs[1,10] = mean(tas[beforeidx])
}
# mean temperature during anthesis
if (is.element(11,vars)){
	inputs[1,11] = mean(tas[duringidx])
}
# mean temperature after anthesis
if (is.element(12,vars)){
	inputs[1,12] = mean(tas[afteridx])
}
# mean precipitation during planting window
if (is.element(13,vars)){
	inputs[1,13] = mean(pr[plantingidx])
}
# mean precipitation before anthesis
if (is.element(14,vars)){
	inputs[1,14] = mean(pr[beforeidx])
}
# mean precipitation during anthesis
if (is.element(15,vars)){
	inputs[1,15] = mean(pr[duringidx])
}
# mean precipitation after anthesis
if (is.element(16,vars)){
	inputs[1,16] = mean(pr[afteridx])
}
# number of planting window days where tasmax exceeds 30C
if (is.element(17,vars)){
	inputs[1,17] = length(which((tasmax[plantingidx])>30))
}
# number of days before anthesis where tasmax exceeds 30C
if (is.element(18,vars)){
	inputs[1,18] = length(which((tasmax[beforeidx])>30))
}
# number of days during anthesis where tasmax exceeds 30C
if (is.element(19,vars)){
	inputs[1,19] = length(which((tasmax[duringidx])>30))
}
# number of days after anthesis where tasmax exceeds 30C
if (is.element(20,vars)){
	inputs[1,20] = length(which((tasmax[afteridx])>30))
}
# number of planting window days where tasmax exceeds 35C
if (is.element(21,vars)){
	inputs[1,21] = length(which((tasmax[plantingidx])>35))
}
# number of days before anthesis where tasmax exceeds 35C
if (is.element(22,vars)){
	inputs[1,22] = length(which((tasmax[beforeidx])>35))
}
# number of days during anthesis where tasmax exceeds 35C
if (is.element(23,vars)){
	inputs[1,23] = length(which((tasmax[duringidx])>35))
}
# number of days after anthesis where tasmax exceeds 35C
if (is.element(24,vars)){
	inputs[1,24] = length(which((tasmax[afteridx])>35))
}
# number of planting window days where tasmin is less than 0C
if (is.element(25,vars)){
	inputs[1,25] = length(which((tasmin[plantingidx])<0))
}
# number of days before anthesis where tasmin is less than 0C
if (is.element(26,vars)){
	inputs[1,26] = length(which((tasmin[beforeidx])<0))
}
# number of days during anthesis where tasmin is less than 0C
if (is.element(27,vars)){
	inputs[1,27] = length(which((tasmin[duringidx])<0))
}
# number of days after anthesis where tasmin is less than 0C
if (is.element(28,vars)){
	inputs[1,28] = length(which((tasmin[afteridx])<0))
}
# number of planting window days where tasmin is less than 5C
if (is.element(29,vars)){
	inputs[1,29] = length(which((tasmin[plantingidx])<5))
}
# number of days before anthesis where tasmin is less than 5C
if (is.element(30,vars)){
	inputs[1,30] = length(which((tasmin[beforeidx])<5))
}
# number of days during anthesis where tasmin is less than 5C
if (is.element(31,vars)){
	inputs[1,31] = length(which((tasmin[duringidx])<5))
}
# number of days after anthesis where tasmin is less than 5C
if (is.element(32,vars)){
	inputs[1,32] = length(which((tasmin[afteridx])<5))
}
# number of planting window days where pr exceeds 1mm
if (is.element(33,vars)){
	inputs[1,33] = length(which((pr[plantingidx])>1))
}
# number of days before anthesis where pr exceeds 1mm
if (is.element(34,vars)){
	inputs[1,34] = length(which((pr[beforeidx])>1))
}
# number of days during anthesis where pr exceeds 1mm
if (is.element(35,vars)){
	inputs[1,35] = length(which((pr[duringidx])>1))
}
# number of days after anthesis where pr exceeds 1mm
if (is.element(36,vars)){
	inputs[1,36] = length(which((pr[afteridx])>1))
}
# maximum consecutive planting window days where pr<0.01mm
if (is.element(37,vars)){
	subpr = pr[plantingidx]
	streak=numeric(0)
	if (subpr[1]<0.01){
		streak[1] = 1
	} else {
		streak[1] = 0
	}
	for (dd in c(2:length(subpr))){
		if (subpr[dd]<0.01){
			streak[dd] = streak[dd-1] + 1
			streak[dd-1] = 0
		} else {
			streak[dd] = 0
		}
	}
	inputs[1,37] = max(streak)
}
# maximum consecutive days before anthesis where pr<0.01mm
if (is.element(38,vars)){
	subpr = pr[beforeidx]
	streak=numeric(0)
	if (subpr[1]<0.01){
		streak[1] = 1
	} else {
		streak[1] = 0
	}
	for (dd in c(2:length(subpr))){
		if (subpr[dd]<0.01){
			streak[dd] = streak[dd-1] + 1
			streak[dd-1] = 0
		} else {
			streak[dd] = 0
		}
	}
	inputs[1,38] = max(streak)
}
# maximum consecutive days during anthesis where pr<0.01mm
if (is.element(39,vars)){
	subpr = pr[duringidx]
	streak=numeric(0)
	if (subpr[1]<0.01){
		streak[1] = 1
	} else {
		streak[1] = 0
	}
	for (dd in c(2:length(subpr))){
		if (subpr[dd]<0.01){
			streak[dd] = streak[dd-1] + 1
			streak[dd-1] = 0
		} else {
			streak[dd] = 0
		}
	}
	inputs[1,39] = max(streak)
}
# maximum consecutive days after anthesis where pr<0.01mm
if (is.element(40,vars)){
	subpr = pr[afteridx]
	streak=numeric(0)
	if (subpr[1]<0.01){
		streak[1] = 1
	} else {
		streak[1] = 0
	}
	for (dd in c(2:length(subpr))){
		if (subpr[dd]<0.01){
			streak[dd] = streak[dd-1] + 1
			streak[dd-1] = 0
		} else {
			streak[dd] = 0
		}
	}
	inputs[1,40] = max(streak)
}

v = inputs[vars]

emuyield = coeffs[1] + (coeffs[2]*v[1]) + (coeffs[3]*v[2]) + (coeffs[4]*v[3]) + (coeffs[5]*v[4]) + (coeffs[6]*v[5]) + (coeffs[7]*v[1]*v[2]) + (coeffs[8]*v[1]*v[3]) + (coeffs[9]*v[1]*v[4]) + (coeffs[10]*v[1]*v[5]) + (coeffs[11]*v[2]*v[3]) + (coeffs[12]*v[2]*v[4]) + (coeffs[13]*v[2]*v[5]) + (coeffs[14]*v[3]*v[4]) + (coeffs[15]*v[3]*v[5]) + (coeffs[16]*v[4]*v[5]) + (coeffs[17]*v[1]*v[1]) + (coeffs[18]*v[2]*v[2]) + (coeffs[19]*v[3]*v[3]) + (coeffs[20]*v[4]*v[4]) + (coeffs[21]*v[5]*v[5])

return(emuyield)

}


