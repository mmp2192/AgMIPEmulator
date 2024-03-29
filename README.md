# AgMIPEmulator

This emulator tool is trained on crop yield simulations from the AgMIP Global Gridded Crop Model Intercomparison 
(GGCMI) Project Phase II, which conducted a sensitivity analysis on the response of yield to varying climate inputs.
(https://gmd.copernicus.org/articles/13/2315/2020/)  
  
We select yield simulations from a temperature range of -1C to +2C and precipitation range of -20% to +20% (changes 
from historical 1980-2010 AgMERRA climate data) to more accurately represent seasonal variability in the historical 
period as well as potential seasonal climate forecasts. Carbon and nitrogen dimensions are kept constant at baseline 
values (360 ppm CO2 and 60 kgN/ha).  
  
For each 0.5-degree pixel in a given region, we fit a quadratic least-squares regression of climate inputs against 
simulated yields, starting with the 40 available climate indices listed below.   
v1 = mean growing season temperature in degrees C  
v2 = mean growing season precipitation in mm/day  
v3 = number of days in growing season when TMAX>30C  
v4 = number of days in growing season when TMAX>35C  
v5 = number of days in growing season when TMIN<0C  
v6 = number of days in growing season when TMIN<5C  
v7 = number of days in growing season when PR>1mm  
v8 = maximum number of consecutive days in growing season when P<0.01mm  
v9 = mean growing season temperature in degrees C during planting window  
v10 = mean growing season temperature in degrees C before anthesis  
v11 = mean growing season temperature in degrees C during anthesis  
v12 = mean growing season temperature in degrees C after anthesis  
v13 = mean growing season precipitation in mm/day during planting window  
v14 = mean growing season precipitation in mm/day before anthesis  
v15 = mean growing season precipitation in mm/day during anthesis  
v16 = mean growing season precipitation in mm/day after anthesis  
v17 = number of days when TMAX>30C during planting window  
v18 = number of days when TMAX>30C before anthesis  
v19 = number of days when TMAX>30C during anthesis  
v20 = number of days when TMAX>30C after anthesis  
v21 = number of days when TMAX>35C during planting window  
v22 = number of days when TMAX>35C before anthesis  
v23 = number of days when TMAX>35C during anthesis  
v24 = number of days when TMAX>35C after anthesis  
v25 = number of days when TMIN<0C during planting window  
v26 = number of days when TMIN<0C before anthesis  
v27 = number of days when TMIN<0C during anthesis  
v28 = number of days when TMIN<0C after anthesis  
v29 = number of days when TMIN<5C during planting window  
v30 = number of days when TMIN<5C before anthesis  
v31 = number of days when TMIN<5C during anthesis  
v32 = number of days when TMIN<5C after anthesis  
v33 = number of days when PR>1mm during planting window  
v34 = number of days when PR>1mm before anthesis  
v35 = number of days when PR>1mm during anthesis  
v36 = number of days when PR>1mm after anthesis  
v37 = maximum number of consecutive days when P<0.01mm during planting window  
v38 = maximum number of consecutive days when P<0.01mm before anthesis  
v39 = maximum number of consecutive days when P<0.01mm during anthesis  
v40 = maximum number of consecutive days when P<0.01mm after anthesis  
  
We reduce the number of input variables to five per pixel by running the regression four times, eliminating half of 
the input variables between rounds. The process of elimination involves calculating the percent of statistically 
significant (p-value < 0.05) terms associated with a given climate input, including cross-terms, out of the total number 
of terms for that variable. We then rank the climate inputs based on this percentage and remove the lower-ranking half.  
  
The five most significant variables for each crop are stored as a list of five indicators (corresponding to the list 
above) per pixel in a .mat file online. Similarly, the 21 coefficients for a 5-variable quadratic linear regression equation 
for each crop are stored as 21 values per pixel in a .mat file. In order to run the emulator using any kind of climate data 
for a given crop and pixel, you will need to be connected to the internet in order to access both the list of indicators 
and the list of coefficients (and climate input data if needed, see below).  
  
To run, you will also need climate input data. This must be on a daily scale for a full year, one vector for each maximum
temperature, minimum temperature, and precipitation variables. If you don't have a particular set of climate data or forecast, 
you can run the emulator in a mode that will use a dataset of combined CHIRTS temperature and CHIRPS precipitation. If so, 
you'll need to input the year of CHIRTS/CHIRPS data you want to use (must be between 1983 and 2016).  
  
In either case, you'll also need to input the latitude and longitude coordinates of the area you're interested (can be 
anywhere across the globe as long as it's not over the ocean).   
  
The scripts in this repository work as follows:  
  
agmip_emulator.R: converts daily climate data into the 40 possible climate input variables, then uses the relevant climate 
input variables and pre-calculated coefficients to calculate emulated yield anomaly for the given pixel from the 1980-2010 
baseline.  
  
Currently this emulator works for rainfed maize only, until we are able to run multiple crops and irrigation types on a 
global scale.  
