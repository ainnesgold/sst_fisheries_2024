##CMIP6 Data
##This code reads in the CMIP6 data and creates an anomaly dataframe

#downloaded from: https://psl.noaa.gov/ipcc/cmip6/timeseries6.html

library(ncdf4)

# Open the NetCDF file
ncfile <- nc_open("CMIP6 trials/noaasst.nc")

# Read a variable
anomaly <- ncvar_get(ncfile, "anomaly")
year <- ncvar_get(ncfile, "year")

anomaly_df <- cbind(year, anomaly)
anomaly_df <- as.data.frame(anomaly_df) %>%
  filter(year > 2020 & year < 2090)

write.csv(anomaly_df, "anomaly_df.csv")