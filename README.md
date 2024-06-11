## Regional_Glaciation_data

The purpose of this repository is to download original data, preprocess them, convert and save them to GeoTIFF and NetCDF formats, then visualize them for a watershed covering western Canada.  

See the [Clarke et al., 2015]( https://www.nature.com/articles/ngeo2407) and data source [RGM_archive]( https://couplet.unbc.ca/data/RGM_archive/) for more information.  

This repository is a work in progress! 

# Scripts

## RGM_convert.py

The purpose of this script is to combine glacier regions produced from the Regional Glaciation Model for Western Canada, which are saved in .mat format, and save data to NetCDF and GeoTIFF formats by considering its original CRS projection information provided in the readme file. As a test case example, I used CanESM-320km/RCP85 data. 
