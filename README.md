# Bharatpur-LULC-Change-GEE
Land use/land cover change detection (2005-2025) for Bharatpur Metropolitan City, Nepal, using Google Earth Engine, Landsat and Random Forest
 Land Use / Land Cover Change in Bharatpur Metropolitan City, Nepal (2005-2025)
 
Land cover change detection for Bharatpur Metropolitan City using **Google Earth Engine (GEE)**, multi-temporal **Landsat** imagery and a **Random Forest** classifier, with a focus on forest gain, loss and persistence.
 
*MSc Sustainability and Environmental Management thesis, Coventry University (2026).*
 
## Land Cover Maps
 
### 2005
<img width="505" height="352" alt="image" src="https://github.com/user-attachments/assets/327b4076-7cc5-4d01-9928-07af0aeaa74b" />

 
### 2015
<img width="468" height="351" alt="image" src="https://github.com/user-attachments/assets/ff926f4c-89fc-4cfd-bdca-5887f6007f87" />

 
### 2025
<img width="464" height="389" alt="image" src="https://github.com/user-attachments/assets/46df1809-073f-4a0b-8163-1efbafd118db" />

 
## Overview
 
Bharatpur sits next to Chitwan National Park and is one of Nepal's fastest-growing cities. This project maps five land cover classes (Forest, Water Body, Built-up, Agriculture, Other Land) for 2005, 2015 and 2025, then measures how the city changed and what land urban growth consumed.
 
**Main findings**
- Built-up area grew by **28.53 km² (+30.0%)** between 2005 and 2025, and the pace nearly doubled in the second decade (+0.92%/yr to +1.91%/yr).
- Urban growth came mostly from **agricultural land** (40.94 km² converted to built-up), not forest (12.19 km²).
- Forest showed a **net loss in 2005-2015 and a net gain in 2015-2025**, so net 20-year change hides a lot of simultaneous loss and gain.
## Data and Methods
 
| Item | Detail |
|------|--------|
| Imagery | Landsat 5 TM (2005), Landsat 8 OLI (2015), Landsat 9 OLI-2 (2025), Collection 2 Level 2 surface reflectance, 30 m |
| Compositing | Median composite, 1 September to 31 December each year (post-monsoon) |
| Cloud masking | Bitwise mask on `QA_PIXEL` and `QA_RADSAT` |
| Features | 6 reflectance bands + NDVI, NDBI, MNDWI, NDMI, SAVI, BSI, IBI (13 variables) |
| Classifier | Random Forest, 200 trees, bag fraction 0.6, min leaf population 2, seed 42 |
| Samples | Digitized per year from the composite and high-resolution Google Earth imagery, 70/30 train/validation split |
| Post-processing | Rule-based water and built-up refinement, focal mode smoothing (1 pixel) |
| Change analysis | Post-classification comparison, transition matrices, forest gain/loss/persistence, annual rate of change |
 
## Results
 
| Class | 2005 (km²) | 2015 (km²) | 2025 (km²) | Net change 2005-2025 |
|-------|-----------|-----------|-----------|----------------------|
| Forest | 175.69 | 157.93 | 164.44 | -11.24 km² (-6.4%) |
| Built-up | 94.97 | 103.67 | 123.51 | +28.53 km² (+30.0%) |
| Agriculture | 130.23 | 126.64 | 113.12 | -17.11 km² (-13.1%) |
| Water Body | 18.02 | 27.05 | 12.65 | -5.36 km² (-29.8%) |
 
Other Land increased by 5.87 km² (+43.4%).
 
**Accuracy**
 
| Year | Overall accuracy | Kappa |
|------|------------------|-------|
| 2005 | 80.59% | 0.743 |
| 2015 | 87.80% | 0.841 |
| 2025 | 93.04% | 0.908 |
 
**Forest change components**
 
| Period | Gross loss | Gross gain | Net |
|--------|-----------|-----------|-----|
| 2005-2015 | 29.88 km² | 14.41 km² | -15.47 km² |
| 2015-2025 | 8.96 km² | 16.24 km² | +7.28 km² |
 

 
## Limitations
 
- Based only on satellite data. No field-based ground truth was collected, so training and validation samples rely on visual interpretation.
- Three different Landsat sensors were used, which can cause small inconsistencies between years. The **Water Body** trend (18.02, 27.05, then 12.65 km²) is probably affected by this and by river flow variation, so treat it with caution.
- 30 m resolution mixes pixels in dense or fragmented areas such as riverine forest and the urban edge.
- The 20-year net forest change from area totals (-11.24 km²) differs from the pixel-based gain/loss result (-5.97 km²) because the endpoint maps were classified independently.
- The study links forest recovery to community forestry only as a possible explanation. No spatial data on community forests was available to test it.
