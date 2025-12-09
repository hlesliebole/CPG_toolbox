 

      Quarterly Averaged Subaerial Cross-Shore Profiles on State Beach MOP Lines 


1. Filename Convention

  StateBeachName_MopID_YYYY_Q.dat

  eg. SanElijo_D0685_2023_3.dat


2. State Beaches and Mop Ranges

  State Beach                Mop Range

  BorderField              D0002 - D0028
  SilverStrand             D0085 - D0157
  TorreyPines              D0536 - D0607
  Cardiff                  D0664 - D0683
  SanElijo                 D0684 - D0706
  Moonlight                D0717 - D0725
  Leucadia                 D0740 - D0757
  SouthCarlsbad            D0762 - D0819
  Carlsbad                 D0825 - D0854


3. File Format

  Column           Variable

   1         - Profile elevation cross-chore location (m) on the MOP transect line relative
              to the Mop back beach point.
   2         - UTM Zone '11 S' Easting (m) of profile elevation location.
   3         - UTM Zone '11 S' Northing (m) of profile elevation location.
   4         - Longitude of profile elevation location.
   5         - Latitude of profile elevation location.
   6         - Profile elevation (m, NAVD88) ; NaN = No Data ; MSL = 0.774m NAVD88


4. Misc Notes

  Profiles are derived from Truck LiDAR and ATV MiniRanger LiDAR surveys only.

  Individual survey profiles are first reduced to monthly averages and then a quarterly
  (JFM, AMJ, JAS, OND) average.

  Profiles are screened for elevation outliers using the MATLAB rmoutliers.m function with 
  the default median absolute deviations (MAD) method.   