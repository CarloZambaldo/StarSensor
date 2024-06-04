# StarSensor
Model of a star sensor using MATLAB

**To run the code write:**
```
stars = load('stars_cat.mat');
stars.dec = dec;
stars.mag = mag;
stars.ra = ra;
[A_eps, A_BN_tilde] = star_sensor(stars, ATT_mat, FOV_in_deg)
```

where:
* ```ATT_mat``` is the attitude matrix
* ```FOV_in_deg``` is the field of view of the startracker in degrees
* ```A_BN_tilde``` is the attitude matrix read by the sensor
