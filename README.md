# StarSensor
Model of a star sensor using MATLAB

**To run the code write:**
stars = load('stars_cat.mat');
stars.dec = dec;
stars.mag = mag;
stars.ra = ra;
[A_eps, A_BN_tilde] = star_sensor(stars, ATT_mat, FOV_in_deg)
