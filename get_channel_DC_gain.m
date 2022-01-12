function h = get_channel_DC_gain(xt,yt,zt,xr,yr,zr)

%%%---Parameters---%%%
distance = sqrt((xr - xt)^2 + (yr - yt)^2 + (zr - zt)^2); % distance between Tx and Rx (m)
height = zt - zr;
phi_irradiance = acosd(height/distance); % angle of irradiance (degree)
phi_incidence = acosd(height/distance); % angle of incidence (degree)

FOV_LED = 60; % transmitter semi-angle at half power (degree)
m = -log(2)/log(cosd(FOV_LED)); % order of Lambertian emission
area_PD = 10^(-4); % Physical area of the PD, 1 cm^2
gain_filter = 0.9; % gain of the optical filter
RI = 1.5; % refraction index of lens
FOV_lens = 72; % half-angle FOV of lens
gain_lens = (RI/sind(FOV_lens))^2; % gain of the imaging lens

%%%---Parameters---%%%
h = 2*(m+1)*area_PD*((cosd(phi_irradiance))^m)*gain_filter*gain_lens*...
    cosd(phi_incidence)/2/pi/(distance^2); % channel dc gain
