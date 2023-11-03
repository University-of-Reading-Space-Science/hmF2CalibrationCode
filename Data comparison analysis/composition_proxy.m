function [proxy, az_corrected_proxy] = composition_proxy(foF2, F107, lat, lon, year, month)

% Estimate composition from noon foF2 values.

proxy = (foF2.^2)./F107;

% Correct for solar zenith angle (assume overhead sun)
path(path,'C:\Users\jk902720\Dropbox\Research\useful matlab files\');

dims = size(year); % assume a 2D array

LT = 12-(lon/15) % calculate local noon, as sunzena requires UT

for i=1:dims(1)
    for j=1:dims(2)
        [day(i,j)] = day_no(year(i,j), month(i,j), 15); % Take the 15th as approximately mid-month
        
        [sunza(i,j),sunaz(i,j)] = sunzena((pi/180)*lat,(pi/180)*lon,year(i,j),day(i,j),LT,0,0)
    end
end

plot(180*sunza./pi)

az_corrected_proxy = proxy./cos(sunza);

end