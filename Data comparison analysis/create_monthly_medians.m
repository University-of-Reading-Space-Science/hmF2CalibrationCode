load kokobunji_corrected_all.mat

% path(path,'C:\Users\jk902720\Dropbox\Research\useful matlab files')

% Need to fill in nighttime foE values equivalent to 0.4MHz as done by
% Jarvis et al 1998. To do this, need solar zenith angle

%% The following section is commented out as solar zenith angle was calculated independenty and stored in directory as 'matched_sza.mat' and loaded where needed
%% Code is repeated here to demonstrate how it was calculated

% Zenith angle calculation (needed for estimating nighttime hmF2) requires UT
% LT = (0:23);
% UT=LT-9;
% UT(UT<0) = UT(UT<0)+24;
% 
% year_loop = 1986:2020;
% 
% for i=1:35
%   for j=1:12
%       LT = (0:23);
%       UT=LT-9;
%       UT(UT<0) = UT(UT<0)+24;
%       for k=1:24
%         day_number = day_no(year_loop(i),j,15); % Day 15 taken as approx middle day of each month
%        [sza(i,j,k), ~] = sunzena((pi/180)*35.7104,(pi/180)*139.4632,year_loop(i),day_number,UT(k),0,0);
%       end
%   end
% end

% save matched_sza sza  

load matched_sza % Calculated from code in loop above so no need for both.

% initialise arrays
median_foF2 = NaN*ones(length(u_year),12,24);
median_foE = NaN*ones(length(u_year),12,24);
median_foF1 = NaN*ones(length(u_year),12,24);
median_M3000F2 = NaN*ones(length(u_year),12,24);
dmedian_foF2 = NaN*ones(length(u_year),12,24);
dmedian_foE = NaN*ones(length(u_year),12,24);
dmedian_foF1 = NaN*ones(length(u_year),12,24);
median_foF1_occurrence = NaN*ones(length(u_year),12,24);
dmedian_M3000F2 = NaN*ones(length(u_year),12,24);

hmF2_BradDud = NaN*ones(length(u_year),12,24);
dhmF2_BradDud = NaN*ones(length(u_year),12,24);

% median_hmF2_BradDud = NaN*ones(length(u_year),12,24);
% dmedian_hmF2_BradDud = NaN*ones(length(u_year),12,24);

median_am_hourly = NaN*ones(length(u_year),12,24);  % Hourly matrices, so will need to interpolate
dmedian_am_hourly = NaN*ones(length(u_year),12,24);
median_f107_hourly = NaN*ones(length(u_year),12,24); % Hourly matrices, so will need to interpolate
dmedian_f107_hourly = NaN*ones(length(u_year),12,24);

year_index_iono = NaN*median_foF2;
month_index_iono = year_index_iono;
hour_index_iono = year_index_iono;
datenum_index_iono = year_index_iono;

minN = 5; % minimum number of counts which create a viable median.

for i=1:length(u_year)
    for j = 1:12
        for k=1:24
          year_index_iono(i,j,k) = u_year(i);
          month_index_iono(i,j,k) = j;
          hour_index_iono(i,j,k) = k-1;
          datenum_index_iono(i,j,k) = datenum(year_index_iono(i,j,k),month_index_iono(i,j,k),15,hour_index_iono(i,j,k),0,0); % For the purposes of plot assume day 15 of each month as central.
          subset = find(year == u_year(i) & month==j & hour == (k-1));
          [median_foF2(i,j,k), dmedian_foF2(i,j,k)] = median_bootstrap(foF2(subset),minN);
          [median_foE(i,j,k), dmedian_foE(i,j,k)] = median_bootstrap(foE(subset),minN);
          % If approximate night-time values should be included, uncomment the three
          % lines below to add an approximate foE of 0.4MHz.
          
          % NB assumes sza values do not vary by year and so uses the first
          % year in the sza array for all years
          if sza(1,j,k) >= (90*pi/180) && isnan(median_foE(i,j,k)) == 1  % If nighttime, set foE to 0.4 MHz as done in Jarvis et al. - should be 100 degrees so this is conservative but this produces a gap in plots when foE < min ionosonde freq (~1MHz)
              median_foE(i,j,k) = 0.4;
          end
          
          [median_M3000F2(i,j,k), dmedian_M3000F2(i,j,k)] = median_bootstrap(M3000F2(subset),minN);
          [median_foF1(i,j,k), dmedian_foF1(i,j,k)] = median_bootstrap(foF1(subset),minN);
         
          % Calculate hmF2 using BRadley-Dudeney formula
          [hmF2_BradDud(i,j,k), dhmF2_BradDud(i,j,k)] = hmF2BradDud(median_foF2(i,j,k),median_foE(i,j,k), median_M3000F2(i,j,k)/10,dmedian_foF2(i,j,k),dmedian_foE(i,j,k), dmedian_M3000F2(i,j,k)/10 );
        end
    end
end

% Then need to extrapolate am values which are 3 hourly

for i=1:length(am_u_year)
    for j=1:12
        for k=1:22
           if length(find(am_hour == (k-1) ))>0
                am_subset = find(am_year == am_u_year(i) & am_month == j & am_hour == (k-1));
               [median_am_hourly(i,j,k), dmedian_am_hourly(i,j,k)] = median_bootstrap(am(am_subset),minN);        
                median_am_hourly(i,j,k+1) = median_am_hourly(i,j,k);
                median_am_hourly(i,j,k+2) = median_am_hourly(i,j,k);
                dmedian_am_hourly(i,j,k+1) = dmedian_am_hourly(i,j,k);
                dmedian_am_hourly(i,j,k+2) = dmedian_am_hourly(i,j,k);
           end
        end
    end
end

% Then need to extrapolate f107 values which are taken at 17:00 and 00:00
% until the 1990s so calculate simple monthly average and store hourly for
% ease of comparison
load f107_all_data

 for i=1:length(f107_u_year)
     for j=1:12
            if length(find(f107_month == j ))>0
                 f107_subset = find(f107_year == f107_u_year(i) & f107_month == j);
               [median_f107_hourly(i,j,1:24), dmedian_f107_hourly(i,j,1:24)] = median_bootstrap(observed_f107(f107_subset),minN);        
           end
     end
 end

load SSN_data.mat

SSN_u_year = unique(SSNyear);
 monthly_SSN = NaN*ones(length(SSN_u_year),12,24);

 for i=1:length(SSN_u_year)
     for j=1:12
            SSN_subset = find(SSNyear == SSN_u_year(i) & SSNmonth == j);
            if length(SSN_subset)>0
              monthly_SSN(i,j,1:24) = SSN_monthly(SSN_subset);
            end
     end
 end
 