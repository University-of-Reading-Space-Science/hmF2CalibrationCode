% Script to determine the monthly median and average F2 layer height from
% the Japan MU radar data. Adapted to take a time of interest (TOI) as an
% initial argument so that the nature of the long-term variation can be
% determined as a time of day as well as seasonally.
%
% CJS April 2023
%


years = 1986:2020;

% Set height range for calculating average profile (to look for changes in gradiant
% of lower edge) and intitialise arrays
profile_fit_height = [100:5:500];
annual_average_profile = NaN*ones(length(profile_fit_height),length(years));
dannual_average_profile = annual_average_profile;

% Set min and max heights from which the F2 peak is both fitted and
% identified from the max value in thei height range.
min_ht = 180;
max_ht = 500;

% Set time of interest (assume local time)
TOI = 0:23;

% TOI = TOI-9; % Convert to UT
% tmp = find(TOI<0);
% TOI(tmp) = TOI(tmp)+24;

%Initialise arrays
mean_all_F2_max_h_TOI = NaN*ones(length(years),12,24);
dmean_all_F2_max_h_TOI = NaN*ones(length(years),12,24);
median_all_h = NaN*ones(length(years),12,24);
mean_all_h = NaN*ones(length(years),12,24);
dmedian_all_h = NaN*ones(length(years),12,24);
dmean_all_h = NaN*ones(length(years),12,24);
geo_mean_h = NaN*ones(length(years),12,24);
dgeo_mean_h = NaN*ones(length(years),12,24);

for hr_loop = 1:length(TOI)
for i=1:length(years)
    for j = 1:12
    filelist_str = [num2str(years(i)) num2str(j,'%02d') '*pwr.nc'];

    % filelist = dir('2001*pwr.nc');
    filelist = dir(filelist_str);

    all_h = [];
    all_t = [];
    all_mean_F2_max_h = [];
    all_median_F2_max_h = [];
    all_meanISRpower = [];
    all_medianISRpower = [];

    for m=1:length(filelist)
        % [mean_F2_max_h, F2_max_t, meanISRpeak, meanISRheight, meanISRpower, medianISRheight, medianISRpower] = find_ISR_F2_peak(filelist(m).name, min_ht, max_ht);
        [mean_F2_max_h, F2_max_t, meanISRpeak] = find_ISR_F2_peak(filelist(m).name, min_ht, max_ht);
        
        all_h = [all_h; meanISRpeak];
        all_t = [all_t; F2_max_t];
        
        % concatenate into an array of all times
        all_mean_F2_max_h = [all_mean_F2_max_h; mean_F2_max_h];
%        all_median_F2_max_h = [all_median_F2_max_h, median_F2_max_h];
    end
   
    datevec_all_t = datevec(all_t);
    % limit to non values and calculate annual averages
    subset_TOI = find(datevec_all_t(:,4) == TOI(hr_loop)); % assuming local time...
%     noon = find(datevec_all_t(:,4) == 21); % assuming GMT time...
    
%     median_all_F2_max_h_noon(i) = nanmedian(all_median_F2_max_h(noon));
%     dmedian_all_F2_max_h_noon(i) = 1.253*nanstd(all_median_F2_max_h(noon))./sqrt(length(years));
    
    mean_all_F2_max_h_TOI(i,j,hr_loop) = nanmean(all_mean_F2_max_h(subset_TOI));
    dmean_all_F2_max_h_TOI(i,j,hr_loop) = nanstd(all_mean_F2_max_h(subset_TOI))./sqrt(length(subset_TOI));
    
%     median_all_h(i,j, hr_loop) = nanmedian(all_h(subset_TOI));
%     dmedian_all_h(i,j, hr_loop) = 1.253*nanstd(all_h(subset_TOI))/sqrt(length(subset_TOI));

% Implicily contains threshold of (usually > 5 values required before a value is calculated) 
% - used to be consistent with ionosonde analysis. CJS April 2023
% This produces no data. Reduced to 3 (introduced as a second variable in
% median_bootstrap) CJS 2023
    [median_all_h(i,j,hr_loop), dmedian_all_h(i,j,hr_loop)] = median_bootstrap(all_h(subset_TOI),3);
    
    geo_mean_h(i, j, hr_loop) = 10.^(nanmean(log10(all_h(subset_TOI))));
    dgeo_mean_h(i,j, hr_loop) = 10.^(nanstd(log10(all_h(subset_TOI)))./sqrt(length(subset_TOI(~isnan(all_h(subset_TOI))))));

    mean_all_h(i,j, hr_loop) = nanmean(all_h(subset_TOI));
    % dmean_all_h(i,j, hr_loop) = nanstd(all_h(subset_TOI))./sqrt(length(~isnan(all_h(subset_TOI))));
    dmean_all_h(i,j, hr_loop) = nanstd(all_h(subset_TOI))./sqrt(length(find(isnan(all_h(subset_TOI))==0))); %less efficient but at least correct!
      
    end
end
end

figure(999),errorbar(years, nanmedian(median_all_h(:,:,13), 2), 1.253*nanstd(median_all_h(:,:,13),0,2)./sqrt(12));
hold on, errorbar(years, nanmedian(mean_all_F2_max_h_TOI(:,:,13),2), 1.253*nanstd(mean_all_F2_max_h_TOI(:,:,13),0,2)./sqrt(12))
hold off

figure(1000),plot(years, nanmedian(mean_all_F2_max_h_TOI(:,:,13),2))
hold
figure(1000),plot(years, nanmean(mean_all_F2_max_h_TOI(:,:,13),2))

figure(1001),plot(1:12, nanmean(mean_all_F2_max_h_TOI(:,:,13),1))
hold
figure(1001),plot(1:12, nanmedian(mean_all_F2_max_h_TOI(:,:,13),1))

figure(1003),plot(reshape(mean_all_F2_max_h_TOI(:,:,13)',1,length(years)*12))
reshaped_mean_all_F2_max_h_noon = reshape(mean_all_F2_max_h_TOI(:,:,13)',1,length(years)*12);
figure(10005),plot(reshaped_mean_all_F2_max_h_noon(6:12:end))
figure(10005),hold on, plot(reshaped_mean_all_F2_max_h_noon(1:12:end))
figure(10005),hold on, plot(reshaped_mean_all_F2_max_h_noon(3:12:end))
figure(10005),hold on, plot(reshaped_mean_all_F2_max_h_noon(9:12:end))

figure(1007),plot(1:12, (mean_all_F2_max_h_TOI(:,:,13)))
