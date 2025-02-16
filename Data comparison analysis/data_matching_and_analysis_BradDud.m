%% Calibrate with ISR data and estimate long-term trends unaccounted for by
%% geomagnetic activity and solar variability

%% load in the ISR monthly medians (calculated using)ISR_F2_loop_monthly_hourly.m
load post_MIST_MU_analysis_allUT_arrays_SNR0_05_mean_error_corrected.mat

%% load in geomagnetic 'am' index
load am_data.mat

%% Now load in the F10.7cm radio flux data
load f107_data.mat

%% Create monthly medians from ionosonde data and am index
create_monthly_medians

% Now need to find common start and end years

ISR_years = years;
start_year = max([f107_u_year(1), am_u_year(1), ISR_years(1), u_year(1), SSN_u_year(1)]);
end_year = min([f107_u_year(end), am_u_year(end), ISR_years(end), u_year(end), SSN_u_year(end)]);
matched_years = start_year:end_year;

%% create matched arrays

sub_f107 = find(f107_u_year >= start_year & f107_u_year <=end_year);
matched_f107 = f107_monthly_median(sub_f107,:,:);

sub_am = find(am_u_year >= start_year & am_u_year <=end_year);
matched_am = median_am_hourly(sub_am,:,:);

sub_SSN = find(SSN_u_year >= start_year & SSN_u_year <=end_year);
matched_SSN = monthly_SSN(sub_SSN,:,:);

sub_iono = find(u_year >= start_year & u_year <=end_year);

matched_BradDud = hmF2_BradDud(sub_iono,:,:); % Calculated from the monthly median data
dmatched_BradDud = dhmF2_BradDud(sub_iono,:,:);

% % Match date indexes too...
year_index = year_index_iono(sub_iono,:,:);
month_index = month_index_iono(sub_iono,:,:);
hour_index = hour_index_iono(sub_iono,:,:);

sub_ISR = find(ISR_years >= start_year & ISR_years <=end_year);
% matched_ISR = median_all_h(sub_ISR,:,:);
% dmatched_ISR = dmedian_all_h(sub_ISR,:,:);
matched_ISR = mean_all_h(sub_ISR,:,:);
dmatched_ISR = dmean_all_h(sub_ISR,:,:);

%% Now can do the analysis

%% Plot everything ionosonde against everything ISR

%% convert to vectors
[m,n,p] = size(matched_ISR);

reshaped_BradDud_all = reshape(matched_BradDud,1,m*n*p);
dreshaped_BradDud_all = reshape(dmatched_BradDud,1,m*n*p);

reshaped_ISR_all = reshape(matched_ISR,1,m*n*p);
dreshaped_ISR_all = reshape(dmatched_ISR,1,m*n*p);

% Identify those points with foE set to 0.4MHz
reshaped_foE_all = reshape(median_foE(sub_iono,:,:),1,[]);

reshaped_year_index = reshape(year_index,1,m*n*p);
reshaped_month_index = reshape(month_index,1,m*n*p);
reshaped_hour_index = reshape(hour_index,1,m*n*p);

reshaped_sza = reshape(sza,1,m*n*p);

%% Throw out NaNs before fit
% bad = find( isnan(reshaped_BradDud_all+reshaped_ISR_all+dreshaped_BradDud_all+dreshaped_ISR_all)==1);
bad = find( isnan(reshaped_BradDud_all+reshaped_ISR_all)==1);

reshaped_BradDud_all(bad) = [];
reshaped_ISR_all(bad) = [];
dreshaped_BradDud_all(bad) = [];
dreshaped_ISR_all(bad) = [];

reshaped_foE_all(bad) = [];
foE_set_to_04 = find(reshaped_foE_all == 0.4);

reshaped_year_index(bad) = [];
reshaped_month_index(bad) = [];
reshaped_hour_index(bad) = [];

reshaped_sza(bad) = [];

% sza_lt_90 = find(reshaped_sza < 90*pi/180); % added in response to ref #2's request for revision of figure 2
% sza_ge90_le100 = find(reshaped_sza >= 90*pi/180 & reshaped_sza <= 100*pi/180);
% sza_gt_100 = find(reshaped_sza > 100*pi/180);

% Robust fit, accounting for outliers
    if length(reshaped_ISR_all)>2
        mdlr = fitlm(reshaped_ISR_all,reshaped_BradDud_all, 'RobustOpts','on');
        fitvalues_iono_corr_vs_ISR = table2array(mdlr.Coefficients);
        b_all = fitvalues_iono_corr_vs_ISR(2,1);
        db_all = fitvalues_iono_corr_vs_ISR(2,2);
        a_all = fitvalues_iono_corr_vs_ISR(1,1);
        da_all = fitvalues_iono_corr_vs_ISR(1,2);
%%
        rsquare_robustfit_all = corr(reshaped_BradDud_all',(a_all+b_all*reshaped_ISR_all)')^2;
%%
    else
        b_all = NaN;
        db_all = NaN;
        a_all = NaN;
        da_all = NaN;
    end

% % Create figure 2 in the paper
% % Print Bradley-Dudeney hmF2 values vs ISR values for all times, along with
% % robust fit (to minimise influence of outliers)
% figure(100)
% errorbar(reshaped_ISR_all, reshaped_BradDud_all, dreshaped_ISR_all, dreshaped_ISR_all, dreshaped_BradDud_all, dreshaped_BradDud_all,'k.','capsize',0,'color',[0.5 0.5 0.5])
% hold on
% plot(reshaped_ISR_all, reshaped_BradDud_all,'k.','markersize',6)
% % Now separate foE = 0.4 MHz points
% % plot(reshaped_ISR_all(foE_set_to_04),reshaped_BradDud_all(foE_set_to_04),'r.') 
% 
% % Now plot night and dawn/dusk values separately
% plot(reshaped_ISR_all(sza_gt_100),reshaped_BradDud_all(sza_gt_100),'c.')
% plot(reshaped_ISR_all(sza_ge90_le100),reshaped_BradDud_all(sza_ge90_le100),'m.')
% 
% 
% axis([50 500 50 500])
% axis('square')
% plot([50 500], b_all*[50 500]+a_all,'r')
% text(60,480,['All y=(' num2str(b_all,'%2.2f') ' \pm ' num2str(db_all,'%2.2f') ')x + (' num2str(a_all,'%2.2f') ' \pm ' num2str(da_all,'%2.2f') ')']);
% xlabel('ISR hmF2','fontsize',14)
% ylabel('Ionosonde hmF2','fontsize',14)
% hold off

%%
sza_lt_90 = find(reshaped_sza < 90*pi/180); % added in response to ref #2's request for revision of figure 2
sza_ge90_le100 = find(reshaped_sza >= 90*pi/180 & reshaped_sza <= 100*pi/180);
sza_gt_100 = find(reshaped_sza > 100*pi/180);

% Create revised figure 2 after referee #2 comments
% Print Bradley-Dudeney hmF2 values vs ISR values for all times, along with
% robust fit (to minimise influence of outliers)
figure(100)
% subplot(221)
tiledlayout(2,2)
nexttile
%sp211 = subplot('position',[0.15 0.55 0.4 0.4]);
errorbar(reshaped_ISR_all, reshaped_BradDud_all, dreshaped_ISR_all, dreshaped_ISR_all, dreshaped_BradDud_all, dreshaped_BradDud_all,'k.','capsize',0,'color',[0.5 0.5 0.5])
hold on
plot(reshaped_ISR_all, reshaped_BradDud_all,'k.','markersize',6)
% Now separate foE = 0.4 MHz points
% plot(reshaped_ISR_all(foE_set_to_04),reshaped_BradDud_all(foE_set_to_04),'r.') 

% Now plot night and dawn/dusk values separately
plot(reshaped_ISR_all(sza_gt_100),reshaped_BradDud_all(sza_gt_100),'c.')
plot(reshaped_ISR_all(sza_ge90_le100),reshaped_BradDud_all(sza_ge90_le100),'m.')

axis([50 500 50 500])
axis('square')
plot([50 500], b_all*[50 500]+a_all,'k.-')
text(60,480,['y=(' num2str(b_all,'%2.2f') '\pm' num2str(db_all,'%2.2f') ')x+(' num2str(a_all,'%2.2f') '\pm' num2str(da_all,'%2.2f') ')'],'fontsize',6);
title('All','fontsize',14)
% xlabel('ISR hmF2','fontsize',14)
ylabel('Ionosonde hmF2','fontsize',14)
hold off

% subplot(222)
%sp222 = subplot('position',[0.55 0.55 0.4 0.4]);
nexttile
% Daytime values only
errorbar(reshaped_ISR_all(sza_lt_90), reshaped_BradDud_all(sza_lt_90), dreshaped_ISR_all(sza_lt_90), dreshaped_ISR_all(sza_lt_90), dreshaped_BradDud_all(sza_lt_90), dreshaped_BradDud_all(sza_lt_90),'k.','capsize',0,'color',[0.5 0.5 0.5])
hold on
plot(reshaped_ISR_all(sza_lt_90), reshaped_BradDud_all(sza_lt_90),'k.','markersize',6)
axis([50 500 50 500])
axis('square')
title('Daytime','fontsize',14)
% xlabel('ISR hmF2','fontsize',14)
% ylabel('Ionosonde hmF2','fontsize',14)
hold off

% Robust fit, accounting for outliers


    if length(sza_lt_90)>2
        mdlr_all_day = fitlm(reshaped_ISR_all(sza_lt_90),reshaped_BradDud_all(sza_lt_90), 'RobustOpts','on');
        fitvalues_iono_corr_vs_ISR_day = table2array(mdlr_all_day.Coefficients);
        b_day = fitvalues_iono_corr_vs_ISR_day(2,1);
        db_day = fitvalues_iono_corr_vs_ISR_day(2,2);
        a_day = fitvalues_iono_corr_vs_ISR_day(1,1);
        da_day = fitvalues_iono_corr_vs_ISR_day(1,2);
%%
        rsquare_robustfit_day = corr(reshaped_BradDud_all(sza_lt_90)',(a_day+b_day*reshaped_ISR_all(sza_lt_90))')^2;
%%

    else
        b_day = NaN;
        db_day = NaN;
        a_day = NaN;
        da_day = NaN;
    end
% subplot(sp222)
    hold on
    plot([50 500], b_day*[50 500]+a_day,'k.-')
    text(60,480,['y=(' num2str(b_day,'%2.2f') '\pm' num2str(db_day,'%2.2f') ')x+(' num2str(a_day,'%2.2f') '\pm' num2str(da_day,'%2.2f') ')'],'fontsize',6);
    % xlabel('ISR hmF2','fontsize',14)
    % ylabel('Ionosonde hmF2','fontsize',14)
    hold off

% Now separate foE = 0.4 MHz points
% plot(reshaped_ISR_all(foE_set_to_04),reshaped_BradDud_all(foE_set_to_04),'r.') 

% Now plot night and dawn/dusk values separately
% subplot(223)
nexttile
% sp223 = subplot('position',[0.15 0.05 0.4 0.4]);
errorbar(reshaped_ISR_all(sza_gt_100), reshaped_BradDud_all(sza_gt_100), dreshaped_ISR_all(sza_gt_100), dreshaped_ISR_all(sza_gt_100), dreshaped_BradDud_all(sza_gt_100), dreshaped_BradDud_all(sza_gt_100),'k.','capsize',0,'color',[0.5 0.5 0.5]);
hold on
plot(reshaped_ISR_all(sza_gt_100), reshaped_BradDud_all(sza_gt_100),'c.','markersize',6)
axis([50 500 50 500])
axis('square')
title('Night-time','fontsize',14)
xlabel('ISR hmF2','fontsize',14)
ylabel('Ionosonde hmF2','fontsize',14)
hold off
% Robust fit, accounting for outliers


    if length(sza_gt_100)>2
        mdlr_night = fitlm(reshaped_ISR_all(sza_gt_100),reshaped_BradDud_all(sza_gt_100), 'RobustOpts','on');
        fitvalues_iono_corr_vs_ISR_night = table2array(mdlr_night.Coefficients);
        b_night = fitvalues_iono_corr_vs_ISR_night(2,1);
        db_night = fitvalues_iono_corr_vs_ISR_night(2,2);
        a_night = fitvalues_iono_corr_vs_ISR_night(1,1);
        da_night = fitvalues_iono_corr_vs_ISR_night(1,2);
%%
        rsquare_robustfit_night = corr(reshaped_BradDud_all(sza_gt_100)',(a_night+b_night*reshaped_ISR_all(sza_gt_100))')^2;
%%
        
    else
        b_night = NaN;
        db_night = NaN;
        a_night = NaN;
        da_night = NaN;
    end
% subplot(sp222)
    hold on
    plot([50 500], b_night*[50 500]+a_night,'k.-')
    text(60,480,['y=(' num2str(b_night,'%2.2f') '\pm' num2str(db_night,'%2.2f') ')x+(' num2str(a_night,'%2.2f') '\pm' num2str(da_night,'%2.2f') ')'],'fontsize',6);
    % xlabel('ISR hmF2','fontsize',14)
    % ylabel('Ionosonde hmF2','fontsize',14)
    hold off

% subplot(224)
% sp224 = subplot('position',[0.55 0.05 0.4 0.4]);
nexttile
errorbar(reshaped_ISR_all(sza_ge90_le100), reshaped_BradDud_all(sza_ge90_le100), dreshaped_ISR_all(sza_ge90_le100), dreshaped_ISR_all(sza_ge90_le100), dreshaped_BradDud_all(sza_ge90_le100), dreshaped_BradDud_all(sza_ge90_le100),'k.','capsize',0,'color',[0.5 0.5 0.5]);
hold on
plot(reshaped_ISR_all(sza_ge90_le100), reshaped_BradDud_all(sza_ge90_le100),'m.','markersize',6)
axis([50 500 50 500])
axis('square')
title('Twilight','fontsize',14)
xlabel('ISR hmF2','fontsize',14)
hold off
% Robust fit, accounting for outliers


    if length(sza_ge90_le100)>2
        mdlr_twilight = fitlm(reshaped_ISR_all(sza_ge90_le100),reshaped_BradDud_all(sza_ge90_le100), 'RobustOpts','on');
        fitvalues_iono_corr_vs_ISR_twilight = table2array(mdlr_twilight.Coefficients);
        b_twilight = fitvalues_iono_corr_vs_ISR_twilight(2,1);
        db_twilight = fitvalues_iono_corr_vs_ISR_twilight(2,2);
        a_twilight = fitvalues_iono_corr_vs_ISR_twilight(1,1);
        da_twilight = fitvalues_iono_corr_vs_ISR_twilight(1,2);
%%
        rsquare_robustfit_twilight = corr(reshaped_BradDud_all(sza_ge90_le100)',(a_twilight+b_twilight*reshaped_ISR_all(sza_ge90_le100))')^2;
%%
    else
        b_twilight = NaN;
        db_twilight = NaN;
        a_twilight = NaN;
        da_twilight = NaN;
    end
% subplot(sp222)
    hold on
    plot([50 500], b_twilight*[50 500]+a_twilight,'k.-')
    text(60,480,['y=(' num2str(b_twilight,'%2.2f') '\pm' num2str(db_twilight,'%2.2f') ')x+(' num2str(a_twilight,'%2.2f') '\pm' num2str(da_twilight,'%2.2f') ')'],'fontsize',6);
    % xlabel('ISR hmF2','fontsize',14)
    % ylabel('Ionosonde hmF2','fontsize',14)
    hold off



% ylabel('Ionosonde hmF2','fontsize',14)

% plot(reshaped_ISR_all(sza_gt_100),reshaped_BradDud_all(sza_gt_100),'c.')
% plot(reshaped_ISR_all(sza_ge90_le100),reshaped_BradDud_all(sza_ge90_le100),'m.')

% axis([50 500 50 500])
% axis('square')
% plot([50 500], b_all*[50 500]+a_all,'r')
% text(60,480,['All y=(' num2str(b_all,'%2.2f') ' \pm ' num2str(db_all,'%2.2f') ')x + (' num2str(a_all,'%2.2f') ' \pm ' num2str(da_all,'%2.2f') ')']);
% xlabel('ISR hmF2','fontsize',14)
% ylabel('Ionosonde hmF2','fontsize',14)
% hold off


% Now repeat fit without the nighttime values where foE is set to 0.4

% daytime_data = find(reshaped_foE_all ~= 0.4 & isnan(reshaped_foE_all)==0);
% daytime_data = find(reshaped_sza < 90*pi/180);



print -djpeg figure02_revised.jpg

%%

%% Plot annual ionosonde against ISR
% Plot not used in paper but it is interesting to see how the annual
% populations contribute to the distribution as a whole
% The fit parameters are subsequently used to reconstruct the annual fit
% error for an example height of 250 km

[m,n,p] = size(matched_ISR);

% Daytime values only used
matched_foE = median_foE(sub_iono,:,:);

for y = 1:length(ISR_years)

    %% convert to vectors
    reshaped_BradDud = reshape(matched_BradDud(y,:,:),1,n*p);
    dreshaped_BradDud = reshape(dmatched_BradDud(y,:,:),1,n*p);
    reshaped_ISR = reshape(matched_ISR(y,:,:),1,n*p);
    dreshaped_ISR = reshape(dmatched_ISR(y,:,:),1,n*p);
    reshaped_foE = reshape(matched_foE(y,:,:),1,n*p);

    %% Throw out NaNs and nighttime values before fit
    bad = find( isnan(reshaped_BradDud+reshaped_ISR)==1 | reshaped_foE == 0.4);
    
    reshaped_BradDud(bad) = [];
    reshaped_ISR(bad) = [];
    dreshaped_BradDud(bad) = [];
    dreshaped_ISR(bad) = [];

    % Robust fit, accounting for outliers
    if length(reshaped_ISR) > 2
        mdlr = fitlm(reshaped_ISR,reshaped_BradDud, 'RobustOpts','on');
        fitvalues_iono_corr_vs_ISR = table2array(mdlr.Coefficients);
        b_yearly(y) = fitvalues_iono_corr_vs_ISR(2,1);
        db_yearly(y) = fitvalues_iono_corr_vs_ISR(2,2);
        a_yearly(y) = fitvalues_iono_corr_vs_ISR(1,1);
        da_yearly(y) = fitvalues_iono_corr_vs_ISR(1,2);
        p_a_yearly(y) = fitvalues_iono_corr_vs_ISR(1,4);
        p_b_yearly(y) = fitvalues_iono_corr_vs_ISR(2,4);
    else
        b_yearly(y) = NaN;
        db_yearly(y) = NaN;
        a_yearly(y) = NaN;
        da_yearly(y) = NaN;
        p_a_yearly(y) = NaN;
        p_b_yearly(y) = NaN;
    end

    % Plot cumulative plot of each annual fit ( not used in paper but this demonstrates the variability in fits parameters over the 35 years of data.
    figure(106)
    errorbar(reshaped_ISR, reshaped_BradDud, dreshaped_ISR, dreshaped_ISR, dreshaped_BradDud, dreshaped_BradDud,'k.','capsize',0,'color',[0.5 0.5 0.5])
    hold on
    plot(reshaped_ISR, reshaped_BradDud,'k.','markersize',6)
    axis([50 500 50 500])
    axis('square')
    plot([50 500], b_yearly(y)*[50 500]+a_yearly(y),'r')
    text(60,450,['y=(' num2str(b_yearly(y)) ' \pm ' num2str(db_yearly(y)) ')x + (' num2str(a_yearly(y)) ' \pm ' num2str(da_yearly(y)),')']);
    xlabel('ISR hmF2','fontsize',14)
    ylabel('Ionosonde hmF2','fontsize',14)
    % hold off
    pause(0.1)
end

% Create composition proxy array from foF2.^2 at local noon (array index 13) scaled by f10.7 index  
composition_proxy = (reshape(median_foF2(30:end,:,13),35,12).^2)./matched_f107;

% reconstruct model values for a known height and see how they compare with
% possible proxies for the innaccuracy
%
annual_hmF2_fit_for_250_km = b_yearly*250+a_yearly;
dmx = (b_yearly*250).*sqrt( (db_yearly./b_yearly).^2);
dannnual_hmF2_fit_for_250_km = sqrt(dmx.^2 + da_yearly.^2);

% Create annual am index
for i=1:length(am_u_year)
    am_subset = find(am_year == am_u_year(i) & isnan(am)==0);
    mean_annual_am(i) = nanmean(am(am_subset));
    dmean_annual_am(i) = nanstd(am(am_subset))./sqrt(length(am_subset));
end

% remove nighttime foE values before calculating annual median.
badfoE = find(median_foE == 0.4);
median_foE(badfoE) = NaN;

% Create annual F2/F1 ratio and F2/E ratio
F2F1_ratio = median_foF2./median_foF1;
F2E_ratio = median_foF2./median_foE;
for i=1:length(u_year)
    F2F1_subset = find(year_index_iono == u_year(i) & isnan(median_foF2)==0 & isnan(median_foF1) == 0);
    mean_annual_F2F1_ratio(i) = nanmean(F2F1_ratio(F2F1_subset));
    dmean_annual_F2F1_ratio(i) = nanstd(F2F1_ratio(F2F1_subset))./sqrt(length(F2F1_subset));

    F2E_subset = find(year_index_iono == u_year(i) & isnan(median_foF2)==0 & isnan(median_foE) == 0);
    % mean_annual_F2E_ratio(i) = nanmean(F2E_ratio(F2E_subset));
    mean_annual_F2E_ratio(i) =   mean(F2E_ratio(F2E_subset),'omitnan');
    dmean_annual_F2E_ratio(i) = nanstd(F2E_ratio(F2E_subset))./sqrt(length(F2E_subset));
end

% am data start in 1959 so 1986 is am_u_year(28) and 2020 is am_u_year(end-2)

%%

% Plot four parameters for visual comparison (figure 8 in paper)
figure(109)
subplot(221)
    ydata = 100*(annual_hmF2_fit_for_250_km-250)/250;
    dydata = ydata.*(dannnual_hmF2_fit_for_250_km./annual_hmF2_fit_for_250_km);
    errorbar(years, ydata,dydata,'k','capsize',0)
    xlabel('Year')
    ylabel('\epsilon hmF2 (%)')
    title('Model error for 250 km')
    text(1981,16,'a)')
subplot(222)
    plot(years, nanmean(composition_proxy,2),'k')
    xlabel('Year')
    ylabel('composition proxy')
    title('Composition Proxy')
    text(1981,0.76,'b)')
subplot(223)
    errorbar(years, mean_annual_am(28:end-2), dmean_annual_am(28:end-2),'k','capsize',0)
    xlabel('Year')
    ylabel('am')
    title('annual mean am index')
    text(1981,36,'c)')
subplot(224)
    yyaxis left
    set(gca,'ycolor',[0 0 0]);
    errorbar(years, mean_annual_F2F1_ratio(30:end), dmean_annual_F2F1_ratio(30:end),'k','capsize',0) 
    ylabel('F2F1 Ratio','color',[0 0 0])
    hold on
    text(1981,1.77,'d)')
    % hold off
    yyaxis right
    set(gca,'ycolor',[0 0 0]);
    % hold on
    errorbar(years, mean_annual_F2E_ratio(30:end), dmean_annual_F2E_ratio(30:end),'k:','capsize',0)
    xlabel('Year')
    ylabel('F2E Ratio','color',[0 0 0])
    title('annual mean F2/E & F2/F1 ratios')
    hold off
print -djpeg fit_err_vs_time_cf_proxies_updated.jpg

%%

% Plot model error and plots of fits between model error and other three parameters
% This plot is not used in the paper but the correlations are quoted
figure(110)
subplot(221)
    errorbar(years, ydata,dydata,'k','capsize',0)
    xlabel('years')
    ylabel('\epsilon hmF2 (%)')
    title('Model error for 250 km vs year')
subplot(222)
    hold off
    errorbar(nanmean(composition_proxy,2), ydata, dydata,'k.', 'capsize',0)
    ylabel('\epsilon hmF2 (%)')
    xlabel('composition proxy')
    title('\delta hmF2 vs Composition Proxy')
    
    [corcoef_compproxy, pval_compproxy] = corrcoef(nanmean(composition_proxy,2), ydata)
    
     % fit line to data
     mdlr_comp = fitlm(nanmean(composition_proxy,2),ydata, 'RobustOpts','on');
     fitvalues_iono_corr_vs_ISR = table2array(mdlr_comp.Coefficients);
     b_comp = fitvalues_iono_corr_vs_ISR(2,1);
     db_comp = fitvalues_iono_corr_vs_ISR(2,2);
     a_comp = fitvalues_iono_corr_vs_ISR(1,1);
     da_comp = fitvalues_iono_corr_vs_ISR(1,2);
    
     hold on
     
     xvals = get(gca,'xlim');
     plot(xvals,b_comp*xvals+a_comp,'r')
     hold off
     
subplot(223)
    errorbar(mean_annual_am(28:end-2), ydata,-dydata,dydata,-dmean_annual_am(28:end-2),dmean_annual_am(28:end-2) ,'k.','capsize',0)
    ylabel('\epsilon hmF2 (%)')
    xlabel('am')
    title('\delta hmF2 vs annual am index')

    [corcoeff_am, pval_am] = corrcoef(mean_annual_am(28:end-2), ydata)
    
     % fit line to data
     mdlr_am = fitlm(mean_annual_am(28:end-2),ydata, 'RobustOpts','on');
     fitvalues_iono_corr_vs_ISR = table2array(mdlr_am.Coefficients);
     b_am = fitvalues_iono_corr_vs_ISR(2,1);
     db_am = fitvalues_iono_corr_vs_ISR(2,2);
     a_am = fitvalues_iono_corr_vs_ISR(1,1);
     da_am = fitvalues_iono_corr_vs_ISR(1,2);
    
     hold on
     
     xvals = get(gca,'xlim');
     plot(xvals,b_am*xvals+a_am,'r')
     hold off

    
subplot(224)
    errorbar(mean_annual_F2F1_ratio(30:end), ydata, -dydata, dydata, -dmean_annual_F2F1_ratio(30:end), dmean_annual_F2F1_ratio(30:end),'k.','capsize',0) 
    xlabel('\epsilon hmF2 (%)')
    xlabel('foF2/foF1 ratio')
    title('\delta hmF2 vs annual mean foF2/foF1')

    [corcoeff_F2F1, pval_F2F1] = corrcoef(mean_annual_F2F1_ratio(30:end), ydata)
    [corcoeff_F2E, pval_F2E] = corrcoef(mean_annual_F2E_ratio(30:end), ydata)
    
     % fit line to F2F1 data
     mdlr_F2F1 = fitlm(mean_annual_F2F1_ratio(30:end), ydata, 'RobustOpts','on');
     fitvalues_iono_corr_vs_ISR = table2array(mdlr_F2F1.Coefficients);
     b_F2F1 = fitvalues_iono_corr_vs_ISR(2,1);
     db_F2F1 = fitvalues_iono_corr_vs_ISR(2,2);
     a_F2F1 = fitvalues_iono_corr_vs_ISR(1,1);
     da_F2F1 = fitvalues_iono_corr_vs_ISR(1,2);
    
     hold on
     
     xvals = get(gca,'xlim');
     plot(xvals,b_F2F1*xvals+a_F2F1,'r')
     hold off
        
print -djpeg fit_err_vs_proxies.jpg