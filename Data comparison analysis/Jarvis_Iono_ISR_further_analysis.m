F2F1_threshold = 1.6; % ratio below which, formula starts to deviate

F2E_threshold = 2.5; % Threshold below which, formula starts to deviate.

F2F1_ratio = median_foF2./median_foF1;
dF2F1_ratio = F2F1_ratio.*sqrt( (dmedian_foF2./median_foF2).^2 + (dmedian_foF1./median_foF1).^2 );

figure(10)
% This creates figure 4 in the paper
diff_hmF2 = (mean_hmF2_ISR - mean_hmF2);
ddiff_hmF2 = sqrt(dmean_hmF2_ISR.^2 + dmean_hmF2.^2);
meanF2F1ratio = nanmean(F2F1_ratio(30:end,:,:),1);
dmeanF2F1ratio = reshape(nanstd(dF2F1_ratio(30:end,:,:),0,1),12,24)./reshape(sqrt(nansum(~isnan(dF2F1_ratio(30:end,:,:)),1)),12,24);
% Need to work out the uncertainty on the median F2/F1 ratio over time...
errorbar(reshape(nanmean(F2F1_ratio(30:end,:,:),1),12,24),diff_hmF2, -ddiff_hmF2, ddiff_hmF2,-dmeanF2F1ratio,dmeanF2F1ratio,'.','capsize',0,'color',[0.5 0.5 0.5])  
xlabel('foF2/foF1 ratio')
ylabel('\delta hmF2 (km)')
title('ISR-Ionosonde hmF2 vs foF2/foF1 ratio')

% print -djpeg ISR_minus_Ionosonde_hmF2_vs_foF1_foF2_ratio.jpg

x_variable = reshape(nanmean(F2F1_ratio(30:end,:,:),1),1,[]);
dx_variable = reshape(nanstd(F2F1_ratio(30:end,:,:),[],1)./sqrt(35),1,[]);
y_variable = reshape((mean_hmF2_ISR - mean_hmF2),1,[]);
dy_variable = reshape(sqrt( dmean_hmF2_ISR.^2 +dmean_hmF2.^2),1,[]);
bad_ratio = find( x_variable <=F2F1_threshold);

mdlr_fit = fitlm(x_variable(bad_ratio),y_variable(bad_ratio), 'RobustOpts','on');
fitvalues_cal = table2array(mdlr_fit.Coefficients);
b_cal = fitvalues_cal(2,1);
db_cal = fitvalues_cal(2,2);
a_cal = fitvalues_cal(1,1);
da_cal = fitvalues_cal(1,2);
p_a_cal = fitvalues_cal(1,4);
p_b_cal = fitvalues_cal(2,4);

% overplot fit
hold on 
plot([min(x_variable(bad_ratio)) max(x_variable(bad_ratio))], b_cal*[min(x_variable(bad_ratio)) max(x_variable(bad_ratio))] +a_cal,'r')
hold off

print -djpeg ISR_minus_Ionosonde_hmF2_vs_foF2_foF1_ratio.jpg

% now correct data 

y_variable_cal = y_variable;
y_variable_cal(bad_ratio) = y_variable(bad_ratio) - (b_cal*x_variable(bad_ratio)+a_cal);
dmx = (b_cal.*x_variable(bad_ratio)).*sqrt( (db_cal./b_cal).^2 + (dx_variable(bad_ratio)./x_variable(bad_ratio)).^2);
dy_variable_cal = dy_variable;
dy_variable_cal(bad_ratio) = sqrt(dmx.^2 + (da_cal./a_cal).^2);
  

% Plot calibrated data (not used in the paper)
figure(15)
% plot(x_variable, y_variable_cal,'o')
% hold on
errorbar(x_variable, y_variable_cal, -dy_variable_cal,dy_variable_cal,-dx_variable,dx_variable,'k.','capsize',0,'color',[0.5 0.5 0.5]);
xlabel('foF2/foF1 ratio')
ylabel('ISR-Ionosonde F2 height, Corrected (km)') 
hold off

% Now calibrate Ionosonde hmF2 values
mean_hmF2_cal = reshape(mean_hmF2,1,[]);
dmean_hmF2_cal = reshape(dmean_hmF2,1,[]);
mean_hmF2_cal(bad_ratio) = mean_hmF2_cal(bad_ratio) + (b_cal*x_variable(bad_ratio)+a_cal); 
% dmean_hmF2_cal(bad_ratio) = mean_hmF2_cal(bad_ratio) + (b_cal*x_variable(bad_ratio)+a_cal); 
dmx = (b_cal.*x_variable(bad_ratio)).*sqrt( (db_cal./b_cal).^2 + (dx_variable(bad_ratio)./x_variable(bad_ratio)).^2);
dmean_hmF2_cal = dy_variable;
dmean_hmF2_cal(bad_ratio) = sqrt(dmx.^2 + (da_cal./a_cal).^2);


% OK, now keep shape of arrays and produce a corrected version of the
% diurnal & seasonal chart.

array_foF2 = reshape(nanmedian(median_foF2(30:end,:,:)),12,24);
array_foF1 = reshape(nanmedian(median_foF1(30:end,:,:)),12,24);
% array_F2F1_ratio = nanmean(array_foF2./array_foF1,1);   % Bug found. This
% was identifying nighttime values to be corrected.
array_F2F1_ratio = array_foF2./array_foF1;
bad_ratio = find( (array_F2F1_ratio) <= F2F1_threshold);

mean_hmF2_corrected = mean_hmF2;
mean_hmF2_corrected(bad_ratio) = mean_hmF2_corrected(bad_ratio) + (b_cal*array_F2F1_ratio(bad_ratio)+a_cal);

dmean_hmF2_corrected = dmean_hmF2;
dmx = (b_cal.*x_variable(bad_ratio)).*sqrt( (db_cal./b_cal).^2 + (dx_variable(bad_ratio)./x_variable(bad_ratio)).^2);
dmean_hmF2_corrected(bad_ratio) = sqrt(dmean_hmF2(bad_ratio).^2 + dmx'.^2);


% Create revised version of figure 3 in the paper (after comment from Ref
% #2)

%%

% Indicate which bins lie in the dawn-dusk zones (sza between 100 and 90)
% need to index a single year (assume variation all the same)
sza_annual = reshape(sza(1,:,:),12,24);
hour_index_annual = reshape(hour_index(1,:,:),12,24);
month_index_annual = reshape(month_index(1,:,:),12,24);

for i=1:12
[tmp_val, tmp_index] = min( abs(sza_annual(i,1:12)-(100*pi/180)),[],'linear');
dawn_val_100(i) = tmp_val(1); 
dawn_index_100(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,1:12)-(90*pi/180)),[],'linear');
dawn_val_90(i) = tmp_val(1); 
dawn_index_90(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,13:24)-(90*pi/180)),[],'linear');
dusk_val_90(i) = tmp_val(1); 
dusk_index_90(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,13:24)-(100*pi/180)),[],'linear');
dusk_val_100(i) = tmp_val(1); 
dusk_index_100(i) = tmp_index(1);
end

figure(3)
tiledlayout(2,2)

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off

set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('Ionosonde ')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('ISR')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR' - mean_hmF2');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[-75 75])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'\deltaMean hmF2 (km)','fontsize',12)
title('ISR-Ionosonde')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR_count');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 35])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Years in ISR mean','fontsize',12)
title('ISR Data count')

print -djpeg figure3_revised.jpg

%%

figure(17)
% Figure 5 in the paper
tmpIm = imagesc(1:12,0:23,mean_hmF2_corrected');

set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('Ionosonde Corrected')

print -djpeg Ionosonde_diurnal_seasonal_hmF2_corrected.jpg

%% Revised figure 5 as suggested by referee #2
figure(5)
tiledlayout(2,2)
nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_corrected');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('Ionosonde Corrected')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('ISR')

nexttile(3)
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR' - mean_hmF2_corrected');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[-75 75])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'\deltaMean hmF2 (km)','fontsize',12)
title('ISR - Ionosonde Corrected')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR_count');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 35])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Years in ISR mean','fontsize',12)
title('ISR Data count')

print -djpeg Figure05_revised.jpg

%%



% Now apply to time sequence (i.e. don't average out over years first)

matched_BradDud_corrected = matched_BradDud;
dmatched_BradDud_corrected = dmatched_BradDud;
bad_ratio = find(F2F1_ratio(30:end,:,:) <= F2F1_threshold);


corr_error = sqrt(   ((b_cal.*F2F1_ratio).*sqrt( (db_cal./b_cal).^2 + (dF2F1_ratio./F2F1_ratio).^2)).^2 + da_cal.^2);

matched_BradDud_corrected(bad_ratio) = matched_BradDud_corrected(bad_ratio) + (b_cal*F2F1_ratio(bad_ratio) + a_cal);
dmatched_BradDud_corrected(bad_ratio) = sqrt(dmatched_BradDud_corrected(bad_ratio).^2 + corr_error(bad_ratio).^2);

%%
figure(12)
% Figure 6 in the paper

% subplot(311)
subplot('position',[0.15 0.66 0.8 0.26])
plot(reshape(year_index_iono+month_index_iono/12,1,[]), reshape(F2F1_ratio,1,[]),'k.')
set(gca,'xlim',[u_year(1) u_year(end)+1])
set(gca,'ylim',[1 2.8])
set(gca,'xticklabel',[])
hold on
plot([u_year(1) u_year(end)+1],[F2F1_threshold F2F1_threshold],'k','LineWidth',1);
hold off

% xlabel('Year','fontsize',14)
ylabel('F2/F1','fontsize',14)

%subplot(312)
subplot('position',[0.15 0.38 0.8 0.26])
% y_in = nanmean(nanmean(F2F1_ratio,3),2);
y_in = mean_annual_F2F1_ratio; % calculated in data_matching_and_analysis_BradDud
dy_in = nanstd(F2F1_ratio,[],3)./sqrt(sum(~isnan(nanmean(F2F1_ratio,3))));
dy_in = sqrt(nansum(dy_in.^2,2));
errorbar(1957:2020, y_in, dy_in,'k','capsize',0,'color',[0.5 0.5 0.5])
set(gca,'xlim',[u_year(1) u_year(end)+1])
set(gca,'xticklabel',[])
set(gca,'ylim',[1 2.8])
hold on
plot([u_year(1) u_year(end)+1],[F2F1_threshold F2F1_threshold],'k','LineWidth',1);
hold off

ylabel('foF2/foF1','fontsize',14)

% subplot(313)
subplot('position',[0.15 0.1 0.8 0.26])
for i=1:length(u_year)
    bad = find(F2F1_ratio(i,:,:) <=F2F1_threshold);
    total_number = find(isnan(F2F1_ratio(i,:,:)) == 0);
    bin_fraction(i) = length(bad)/length(total_number);
end

ydata_all(1,:) = bin_fraction;
ydata_all(2,:) = bin_fraction;
xdata_all(1,:) = (u_year(1):u_year(end));
xdata_all(2,:) = (u_year(1):u_year(end))+1;

plot(reshape(xdata_all,1,[]),100*reshape(ydata_all,1,[]),'k')

set(gca,'xlim',[u_year(1) u_year(end)+1],'ylim',[0 100])
ylabel(['% F2/F1 \leq ' num2str(F2F1_threshold)],'fontsize',14)
xlabel('year','fontsize',14)
xlabel('Year','fontsize',14)

orient('tall')
print -djpeg F2F1_vs_year_3panel.jpg
%%
figure(121)
% additional figure to test F2/E ratio change with time (similar to figure
% 12 - figure 6 in the original submitted paper)

% subplot(311)
subplot('position',[0.15 0.66 0.8 0.26])
plot(reshape(year_index_iono+month_index_iono/12,1,[]), reshape(F2E_ratio,1,[]),'k.')
hold on
plot([u_year(1) u_year(end)+1],[F2E_threshold F2E_threshold],'k','LineWidth',1);
hold off
set(gca,'xlim',[u_year(1) u_year(end)+1])
set(gca,'ylim',[1 5.0])
set(gca,'xticklabel',[])
% xlabel('Year','fontsize',14)
ylabel('F2/E','fontsize',14)

%subplot(312)
subplot('position',[0.15 0.38 0.8 0.26])
% y_in = nanmean(nanmean(F2F1_ratio,3),2);
y_in = mean_annual_F2E_ratio; % calculated in data_matching_and_analysis_BradDud
dy_in = nanstd(F2E_ratio,[],3)./sqrt(sum(~isnan(nanmean(F2E_ratio,3))));
dy_in = sqrt(nansum(dy_in.^2,2));
errorbar(1957:2020, y_in, dy_in,'k','capsize',0,'color',[0.5 0.5 0.5])
set(gca,'xlim',[u_year(1) u_year(end)+1])
set(gca,'xticklabel',[])
set(gca,'ylim',[1 5.0])

hold on
plot([u_year(1) u_year(end)+1],[F2E_threshold F2E_threshold],'k','LineWidth',1);
hold off

ylabel('foF2/foE','fontsize',14)

% subplot(313)
subplot('position',[0.15 0.1 0.8 0.26])
for i=1:length(u_year)
    bad = find(F2E_ratio(i,:,:) <=F2E_threshold);
    total_number = find(isnan(F2E_ratio(i,:,:)) == 0);
    bin_fraction(i) = length(bad)/length(total_number);
end

ydata_all(1,:) = bin_fraction;
ydata_all(2,:) = bin_fraction;
xdata_all(1,:) = (u_year(1):u_year(end));
xdata_all(2,:) = (u_year(1):u_year(end))+1;

plot(reshape(xdata_all,1,[]),100*reshape(ydata_all,1,[]),'k')
set(gca,'xlim',[u_year(1) u_year(end)+1],'ylim',[0 100])
ylabel(['% F2/E \leq ' num2str(F2E_threshold)],'fontsize',14)
xlabel('year','fontsize',14)
xlabel('Year','fontsize',14)

orient('tall')
print -djpeg F2E_vs_year_3panel.jpg

%%
figure(22)
% Uncorrected daytime data only - not used in the paper
daytime_subset = find(sza <(90*pi/180));
y_variable = matched_BradDud_corrected;
dy_variable = dmatched_BradDud_corrected;
hold off
errorbar(reshape(matched_ISR(daytime_subset),1,[]), reshape(y_variable(daytime_subset),1,[]),-reshape(dy_variable(daytime_subset),1,[]),reshape(dy_variable(daytime_subset),1,[]),-reshape(dmatched_ISR(daytime_subset),1,[]),reshape(dmatched_ISR(daytime_subset),1,[]),'k.','capsize',0,'color',[0.5 0.5 0.5])
hold on
plot(reshape(matched_ISR(daytime_subset),1,[]), reshape(y_variable(daytime_subset),1,[]),'k.')
hold off
mdlr_fit_day = fitlm(reshape(matched_ISR(daytime_subset),1,[]),reshape(y_variable(daytime_subset),1,[]), 'RobustOpts','on');
fitvalues_cal_day = table2array(mdlr_fit_day.Coefficients);
b_cal_day = fitvalues_cal_day(2,1);
db_cal_day = fitvalues_cal_day(2,2);
a_cal_day = fitvalues_cal_day(1,1);
da_cal_day = fitvalues_cal_day(1,2);
p_a_cal_day = fitvalues_cal_day(1,4);
p_b_cal_day = fitvalues_cal_day(2,4);
hold on
axis([0 600 0 600])
axis('square')
plot([100 500],[100 500]*b_cal_day + a_cal_day,'r')
text(60,550,['Day y=(' num2str(b_cal_day,'%2.2f') ' \pm ' num2str(db_cal_day,'%2.2f') ')x + (' num2str(a_cal_day,'%2.2f') ' \pm ' num2str(da_cal_day,'%2.2f') ')']);
hold off
xlabel('ISR hmF2 (km)','fontsize',14)
ylabel('Corrected daytime Ionosonde hmF2 (km)','fontsize',14)

print -djpeg Ionosonde_corrected_vs_ISR_hmF2_daytime.jpg

figure(23)
% Figure 7 in the paper
daytime_subset = find(sza <(90*pi/180));
y_variable = matched_BradDud;
dy_variable = dmatched_BradDud;
hold off
% plot(reshape(matched_ISR,1,[]), reshape(hmF2_BradDud_corrected(30:end,:,:),1,[]),'.')
errorbar(reshape(matched_ISR(daytime_subset),1,[]), reshape(y_variable(daytime_subset),1,[]),-reshape(dy_variable(daytime_subset),1,[]),reshape(dy_variable(daytime_subset),1,[]),-reshape(dmatched_ISR(daytime_subset),1,[]),reshape(dmatched_ISR(daytime_subset),1,[]),'k.','capsize',0,'color',[0.5 0.5 0.5])
hold on
plot(reshape(matched_ISR(daytime_subset),1,[]), reshape(y_variable(daytime_subset),1,[]),'k.')
hold off
mdlr_fit_day = fitlm(reshape(matched_ISR(daytime_subset),1,[]),reshape(y_variable(daytime_subset),1,[]), 'RobustOpts','on');
fitvalues_cal_day = table2array(mdlr_fit_day.Coefficients);
b_day = fitvalues_cal_day(2,1);
db_day = fitvalues_cal_day(2,2);
a_day = fitvalues_cal_day(1,1);
da_day = fitvalues_cal_day(1,2);
p_a_day = fitvalues_cal_day(1,4);
p_b_day = fitvalues_cal_day(2,4);
hold on
axis([0 600 0 600])
axis('square')
plot([100 500],[100 500]*b_day + a_day,'r')
text(60,550,['Day y=(' num2str(b_day,'%2.2f') ' \pm ' num2str(db_day,'%2.2f') ')x + (' num2str(a_day,'%2.2f') ' \pm ' num2str(da_day,'%2.2f') ')']);
hold off
xlabel('ISR hmF2 (km)','fontsize',14)
ylabel('Daytime Ionosonde hmF2 (km)','fontsize',14)

print -djpeg Ionosonde_vs_ISR_hmF2_daytime.jpg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Further plots to answer referees comments

%%

% calculate ratio with night values
median_foE_with_night = matched_foE;
dmedian_foE_with_night = dmedian_foE(30:end,:,:);
matched_foF2 = median_foF2(30:end,:,:);
dmatched_foF2 = dmedian_foF2(30:end,:,:);
% set error of nighttime values to +/- 0.4 MHz
baderrors = find(isnan(dmedian_foE_with_night) == 1);
dmedian_foE_with_night(baderrors) = 0.4;
% night_values = find(sza < 100*pi/180);
matched_F2E_ratio_with_night = matched_foF2./median_foE_with_night;
dmatched_F2E_ratio_with_night = matched_F2E_ratio_with_night.*sqrt( (dmatched_foF2./matched_foF2).^2 + (dmedian_foE_with_night./median_foE_with_night).^2 );

% calculate ratio ignoring nighttime foE=0.4;
F2E_ratio = median_foF2./median_foE;
bad = find(median_foE==0.4); % identify artificial night-time values
F2E_ratio(bad) = NaN;
dF2E_ratio = F2E_ratio.*sqrt( (dmedian_foF2./median_foF2).^2 + (dmedian_foE./median_foE).^2 );
matched_F2E_ratio = F2E_ratio(30:end,:,:);
dmatched_F2E_ratio = dF2E_ratio(30:end,:,:);

%%

figure(10421)
hold off
% This creates equivalent to figure 4 in the paper but uses F2E to sort
% the differences between ISR and ionosonde derived hmF2 values
% This version includes all the night-time estimates of foE=0.4 MHz
diff_hmF2 = (mean_hmF2_ISR - mean_hmF2);

ddiff_hmF2 = sqrt(dmean_hmF2_ISR.^2 + dmean_hmF2.^2);
meanF2Eratio_with_night = nanmean(matched_F2E_ratio_with_night,1);
dmeanF2Eratio_with_night = reshape(nanstd(dmatched_F2E_ratio_with_night,0,1),12,24)./reshape(sqrt(nansum(~isnan(dmatched_F2E_ratio_with_night),1)),12,24);
% Need to work out the uncertainty on the median F2/F1 ratio over time...
errorbar(reshape(nanmean(matched_F2E_ratio_with_night,1),12,24),diff_hmF2, -ddiff_hmF2, ddiff_hmF2,-dmeanF2Eratio_with_night,dmeanF2Eratio_with_night,'.','capsize',0,'color',[0.5 0.5 0.5])  

xlabel('foF2/foE ratio (with night-time estimates)')
ylabel('\delta hmF2 (km)')
title('ISR-Ionosonde hmF2 vs foF2/foE ratio with night-time estimates')

% F2E_threshold = 2.5;
% Fit_order = 3;
% x2_variable = reshape(nanmean(matched_F2E_ratio,1),1,[]);
% dx2_variable = reshape(nanstd(matched_F2E_ratio,[],1)./sqrt(35),1,[]);
x2_variable = reshape(meanF2Eratio_with_night,1,[]);
dx2_variable = reshape(dmeanF2Eratio_with_night,1,[]);
y2_variable = reshape(diff_hmF2,1,[]);
dy2_variable = reshape(ddiff_hmF2,1,[]);
bad_F2Eratio = find( x2_variable <= F2E_threshold & isreal(x2_variable) & isreal(y2_variable));

mdlr_fit_F2E_with_night = fitlm(x2_variable(bad_F2Eratio),y2_variable(bad_F2Eratio), 'RobustOpts','on');
fitvalues_cal_F2E_with_night = table2array(mdlr_fit_F2E_with_night.Coefficients);
b_cal_F2E_with_night = fitvalues_cal_F2E_with_night(2,1);
db_cal_F2E_with_night = fitvalues_cal_F2E_with_night(2,2);
a_cal_F2E_with_night = fitvalues_cal_F2E_with_night(1,1);
da_cal_F2E_with_night = fitvalues_cal_F2E_with_night(1,2);
p_a_cal_F2E_with_night = fitvalues_cal_F2E_with_night(1,4);
p_b_cal_F2E_with_night = fitvalues_cal_F2E_with_night(2,4);

% [sorted_x2_variable, x2_sort_index] = sort(x2_variable(bad_F2Eratio));
% [F2E_coeffs, F2E_S] = polyfit(x2_variable(bad_F2Eratio(x2_sort_index)), y2_variable(bad_F2Eratio(x2_sort_index)),Fit_order);
% [fitted_y2, dfitted_y2] = polyval(F2E_coeffs, x2_variable(bad_F2Eratio(x2_sort_index)),F2E_S);

% overplot fit
hold on 
% plot(x2_variable(bad_F2Eratio(x2_sort_index)),fitted_y2,'r')
plot([1.7 2.5], b_cal_F2E_with_night*[1.7 2.5]+a_cal_F2E_with_night,'r')
hold off


%%
% Above shows that large variable values of hmF2 are calculated from
% night-time values of foE assumed to be 0.4 MHz;
% Ratios below 2.5 are unaffected by these large values.
% if we correct for the lower values only, what does it do to the seasonal_diurnal plot?

%%

figure(1042)
hold off
% This creates equivalent to figure 4 in the paper but uses F2E to sort
% the differences between ISR and ionosonde derived hmF2 values
diff_hmF2 = (mean_hmF2_ISR - mean_hmF2);

ddiff_hmF2 = sqrt(dmean_hmF2_ISR.^2 + dmean_hmF2.^2);
meanF2Eratio = nanmean(matched_F2E_ratio,1);
dmeanF2Eratio = reshape(nanstd(dmatched_F2E_ratio,0,1),12,24)./reshape(sqrt(nansum(~isnan(dmatched_F2E_ratio),1)),12,24);
% Need to work out the uncertainty on the median F2/F1 ratio over time...
errorbar(reshape(nanmean(matched_F2E_ratio,1),12,24),diff_hmF2, -ddiff_hmF2, ddiff_hmF2,-dmeanF2Eratio,dmeanF2Eratio,'.','capsize',0,'color',[0.5 0.5 0.5])  

xlabel('foF2/foE ratio')
ylabel('\delta hmF2 (km)')
title('ISR-Ionosonde hmF2 vs foF2/foE ratio')

% F2E_threshold = 2.5;
% Fit_order = 3;
% x2_variable = reshape(nanmean(matched_F2E_ratio,1),1,[]);
% dx2_variable = reshape(nanstd(matched_F2E_ratio,[],1)./sqrt(35),1,[]);
x2_variable = reshape(meanF2Eratio,1,[]);
dx2_variable = reshape(dmeanF2Eratio,1,[]);
y2_variable = reshape(diff_hmF2,1,[]);
dy2_variable = reshape(ddiff_hmF2,1,[]);
bad_F2Eratio = find( x2_variable <= F2E_threshold & isreal(x2_variable) & isreal(y2_variable));

mdlr_fit_F2E = fitlm(x2_variable(bad_F2Eratio),y2_variable(bad_F2Eratio), 'RobustOpts','on');
fitvalues_cal_F2E = table2array(mdlr_fit_F2E.Coefficients);
b_cal_F2E = fitvalues_cal_F2E(2,1);
db_cal_F2E = fitvalues_cal_F2E(2,2);
a_cal_F2E = fitvalues_cal_F2E(1,1);
da_cal_F2E = fitvalues_cal_F2E(1,2);
p_a_cal_F2E = fitvalues_cal_F2E(1,4);
p_b_cal_F2E = fitvalues_cal_F2E(2,4);

% [sorted_x2_variable, x2_sort_index] = sort(x2_variable(bad_F2Eratio));
% [F2E_coeffs, F2E_S] = polyfit(x2_variable(bad_F2Eratio(x2_sort_index)), y2_variable(bad_F2Eratio(x2_sort_index)),Fit_order);
% [fitted_y2, dfitted_y2] = polyval(F2E_coeffs, x2_variable(bad_F2Eratio(x2_sort_index)),F2E_S);

% overplot fit
hold on 
% plot(x2_variable(bad_F2Eratio(x2_sort_index)),fitted_y2,'r')
plot([1.7 2.5], b_cal_F2E*[1.7 2.5]+a_cal_F2E,'r')
hold off

%%



%%

figure(1042),print -djpeg delta_hmF2_vs_F2E_ratio.jpg

% So, there is a bias due to the foF2/foE ratio too. removing this adds lots of noise.
% If the foF1 effect is corrected for first, is there still an foE effect?

figure(1043)
hold off
% This creates equivalent to figure 4 in the paper but uses F2E to sort
% the differences between ISR and ionosonde derived hmF2 values
diff_hmF2_F2F1corrected = (mean_hmF2_ISR - mean_hmF2_corrected);

ddiff_hmF2_F2F1corrected = sqrt(dmean_hmF2_ISR.^2 + dmean_hmF2_corrected.^2);
meanF2Eratio = nanmean(matched_F2E_ratio,1);
dmeanF2Eratio = reshape(nanstd(dmatched_F2E_ratio,0,1),12,24)./reshape(sqrt(nansum(~isnan(dmatched_F2E_ratio),1)),12,24);
% Need to work out the uncertainty on the median F2/F1 ratio over time...
errorbar(reshape(nanmean(matched_F2E_ratio,1),12,24),diff_hmF2_F2F1corrected, -ddiff_hmF2_F2F1corrected, ddiff_hmF2_F2F1corrected,-dmeanF2Eratio,dmeanF2Eratio,'.','capsize',0,'color',[0.5 0.5 0.5])  

xlabel('foF2/foE ratio')
ylabel('\delta hmF2 (F2F1 corrected), km')
title('ISR-Ionosonde hmF2 (F2F1 corrected) vs foF2/foE ratio')

x2b_variable = reshape(meanF2Eratio,1,[]);
dx2b_variable = reshape(dmeanF2Eratio,1,[]);
y2b_variable = reshape(diff_hmF2_F2F1corrected,1,[]);
dy2b_variable = reshape(ddiff_hmF2_F2F1corrected,1,[]);
bad_F2Eratio = find( x2b_variable <= F2E_threshold & isreal(x2b_variable) & isreal(y2b_variable));

mdlr_fit_corr_F2E = fitlm(x2b_variable(bad_F2Eratio),y2b_variable(bad_F2Eratio), 'RobustOpts','on');
fitvalues_cal_corr_F2E = table2array(mdlr_fit_corr_F2E.Coefficients);
b_cal_corr_F2E = fitvalues_cal_corr_F2E(2,1);
db_cal_corr_F2E = fitvalues_cal_corr_F2E(2,2);
a_cal_corr_F2E = fitvalues_cal_corr_F2E(1,1);
da_cal_corr_F2E = fitvalues_cal_corr_F2E(1,2);
p_a_cal_corr_F2E = fitvalues_cal_corr_F2E(1,4);
p_b_cal_corr_F2E = fitvalues_cal_corr_F2E(2,4);

% [sorted_x2_variable, x2_sort_index] = sort(x2_variable(bad_F2Eratio));
% [F2E_coeffs, F2E_S] = polyfit(x2_variable(bad_F2Eratio(x2_sort_index)), y2_variable(bad_F2Eratio(x2_sort_index)),Fit_order);
% [fitted_y2, dfitted_y2] = polyval(F2E_coeffs, x2_variable(bad_F2Eratio(x2_sort_index)),F2E_S);

% overplot fit
hold on 
% plot(x2_variable(bad_F2Eratio(x2_sort_index)),fitted_y2,'r')
plot([1.7 2.5], b_cal_corr_F2E*[1.7 2.5]+a_cal_corr_F2E,'r')
hold off


% Does this correction remove the F2E bias?

F2Ecorrection_factors = zeros(size(x2_variable));
dF2Ecorrection_factors = zeros(size(x2_variable));
bad_F2Eratio = find(x2_variable <= F2E_threshold & x2_variable > 1.7); 

F2Ecorrection_factors(bad_F2Eratio) = b_cal_F2E.*(x2_variable(bad_F2Eratio)) + a_cal_F2E;
dF2Ecorrection_factors(bad_F2Eratio) = sqrt((db_cal_F2E.*x2_variable(bad_F2Eratio)).^2 + da_cal_F2E.^2);

y2_variable_F2Ecorrected = y2_variable - F2Ecorrection_factors;
dy2_variable_F2Ecorrected = sqrt(dy2_variable.^2 + dF2Ecorrection_factors.^2);

% replot corrected values to see if it has removed the bias.
figure(2345)
errorbar(x2_variable,y2_variable_F2Ecorrected, -dy2_variable_F2Ecorrected, dy2_variable_F2Ecorrected,-dx2_variable,dx2_variable,'.','capsize',0,'color',[0.5 0.5 0.5])  

xlabel('foF2/foE ratio')
ylabel('\delta hmF2 (km)')
title('ISR-Ionosonde (F2E corrected) hmF2 vs foF2/foE ratio')

%%
% How does this affect the seasonal/dirunal plot?
% Create revised version of figure 3 for F2 E correction

% First need to apply the corrections to the seasonal/dirunal grid of hmF2
% values.

% OK, now keep shape of arrays and produce a corrected version of the
% diurnal & seasonal chart.

array_foF2 = reshape(nanmedian(median_foF2(30:end,:,:)),12,24);
array_foE = reshape(nanmedian(median_foE_with_night),12,24);
% array_F2F1_ratio = nanmean(array_foF2./array_foF1,1);   % Bug found. This
% was identifying nighttime values to be corrected.
array_F2E_ratio = array_foF2./array_foE;
bad_ratio = find( (array_F2E_ratio) <= F2E_threshold);

mean_hmF2_corrected_F2E = mean_hmF2;
mean_hmF2_corrected_F2E(bad_F2Eratio) = mean_hmF2_corrected_F2E(bad_F2Eratio) + (b_cal_F2E*array_F2E_ratio(bad_F2Eratio)+a_cal_F2E);

dmean_hmF2_corrected_F2E = dmean_hmF2;
dmx = (b_cal.*x_variable(bad_F2Eratio)).*sqrt( (db_cal./b_cal).^2 + (dx_variable(bad_F2Eratio)./x_variable(bad_F2Eratio)).^2);
dmean_hmF2_corrected_F2E(bad_F2Eratio) = sqrt(dmean_hmF2(bad_F2Eratio)'.^2 + dmx'.^2);



figure(3000)
tiledlayout(2,2)

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_corrected_F2E');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('Ionosonde F2E')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('ISR')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR' - mean_hmF2_corrected_F2E');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[-75 75])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'\deltaMean hmF2 (km)','fontsize',12)
title('ISR-Ionosonde F2E')

nexttile
tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR_count');
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 35])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Years in ISR mean','fontsize',12)
title('ISR Data count')

print -djpeg figure3_equivalent_for_F2E.jpg

%%


%%


% Now apply to time sequence (i.e. don't average out over years first)

% matched_BradDud_corrected = matched_BradDud;
% dmatched_BradDud_corrected = dmatched_BradDud;
% matched_F2E_ratio = F2E_ratio(30:end,:,:);
bad_F2Eratio_all = find(matched_F2E_ratio <= F2E_threshold & matched_F2E_ratio > 1.7); 
% Lower limit included here so that corrections are not applied to data where no hmF2 should be calculated according to the Bradley Dudeney formula

F2Ecorrection_factors = zeros(size(matched_BradDud));
dF2Ecorrection_factors = zeros(size(matched_BradDud));
% for cc = 1:length(bad_F2Eratio_all)
  % [F2Ecorrection_factors(bad_F2Eratio_all),dF2Ecorrection_factors(bad_F2Eratio_all)] = polyval(F2E_coeffs,matched_F2E_ratio(bad_F2Eratio_all),F2E_S);
  F2Ecorrection_factors(bad_F2Eratio_all) = b_cal_F2E.*(matched_F2E_ratio(bad_F2Eratio_all)) + a_cal_F2E;
  dF2Ecorrection_factors(bad_F2Eratio_all) = sqrt((db_cal_F2E.*matched_F2E_ratio(bad_F2Eratio_all)).^2 + da_cal_F2E.^2);

  % polyval(F2E_coeffs,matched_F2E_ratio(bad_F2Eratio_all),F2E_S);
% end
matched_BradDud_F2Ecorrected = matched_BradDud + F2Ecorrection_factors;
%dmatched_BradDud_corrected(bad_F2Eratio) = sqrt(dmatched_BradDud_corrected(bad_F2Eratio).^2 + dF2Ecorrection_factors.^2);
dmatched_BradDud_F2Ecorrected = sqrt(dmatched_BradDud.^2 + dF2Ecorrection_factors.^2);

% Now plot data corrected by F2E correction
day_foE = find(matched_foE ~=0.4);
figure(999)
hold off
errorbar(reshape(matched_ISR(day_foE),1,[]), reshape(matched_BradDud_F2Ecorrected(day_foE),1,[]),-reshape(dmatched_ISR(day_foE),1,[]), reshape(dmatched_ISR(day_foE),1,[]), -reshape(dmatched_BradDud_F2Ecorrected(day_foE),1,[]), reshape(dmatched_BradDud_F2Ecorrected(day_foE),1,[]),'.','capsize',0,'color',[0.5 0.5 0.5])
hold on
plot(reshape(matched_ISR(day_foE),1,[]), reshape(matched_BradDud_F2Ecorrected(day_foE),1,[]),'k.')
axis([0 600 0 600])
xlabel('ISR hmF2, km')
ylabel('Ionosonde hmF2 (F2/E corrected), km')

x3_variable = reshape(matched_ISR(day_foE),1,[]);
y3_variable = reshape(matched_BradDud_F2Ecorrected(day_foE),1,[]);
mdlr_fit3 = fitlm(x3_variable,y3_variable, 'RobustOpts','on');
fitvalues_cal3 = table2array(mdlr_fit3.Coefficients);
b_cal3 = fitvalues_cal3(2,1);
db_cal3 = fitvalues_cal3(2,2);
a_cal3 = fitvalues_cal3(1,1);
da_cal3 = fitvalues_cal3(1,2);
p_a_cal3 = fitvalues_cal3(1,4);
p_b_cal3 = fitvalues_cal3(2,4);

hold on
plot([100 500],[100 500]*b_cal3 + a_cal3,'r')
text(60,550,['F2E corrected y=(' num2str(b_cal3,'%2.2f') ' \pm ' num2str(db_cal3,'%2.2f') ')x + (' num2str(a_cal3,'%2.2f') ' \pm ' num2str(da_cal3,'%2.2f') ')']);

print -djpeg ionosonde_F2E_corrected_vs_ISR_hmf2.jpg

% % Does this still leave an effect due to foF1?
% % repeat with corrected values...
% 

y5_variable = y2_variable_F2Ecorrected;
dy5_variable = dy2_variable_F2Ecorrected;
x5_variable = reshape(meanF2F1ratio,1,[]);
dx5_variable = reshape(dmeanF2F1ratio,1,[]);

% % now plot the new figure...
figure(99999)
errorbar(x5_variable,y5_variable,-dy5_variable,dy5_variable,-dx5_variable,dx5_variable,'.','capsize',0,'color',[0.5 0.5 0.5])  
xlabel('foF2/foF1 ratio')
ylabel('\delta hmF2 F2E corrected (km)')
title('ISR-Ionosonde hmF2 (corrected for F2/E bias) vs foF2/foF1 ratio')

%%
% Generate seasonal/diurnal plots of F2E and F2F1 ratios

figure(10023)
subplot(211)
imagesc(array_F2E_ratio')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
colorbar
set(gca,'clim',[0 3])
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean foF2/foE','fontsize',14)

subplot(212)
imagesc(array_F2F1_ratio')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 3])
colorbar
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean foF2/foF1','fontsize',14)

orient tall

print -djpeg seasonal_diurnal_F2E_F2F1_ratios.jpg

%%
% Create seasonal, diurnal plots of individual parameters

figure(10024)
tiledlayout(2,2)

nexttile
imagesc(reshape(nanmean(median_foE(30:end,:,:)),12,24)')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 5])
cb = colorbar;
ylabel(cb,'MHz','fontsize',14);
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean foE','fontsize',14)

nexttile
imagesc(reshape(nanmean(median_foF1(30:end,:,:)),12,24)')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 10])
cb = colorbar;
ylabel(cb,'MHz','fontsize',14);
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean foF1','fontsize',14)

nexttile
imagesc(reshape(nanmean(median_foF2(30:end,:,:)),12,24)')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
set(gca,'clim',[0 15])
cb = colorbar;
ylabel(cb,'MHz','fontsize',14);
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean foF2','fontsize',14)

nexttile
imagesc(reshape(nanmean(median_M3000F2(30:end,:,:)),12,24)')
hold on
plot(1:12,hour_index_annual(1,12+dusk_index_100),'w','linewidth',1.5)
plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_90),'w:','linewidth',1.5)
plot(1:12,hour_index_annual(1,dawn_index_100),'w','linewidth',1.5)
hold off
set(gca,'YDir','normal')
%set(gca,'clim',[0 10])
colorbar
xlabel('Month','fontsize',14)
ylabel('Hour','fontsize',14)
title('Mean M(3000)F2','fontsize',14)

print -djpeg seasonal_diurnal_parameters.jpeg

%%