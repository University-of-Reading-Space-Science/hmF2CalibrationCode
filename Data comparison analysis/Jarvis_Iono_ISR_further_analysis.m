F2F1_threshold = 1.6; % ratio below which, formula starts to deviate

F2F1_ratio = median_foF2./median_foF1;
dF2F1_ratio = F2F1_ratio.*sqrt( (dmedian_foF2./median_foF2).^2 + (dmedian_foF1./median_foF1).^2 );

figure(10)
% This creates figure 4 in the paper
diff_hmF2 = (mean_hmF2_ISR - mean_hmF2);
ddiff_hmF2 = sqrt(dmean_hmF2_ISR.^2 + dmean_hmF2.^2);
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

array_foF2 = nanmedian(median_foF2(30:end,:,:));
array_foF1 = nanmedian(median_foF1(30:end,:,:));
array_F2F1_ratio = nanmean(array_foF2./array_foF1,1);
bad_ratio = find( (array_F2F1_ratio) <= F2F1_threshold);

mean_hmF2_corrected = mean_hmF2;
mean_hmF2_corrected(bad_ratio) = mean_hmF2_corrected(bad_ratio) + (b_cal*array_F2F1_ratio(bad_ratio)+a_cal);

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


% Now apply to time sequence (i.e. don't average out over years first)

matched_BradDud_corrected = matched_BradDud;
dmatched_BradDud_corrected = dmatched_BradDud;
bad_ratio = find(F2F1_ratio(30:end,:,:) <= F2F1_threshold);


corr_error = sqrt(   ((b_cal.*F2F1_ratio).*sqrt( (db_cal./b_cal).^2 + (dF2F1_ratio./F2F1_ratio).^2)).^2 + da_cal.^2);

matched_BradDud_corrected(bad_ratio) = matched_BradDud_corrected(bad_ratio) + (b_cal*F2F1_ratio(bad_ratio) + a_cal);
dmatched_BradDud_corrected(bad_ratio) = sqrt(dmatched_BradDud_corrected(bad_ratio).^2 + corr_error(bad_ratio).^2);


figure(12)
% Figure 6 in the paper

% subplot(311)
subplot('position',[0.15 0.66 0.8 0.26])
plot(reshape(year_index_iono+month_index_iono/12,1,[]), reshape(F2F1_ratio,1,[]),'k.')
set(gca,'xlim',[u_year(1) u_year(end)+1])
set(gca,'ylim',[1 2.8])
set(gca,'xticklabel',[])
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
ylabel(['F2/F1 \leq ' num2str(F2F1_threshold) '(%)'],'fontsize',14)
xlabel('year','fontsize',14)
xlabel('Year','fontsize',14)

orient('tall')
print -djpeg F2F1_vs_year_3panel.jpg


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
