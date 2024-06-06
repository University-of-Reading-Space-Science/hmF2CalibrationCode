load new_Slough_F2EF1.mat
load Stanley_F2EF1.mat

% Adjust for mis-digitised data points
bad = find(Slough_foF1 < 100);
Slough_foF1(bad) = Slough_foF1(bad)*10;
bad = find(Stanley_foF1 < 100);
Stanley_foF1(bad) = Stanley_foF1(bad)*10;

Slough_datevec = datevec(Slough_datenumber);
Slough_U_year = unique(Slough_datevec(:,1));

for i=1:length(Slough_U_year)
    for j = 1:12
        for k=1:24
            subset = find(Slough_datevec(:,1) == Slough_U_year(i) & Slough_datevec(:,2) == j & Slough_datevec(:,4) == (k-1));
            [Slough_median_foF2(i,j,k), Slough_dmedian_foF2(i,j,k)] = median_bootstrap(Slough_foF2(subset),5);
            [Slough_median_foF1(i,j,k), Slough_dmedian_foF1(i,j,k)] = median_bootstrap(Slough_foF1(subset),5);
            [Slough_median_foE(i,j,k), Slough_dmedian_foE(i,j,k)] = median_bootstrap(Slough_foE(subset),5);
         end
    end
end

Stanley_datevec = datevec(Stanley_datenumber);
Stanley_U_year = unique(Stanley_datevec(:,1));

for i=1:length(Stanley_U_year)
    for j = 1:12
        for k=1:24
            subset = find(Stanley_datevec(:,1) == Stanley_U_year(i) & Stanley_datevec(:,2) == j & Stanley_datevec(:,4) == (k-1));
            [Stanley_median_foF2(i,j,k), Stanley_dmedian_foF2(i,j,k)] = median_bootstrap(Stanley_foF2(subset),5);
            [Stanley_median_foF1(i,j,k), Stanley_dmedian_foF1(i,j,k)] = median_bootstrap(Stanley_foF1(subset),5);
            [Stanley_median_foE(i,j,k), Stanley_dmedian_foE(i,j,k)] = median_bootstrap(Stanley_foE(subset),5);
        end
    end
end

% Convert to MHz
Stanley_median_foF2 = Stanley_median_foF2/10;
Stanley_median_foF1 = Stanley_median_foF1/100;
Stanley_median_foE = Stanley_median_foE/100;

Slough_median_foF2 = Slough_median_foF2/10;
Slough_median_foF1 = Slough_median_foF1/100;
Slough_median_foE = Slough_median_foE/100;


Stanley_F2F1_ratio = (Stanley_median_foF2./Stanley_median_foF1);
Stanley_F2E_ratio = (Stanley_median_foF2./Stanley_median_foE);

Slough_F2F1_ratio =  (Slough_median_foF2./Slough_median_foF1);
Slough_F2E_ratio =  (Slough_median_foF2./Slough_median_foE);


for i=1:length(Stanley_U_year)
      Stanley_annual_mean_F2F1_ratio(i) = nanmean(reshape(Stanley_F2F1_ratio(i,:,:),1,[]));
      dStanley_annual_mean_F2F1_ratio(i) = nanstd(reshape(Stanley_F2F1_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Stanley_F2F1_ratio(i,:,:))==1)));
      Stanley_annual_mean_F2E_ratio(i) = nanmean(reshape(Stanley_F2E_ratio(i,:,:),1,[]));
      dStanley_annual_mean_F2E_ratio(i) = nanstd(reshape(Stanley_F2E_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Stanley_F2E_ratio(i,:,:))==1)));

end

for i=1:length(Slough_U_year)
   Slough_annual_mean_F2F1_ratio(i) = nanmean(reshape(Slough_F2F1_ratio(i,:,:),1,[]));
   dSlough_annual_mean_F2F1_ratio(i) = nanstd(reshape(Slough_F2F1_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Slough_F2F1_ratio(i,:,:))==1)));
   Slough_annual_mean_F2E_ratio(i) = nanmean(reshape(Slough_F2E_ratio(i,:,:),1,[]));
   dSlough_annual_mean_F2E_ratio(i) = nanstd(reshape(Slough_F2E_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Slough_F2E_ratio(i,:,:))==1)));

end

% % Figure 9 in the paper
% figure(2),subplot(211),errorbar(Slough_U_year, Slough_annual_mean_F2F1_ratio,dSlough_annual_mean_F2F1_ratio,'k','capsize',0)
% set(gca,'ylim',[1.2 2.2],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
% subplot(211),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[1.6 1.6],'k-.'),hold off
% ylabel('foF2/foF1','fontsize',14);
% title('Slough/Chilton','fontsize',14)
% figure(2),subplot(212),errorbar(Stanley_U_year, Stanley_annual_mean_F2F1_ratio,dStanley_annual_mean_F2F1_ratio, 'k', 'capsize',0)
% set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)],'ylim',[1.2 2.2]);
% subplot(212),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[1.6 1.6],'k-.'),hold off
% title('Stanley','fontsize',14)
% ylabel('foF2/foF1','fontsize',14);
% xlabel('Year','fontsize',14);
% print -djpeg Slough_Stanley_annual_F2F1_ratio.jpg
% %% 
% 
% 
% % Figure ?? in the paper - after review
% figure(3),subplot(211),errorbar(Slough_U_year, Slough_annual_mean_F2E_ratio,dSlough_annual_mean_F2E_ratio,'k','capsize',0)
% set(gca,'ylim',[1.2 4.0],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
% subplot(211),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[2.5 2.5],'k-.'),hold off
% ylabel('foF2/foE','fontsize',14);
% title('Slough/Chilton','fontsize',14)
% figure(3),subplot(212),errorbar(Stanley_U_year, Stanley_annual_mean_F2E_ratio,dStanley_annual_mean_F2E_ratio, 'k', 'capsize',0)
% set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)],'ylim',[1.2 4.0]);
% subplot(212),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[2.5 2.5],'k-.'),hold off
% title('Stanley','fontsize',14)
% ylabel('foF2/foE','fontsize',14);
% xlabel('Year','fontsize',14);
% print -djpeg Slough_Stanley_annual_F2E_ratio.jpg



%%
% portions of 6 panel plot showing ratios for two different stations
% Figure 13 in the revised paper

F2F1_threshold = 1.6;
F2E_threshold = 2.5;

figure(13)

subplot(421)
errorbar(Slough_U_year, Slough_annual_mean_F2E_ratio,dSlough_annual_mean_F2E_ratio,'k','capsize',0)
set(gca,'ylim',[1.2 4.0],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
hold on 
plot([Slough_U_year(1) Slough_U_year(end)],[F2E_threshold F2E_threshold],'k-.'),hold off
ylabel('foF2/foE','fontsize',14);
title('Slough/Chilton','fontsize',14)

% subplot('position',[0.15 0.1 0.8 0.26])
subplot(423)
for i=1:length(Slough_U_year)
    Slough_bad = find(Slough_F2E_ratio(i,:,:) <=F2E_threshold);
    Slough_total_number = find(isnan(Slough_F2E_ratio(i,:,:)) == 0);
    Slough_bin_fraction(i) = length(Slough_bad)/length(Slough_total_number);
end

Slough_ydata_all(1,:) = Slough_bin_fraction;
Slough_ydata_all(2,:) = Slough_bin_fraction;
Slough_xdata_all(1,:) = (Slough_U_year(1):Slough_U_year(end));
Slough_xdata_all(2,:) = (Slough_U_year(1):Slough_U_year(end))+1;

plot(reshape(Slough_xdata_all,1,[]),100*reshape(Slough_ydata_all,1,[]),'k')

set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)+1],'ylim',[0 100])
ylabel(['% F2/E \leq ' num2str(F2E_threshold)],'fontsize',14)
% xlabel('year','fontsize',14)
% xlabel('Year','fontsize',14)

subplot(425)
errorbar(Slough_U_year, Slough_annual_mean_F2F1_ratio,dSlough_annual_mean_F2F1_ratio,'k','capsize',0)
set(gca,'ylim',[1.2 4.0],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
hold on 
plot([Slough_U_year(1) Slough_U_year(end)],[F2F1_threshold F2F1_threshold],'k-.'),hold off
ylabel('foF2/foF1','fontsize',14);

% subplot('position',[0.15 0.1 0.8 0.26])
subplot(427)
for i=1:length(Slough_U_year)
    Slough_bad = find(Slough_F2F1_ratio(i,:,:) <=F2F1_threshold);
    Slough_total_number = find(isnan(Slough_F2F1_ratio(i,:,:)) == 0);
    Slough_bin_fraction(i) = length(Slough_bad)/length(Slough_total_number);
end

Slough_ydata_all(1,:) = Slough_bin_fraction;
Slough_ydata_all(2,:) = Slough_bin_fraction;
Slough_xdata_all(1,:) = (Slough_U_year(1):Slough_U_year(end));
Slough_xdata_all(2,:) = (Slough_U_year(1):Slough_U_year(end))+1;

plot(reshape(Slough_xdata_all,1,[]),100*reshape(Slough_ydata_all,1,[]),'k')

set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)+1],'ylim',[0 100])
ylabel(['% F2/F1 \leq ' num2str(F2F1_threshold)],'fontsize',14)
xlabel('year','fontsize',14)

% Now the same for Stanley data...
subplot(422)
errorbar(Stanley_U_year, Stanley_annual_mean_F2E_ratio,dStanley_annual_mean_F2E_ratio,'k','capsize',0)
set(gca,'ylim',[1.2 4.0],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
hold on 
plot([Slough_U_year(1) Slough_U_year(end)],[F2E_threshold F2E_threshold],'k-.'),hold off
ylabel('foF2/foE','fontsize',14);
title('Stanley','fontsize',14)

% subplot('position',[0.15 0.1 0.8 0.26])
subplot(424)
for i=1:length(Stanley_U_year)
    Stanley_bad = find(Stanley_F2E_ratio(i,:,:) <=F2E_threshold);
    Stanley_total_number = find(isnan(Stanley_F2E_ratio(i,:,:)) == 0);
    Stanley_bin_fraction(i) = length(Stanley_bad)/length(Stanley_total_number);
end

Stanley_ydata_all(1,:) = Stanley_bin_fraction;
Stanley_ydata_all(2,:) = Stanley_bin_fraction;
Stanley_xdata_all(1,:) = (Stanley_U_year(1):Stanley_U_year(end));
Stanley_xdata_all(2,:) = (Stanley_U_year(1):Stanley_U_year(end))+1;

plot(reshape(Stanley_xdata_all,1,[]),100*reshape(Stanley_ydata_all,1,[]),'k')

set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)+1],'ylim',[0 100])
ylabel(['% F2/E \leq ' num2str(F2E_threshold)],'fontsize',14)
% xlabel('year','fontsize',14)
% xlabel('Year','fontsize',14)

subplot(426)
errorbar(Stanley_U_year, Stanley_annual_mean_F2F1_ratio,dStanley_annual_mean_F2F1_ratio,'k','capsize',0)
set(gca,'ylim',[1.2 4.0],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
hold on 
plot([Stanley_U_year(1) Stanley_U_year(end)],[F2F1_threshold F2F1_threshold],'k-.'),hold off
ylabel('foF2/foF1','fontsize',14);

subplot(428)
for i=1:length(Stanley_U_year)
    Stanley_bad = find(Stanley_F2F1_ratio(i,:,:) <=F2F1_threshold);
    Stanley_total_number = find(isnan(Stanley_F2F1_ratio(i,:,:)) == 0);
    Stanley_bin_fraction(i) = length(Stanley_bad)/length(Stanley_total_number);
end

Stanley_ydata_all(1,:) = Stanley_bin_fraction;
Stanley_ydata_all(2,:) = Stanley_bin_fraction;
Stanley_xdata_all(1,:) = (Stanley_U_year(1):Stanley_U_year(end));
Stanley_xdata_all(2,:) = (Stanley_U_year(1):Stanley_U_year(end))+1;

plot(reshape(Stanley_xdata_all,1,[]),100*reshape(Stanley_ydata_all,1,[]),'k')

set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)+1],'ylim',[0 100])
ylabel(['% F2/F1 \leq ' num2str(F2F1_threshold)],'fontsize',14)
xlabel('year','fontsize',14)

orient('landscape')

print -djpeg Slough_Stanley_annual_F2E_F2F1_ratios.jpg

%%



% figure(4)
% subplot(311),errorbar(Stanley_U_year, Stanley_annual_mean_F2E_ratio,dStanley_annual_mean_F2E_ratio, 'k', 'capsize',0)
% set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)],'ylim',[1.2 4.0]);
% subplot(212),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[2.5 2.5],'k-.'),hold off
% title('Stanley','fontsize',14)
% ylabel('foF2/foE','fontsize',14);
% xlabel('Year','fontsize',14);
% 
% subplot('position',[0.15 0.1 0.8 0.26])
% for i=1:length(u_year)
%     bad = find(F2F1_ratio(i,:,:) <=F2F1_threshold);
%     total_number = find(isnan(F2F1_ratio(i,:,:)) == 0);
%     bin_fraction(i) = length(bad)/length(total_number);
% end
% 
% ydata_all(1,:) = bin_fraction;
% ydata_all(2,:) = bin_fraction;
% xdata_all(1,:) = (u_year(1):u_year(end));
% xdata_all(2,:) = (u_year(1):u_year(end))+1;
% 
% plot(reshape(xdata_all,1,[]),100*reshape(ydata_all,1,[]),'k')
% 
% set(gca,'xlim',[u_year(1) u_year(end)+1],'ylim',[0 100])
% ylabel(['% F2/F1 \leq ' num2str(F2F1_threshold)],'fontsize',14)
% xlabel('year','fontsize',14)
% xlabel('Year','fontsize',14)
% 
% orient('tall')