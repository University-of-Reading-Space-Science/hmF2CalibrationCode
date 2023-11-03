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
        end
    end
end

Stanley_median_foF2 = Stanley_median_foF2/10;
Stanley_median_foF1 = Stanley_median_foF1/100;

Slough_median_foF2 = Slough_median_foF2/10;
Slough_median_foF1 = Slough_median_foF1/100;


Stanley_F2F1_ratio = (Stanley_median_foF2./Stanley_median_foF1);
Slough_F2F1_ratio =  (Slough_median_foF2./Slough_median_foF1);

for i=1:length(Stanley_U_year)
      Stanley_annual_mean_F2F1_ratio(i) = nanmean(reshape(Stanley_F2F1_ratio(i,:,:),1,[]));
      dStanley_annual_mean_F2F1_ratio(i) = nanstd(reshape(Stanley_F2F1_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Stanley_F2F1_ratio(i,:,:))==1)));
end

for i=1:length(Slough_U_year)
   Slough_annual_mean_F2F1_ratio(i) = nanmean(reshape(Slough_F2F1_ratio(i,:,:),1,[]));
   dSlough_annual_mean_F2F1_ratio(i) = nanstd(reshape(Slough_F2F1_ratio(i,:,:),1,[]))./sqrt(length(find(isnan(Slough_F2F1_ratio(i,:,:))==1)));
end

% Figure 9 in the paper
figure(2),subplot(211),errorbar(Slough_U_year, Slough_annual_mean_F2F1_ratio,dSlough_annual_mean_F2F1_ratio,'k','capsize',0)
set(gca,'ylim',[1.2 2.2],'xlim',[Slough_U_year(1) Slough_U_year(end)]);
subplot(211),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[1.6 1.6],'k-.'),hold off
ylabel('foF2/foF1','fontsize',14);
title('Slough/Chilton','fontsize',14)
figure(2),subplot(212),errorbar(Stanley_U_year, Stanley_annual_mean_F2F1_ratio,dStanley_annual_mean_F2F1_ratio, 'k', 'capsize',0)
set(gca,'xlim',[Slough_U_year(1) Slough_U_year(end)],'ylim',[1.2 2.2]);
subplot(212),hold on, plot([Slough_U_year(1) Slough_U_year(end)],[1.6 1.6],'k-.'),hold off
title('Stanley','fontsize',14)
ylabel('foF2/foF1','fontsize',14);
xlabel('Year','fontsize',14);
print -djpeg Slough_Stanley_annual_F2F1_ratio.jpg