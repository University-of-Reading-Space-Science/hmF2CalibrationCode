mean_hmF2_ISR = NaN*ones(12,24);
dmean_hmF2_ISR = NaN*ones(12,24);
mean_hmF2_ISR_count = NaN*ones(12,24);

for j=1:12
  for k=1:24

    z_in = reshape(matched_ISR(:,j,k),1,[])';

    fit_date = datenum(year_index,month_index, 15*ones(size(year_index)),zeros(size(year_index)),zeros(size(year_index)),zeros(size(year_index)));
    reshape_fit_date = reshape(fit_date,1,[]);

    % estimate average height of subset over all years
    mean_hmF2_ISR(j,k) = nanmean(z_in);
    dmean_hmF2_ISR(j,k) = nanstd(z_in)./sqrt((sum(~isnan(z_in))));  
    mean_hmF2_ISR_count(j,k) = sum(~isnan(z_in));
  end
end


figure(4)
% This creates figure 3b in the paper

tmpIm = imagesc(1:12,0:23,mean_hmF2_ISR');
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('ISR')

print -djpeg ISR_diurnal_seasonal_mean_hmF2.jpg



