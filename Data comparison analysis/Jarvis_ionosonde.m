mean_hmF2 = NaN*ones(12,24);
dmean_hmF2 = NaN*ones(12,24);
mean_hmF2_count = NaN*ones(12,24);

for j=1:12
  for k=1:24

    z_in = reshape(matched_BradDud(:,j,k),1,[])';

    fit_date = datenum(year_index,month_index, 15*ones(size(year_index)),zeros(size(year_index)),zeros(size(year_index)),zeros(size(year_index)));
    reshape_fit_date = reshape(fit_date,1,[]);
    good = find(isnan(z_in)==0);

    if length(good)>2
        % estimate average height of subset
        mean_hmF2(j,k) = nanmean(z_in);
        dmean_hmF2(j,k) = nanstd(z_in)./sqrt((sum(~isnan(z_in))));
        mean_hmF2_count(j,k) = sum(~isnan(z_in));
    end
  end
end


figure(5)
% This creates figure 3a in the paper

tmpIm = imagesc(1:12,0:23,mean_hmF2');
set(gca,'YDir','normal')
set(gca,'clim',[200 400])
xlabel('Month')
ylabel('Hour (LT)')
cb = colorbar;
ylabel(cb,'Mean hmF2 (km)','fontsize',12)
title('Ionosonde')

% Indicate which bins lie in the dawn-dusk zones (sza between 100 and 90)
% need to index a single year (assume variation all the same)
sza_annual = reshape(sza(1,:,:),12,24);
hour_index_annual = reshape(hour_index(1,:,:),12,24);
month_index_annual = reshape(month_index(1,:,:),12,24);

for i=1:12
[tmp_val, tmp_index] = min( abs(sza_annual(i,1:12)-(100*pi/180)),[],'linear');
dawn_val_100(i) = tmp_val(1); dawn_index_100(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,1:12)-(90*pi/180)),[],'linear');
dawn_val_90(i) = tmp_val(1); dawn_index_90(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,13:24)-(90*pi/180)),[],'linear');
dusk_val_90(i) = tmp_val(1); dusk_index_90(i) = tmp_index(1);
[tmp_val, tmp_index] = min( abs(sza_annual(i,13:24)-(100*pi/180)),[],'linear');
dusk_val_100(i) = tmp_val(1); dusk_index_100(i) = tmp_index(1);
end
hold on
figure(5),plot(1:12,hour_index_annual(1,12+dusk_index_100),'w')
figure(5),plot(1:12,hour_index_annual(1,12+dusk_index_90),'w:')
figure(5),plot(1:12,hour_index_annual(1,dawn_index_90),'w:')
figure(5),plot(1:12,hour_index_annual(1,dawn_index_100),'w')
hold off

print -djpeg Ionosonde_diurnal_seasonal_mean_hmF2.jpg
