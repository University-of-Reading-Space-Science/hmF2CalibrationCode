bad_foE = find(median_foE == 0.4);
median_foE(bad_foE) = NaN;

scaling_factors = [0.5:0.1:2.5];

for scaling_loop = 1:length(scaling_factors)
for i=1:length(u_year)
    for j = 1:12
        for k=1:24
          % year_index_iono(i,j,k) = u_year(i);
          % month_index_iono(i,j,k) = j;
          % hour_index_iono(i,j,k) = k-1;
          % datenum_index_iono(i,j,k) = datenum(year_index_iono(i,j,k),month_index_iono(i,j,k),15,hour_index_iono(i,j,k),0,0); % For the purposes of plot assume day 15 of each month as central.
          % subset = find(year == u_year(i) & month==j & hour == (k-1));
          % [median_foF2(i,j,k), dmedian_foF2(i,j,k)] = median_bootstrap(foF2(subset),minN);
          % [median_foE(i,j,k), dmedian_foE(i,j,k)] = median_bootstrap(foE(subset),minN);
          % [median_foF1(i,j,k), dmedian_foF1(i,j,k)] = median_bootstrap(foF1(subset),minN);        
          % [median_M3000F2(i,j,k), dmedian_M3000F2(i,j,k)] = median_bootstrap(M3000F2(subset),minN);
          
         
          % Calculate hmF2 using Bradley-Dudeney formula - but with x
          % scaled by a particular factor.
          [hmF2_BradDud_scaled(i,j,k), dhmF2_BradDud_scaled(i,j,k)] = hmF2BradDud_scaled(median_foF2(i,j,k),median_foE(i,j,k), median_M3000F2(i,j,k)/10,dmedian_foF2(i,j,k),dmedian_foE(i,j,k), dmedian_M3000F2(i,j,k)/10, scaling_factors(scaling_loop) );
        end
    end
end

% now plot point versus ISR values
figure(12345)
hold off
plot(reshape(matched_ISR,1,[]), reshape(hmF2_BradDud_scaled(30:end,:,:),1,[]),'.')

% Robust fit, accounting for outliers
mdlr = fitlm(reshape(matched_ISR,1,[]), reshape(hmF2_BradDud_scaled(30:end,:,:),1,[]), 'RobustOpts','on');
fitvalues_iono_corr_vs_ISR = table2array(mdlr.Coefficients);
b_all(scaling_loop) = fitvalues_iono_corr_vs_ISR(2,1);
db_all(scaling_loop) = fitvalues_iono_corr_vs_ISR(2,2);
a_all(scaling_loop) = fitvalues_iono_corr_vs_ISR(1,1);
da_all(scaling_loop) = fitvalues_iono_corr_vs_ISR(1,2);

hold on
plot([0 500], b_all(scaling_loop)*[0 500]+a_all(scaling_loop))

pause(0.1)
hold off
end