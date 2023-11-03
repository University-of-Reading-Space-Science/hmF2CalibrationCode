close 700
figure(700)
subplot(211)
plot(reshape(datenum_index_iono,1,[]),reshape(hmF2_BradDud,1,[]),'k.')
axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
set(gca,'xtick',[datenum(1955:5:2025,1,1)])
set(gca,'xticklabel',[1955:5:2025])
% set(gca,'xticklabel',[])
figure(700)
title('Kokubunji ionosonde 1957-2020','fontsize',14)
subplot(211)
ylabel('hmF2 (km)','fontsize',14')

subplot(212)
plot(reshape(datenum_index_iono(30:end,:,:),1,[]),reshape(matched_ISR,1,[]),'k.')
axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
set(gca,'xtick',[datenum(1955:5:2025,1,1)])
set(gca,'xticklabel',[1955:5:2025])
title('MU radar 1986-2020','fontsize',14)
ylabel('hmF2 (km)','fontsize',14')
xlabel('Year','fontsize',14)

orient('landscape')

print -djpeg figure01.jpg