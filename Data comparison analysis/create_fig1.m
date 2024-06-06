% close 700
% figure(700)
% subplot(211)
% plot(reshape(datenum_index_iono,1,[]),reshape(hmF2_BradDud,1,[]),'k.')
% axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
% set(gca,'xtick',[datenum(1955:5:2025,1,1)])
% set(gca,'xticklabel',[1955:5:2025])
% % set(gca,'xticklabel',[])
% figure(700)
% title('Kokubunji ionosonde 1957-2020','fontsize',14)
% subplot(211)
% ylabel('hmF2 (km)','fontsize',14')
% 
% subplot(212)
% plot(reshape(datenum_index_iono(30:end,:,:),1,[]),reshape(matched_ISR,1,[]),'k.')
% axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
% set(gca,'xtick',[datenum(1955:5:2025,1,1)])
% set(gca,'xticklabel',[1955:5:2025])
% title('MU radar 1986-2020','fontsize',14)
% ylabel('hmF2 (km)','fontsize',14')
% xlabel('Year','fontsize',14)
% 
% orient('landscape')
% 
% print -djpeg figure01.jpg


% Plot revised figure 1 (after ref #2 comment)

figure(700)

subplot(311)
plot(reshape(permute(datenum_index_iono,[3 2 1]),1,[]), reshape(permute(hmF2_BradDud,[3 2 1]),1,[]),'k.') %,'markersize',5)
axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
set(gca,'xtick',[datenum(1955:5:2025,1,1)])
set(gca,'xticklabel',[1955:5:2025])
ylabel('hmF2 (km)','fontsize',14')
% figure(700)
title('Kokubunji ionosonde 1957-2020','fontsize',14)

subplot(312)
hold off
plot(reshape(permute(datenum_index_iono(30:end,:,:),[3 2 1]),1,[]), reshape(permute(median_all_h,[3 2 1]),1,[]),'k.')
axis([datenum(1955,1,1) datenum(2025,1,1), 150 500])
set(gca,'xtick',[datenum(1955:5:2025,1,1)])
set(gca,'xticklabel',[1955:5:2025])
hold on
plot([datenum(1955,1,1) datenum(2025,1,1)], [180 180],'k--')
ylabel('hmF2 (km)','fontsize',14')
title('MU radar 1986-2020','fontsize',14)
hold off

subplot(313),plot(reshape(permute(datenum_index_iono(30:end,:,:),[3 2 1]),1,[]), reshape(permute(hmF2_BradDud(30:end,:,:),[3 2 1]),1,[]) - reshape(permute(median_all_h,[3 2 1]),1,[]),'k.','markersize',5)
axis([datenum(1955,1,1) datenum(2025,1,1), -250 250])
set(gca,'xtick',[datenum(1955:5:2025,1,1)])
set(gca,'xticklabel',[1955:5:2025])
ylabel('\deltahmF2 (km)','fontsize',14')
title('MU radar - Ionosonde 1986-2020','fontsize',14)

orient('landscape')

print -djpeg figure01_revised.jpg

