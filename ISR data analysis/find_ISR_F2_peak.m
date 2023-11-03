function [F2_max_h, F2_max_t, mean_max_ISRheight] = find_ISR_F2_peak(filename, min_ht, max_ht) 

% Returns both maximum (mean_max_ISRheight) and fitted (F2_max_h) F2 peaks
% for all times in each datafile.
% currently averages over 4 scan positions before identifying peak in order
% to reduce the noise in the profile.
% 
% CJS April 2023
% Fitting suppressed as not useful

% Threshold SNR is set to 0.05. Data with SNR below this produce noisy
% profles and are discarded.

ISRheight = ncread(filename,'height');
ISRpowerdB = ncread(filename,'pwr');
ze = ncread(filename,'ze');

% Correct to received power
% Ps = (C*Rmax.^2)./((1+Tr)*R.^2)*N;
% Where Rmax is the height of the F2 peak
% Where C contains the transmitter power ~ 1MW
% Power is stored in dB (arbitrary units).
% C = Pmax/Nmax where Pmax and Nmax are determined at the F2 peak.

ISRpower = (10.^(ISRpowerdB/10)); % Convert from dB (arbitrary units)

% average over all four beams (third dimension)
ISRpower = mean(ISRpower,3);

% convert to range
[~, timeloop, posloop] = size(ISRpower); % Takes third dimension as '1' if averaging over positions is done beforehand.

ISRrange = NaN*(ones(size(ISRheight)));
ISRnoise = NaN*ones(timeloop,posloop);
ISRSNR = NaN*ones(size(ISRpower));
for j=1:posloop
   ISRrange(:,j) = ht2rng(ISRheight(:,j),ze(j)); 
   for p=1:timeloop
       % sub = find(ISRrange(:,j) >= 700);
       % Estimate noise power from gates above 700 km. 
       ISRnoise(p,j) = mean(ISRpower((ISRrange(:,j) >= 700),p,j));
       ISRpower(:,p,j) = (ISRpower(:,p,j) - ISRnoise(p,j))./((ISRrange(39,j).^2)./(ISRrange(:,j).^2));
       ISRSNR(:,p,j) = ISRpower(:,p,j)./ISRnoise(p,j);
       % Need to find appropriate threshold and set it to avoid bad data
       bad = find(isnan(ISRSNR(:,p,j))==1 | ISRSNR(:,p,j)<0.05);  %set signal to noise threshold to 0.05
       ISRpower(bad,p,j) = NaN;
       %ISRpower(ISRSNR<0.01) = NaN;
   end
end

start_time = ncread(filename,'stime');
obsdatestr = num2str(ncread(filename,'obsdate'));

obsyear = str2double(obsdatestr(1:4));
obsmonth = str2double(obsdatestr(5:6));
obsday = str2double(obsdatestr(7:8));

obsdatenum = datenum(obsyear, obsmonth, obsday);

% Initialise arrays
F2_max_h = NaN*ones(timeloop,posloop);
F2_max_t = NaN*ones(timeloop,posloop);
%max_meanISRheight = NaN*ones(1,timeloop);
max_ISRheight = NaN*ones(timeloop,posloop);
% ISRpower = NaN*ones(size(ISRpower));


for i=1:timeloop
    for k = 1:posloop

%       if exist('plotflag','var') == 1
%         figure(1)
%         hold off
%         plot(ISRpower(:,i,k),ISRheight(:,k),'k')
%         title(num2str(datevec(obsdatenum + double(start_time(i))/(3600*24))))
%         hold off
%         pause
%       end
        
    F2_sub = find(ISRheight(:,k) > min_ht & ISRheight(:,k) < max_ht);   
    [val_max, pos1] = nanmax(ISRpower(F2_sub,i,k));
    
   if ((pos1 ~= 1) && (pos1~= length(F2_sub)))
      max_ISRheight(i,k) = ISRheight(F2_sub(pos1),k);
   end
%     good = find( isnan(ISRpower(F2_sub,i,k)) == 0); 
%     if length(good)
%         [coeffs3] = polyfit(ISRheight(F2_sub(good),k),ISRpower(F2_sub(good),i,k),2);
%         [F2_fit3] = polyval(coeffs3, ISRheight(F2_sub(good)),k);
%         [val,pos2] = nanmax(F2_fit3);
%        if ((pos2 ~= 1) && (pos2~= length(F2_sub(good))))
%             F2_max_h(i,k) = ISRheight(F2_sub(good(pos2)),k);
%        end
%     end
    F2_max_t(i,k) = obsdatenum + double(start_time(i))/(3600*24);

        if exist('plotflag','var') == 1
             figure(1)
             plot(ISRpower(:,i,k),ISRheight(:,k))
             title(num2str(datevec(obsdatenum + double(start_time(i))/(3600*24))))
             hold on
             % plot([ISRnoise(i,k),ISRnoise(i,k)],[ISRheight(1,k),ISRheight(end,k)],'r:')
%              if length(good)
%                 plot(F2_fit3,ISRheight(F2_sub(good),k),'r')
%                 plot(val, ISRheight(F2_sub(good(pos2)),k),'b*')
%              end
             plot(val_max,max_ISRheight(i,k),'ro')
             hold off

             pause
             
             figure(2), plot(ISRSNR(:,i,k),ISRheight(:,k))

             pause
        end
    end

mean_max_ISRheight = nanmean(max_ISRheight,2); 
% dmean_max_ISRheight = nanstd(max_ISRheight,2)/sqrt(posloop); % Irrelevant if already averaged over scan positions
F2_max_h = nanmean(F2_max_h,2);
% dF2_max_h = nanstd(F2_max_h,2)/sqrt(timeloop); % Irrelevant if already averaged over scan positions
F2_max_t = nanmean(F2_max_t,2);
% dF2_max_t = nanstd(F2_max_t,2)/sqrt(timeloop); % Irrelevant if already averaged over scan positions

end  