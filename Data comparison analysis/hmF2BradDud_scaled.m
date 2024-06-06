function [hmF2_BradDud, dhmF2_BradDud] = hmF2BradDud_scaled(foF2,foE, M3000F2,dfoF2, dfoE, dM3000F2, scaling_factor)
% function [hmF2_BradDud, dhmF2_BradDud] = hmF2Dud(foF2,foE, M3000F2, dfoF2, dfoE, dM3000F2)
%
%   A function to calculate the true height of the F2 peak from an
%   empirical formula derived from common ionospheric parameters
%   Ref Bradley & Dudeney, 1973
%
%   Inputs:
%   foF2 - the peak frequency of the F-region (in MHz)
%   foE - the peak frequency of the E-region (in MHz)
%   M3000F2 -  the M3000 index
%   Optional inpus:
%   dfoF2 - uncertainty in the peak frequency of the F-region (in MHz)
%   dfoE - uncertainty in the peak frequency of the E-region (in MHz)
%   dM3000F2 -  uncertainty in the M3000 index
%
%   The function returns hmF2 calculated from the formula
%
%   Christopher Scott, October 2023.

% Bradley and Dudeney 1973 (for X > 1.7)
x = scaling_factor*(foF2./foE);
bad = find( (x <1.7)); % < 1.05*(355/1890+1.4)) );
x(bad) = NaN;

b = (2.5*x-3).^(-2.35) - 1.6;
a = 1890-355./(x-1.4);
hmF2_BradDud = a.*(M3000F2).^b;

if exist('dfoF2')==1
  dx = x.*sqrt( (dfoF2./foF2).^2 + (dfoE./foE).^2 );
  dx(bad) = NaN;
  dbdx = -2.5*2.35*(2.5*x-3).^-3.35;
  db = sqrt( (dx.^2).*(dbdx.^2) );
  dadx = 355*(x-1.4);
  da = sqrt( (dx.^2).*(dadx).^2 );
  dDda = M3000F2.^b;                      % differential of formula wrt a
  dDdb = a.*(M3000F2.^b).*log10(M3000F2);   % differential of formula wrt b
  dDdM = a.*b.*(M3000F2).^(b-1);          % differential of formula wrt M3000F2

  dhmF2_BradDud = sqrt( (da.*da).*(dDda).^2 + (db.*db).*(dDdb).^2 + (dM3000F2.^2).*dDdM.^2);
end

end