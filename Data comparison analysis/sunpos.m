function [sunlon,sunra,sundec,gst] = sunpos(year,day,hour,min,sec)

% SUNPOS:   Calculates the position of the Sun as seen from the Earth.
%
%          Allowed Calls:-
%            SUNPOS(year,day,hour,min,sec).
%
%          Given Arguments:-
%            year: integral part of the year A.D..
%            day:  integral part of the day number in the year.
%            hour: integral part of the hour U.T..
%            min:  integral part of the minutes U.T..
%            sec:  seconds U.T..
%
%          Returned Arguments:-
%            sunlon: longitude of the Sun along the ecliptic in rad..
%            sunra:  right ascension of the Sun in rad..
%            sundec: declination of the Sun in rad..
%            gst:    Greenwich mean sidereal time in rad..
%
%          Input Parameters: (none).
%
%          Output Parameters: (none).
%
%          Notes:-
%            a) "day" for January 1 = 1.
%            b) any one given parameter may be a vector.
%            c) if so: returned parameters will be vectors of equal length.
%
%          Warnings:-
%            a) "year" must lie between 1901 and 2099 inclusive.
%
%          Defaults:-
%            year: none - must be given.
%            day:  none - must be given.
%            hour: none - must be given.
%            min:  none - must be given.
%            sec:  none - must be given.
%
%          Abends:-
%            a) Error if incorrect number of arguments given.
%            b) Error if "year" lies outside the above limits.
%
%          Source:-
%            C.T.RUSSELL:
%              "Geophysical Coordinate Transformations",
%                Cosmic Electrodyn., Vol.2, 184-196, (1971).
%
%          This function uses ATAN, ATAN2, COS, ERROR, FIX, LENGTH, NARGIN,
%           NARGOUT, PI, REM, SIN, SQRT.
%
%          K.S.C. Freeman (15 August 1995) after RUSSELL(1971).

%---Check number of arguments.
%     Proceed only if:-
%       No. Given:    5.
%       No. Returned: 4.
if nargin == 5 & nargout == 4

%---Check that "year" is in the correct range.
%     Proceed only if 1901 < "year" < 2099.
  if year >= 1901 & year <= 2099

%---Set multiples of pi and degree-to-radian conversion.
    twopi = 2.0*pi;
    degtrad = pi/180.0;

%---Calculate fraction of day.
    fracday = (3600*hour + 60*min + sec)/86400;

%---Calculate number of days since 1900.0 and fraction of Julian century.
    day1900 = 365*(year-1900) + fix((year-1901)/4.0) + day - 0.5 + fracday;
    fracjc = day1900/36525.0;

%---Determine Greenwich mean sidereal time.
    gstpar = 279.690983 + 0.9856473354*day1900 + 360.0*fracday + 180.0;
    gst = rem(gstpar,360.0)*degtrad;

%---Obtain geometric mean longitude of the Sun.
    sungml = rem(279.696678 + 0.9856473354*day1900,360.0);

%---Compute mean anomaly of the Sun.
    g = rem(358.475845 + 0.985600267*day1900,360.0)*degtrad;

%---Calculate (true) longitude of the Sun.
    sunljc = 1.91946 - 0.004789*fracjc;
    sunlon = (sungml + sunljc.*sin(g) + 0.020094*sin(2.0*g))*degtrad;

%---Determine number of time points required.
    ntimes = length(sunlon);

%---Correct "sunlon" to lie between 0 and 2pi.
    for i = 1:ntimes
      if sunlon(i) < 0.0
        sunlon(i) = sunlon(i) + twopi;
      elseif sunlon(i) > twopi
        sunlon(i) = sunlon(i) - twopi;
      end
    end

%---Obtain obliquity of the ecliptic.
    obliq = (23.45229 - 0.0130125*fracjc)*degtrad;
    sinobl = sin(obliq);

%---Compute declination of the Sun.
    sindec = sinobl.*sin(sunlon - 9.924e-5);
    cosdec = sqrt(1.0 - sindec.^2);
    tandec = sindec./cosdec;
    sundec = atan(tandec);

%---Calculate right ascension of the Sun.
    sunra = pi - atan2(cos(obliq).*tandec./sinobl,-cos(sunlon)./cosdec);

%---Error if "year" out of range.
  else
    error('Year out of range - should be between 1901 and 2099')

%---End block over range of "year".
  end

%---Error if incorrect number of arguments.
else
  error('Incorrect number of arguments')

%---End block over number of arguments.
end
