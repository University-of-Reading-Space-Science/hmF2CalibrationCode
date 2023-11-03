function [sunza,sunaz] = sunzena(lat,lon,year,day,hour,min,sec)

% SUNZENA:   Calculates the solar zenith angle and azimuth.
%
%          Allowed Calls:-
%            SUNZENA(lat,lon,year,day,hour,min,sec).
%
%          Given Arguments:-
%            lat:  latitude of point on the surface of the Earth in rad..
%            lon:  longitude of point on the surface of the Earth in rad..
%            year: integral part of the year A.D..
%            day:  integral part of the day number in the year.
%            hour: integral part of the hour U.T..
%            min:  integral part of the minutes U.T..
%            sec:  seconds U.T..
%
%          Returned Arguments:-
%            sunza: solar zenith angle from above point.
%            sunaz:   solar azimuth from above point.
%
%          Input Parameters: (none).
%
%          Output Parameters: (none).
%
%          Notes:-
%            a) "day" for January 1 = 1.
%
%          Warnings:-
%            a) all given parameters must be scalars.
%
%          Defaults:-
%            lat:  none - must be given.
%            lon:  none - must be given.
%            year: none - must be given.
%            day:  none - must be given.
%            hour: none - must be given.
%            min:  none - must be given.
%            sec:  none - must be given.
%
%          Abends:-
%            a) Error if incorrect number of arguments given.
%
%          Sources:-
%            C.A.MURRAY:
%              "Vectorial Astrometry",
%                Adam Hilger, Bristol, (1983).
%
%          This function uses ATAN2, COS, ERROR, NARGIN, NARGOUT, SIN,
%           SQRT, SUNPOS.
%
%          K.S.C. Freeman (8 September 1995).

%---Check number of arguments.
%     Proceed only if:-
%       No. Given:    7. 
%       No. Returned: 2.
if nargin == 7 & nargout == 2

%---Set up latitude rotation matrix.
  sinlat = sin(lat);
  coslat = cos(lat);
  latrotn = [0,-sinlat,coslat;1,0,0;0,coslat,sinlat];

%---Determine the position of the Sun.
  [sunlon,sunra,sundec,gst] = sunpos(year,day,hour,min,sec);

%---Set up solar declination and right ascension vector.
  sinsdec = sin(sundec);
  cossdec = cos(sundec);
  sinsra = sin(sunra);
  cossra = cos(sunra);
  sradec = [cossdec*cossra;cossdec*sinsra;sinsdec];

%---Calculate local mean sidereal time. 
  lst = gst + lon;

%---Set up local mean sidereal time rotation matrix.
  sinlst = sin(lst);
  coslst = cos(lst);
  lstrotn = [coslst,-sinlst,0;sinlst,coslst,0;0,0,1];

%---Compute transformation matrix.
  tranmat = (lstrotn*latrotn)';

%---Obtain solar zenith angle and azimuth vector.
  szenaz = tranmat*sradec;

%---Determine solar zenith angle and azimuth.
  sunza = atan2(sqrt(szenaz(1)^2 + szenaz(2)^2),szenaz(3));
  sunaz = atan2(szenaz(1),szenaz(2));

%---Error if incorrect number of arguments.
else
 error('Incorrect number of arguments')

%---End block over number of arguments.
end
