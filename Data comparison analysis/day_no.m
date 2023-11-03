function date = day_no(year, month, day)
% FUNCTION DATE = day_no(YEAR, MONTH, DAY) outputs the date as a
% day number. It accounts for leap years.
% Made millennium compliant 5.11.2001 CJD.

ord_year = [31 28 31 30 31 30 31 31 30 31 30 31];
leap_year =[31 29 31 30 31 30 31 31 30 31 30 31];

if ((fix(year)/4 -fix(year/4) == 0) & (fix(year)/100 - fix(year/100)) ~= 0)...
    | (fix(year)/400 -fix(year/400)) == 0
 date = sum(leap_year(1:month-1)') + day;
else
  date = sum(ord_year(1:month-1)) + day;
end
