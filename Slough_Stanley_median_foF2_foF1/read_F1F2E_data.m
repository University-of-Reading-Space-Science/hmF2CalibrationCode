% function [datenumber, foF2, foF1] = read_F1F2E_data(infile)
function [datenumber, foF2, foF1, foE] = read_F1F2E_data(infile)

% infile = 'Slough_E_F2_F1.txt'

i=1;
fid = fopen(infile);
line = fgets(fid);

while(~feof(fid))

    if strcmp(line(1), 'T') == 1
        UT_offset = str2num(line(20:22));
        line = fgets(fid);
        % line = fgets(fid);

    else

        year(i) = str2num(line(1:4));
        month(i) = str2num(line(5:6));
        day(i) = str2num(line(7:8));
        hour(i) = str2num(line(10:11));
        minute(i) = str2num(line(12:13));
        second(i) = str2num(line(14:15));
        
        datenumber(i) = datenum(year(i), month(i), day(i), hour(i)-UT_offset, minute(i), second(i));
        
        foF2(i) = str2double(line(17:19));

        if length(line) > 23
           foF1(i) = str2double(line(23:25));
        end

        if length(line) > 29
            foE(i) = str2double(line(29:31));
        end
        i=i+1;
    end

    line = fgets(fid);
end

bad = find(foF1 == 0);
foF1(bad) = NaN;

fclose(fid)