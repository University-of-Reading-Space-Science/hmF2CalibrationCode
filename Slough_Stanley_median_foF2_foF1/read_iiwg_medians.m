function [year, month, CC, UT_offset, medians, counts, range, upper_q, lower_q, upper_d, lower_d] = read_iiwg_medians(infile);
% function [year, month, CC, UT_offset, medians, counts, range, upper_q, lower_q, upper_d, lower_d] = read_iiwg_medians(infile);
% function to read monthly median data files of foF2 obtained from
% www.ukssdc.ac.uk and correct for local time offsets
% CJS November 2013

% Bug fixed 05/02/2014 if index < 1 now results in index = 24 + index
% rather than 24-index which was giving values > 24 hours in a day
% CJS 05/02/2014

% fid=fopen('RL052_193001_201212_medians.txt');

fid=fopen(infile);
%# Get file size.
i=1;
while(~feof(fid))
    fgets(fid);
    i=1+1;
end

filesize=i;
frewind(fid);

i=1;
line=fgets(fid);

medians=NaN*ones(24,round(filesize/9));
counts=NaN*ones(24,round(filesize/9));
range=NaN*ones(24,round(filesize/9));
upper_q=NaN*ones(24,round(filesize/9));
lower_q = NaN*ones(24,round(filesize/9));
upper_d = NaN*ones(24,round(filesize/9));
lower_d = NaN*ones(24,round(filesize/9));
year = NaN*ones(1,round(filesize/9));
month = NaN*ones(1,round(filesize/9));
CC = NaN*ones(1,round(filesize/9));

while(~feof(fid))
if strcmp(line(1),'T')==1
     UT_offset(i) = str2double(line(20:22));
     line=fgets(fid);
elseif strcmp(line(1),'Y')==1
    % do nothing
    line=fgets(fid);
else 
    year(i) = str2double(line(1:4));
    month(i) = str2double(line(5:6));
    CC(i) = str2double(line(8:9));
    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end
       medians(index,i) = str2double(line(13+(j-1)*6:15+(j-1)*6));
    end
    line=fgets(fid);
    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end
       counts(index,i) = str2double(line(15+(j-1)*6:17+(j-1)*6));
    end
    line=fgets(fid);
    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end
        range(index,i) = str2double(line(15+(j-1)*6:17+(j-1)*6));
    end
    line=fgets(fid);
    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end
        upper_q(index,i) = str2double(line(13+(j-1)*6:15+(j-1)*6));
    end
    
    line=fgets(fid);

    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end
        lower_q(index,i) = str2double(line(13+(j-1)*6:15+(j-1)*6));
    end
    
    line = fgets(fid);

    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end

       upper_d(index,i) = str2double(line(13+(j-1)*6:15+(j-1)*6));
    end
    
    line=fgets(fid);
    
    for j=1:24
       index = j-UT_offset(i);
       if index < 1
           index = 24+index;
       elseif index > 24
           index = index-24;
       end

        lower_d(index,i) = str2double(line(13+(j-1)*6:15+(j-1)*6));
    end
    
  line = fgets(fid);  
  
  i=i+1;
end
end

end