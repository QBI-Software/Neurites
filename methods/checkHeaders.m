function rtn = checkHeaders(csvfile,hdrs)
    rtn = 0;           
    fid = fopen(csvfile);
    C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',', ...
                 'treatAsEmpty', {'NA', 'na'}, ...
                 'commentStyle', '//');
    fclose(fid);
    h = strsplit(hdrs,',');
    for i=1: length(C)
        v = strjoin(C{1,i}(1));
        v = strrep(v,' ','');
        if (isempty(strfind(hdrs,v)) || isempty(strmatch(v,h,'exact')))
            csvfile
            rtn = v
            break
        end
    end
end

