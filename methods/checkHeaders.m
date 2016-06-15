function rtn = checkHeaders(csvfile,hdrs)
    rtn = 0;           
    fid = fopen(csvfile);
    C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',', ...
                 'treatAsEmpty', {'NA', 'na'}, ...
                 'commentStyle', '//');
    fclose(fid);
    for i=1: length(C)
        v = strjoin(C{1,i}(1));
        v = strrep(v,' ','');
        if (size(strfind(hdrs,v)) == 0)
            csvfile
            rtn = v
            break
        end
    end
end

