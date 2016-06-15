function clearCSVdata(programtype)
    hCf0 = findobj('Tag', 'Menu_File_loadcsv');
    set(hCf0,'UserData',{});
    if strcmp(programtype, 'synapse') > 0
        hCf1 = findobj('Tag', 'editAnalysisfiles1');
        hCf2 = findobj('Tag', 'editAnalysisfiles2');
        set(hCf1,'String','');
        set(hCf2,'String','');
    end
end

