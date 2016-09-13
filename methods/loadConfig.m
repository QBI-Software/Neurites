function loadConfig(M, programtype,clearlabels)
% Loads or clears config values to GUI
% M is a structure containing field values or empty string for clearing
% programtype indicates which program using this eg 'annulus' or 'synapse'
% clearlabels flag will also clear fields with configurable labels
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    initval = '';
    if ~isempty(M)
        set(hScale,'String',num2str(M.Scale));
        set(hSX,'String',num2str(M.Shiftx));
        set(hSY,'String',num2str(M.Shifty));
    else
        set(hScale,'String',initval);
        set(hSX,'String',initval);
        set(hSY,'String',initval);
    end
    
    if (strcmp(programtype,'annulus') > 0)
        hCx = findobj('Tag','editCentroidX');
        hCy = findobj('Tag','editCentroidY');
        hOD = findobj('Tag','editOD');
        hID = findobj('Tag','editID');
        hW = findobj('Tag','editWidth');
        hH = findobj('Tag','editHeight');
        if ~isempty(M)
            set(hCx,'String',num2str(M.CentroidX));
            set(hCy,'String',num2str(M.CentroidY));
            set(hOD,'String',num2str(M.AnnulusOD));
            set(hID,'String',num2str(M.AnnulusID));
            set(hFit,'String',num2str(M.Fit));
            set(hW,'String',num2str(M.Width));
            set(hH,'String',num2str(M.Height));
        else %clear allfields
            set(hScale,'String',initval);
            set(hSX,'String',initval);
            set(hSY,'String',initval);
            set(hCx,'String',initval);
            set(hCy,'String',initval);
            set(hOD,'String',initval);
            set(hID,'String',initval);
            set(hFit,'String',initval);
            set(hW,'String',initval);
            set(hH,'String',initval);
        end
    else
        hCell1 = findobj('Tag','editCell1');
        hCell2 = findobj('Tag','editCell2');
        hFit = findobj('Tag','editFit');
        hMin = findobj('Tag','editMinlength');
        hTol = findobj('Tag','editTolerance');
        if  ~isempty(M)
            set(hCell1,'String',M.Cell1);
            set(hCell2,'String',M.Cell2);
            set(hScale,'String',num2str(M.Scale));
            set(hFit,'String',num2str(M.Fit));
            set(hSX,'String',num2str(M.Shiftx));
            set(hSY,'String',num2str(M.Shifty));
            set(hMin,'String',num2str(M.Minlength));
            set(hTol,'String',num2str(M.Tolerance));
            set(hCell1,'ForegroundColor','k');
            set(hCell2,'ForegroundColor','k');
        else
            if (clearlabels)
                set(hCell1,'String',initval);
                set(hCell2,'String',initval);
            end
            set(hScale,'String',initval);
            set(hFit,'String',initval);
            set(hSX,'String',initval);
            set(hSY,'String',initval);
            set(hMin,'String',initval);
            set(hTol,'String',initval);
        end
    end
    
    
end

