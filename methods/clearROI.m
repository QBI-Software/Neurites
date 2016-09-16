function completed = clearROI( handles )
%clearROI :   Removes any ROIs on axes
completed =0;
axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
end	
completed =1;
end

