function [ labelout ] = cleanlabel( label )
%Check label is valid for table headers
    labelout = strrep(label,' ','');
    labelout = strrep(labelout,'(','_');
    labelout = strrep(labelout,')','_');

end

