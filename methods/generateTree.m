function tree = generateTree( csvdata )
%Generate dendrogram style tree from csvdata
%   Uses Matlab Phylogenetic tree: phytree 
% 1. Extract branch data per tree/branch order from csv
% 2. Determine orientation of branch as key (degrees)
% 3. List branches in order of order and location
CSV = readtable(csvdata);
t = unique(CSV.Tree);
u = unique(CSV.Order);
csvidx = find(CSV.Tree ==1); %indices of matching rows
branches = {}; %cell(t,u);
blength = 0;
for i=1:length(csvidx)
    order = CSV.Order(i);
    btype = CSV.PointType(i);
    blength = blength + CSV.Length__m_(i);
    if (strcmp(btype, 'CP')==0)
        if (length(branches) >= order)
            branches{order} = cat(1,branches{order},blength);
        else
            branches{order} = [blength];
        end
        blength = 0;
    end
        
end
%load tree vars
branchdata = [];%[1 2;3 4];
branchdistances = [];% [290.7;125;324.2;139.7;0]
for j=1:length(u)
    branchdata = cat(2, branchdata,[j j+1]);
    j = j+1;
end


tree = phytree(branchdata,branchdistances)
names = get(tree,'LeafNames')

%Show plot
view(tree);
%h = plot(tree,'orient','top');
xlabel('Distance from soma (um)')
%set(h.terminalNodeLabels,'Rotation',65)

end
