function outcell = Cell_diff(cell1,cell2)

outcell = cell(size(cell1));

for i = 1:length(cell1)
    
    outcell{i} = cell1{i}-cell2{i};
    
end
