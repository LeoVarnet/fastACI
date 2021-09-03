function inoutcell = Cell_mult(inoutcell,num)

b_is_num = ~iscell(num);

for i = 1:length(inoutcell)
    
    if b_is_num == 1
        inoutcell{i} = inoutcell{i}*num;
    else
        inoutcell{i} = inoutcell{i}.*num{i};
    end
    
end
