%-----(Elemental Node - Mapping table)-----------%
function  [element_nodes]= Mapping_table_fun(nel)

element_nodes= zeros(nel,2);
for i = 1:nel
    element_nodes(i,:) = [i, i+1] ;
end
end

