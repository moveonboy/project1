function [issubsetstruct,num] = isSubsetStructure(obj,obj2)
%isSubsetStructure returns true if one glyanstruct(obj) contains another
%  structure(obj2) and the number of obj2 structure that obj structure
%  contains.
% 
% isSubsetStructure(obj,obj2) compares the two structures, and check if
%    obj2 is fragmentation of obj.
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 9/23/14

subsetResidues = obj2.getAllResidues;
objResidues    = obj.getAllResidues;

num = 0;
issubsetstruct = 0;
for i = 1 : length(subsetResidues)
    ithResidue = subsetResidues{1,i};
    if(isequal(ithResidue.residueType.name,'freeEnd'))
        continue
    end
    
    ithParent  = ithResidue.getParent;
    if(isequal(ithParent.residueType.name,'freeEnd'))
        compareResidue = ithResidue;
        break
    end
end

for j = 1 : length(objResidues)
    jthresidue = objResidues{1,j};
    iscontainStructure = equalsResidue(jthresidue,compareResidue);
    if(iscontainStructure)
        num = num+1;
    end
end
if(num~=0)
    issubsetstruct = 1;
end
