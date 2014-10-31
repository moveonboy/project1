function iscontainStructure = SubtreeEquals(Residue1,Residue2)
%SubtreeEquals compare two structures residue by residue to check if one 
% contains another one.
%
% iscontainStructure = SubtreeEquals(Residue1,Residue2) returns 1(true) if  
% the structure A is containing structure B. residue1 is the root of
% structure A, residue2 is the root of structure B.
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 10/15/14

loci = '';
counter = 0;
iscontainStructure = 1;
isTypeEqual      = typeEquals(Residue1,Residue2);
if(~isTypeEqual)
    iscontainStructure = 0;
    return
end

if(~isempty(Residue2.linkageChildren))
    isChildrenmatch  = (length(Residue1.linkageChildren)==length(Residue2.linkageChildren));
    if(~isChildrenmatch)
        iscontainStructure = 0;
        return
    end
else
    iscontainStructure = 1;
    return
end

for i = 1 : length(Residue1.linkageChildren)
    ithchild = Residue1.linkageChildren(i).child;
    isnumRequireapproved = 0;
    for j = 1 : length(Residue2.linkageChildren)
        if(~isempty(strfind(loci,num2str(j))))
            continue
        end
        jthchild = Residue2.linkageChildren(j).child;
        iscontainStructure = SubtreeEquals(ithchild,jthchild);
        if(iscontainStructure)
           isnumRequireapproved =  isnumRequireapproved+1;
           loci(counter+1) = num2str(j);
           counter = counter+1;
           break
        end
    end
    
    if(isnumRequireapproved == 0)
        iscontainStructure = 0;
        break
    else
        iscontainStructure = 1;
    end
        
end
end

function isTypeEqual = typeEquals(obj1,obj2)
isSatisfyRequire = 0;
if(isequal(obj1.residueType,obj2.residueType))
    isSatisfyRequire = isSatisfyRequire+1;
end
if(isequal(obj1.anomer,obj2.anomer))
    isSatisfyRequire = isSatisfyRequire+1;
end
if(isequal(obj1.stereoConfig,obj2.stereoConfig))
    isSatisfyRequire = isSatisfyRequire+1;
end
if(isequal(obj1.ringType,obj2.ringType))
    isSatisfyRequire = isSatisfyRequire+1;
end
if(~isempty(obj1.linkageParent))
    if(isequal(obj1.linkageParent.bonds,obj2.linkageParent.bonds))
        isSatisfyRequire = isSatisfyRequire+1;
    end
    
    if(isSatisfyRequire == 5)
        isTypeEqual = 1;
    else
        isTypeEqual = 0;
    end
else
    if(isSatisfyRequire == 4)
        isTypeEqual = 1;
    else
        isTypeEqual = 0;
    end
end
end