function [issubsetstruct,num] = containsubset(obj,obj2)
subsetResidues = obj2.getAllResidues;
nonresidues    = obj2.getNonRedEndResidue;
branchnum      = length(nonresidues);
objResidues    = obj.getAllResidues;

num = 0;
issubsetstruct = 0;
for i = 1 : length(subsetResidues)
    CompareResidue = subsetResidues{1,i};
    
    if(isequal(CompareResidue.residueType.name,'freeEnd'))
        continue
    end
    
    parentResidue  = CompareResidue.getParent;
    if(~isequal(parentResidue.residueType.name,'freeEnd'))
        continue
    end
    for j = 1 : length(objResidues)
        ObjResidue = objResidues{1,j};
        isbranchvalid  = 0;
        if(isequal(ObjResidue.residueType.name,CompareResidue.residueType.name))&&...
                (length(CompareResidue.linkageChildren)==length(ObjResidue.linkageChildren))
            iscontinue = 1;
            for k = 1 :branchnum
                objResidue     = ObjResidue;
                compareResidue = CompareResidue;
                while(iscontinue)
                    compareChildren = compareResidue.getChildren;
                    objChildren     = objResidue.getChildren;
                    isvalid         = 0;
                    for ii = 1 : length(objChildren)
                        iithresidue = objChildren(ii);
                        for jj = 1 : length(compareChildren)
                            jjthresidue = compareChildren(jj);
                            if(isequal(iithresidue.residueType.name,jjthresidue.residueType.name))&&...
                                    (isequal(iithresidue.linkageParent.bonds,jjthresidue.linkageParent.bonds))
                                isvalid = isvalid+1;
                                break
                            end
                        end
                    end
                    
                    if(isvalid == length(compareChildren))
                        if(length(compareChildren)==1)
                            compareResidue = compareChildren;
                            objResidue     = objChildren;
                        else
                            compareResidue = compareChildren(k);
                            for kk = 1 : length(objChildren)
                                kkthresidue = objChildren(kk);
                                if(isequal(kkthresidue.residueType.name,compareResidue.residueType.name))
                                    objResidue = kkthresidue;
                                end
                            end
                        end
                        
                        if(isempty(compareResidue.linkageChildren))
                            isbranchvalid = 1;
                            break
                        else
                            if(length(compareResidue.linkageChildren)==length(objResidue.linkageChildren))
                                iscontinue = 1;
                            else
                                iscontinue = 0;
                                break
                            end
                        end
                    else
                        iscontinue = 0;
                        break
                    end
                end
            end
        end
        
        if(isbranchvalid==branchnum)
            num = num+1;
        end
    end
    if(num~=0)
        issubsetstruct = 1;
    end 
end