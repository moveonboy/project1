function isSubsetStruct = checksubset(glycanStruct1,glycanStruct2)
% isSubsetStruct returns 1(true) if glycanStruct2 is subset structure of
% glycanStruct1. Otherwise it returns 0(false);

if(isempty(glycanStruct1.bracket))
    isSubsetStruct = glycanStruct1.contains(glycanStruct2);
else
    % step1 : separate glycanStruct1 into two parts, one is core-structure and another one
    % is bracket-structure.
    core_structure         = glycanStruct1.clone;
    core_structure.bracket = [];
    
    bracket_structure      = glycanStruct1.clone;
    bracket_structure.root = []; 
    bracketresidues        = bracket_structure.getAllResidues;
    
    coreterResidue = core_structure.getNonRedEndResidue;
    core_depth     = [];
    for i = 1 : length(coreterResidue)
        ithresidue = coreterResidue{1,i};
        if(isequal(ithresidue.residueType.name,'Fuc'))
            continue
        end
        core_depth = getBranchDepth(core_structure,ithresidue);
        break
    end
    
    % step2 : separate glycanStruct2 into two parts based on the
    % core_structure and bracket_structure from glycanStruct1.
    
    compareobj = glycanStruct2.clone;
    objnreResidue = compareobj.getNonRedEndResidue;
    counter           = 0;
    objdepth      = [];
    fucindex          = [];
    for i = 1 : length(objnreResidue)
        ithresidue = objnreResidue{1,i};
        if(isequal(ithresidue.residueType.name,'Fuc'))
            fucindex = i;
            continue
        end
        objdepth(counter+1) = getBranchDepth(compareobj,ithresidue);
        counter = counter + 1;
    end
    objnreResidue(1,fucindex) = '';
    
    if(core_depth>=max(objdepth))
        isSubsetStruct = core_structure.contains(compareobj);
    else
        counter = 0;
        for i = 1 : length(objnreResidue)
            ithresidue = objnreResidue{1,i};   
            branchdepthdiff = objdepth(i)-core_depth;
            if(branchdepthdiff==0)
                continue
            end
            for j = 1 : abs(branchdepthdiff)-1
                parentresidue = ithresidue.getLinkageParent.getParent;
                ithresidue = parentresidue;
            end
            residuestobracket(counter+1,1) = ithresidue;
            counter                        = counter+1;
        end
        compareobj.bracketResidue(residuestobracket);
        
        core2_struct = compareobj.clone;
        core2_struct.bracket = [];
        
        bracket2_struct      = compareobj.clone;
        bracket2_struct.root = [];
        bracket2residues     = bracket2_struct.getAllResidues;
        
        iscoresubset = core_structure.contains(core2_struct);
        isbracketsubset = bracket_structure.contains(bracket2_struct);

        isSubsetStruct = (iscoresubset)&&(isbracketsubset);
    end
end
end