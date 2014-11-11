function [numSubstr,newsubstrSpecies] = reorgnizestructure(substrSpecies,enzObj)
listofspecis         = CellArrayList;
listofspecis_bracket = CellArrayList;
newsubstrSpecies     = CellArrayList;
resfuncgroupname     = enzObj.resfuncgroup.name;
for i = 1: length(substrSpecies)
    ithSpecies     = substrSpecies.get(i);
    glycanObj      = ithSpecies.glycanStruct;
    commomresidues = getAllResidues(glycanObj,'root');
    counter        = 0;
    for j = 1 : length(commomresidues)
        jthresidue = commomresidues{1,j};
        if(~isempty(jthresidue.linkageChildren))
            continue
        end
        
        if(isequal(jthresidue.residueType.name,'Fuc'))
            continue
        end
        
        if(isequal(jthresidue.residueType.name,'GlcNAc'))&&...
                (isequal(jthresidue.getParent.residueType.name,'Man'))&&...
                (isequal(jthresidue.getParent.anomer.symbol,'b'))
            continue
        end
        
        branchdepth(counter+1) = glycanObj.getBranchDepth(jthresidue);
        counter                = counter+1;
    end
    
    if(isempty(enzObj.isTerminalTarget))
        enzObj.isTerminalTarget = true;
    end
    
    if(enzObj.isTerminalTarget)
        if(~isequal(max(branchdepth),min(branchdepth)))
            listofspecis.add(ithSpecies)
        else
            listofspecis_bracket.add(ithSpecies)
        end
    else
        numoftargetresidue = 0;
        for j = 1 : length(commomresidues)
            jthresidue = commomresidues{1,j};
            if(~isequal(enzObj.resfuncgroup.name,'Fuc'))
                if(isequal(jthresidue.residueType.name,resfuncgroupname))
                    numoftargetresidue = numoftargetresidue+1;
                end
            else
                if(isequal(jthresidue.residueType.name,'Fuc'))&&...
                        (~isequal(jthresidue.linkageParent.bonds.posParent,'6'))
                    numoftargetresidue = numoftargetresidue+1;
                end 
            end
        end
        if(numoftargetresidue<length(branchdepth))&&(numoftargetresidue>0)
            listofspecis.add(ithSpecies)
        else
            listofspecis_bracket.add(ithSpecies)
        end 
    end
end

if(~isempty(listofspecis_bracket))
    species = listofspecis_bracket.get(1);
    if(isempty(species.glycanStruct.bracket.linkageChildren))
       species.glycanStruct.bracket = [];
    end
    newsubstrSpecies.add(listofspecis_bracket.get(1));
end

if(~isempty(listofspecis))
    glycanSpecies  = listofspecis.get(1);
    bracketresidue = getAllResidues(glycanSpecies.glycanStruct,'bracket');
    [corespecies,corebracket] = getstruct(glycanSpecies,resfuncgroupname);
    newglycanSpecies = combinestruct(corespecies,corebracket,bracketresidue,enzObj);
    newsubstrSpecies.add(newglycanSpecies);
end
numSubstr = length(newsubstrSpecies);
end

function [corespecies,corebracket] = getstruct(Species,resfuncgroupname)
% create core structure
corestruct               = Species.glycanStruct.clone;
corestruct.bracket       = [];
corestruct.resetjava;
nonreResidues            = corestruct.getNonRedEndResidue;
counter = 0;
for i = 1 : length(nonreResidues)
    ithresidue = nonreResidues{1,i};
    if(~isequal(resfuncgroupname,'Fuc'))
        if(isequal(ithresidue.residueType.name,resfuncgroupname))
            corebracket(counter+1,1) = ithresidue;
            counter = counter+1;
        end
    else
        if(isequal(ithresidue.residueType.name,resfuncgroupname))&&...
                (~isequal(ithresidue.linkageParent.bonds.posParent,'6'))
            corebracket(counter+1,1) = ithresidue;
            counter = counter+1;
        end
    end
end
corespecies              = Species.clone;
corespecies.glycanStruct = corestruct;
corespecies.name         = corestruct.name;
end

function newglycanSpecies = combinestruct(corespecies,corebracket,bracketresidue,enzObj)
counter = 0;
for i = 1 : length(bracketresidue)
    ithresidue = bracketresidue{1,i};
    if(isequal(ithresidue.residueType.name,'#bracket'))
        continue;
    end
    parent     = ithresidue.getParent;
    if(isequal(parent.residueType.name,'#bracket'))
        ithresidue.linkageParent = [];
        addresidue(counter+1,1) = ithresidue;
        counter = counter+1;
    end
end

counter = 0;
loci    = [];
branchedresidue = '';
for i = 1 : length(addresidue)
    ithaddresidue = addresidue(i);
    if(isequal(ithaddresidue.residueType.name,'Fuc'))
        branchedresidue{counter+1,1} = ithaddresidue;
        loci(counter+1) = i;
        counter = counter+1;
    end
end
addresidue(loci) = '';

if(enzObj.isTerminalTarget)
    for i = 1 : length(corebracket)
        ithresidue    = corebracket(i);
        try
            ithaddresidue = addresidue(i);
        catch
            break
        end
        isResidueAdded = addResidue(corespecies.glycanStruct,ithresidue,ithaddresidue);
    end
    corespecies.glycanStruct.bracketResidue(corebracket);
    if(~isempty(branchedresidue))
        bracketresidues = getAllResidues(corespecies.glycanStruct,'bracket');
        for i = 1 : length(bracketresidues)
            ithresidue = bracketresidues{i};
            if(isequal(ithresidue.residueType.name,'#bracket'))
                for j = 1 : length(branchedresidue)
                    jthresidue = branchedresidue{j};
                    isResidueAdded = addResidue(corespecies.glycanStruct,ithresidue,jthresidue);
                end
                break
            end
        end
    end
else
    corespecies.glycanStruct.bracketResidue(corebracket);
    corebracketresidue = getAllResidues(corespecies.glycanStruct,'bracket');
    for i = 1 : length(corebracketresidue)
        ithresidue = corebracketresidue{i};
        if(isequal(ithresidue.residueType.name,'#bracket'))
            for j =1 : length(bracketresidue)
                jthresidue = bracketresidue{j};
                if(isequal(jthresidue.residueType.name,'#bracket'))
                    continue
                end
                isResidueAdded = addResidue(corespecies.glycanStruct,ithresidue,jthresidue);
            end
        end
    end
end
corespecies.glycanStruct.resetjava;


if(isResidueAdded)
    newglycanSpecies      = corespecies.clone;
    newglycanSpecies.name = newglycanSpecies.glycanStruct.name;
end
end