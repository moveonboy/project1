function newsubstrSpecies = reorgnizestructure(substrSpecies)
%reorgnizestructure_2ndEdition: optimize bracket isomer species into one new bracket form
%
% Syntax:
%   newsubstrSpecies = reorgnizestructure_2ndEdition(substrSpecies)
%
% Input:
%   substrSpecies: list of isomers
%
% Output:
%    newsubstrSpecies:  a species represented in new bracket format
%
% Example:
%  Example 1:
%       tetrabracket1=GlycanSpecies(glycanMLread('tetrabrackettest1.glycoct_xml'));
%       tetrabracket2=GlycanSpecies(glycanMLread('tetrabrackettest2.glycoct_xml'));
%       tetrabracket3=GlycanSpecies(glycanMLread('tetrabrackettest3.glycoct_xml'));
%       tetrabracket4=GlycanSpecies(glycanMLread('tetrabrackettest4.glycoct_xml'));
%       substrSpecies = CellArrayList;
%       substrSpecies.add(tetrabracket1);
%       substrSpecies.add(tetrabracket2);
%       substrSpecies.add(tetrabracket3);
%       substrSpecies.add(tetrabracket4);
%       newsubstrSpecies = reorgnizestructure_2ndEdition(substrSpecies)
%       glycanViewer(newsubstrSpecies.glycanStruct);
%
%
%See also merge2bracket.

% Author:Gang Liu and Yusen Zou
% Date: 11/25/14
speiciesnum  = length(substrSpecies);
newObjspecies = CellArrayList;
for i = 1 : speiciesnum
    ithsubstrObj    = substrSpecies.get(i).glycanStruct.clone;
    branches        = getBranch(ithsubstrObj);
    branchdepth     = [];
    residueterminal = '';
    count = 0;
    for j = 1 : length(branches)
        branchdepth                = [branchdepth branches(1,j).depth];
        residueterminal{count+1} = branches(1,j).residueterminal;
        count = count+1;
    end
    [~,minpos] = min(branchdepth);
    residueterminal(minpos) = '';
    pos    = 0;
    addpos = 0;
    for j = 1 : length(ithsubstrObj.bracket.linkageChildren)
        pos = pos+1;
        jthresidue = ithsubstrObj.bracket.linkageChildren(pos,1).child;
        if(~strcmp(jthresidue.residueType.name,'Fuc'))
            addpos = addpos+1;
            addtoresidue = residueterminal{1,addpos};
            ithsubstrObj.removeResidueObj(jthresidue);
            isResidueAdded = ithsubstrObj.addResidue(addtoresidue,jthresidue);
            pos = pos-1;
        end
    end
    ithsubstrObj.resetjava;
    newObjspecies.add(GlycanSpecies(ithsubstrObj));
end
newsubstrSpecies = merge2bracket(newObjspecies);
end


