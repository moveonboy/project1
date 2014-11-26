function bracketspecies = merge2bracket(isomerspecies)
%MERGE2BRACKET: cluster isomer species to species in bracket form
%
% Syntax:
%   specicesbracket = merge2bracket(isomerspecies)
%
% Input:
%   isomerspecies: list of isomers
%
% Output:
%    specicesbracket: cluster of species represented in bracket format
%
% Example:
%  Example 1:
%       m7{1}=GlycanSpecies(glycanMLread('m7a.glycoct_xml'));
%       m7{2}=GlycanSpecies(glycanMLread('m7b.glycoct_xml'));
%       m7{3}=GlycanSpecies(glycanMLread('m7c.glycoct_xml'));
%       m7{4}=GlycanSpecies(glycanMLread('m7d.glycoct_xml'));
%       m7bracket=merge2bracket(m7);
%       glycanViewer(m7bracket);
%
%  Example 2:
%        tetra{1}=GlycanSpecies(glycanMLread('tetra_test1.glycoct_xml'));
%        tetra{2}=GlycanSpecies(glycanMLread('tetra_test2.glycoct_xml'));
%        tetra{3}=GlycanSpecies(glycanMLread('tetra_test3.glycoct_xml'));
%        tetra{4}=GlycanSpecies(glycanMLread('tetra_test4.glycoct_xml'));
%        for i=1: length(tetra)
%           glycanViewer(tetra{i})
%        end
%        tetrabracket = merge2bracket(tetra);
%        glycanViewer(tetrabracket);
%
%
%See also groupbycomp.

% Author:Gang Liu and Yusen Zou
% Date: 11/16/14

if(iscell(isomerspecies))
    inputtype=1;
elseif(isa(isomerspecies,'GlycanSpecies') && length(isomerspecies)>1)
    inputtype=2;
elseif(isa(isomerspecies,'CellArrayList') && length(isomerspecies)>1)
    inputtype=3;    
else
    error('MATLAB:GLYCOPAT:WRONGINPUTTYPE','WRONG INPUT TYPE');
end

bracketspecies = [];
for j = 1 : length(isomerspecies)
    if(inputtype==2)
        ithisomerspecies = isomerspecies(j,1);
    elseif(inputtype==1)
        ithisomerspecies = isomerspecies{j};
    elseif(inputtype==3)
        ithisomerspecies = isomerspecies.get(j);    
    end
    bracketspecies = addbracket(bracketspecies,ithisomerspecies);
    if(isempty(bracketspecies))
        return
    end
end

end

function bracketspecies = addbracket(glycanspecies1withbracket,...
    glycanspeciestobracket)
if(isempty(glycanspecies1withbracket))
    bracketspecies = glycanspeciestobracket.clone;
    return
else
    glycanbracketObj = glycanspecies1withbracket.glycanStruct.clone;
    glycanbracketObj = subcommonstructplusbracket(glycanbracketObj,...
        glycanbracketObj.root,...
        glycanspeciestobracket.glycanStruct.root);
    bracketspecies = GlycanSpecies(glycanbracketObj);
end

end

function glycanbracketObj = subcommonstructplusbracket(glycanbracketObj,...
    residue1,residue2)

if(~residue1.typequal(residue2)) % is residue types of two residues are different
    glycanbracketObj.bracketResidue(residue1);
    glycanbracketObj.resetjava;
else
    numchild_residue1 = length(residue1.linkageChildren);
    if(numchild_residue1==0)
        return;
    end
    
    removeindex =[];
    for i = 1 : length(residue1.linkageChildren)
        equalnum = 0;
        ithchild1 = residue1.linkageChildren(i).child;
        for j = 1 : length(residue2.linkageChildren)
            ithchild2 = residue2.linkageChildren(j).child;
            if(ithchild1.typequal(ithchild2))
                if(isequal(ithchild1.linkageParent.bonds,ithchild2.linkageParent.bonds))
                    equalnum=equalnum+1;
                    glycanbracketObj = subcommonstructplusbracket(glycanbracketObj,...
                        ithchild1,ithchild2);
                    break
                end
            end
        end
        
        if(equalnum==0)
            removeindex = [removeindex;i];
        end
    end
    
    % remove any unmatched
    removepos = 0;
    for i =  1: length(removeindex)
        childpos = removeindex(i) - removepos;
        residuetoremove = residue1.linkageChildren(childpos).child;
        glycanbracketObj.bracketResidue(residuetoremove);
        glycanbracketObj.resetjava;
        removepos = removepos+1;
    end
end

end

