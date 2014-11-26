function obj = bracketAddResidue(obj,residuetobracket)
%BRACKETADDRESIDUE: Add a residue to the bracket
%
% Syntax: 
%    glycanObj.bracketAddResidue(residuetobracket)
%
% Input: 
%    glycanObj: an object of glycan structure
%    residuetobracket: residue for addition to the bracket
%
% Example:
%    
% 
% See also ADDRESIDUES

%Author: Gang Liu
%Date Lastly Updated: 11/16/2014

for j = 1 : length(residuetobracket)
    jthresiduetobracket = residuetobracket(j,1);
%     if(~isa(jthresiduetobracket,'GlycanResidue'))
%         errorReport(mfilename,'IncorrectInputType');
%     end
    % change parent linkage to the residue to be removed
    linkageParent =  jthresiduetobracket.getLinkageParent;
    parentResidue =  linkageParent.getParent;
    
    linkageChildren = parentResidue.getLinkageChildren;
    nChildren       = length(linkageChildren);
    
    i=1;
    isFindChild=false;
    terminalResiduePos =-1;
    while((i<=nChildren) &&(~isFindChild))
        ithchild = linkageChildren(i,1).getChild;
        if(ithchild==jthresiduetobracket)
            isFindChild=true;
            terminalResiduePos = i;
        end
        i=i+1;
    end
    if(terminalResiduePos==-1)
        error('MATLAB:GNAT:ERRORRESIDUE','RESIDUE INFO INCORRECT');
    end
    parentResidue.unsetLinkageChildren(terminalResiduePos);
end

if(isempty(glycanbracketObj.bracket))
   obj.addBracket;
   obj.addResiduesToBracket(residuetobracket);
else
   obj.addResiduesToBracket(residuetobracket); 
end
    

% update name
obj.glycanjava = obj.structMat2Java;
obj.name = char(obj.glycanjava.toStringOrdered(0));
end
        
