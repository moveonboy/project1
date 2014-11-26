function isStructEqual = equalStruct(obj,obj2)
%equalGlycanStruct returns true if two glycans have the same structure
% equalGlycanStruct(GLYCANSTRUCTobj) compares the structures of two
%  glycans
% 
% Syntax: glycanStructObj1.equalGlycanStruct(glycanStructObj2);
%         equalGlycanStruct(glycanStructObj1,glycanStructObj2);
%
% Example:
%       load('GLYCAN1.mat')
%       load('GLYCAN2.mat')
%       glycanstruct1=glycan1.glycanStruct;
%       glycanstruct2=glycan2.glycanStruct
%       glycanViewer(glycanstruct1);
%       glycanViewer(glycanstruct2);
%       result1 = glycanstruct1.equalGlycanStruct(glycanstruct2);
%       result2 = glycanstruct1.equalStruct(glycanstruct2);
% 
% 
% Note: this comparison solved the bug in equalsStruct used in glycanbuilder
%  library.
%
% See also GlycanStruct

% Author: Gang Liu and Yushen Zhou
% Date Lastly Updated: 10/6/2014
if(isempty(obj2))
    isStructEqual = false;
    return
end

if(isempty(obj.root))
    isStructEqual = false;
    return
end

if(~obj.root.subtreequal(obj2.root))
    isStructEqual = false;
    return
end

if(isempty(obj.bracket))
    isStructEqual = isempty(obj2.bracket);
elseif(isempty(obj2.bracket))
    isStructEqual = false;
else
    isStructEqual = obj.bracket.subtreequal(obj2.bracket);
end

end

