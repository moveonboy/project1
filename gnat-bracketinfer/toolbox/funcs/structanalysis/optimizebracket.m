function glystructObj=optimizebracket(varargin)
%OPTIMIZEBRACKET: optimize the bracket position based on glycan structure
%  
% Syntax: 
%    optimizedglycanstruct = optimizebracket(GlycanStructObj)
%    optimizedglycanstruct = optimizebracket(GlycanSpeciesObj)
%    
% Input: 
%    GlycanStructObj: an object of GlycanStruct class
%    GlycanSpeciesObj: an object of GlycanSpecies class
%
% Output:
%    optimizedglycanstruct: an object of GlycanStruct class
%
% Examples:
%    Example 1:
%       m8Struct    = glycanMLread('m8bracket.glycoct_xml');
%       m8newStruct = optimizebracket(m8Struct)
%       glycanViewer(m8Struct)
%       glycanViewer(m8newStruct)
%
%    Example 2:
%
%See also bracket2Exact.

% Author: Yushen Zhou && Gang Liu
% Date: 11/23/2014

if(length(varargin)==1)
   if(isa(varargin{1},'GlycanSpecies'))
       glystructObj     = varargin{1}.glycanStruct.clone;
   elseif(isa(varargin{1},'GlycanStruct'))
       glystructObj     = varargin{1}.clone;
   else
       error('MATLAB:GNAT:WRONGINPUTTYPE','INCORRECT INPUT TYPE');
   end  
else
   error('MATLAB:GNAT:WRONGINPUTNUMBER','INCORRECT NUMBER OF INPUTS');
end

rootglycanbranches    = glystructObj.getBranch('root');
bracketglycanbranches = glystructObj.getBranch('bracket');
branchnum             = length(rootglycanbranches);
bracketbranchnum      = length(bracketglycanbranches);

branchdepth = zeros(length(rootglycanbranches),1);
for i = 1 : length(rootglycanbranches)
    branchdepth(i)=rootglycanbranches(i).depth;
    nreresidue{1,i}= rootglycanbranches(i).residueterminal;
end

if(min(branchdepth)==max(branchdepth))
    isrootbranchequallength = 1;
else
    isrootbranchequallength = 0;
end

bracketbranchequal = 1;
bracketbranchstr = bracketglycanbranches(1).residuestrwolink;
for i = 1 : length(bracketglycanbranches)
   bracketbranchequal = bracketbranchequal && (...
     fuzzymatch(bracketglycanbranches(i).residuestrwolink,bracketbranchstr));
   bracketresidue{1,i}= bracketglycanbranches(i).residueterminal;
end

if(~bracketbranchequal)||(branchnum~=bracketbranchnum)
   return
end

if((branchnum==bracketbranchnum) && bracketbranchequal)
    if(isrootbranchequallength) 
     % since the bracket is not required, remove the bracket
        for i = 1 : length(nreresidue)
            glystructObj.addResidue(nreresidue{1,i},...
                bracketresidue{1,i});
        end
        glystructObj.removeBracket;
    else  %
        maxrootbranchdepth = max(branchdepth);
        for i = 1 : length(branchdepth)
            if(branchdepth(i)~=maxrootbranchdepth)
                bracketresiduetoremove = bracketresidue{1,i};
                glystructObj.moveResiduesfrombracket2root(...
             bracketresiduetoremove,maxrootbranchdepth-branchdepth(i),...
             nreresidue{1,i})                
            end
        end
    end
end

end

function ismatch=fuzzymatch(str1,str2)
if(length(str1)~=length(str2))
    ismatch = 0;
    return
end

if(isempty(strfind(str1,'?')))&&(isempty(strfind(str2,'?')))
    ismatch  = strcmp(str1,str2);
else
    str2(strfind(str1,'?')) = '?';
    str1(strfind(str2,'?')) = '?';
    ismatch  = strcmp(str1,str2);
end

end