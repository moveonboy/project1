function glycantype=findglycantype(varargin)
%GROUPBYCOMP: group glycan species by their composition
% 
% Syntax:
%     groupspecies = findglycantype(glycanstruct)
%     groupspecies = findglycantype(glycanspecies)
%
% Input:
%     glycanstruct: GlyanStruct object
%     glycanspecies: GlycanSpecies object
%   
% Output:
%   glycantype: numerical value to represent type of the glycan
%      0: high mannose; 1: complex N-linked glycan; 2: hybrid; 3 to be
%      added 
%
% Example:
%     
%
% See also merge2bracket.

% Author: Gang Liu
% Date: 11/16/14
narginchk(1,1);
if(isa(varargin{1},'GlycanStruct'))
  glycanstructobj = varargin{1};  
elseif(isa(varargin{1},'GlycanSpecies'))
  glycanstructobj = varargin{1}.glycanStruct;  
else
   error('MATLAB:GNAT:ERRORINPUTTYPE','WRONG INPUT TYPE'); 
end

if(glycanstructobj.ishighmannose)
    glycantype = 0;
elseif(glycanstructobj.iscomplex)
    glycantype = 1;
elseif(glycanstructobj.ishybrid)
    glycantype = 2;
end

end