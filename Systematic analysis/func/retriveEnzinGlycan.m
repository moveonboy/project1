function glycanEnz = retriveEnzinGlycan(glycanSpeceis,varargin)
%retriveEnzinGlycan returns the infromation of enzymes acted in the
% glycanspcecies
% 
%
%glycanEnz = retriveEnzinGlycan(glycanSpeceis,varargin) loading glycan enzyme
% database and using one certain glycan as input to retrive the enzyme
% information.
% 
%Example:
%      m3gn = GlycanSpecies(glycanMLread('m3gn.glycoct.xml')
%      glycanEnz = retriveEnzinGlycan(m3gn)
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 10/27/14

% narginck(0,2)

if(isempty(varargin))
    load('HL60GlyEnzDB')
else
    glycanEnzDB = varargin{1};
    load(glycanEnzDB)
end

glycanname = glycanSpeceis.glycanStruct.name;
glycanEnz  = GlyEnz(glycanname);
end