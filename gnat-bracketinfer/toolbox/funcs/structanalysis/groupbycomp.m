function listofgroupspecies = groupbycomp(listofSpecies)
%GROUPBYCOMP: group glycan species depending on their compositions
% 
% Syntax:
%      groupspecies = groupbycomp(listofSpecies)
%
% Input:
%      listofSpecies: list of species, variable type: Cell
% 
% Output:
%       listofgroupspecies: list of grouped species, each group containing 
%          different number of isomers, variable type: Cell 
%   
% Example:
%       listofSpecies{1} = GlycanSpecies(glycanMLread('M6a.glycoct_xml'));
%       listofSpecies{2} = GlycanSpecies(glycanMLread('M6b.glycoct_xml'));
%       listofSpecies{3} = GlycanSpecies(glycanMLread('M6c.glycoct_xml'));
%       listofSpecies{4} = GlycanSpecies(glycanMLread('M7a.glycoct_xml'));
%       listofSpecies{5} = GlycanSpecies(glycanMLread('M7b.glycoct_xml'));
%       listofSpecies{6} = GlycanSpecies(glycanMLread('M7c.glycoct_xml'));
%       listofSpecies{7} = GlycanSpecies(glycanMLread('M9.glycoct_xml'));
%       listofgroupspecies = groupbycomp(listofSpecies);
%
% See also merge2bracket.

% Author: Gang Liu,Yusen Zou
% Date: 11/16/14

numspecies = length(listofSpecies);
glycancompstringlist = cell(numspecies,1);
for i = 1 : numspecies
  glycancompstringlist{i} = listofSpecies{i}.glycanStruct.getCompositionString;    
end

uniqueglycomp            = unique(glycancompstringlist);
numgroupspecies          = length(uniqueglycomp);
listofgroupspecies       = cell(numgroupspecies,1);
for i = 1: numgroupspecies
   ithglycomp            = uniqueglycomp{i};
   listfound             = strcmp(ithglycomp,glycancompstringlist);
   listofgroupspecies{i} = listofSpecies(listfound); 
end

end

