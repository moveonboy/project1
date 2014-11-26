function branches = getBranch(obj,varargin)
%getBranches get information about the branch
%
% See also getDepth

if(length(varargin)==1)
   rootorbracket=varargin{1};
elseif(isempty(varargin))
   rootorbracket = 'root';
else
   error('MATLAB:GNAT:WRONGNUMBERPINT','INCORRECT NUMBER OF INPUTS');
end

branches = cell(0,1);
nonredendresidues_species1 =  obj.getNonRedEndResidue(rootorbracket);
for i = 1 : length(nonredendresidues_species1)
    ithresidue = nonredendresidues_species1{1,i};
    if(strcmp(ithresidue.residueType.name,'Fuc'))
       continue; 
    end
    
    branches(end+1).residueterminal  = nonredendresidues_species1{i}; 
    branches(end).depth            = obj.getBranchDepth(nonredendresidues_species1{1,i});
    branches(end).residuestrwolink = obj.getresiduetoroot(ithresidue,'nolinkage');
end

end

