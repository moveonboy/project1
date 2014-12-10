function varargout = MScomparison(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,varargin)
%MScomparison: compare the two MS data based on the rules user defined
% themselves
%  
% Syntax: 
%    [uniqueglycan,relativeabundance] = optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2)
%    [uniqueglycan,relativeabundance,specificglycanabundance]  =...
%               optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,GlycanResidue,option)
%    [uniqueglycan,relativeabundance,specificresidueabundance] =...
%               optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,GlycanResidue,option)
%    [uniqueglycan,relativeabundance,specificstructabundance]  = ...
%               optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,Glycanstruct)
%    [uniqueglycan,relativeabundance,specificglycanabundance,...
%               specificresidueabundance] = optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,GlycanResidue,option)
%    [uniqueglycan,relativeabundance,specificglycanabundance,specificresidueabundance,...
%               specificstructabundance] = optimizebracket(MSdata1,MSdata2,pathwayfilename1,pathwayfilename2,GlycanResidue,Glycanstruct)
%
% Input: 
%    MSdata1: an object of GlycanStruct class
%    MSdata2: an object of GlycanSpecies class
%    option : 0, 1 or 2. If equals to 0, calculating the relativeabundance
%             of glycans with specific GlycanResidues; Else if equals to 1, calculating 
%             the relativeabundance of specific GlycanResidues; Else if
%             euqals to 2, calculating the both relativeabundance of glycans
%             and specific GlycanResidue;
%    GlycanResidue: residues that the glycan requiered to has.
%    Glycanstruct : structures that the glycan requiered to has.
%    
%
% Output:
%    uniqueglycan: an object of GlycanStruct class
%    relativeabundance:
%    specificglycanabundance:
%    specificresidueabundance:
%    specificstructabundance:
%
% Example 1:
%       
%       

%cAuthor: Yusen Zhou
% Date Lastly Updated: 12/9/14
narginchk(4,6);

[glycanDB1,glycanDBstrcomp1] = createLocalDataBase(pathwayfilename1);
[glycanDB2,glycanDBstrcomp2] = createLocalDataBase(pathwayfilename2);
glycanlist1 = msrelativeabundance(MSdata1,glycanDBstrcomp1);
glycanlist2 = msrelativeabundance(MSdata2,glycanDBstrcomp2);
glycaninMS1 = listglycansbyMS(glycanDB1,glycanlist1);
glycaninMS2 = listglycansbyMS(glycanDB2,glycanlist2);


uniqueglycan      = FindDiffGlycan(glycanlist1,glycanlist2);
relativeabundance = struct('MSdata1',[],'MSdata2',[]);
relativeabundance.MSdata1 = glycanlist1;
relativeabundance.MSdata2 = glycanlist2;
varargout{1} = uniqueglycan;
varargout{2} = relativeabundance;
if(nargin==5)
    if(isa(varargin{1},'GlycanStruct'))
        Glycanstruct = varargin{1};
    else
        error('Wrong Input Type');
    end
    specificstructabundance = calculatebystruct(glycaninMS1,glycaninMS2,Glycanstruct);% need to write another function to list species with glycan structure.
    varargout{3} = specificstructabundance;
elseif(nargin==6)
    if(isa(varargin{1},'GlycanResidue'))&&(isnumeric(varargin{2}))
        GlycanResidue = varargin{1};
        option        = varargin{2};
        specificglycanabundance  = calculateglycanbyresidue(glycanlist1,glycanlist2,GlycanResidue);
        specificresidueabundance = calculateresidue(glycanlist1,glycanlist2,GlycanResidue);
        if(option==0)
           varargout{3} = specificglycanabundance;
        elseif(option==1)
           varargout{3} = specificresidueabundance; 
        elseif(option==2)
           varargout{3} = specificglycanabundance;
           varargout{4} = specificresidueabundance; 
        end
    elseif(isa(varargin{1},'GlycanResidue'))&&(isa(varargin{2},'GlycanStruct'))
        GlycanResidue = varargin{1};
        GlycanStruct  = varargin{2};
        specificglycanabundance  = calculateglycanbyresidue(glycanlist1,glycanlist2,GlycanResidue);
        specificresidueabundance = calculateresidue(glycanlist1,glycanlist2,GlycanResidue);
        specificstructabundance  = calculatebystruct(glycaninMS1,glycaninMS2,GlycanStruct);% need to write another function to list species with glycan structure.
        varargout{3} = specificglycanabundance;
        varargout{4} = specificresidueabundance;
        varargout{5} = specificstructabundance;
    else
        error('Wrong Input Type');
    end
end
end

function specificglycanabundance = calculateglycanbyresidue(glycanlist1,glycanlist2,GlycanResidue)
specificglycanabundance = struct('MSdata1',[],'MSdata2',[],'containresidue',GlycanResidue);
glycanlist1abundance = addrequiredglycans(glycanlist1,GlycanResidue);
glycanlist2abundance = addrequiredglycans(glycanlist2,GlycanResidue);
specificglycanabundance.MSdata1 = glycanlist1abundance;
specificglycanabundance.MSdata2 = glycanlist2abundance;
end

function abundance = addrequiredglycans(glycanlist,GlycanResidue)
abundance    = 0;
for i = 1 : length(glycanlist)
    ithglycanstruct = glycanlist(1);
    if(~isempty(strfind(ithglycanstruct.glycancoms,GlycanResidue.name)))
        abundance = abundance+ithglycanstruct.relativeabundance;
    end
end
end

function specificresidueabundance = calculateresidue(glycanlist1,glycanlist2,GlycanResidue)
specificresidueabundance = struct('MSdata1',[],'MSdata2',[],'requiredresidue',GlycanResidue);
glycanlist1abundance = addrequiredresidues(glycanlist1,GlycanResidue);
glycanlist2abundance = addrequiredresidues(glycanlist2,GlycanResidue);
specificresidueabundance.MSdata1 = glycanlist1abundance;
specificresidueabundance.MSdata2 = glycanlist2abundance;
end

function abundance = addrequiredresidues(glycanlist,GlycanResidue)
abundance    = 0;
for i = 1 : length(glycanlist)
    ithglycanstruct = glycanlist(1);
    if(~isempty(strfind(ithglycanstruct.glycancoms,GlycanResidue.name)))
        numofrequiredresidues = length(strfind(ithglycanstruct.glycancoms,GlycanResidue.name));
        abundance = abundance+ithglycanstruct.relativeabundance*numofrequiredresidues;
    end
end
end

function glycaninMS = listglycansbyMS(glycanDB1,glycanlist1)
end

function specificstructabundance = calculatebystruct(glycaninMS1,glycaninMS2,GlycanStruct)
end