function varargout = inferGlyConnPath_bracket(glycanArray,enzObjArray,varargin)
%inferGlyConnPath infer the pathway to connect all glycans based on enzymes.
%
%  [isPathwayFormed,pathway] = inferGlyConnPath(glycanArray,enzArray) infers
%    the pathway connecting all glycans. If no pathway is
%    formed, isPathwayFormed is false.
%
%  [isPathwayFormed,pathway] = inferGlyConnPath(glycanArray,enzArray,
%    'iterativeplot',isdisplay) displays the network at each iterative
%    step if isdisplay is set as true.
%
%    Example: 
%       enzArray = CellArrayList;
%       mgat2     = GTEnz.loadmat('mgat2.mat'); enzArray.add(mgat2);
%       mgat3     = GTEnz.loadmat('mgat3.mat'); enzArray.add(mgat3);
%       mgat4     = GTEnz.loadmat('mgat4.mat'); enzArray.add(mgat4);
%       mgat5     = GTEnz.loadmat('mgat5.mat'); enzArray.add(mgat5);
%       galt          = GTEnz.loadmat('galt.mat'); enzArray.add(galt);
%       glycanArray = CellArrayList; 
%       tetraprimglycan    = GlycanSpecies(glycanMLread('tetriprime_gal_nlinkedglycan.glycoct_xml')) ;
%       m3gngn       = GlycanSpecies(glycanMLread('m3gngn.glycoct_xml')) ;
%       glycanArray.add(tetraprimglycan);
%       glycanArray.add(m3gngn);
%       fprintf(1,'Input of glycan product structure is \n');
%       glycanViewer(tetraprimglycan.glycanStruct);
%       glycanViewer(m3gngn.glycanStruct);
%       [isPath,nlinkedpath]=inferGlyConnPath(glycanArray, enzArray,'iterativedisp',false);
%       fprintf(1,'Inferred network is show below:\n'); 
%       glycanPathViewer(nlinkedpath);
%
% See also inferGlyForwPath,inferGlyRevrPath.

% Author: Gang Liu
% Date Last Updated: 9/18/13

narginchk(2,4);

if(nargin==3)
    error('MATLAB:GNAT:InputNumberWrong','Wrong Input number');
end

if(~isa(glycanArray,'CellArrayList'))
    error('MATLAB:GNAT:InputTypeWrong','Wrong Input Type');
end

if(~isa(enzObjArray,'CellArrayList'))
    error('MATLAB:GNAT:InputTypeWrong','Wrong Input Type');
end

iterativedisplay = false;

if(nargin==4)
    optionformat = varargin{1};
    if(strcmp(optionformat,'iterativedisp'))
        iterativedisplay = varargin{2};
    else
        error('MATLAB:GNAT:InputTypeWrong','Wrong Input Type');
    end
end

nargoutchk(1,2);

glypath         = Pathway;
glycanPairArray = identifyGlycanRole(glycanArray,enzObjArray);

for i = 1 : length(glycanPairArray)
    fprintf(1,'the round: %i\n',i);
    newglypath = Pathway;
    
    theithpair = glycanPairArray.get(i);
    if(checkjglycanPairInPath(glypath,theithpair))
        continue;
    end
    
    newglypath     = buildconnpath(newglypath,theithpair,enzObjArray,iterativedisplay);
    numSpecieAfter = newglypath.getNSpecies;
    
    if(numSpecieAfter>0)
        removeIsolatedSpecies(newglypath);
        % glycanPathViewer(newglypath);
        try
            newglypath2=pathfinding(newglypath,theithpair);
        catch err
            if(strcmp(err.identifier,'MATLAB:GNAT:WrongInput'))
                newglypath2 = [];
            else
                rethrow(err);
            end
        end
        
        %[isvalidpath, newglypath] = removendisnodes(newglypath,glycanPairArray.get(i));
        if(~isempty(newglypath2))
            %glycanPathViewer(newglypath2);
            glypath.addGlyPathByStruct(newglypath2);
            fprintf(1,'the number of total species in the pathway: ');
            disp(num2str(glypath.theSpecies.length));
            fprintf(1,'the number of total reactions in the pathway: ');
            disp(num2str(glypath.theRxns.length)); 
            if(iterativedisplay)
                glycanPathViewer(glypath)
                glycanPathViewer(newglypath2)
            end
        end
        
    end
end

isvalidpathformed = (glypath.theSpecies.length>=2);

if(nargout==1)
    varargout{1} = isvalidpathformed;
elseif(nargout==2)
    varargout{1} = isvalidpathformed;
    varargout{2} = glypath;
end

glypath.clearjava;
clear JAVA
glypath.resetjava;
end

function isglycanPairInPath = checkjglycanPairInPath(glycanpath,glycanpair)
reac = glycanpair.react;
prod = glycanpair.prod;

isglycanPairInPath= glycanpath.isStructinPath(reac)....
    && glycanpath.isStructinPath(prod) ;
end

function glycanpath = buildconnpath(glycanpath, glycanpair,enzObjArray,iterativedisplay)
Yprod    = glycanpair.prod;
XSubstr  = glycanpair.react;

if(glycanpath.isStructinPath(Yprod))
    yprodspecies = glycanpath.theSpecies.get(...
        glycanpath.findsameStructGlycan(Yprod));
    if((yprodspecies.listOfReacRxns.length>0) && ...
            (yprodspecies.listOfProdRxns.length>0))
        return;
    end
end

% Check if Yprod is bracketspeices.
allresiduesforY = Yprod.glycanStruct.getAllResidues;
isbracket       = 0;
for i = 1 : length(allresiduesforY)
    ithresidue = allresiduesforY{i};
    if(isequal(ithresidue.residueType.name,'#bracket'))
        isbracket = 1;
    end
end

% isProdadded = 0;
for i = 1 : length(enzObjArray)
    if(~isbracket)
        [numSubstr,substrSpecies,rxns] = inferGlySubstr(Yprod,enzObjArray.get(i));
    elseif(isbracket)
        if(isa(enzObjArray.get(i),'GHEnz'))
            continue
        end
        [numSubstr,substrSpecies,rxns] = inferBracketGlySubstr(Yprod,enzObjArray.get(i),enzObjArray);
    end
    for j = 1 : numSubstr
        thejthrxn = rxns.get(j);
        substr = substrSpecies.get(j);
        if(substr.glycanStruct.equalStruct(XSubstr.glycanStruct))
            glycanpath.addRxnByStruct(thejthrxn);
        else
            [isPossibleYisXProd,glycanpair]=findYisXProd (substr,XSubstr,enzObjArray);
            if(isPossibleYisXProd)
                glycanpath.addRxnByStruct(thejthrxn);
                glycanpath = buildconnpath(glycanpath,glycanpair,enzObjArray,iterativedisplay);
                if((iterativedisplay)&&(glycanpath.getNSpecies>0))
                      glycanPathViewer(glycanpath);
                end                  
            end
        end
    end
end

end

function glycanPairArray=identifyGlycanRole(glycanArray,enzObjArray)
glycanPairArray = CellArrayList;
for i = 1  :  length(glycanArray)-1
    glycanX = glycanArray.get(i);
    for j  = i+1 :  length(glycanArray)
        glycanY= glycanArray.get(j);
        [isPossibleYisXProd,glycanpair]= findYisXProd (glycanY,glycanX, enzObjArray);
        if(isPossibleYisXProd)
            glycanPairArray.add(glycanpair);
        end
    end
end
end

function [isPossibleProd,glycanpair]= findYisXProd (glycan1,glycan2, enzObjArray)
if(isequal(glycan1,glycan2))
    isPossibleProd = false;
    glycanpair =[];
    return
end

glycan1struct = glycan1.glycanStruct;
glycan2struct = glycan2.glycanStruct;
isPossibleYisXProd = glycan1struct.isProdOfX(glycan2struct,enzObjArray);
isPossibleXisYProd = glycan2struct.isProdOfX(glycan1struct,enzObjArray);

if(isPossibleYisXProd)
    react = glycan2;
    prod = glycan1;
elseif(isPossibleXisYProd)
    react = glycan1;
    prod = glycan2;
end

isPossibleProd =isPossibleYisXProd||isPossibleXisYProd;
if(isPossibleProd)
    glycanpair.react  =  react;
    glycanpair.prod  =  prod;
else
    glycanpair=[];
end

end

