function varargout = inferGHGlySubstr(glycanInputObj,enzObj,varargin)
%inferGHGlySubstr: Infer the substrate based on the glycosidase and its product.
%
% Syntax: 
%       substrSpecies = inferGHGlySubstr(prodObj,ghenz)
%       substrSpecies = inferGHGlySubstr(prodObj,ghenz,usebracket)
%       [numsubstrpecies,substrSpecies] = inferGHGlySubstr(prodObj,ghenz,usebracket)
%       [numsubstrpecies,substrSpecies,rxns] = inferGHGlySubstr(prodObj,ghenz,usebracket)
%       [numsubstrpecies,substrSpecies,rxns,pathway] = inferGHGlySubstr(prodObj,ghenz,usebracket)
%   
% Input:
%       prodObj: the glycan structure, either a GlycanStruct or GlycanSpecies object
%       ghenz: the glycosidase, a GHEnz object
%       usebracket: 1 for true, 0 for false
%    
% Output:
%       numsubstrpecies: the number of substrates
%       substrSpecies: the inferred substrates
%       rxns: the inferred reactions
%       pathway: the inferred pathway 
%
% Examples:
%      Example 1:
%            mani  = GHEnz.loadmat('mani.mat');
%            m5species  = glycanMLread('M5.glycoct_xml');
%            m6species = inferGHGlySubstr(m5species,mani,0);
%            options  = displayset('showmass',true,'showLinkage',true,...
%                             'showRedEnd',true);
%            for i = 1: m6species.length
%                glycanViewer(m6species.get(i).glycanStruct,options);
%            end
%
%      Example 2:
%            mani       = GHEnz.loadmat('mani.mat');
%            m5species  = glycanMLread('M5.glycoct_xml');
%            m6species  = inferGHGlySubstr(m5species,mani,1);
%            options    = displayset('showmass',true,'showLinkage',true,...
%                             'showRedEnd',true);
%            for i = 1: m6species.length
%                glycanViewer(m6species.get(i).glycanStruct,options);
%            end
%
%      Example 3:
%            mani       = GHEnz.loadmat('mani.mat');
%            m5species  = glycanMLread('M5.glycoct_xml');
%            [nsubstr,m6species,m5rxns] = inferGHGlySubstr(m5species,mani,1);
%            options  = displayset('showmass',true,'showLinkage',true,...
%                             'showRedEnd',true);
%            for i = 1: nsubstr
%               glycanRxnViewer(m5rxns.get(i));
%            end
%
%       Example 4:
%            mani       = GHEnz.loadmat('mani.mat');
%            m5species  = glycanMLread('M5.glycoct_xml');
%            [nsubstr, m7species,m6rxns,m6pathway] = inferGHGlySubstr(m5species,mani,1);
%            glycanPathViewer(m6pathway);
%
%      Example 5:
%            mani  = GHEnz.loadmat('mani.mat');
%            m6species  = glycanMLread('m6bracket.glycoct_xml');
%            [nsubstr, m7species,m7rxns,m7pathway] = inferGHGlySubstr(m6species,mani,1);
%            glycanPathViewer(m7pathway);
%
%      Example 6:       
%            mani  = GHEnz.loadmat('mani.mat');
%            m7species  = glycanMLread('m7bracket.glycoct_xml');
%            [nsubstr, m8species,m8rxns,m8pathway] = inferGHGlySubstr(m7species,mani,1);
%            glycanPathViewer(m8pathway);
%            
%      Example 7:
%            mani        = GHEnz.loadmat('mani.mat');
%            m8species  = glycanMLread('m8bracket.glycoct_xml');
%            [nsubstr, m9species,m9rxns,m9pathway] = inferGHGlySubstr(m8species,mani,1);
%            glycanPathViewer(m9pathway);
%
% See also inferGlyProd, inferGlyRevrPath.

% Author: Yusen Zhou and Gang Liu
% Date Last Updated: 11/23/14

if(isa(glycanInputObj,'GlycanStruct'))
    glycanObj = glycanInputObj;
elseif(isa(glycanInputObj,'GlycanSpecies'))
    glycanObj = glycanInputObj.glycanStruct;
else
    error('MATLAB:GNAT:ERRORINPUTTYPE','INCORRECT INPUT TYPE');
end

if(isempty(varargin))
    usebracket = 0;
elseif(length(varargin)==1)
    if(varargin{1})
       usebracket=1;
    else
       usebracket=0;
    end
else
    error('MATLAB:GNAT:ERRORNUMBERINPUT','WRONG INPUT NUMBER')
end

substrSpecies = CellArrayList;
if(nargout==1)
    varargout{1} = CellArrayList;
elseif(nargout==2)
    varargout{1} = 0;
    varargout{2} = CellArrayList;
elseif(nargout==3)
    varargout{1} = 0;
    varargout{2} = substrSpecies;
    varargout{3} = CellArrayList;
elseif(nargout==4)
    varargout{1} = 0;
    varargout{2} = CellArrayList;
    varargout{3} = CellArrayList;
    varargout{4} = Pathway;
end

%rule 1: N-, O-, glycolipids, wt/wo brackets
if(~isempty(enzObj.glycanTypeSpec))
  if(~strcmpi(glycanObj.glycanTypeSpec,...
     enzObj.glycanTypeSpec))
     return
  end
end

if(~isempty(glycanObj.bracket))
    prodbracket   = 1;
    substrbracket = 1;
elseif(usebracket)
    prodbracket   = 0; 
    substrbracket = 1; % inference
else
    prodbracket   = 0;
    substrbracket = 0; % allow multiple substrates 
end

if(prodbracket)
    optimizeglycanobj = glycanObj.clone;
    glycanObj = optimizebracket(optimizeglycanobj);
end

% rule 1: subtrateNAresidue
if(~isempty(enzObj.substNAResidue))
    % prodMinStruct does not have bracket
    if(sum(strcmpi(enzObj.substNAResidue.name,...
            fieldnames(glycanObj.getComposition)))~=0)
        return;
    end
end

allResidues = glycanObj.getAllResidues;  % obtain all residues including the ones in the bracket
for i = 1 : length(allResidues)
    nreResidue = allResidues{1,i};
    if(strcmpi(nreResidue.residueType.name,'freeEnd'))
        continue;
    end
    
    if(prodbracket)
        bracketresidue    = glycanObj.getAllResidues('bracket');
        inBracket = 0;
        for j = 1 : length(bracketresidue)
            jthnreResidue = bracketresidue{1,j};
            if(isequal(jthnreResidue.IDres,nreResidue.IDres))
                inBracket = 1;
                break;
            end
        end
        
        if(inBracket)
            if(strcmp(nreResidue.residueType.name,'#bracket'))
                continue
            end
            substrSpecies = ProceedByRule_withbracket(nreResidue,enzObj,...
                glycanObj,i,substrSpecies);
        else
            substrSpecies = ProceedByRule(nreResidue,enzObj,...
                glycanObj,i,substrSpecies);
        end
        
        if(isempty(substrSpecies))
           continue
        end
    else
        substrSpecies = ProceedByRule(nreResidue,enzObj,...
               glycanObj,i,substrSpecies);
    end
end

if(substrbracket &&...
        substrSpecies.length>1)
    substrbracketspecies = merge2bracket(substrSpecies);
    for i = 1 : length(substrbracketspecies)
        substrbracketspecies(i) = optimizebracket(substrbracketspecies(i));
    end
    substrSpecies = CellArrayList;
    for i = 1 : length(substrbracketspecies)
        substrSpecies.add(substrbracketspecies(i));
    end
end

if(nargout>=3)
    rxns = CellArrayList;
    prodSpecies = GlycanSpecies(glycanObj);
    for i = 1 : length(substrSpecies)
        rxns.add(Rxn(substrSpecies.get(i),prodSpecies,enzObj));
    end
end

if(nargout==1)
    varargout{1} = substrSpecies;
elseif(nargout==2)
    varargout{1} = substrSpecies.length;
    varargout{2} = substrSpecies;
elseif(nargout==3)
    varargout{1} = substrSpecies.length;
    varargout{2} = substrSpecies;
    varargout{3} = rxns;
elseif(nargout==4)
    path=Pathway;
    if(substrSpecies.length>=1)
        path.addGlycans(substrSpecies);
        path.addGlycan(GlycanSpecies(glycanObj));
        path.addRxns(rxns);
    end
    varargout{1}=substrSpecies.length;
    varargout{2}=substrSpecies;
    varargout{3}=rxns;
    varargout{4}=path;
end
end

function ismatch=fuzzymatch(str1,str2)
if(length(str1)~=length(str2))
    ismatch = 0;
    return
end

if(isempty(strfind(str1,'?')))&&(isempty(strfind(str2,'?')))
    ismatch  = strcmp(str1,str2);
elseif(~isempty(strfind(str1,'?')))
    str2(strfind(str1,'?')) = '?';
    ismatch  = strcmp(str1,str2);
elseif(~isempty(strfind(str2,'?')))
    str1(strfind(str2,'?')) = '?';
    ismatch  = strcmp(str2,str1);
end

end

function substrSpecies = ProceedByRule_withbracket(nreResidue,enzObj,...
    glycanObj,i,substrSpecies)
% rule 2: residue next to the functional group
if(~isempty(enzObj.resAtt2FG))
    if(~fuzzymatch(enzObj.dispAttachResLink,...
            nreResidue.dispResidueInfo))
        return
    end
end

% rule 3: the product should not be greater than the specified structure
if(~isempty(enzObj.prodMinStruct))
    % prodMinStruct does not have bracket
    if(~glycanObj.bracketContain(enzObj.prodMinStruct))
        return
    end
end

% rule 4: the product should not work on the children of specificed structure
if(~isempty(enzObj.prodMaxStruct))
    % prodMaxStruct does not have bracket, glycanObj can have
    % bracket
    if(~enzObj.prodMaxStruct.bracketContain(glycanObj))
        return
    end
end

% rule 5: the enzyme should remove the residue from the target branch
% containing specified structure
if(~isempty(enzObj.targetbranchcontain))
    % should not
    chartargetbranch = regexprep(enzObj.targetbranchcontain.name,'freeEnd','');
    charresidue      = glycanObj.getresiduetoroot(nreResidue.getParent);
    charresidue      = regexprep(charresidue,'freeEnd','');
    isValidRxnPos    = (~isempty(strfind(chartargetbranch,charresidue)));
    if(~isValidRxnPos)
        return;
    end
end

nBondChoice = length(enzObj.linkFG.bond);
glycanSubstStructArray = CellArrayList;
isOneResiduedeleted = false;

if(nBondChoice>=1)
    for j = 1 : nBondChoice
        bond                 = enzObj.linkFG.bond(j,1);
        glycanSubstrObj      = glycanObj.clone;
        nreResidues          = glycanSubstrObj.getAllResidues;
        nreResidue           = nreResidues{1,i};
        isResidueAdded       = glycanSubstrObj.addResidue(nreResidue,...
            enzObj.resfuncgroup,enzObj.linkFG.anomer,...
            bond);
        if(isResidueAdded)
            isOneResiduedeleted = true;
            glycanSubstStructArray.add(glycanSubstrObj);
        end
    end
end

if(~isOneResiduedeleted)
    return;
end

for k = 1: glycanSubstStructArray.length
    glycanSubstrObj = glycanSubstStructArray.get(k);
    
    % rule 6, the enzyme should not work on the parent of specified structure
    if(~isempty(enzObj.substMaxStruct))  % if the rule exists
        if(~enzObj.substMaxStruct.bracketContain(glycanSubstrObj))
            continue;
        end
    end
    
    % rule 7, the enzyme should not work on the children of specified structure
    if(~isempty(enzObj.substMinStruct))  % if the rule exists
        if(~glycanSubstrObj.bracketContain(enzObj.substMinStruct))
            continue;
        end
    end
    
    substr = GlycanSpecies(glycanSubstrObj);
    substrSpecies.add(substr);
end

end

function substrSpecies = ProceedByRule(nreResidue,enzObj,...
    glycanObj,i,substrSpecies) 
% rule 2: residue next to the functional group
if(~isempty(enzObj.resAtt2FG))
    if(~fuzzymatch(enzObj.dispAttachResLink,...
            nreResidue.dispResidueInfo))
        return
    end
end

% rule 3: the product should not be greater than the specified structure
if(~isempty(enzObj.prodMinStruct))
    % prodMinStruct does not have bracket
    if(~glycanObj.bracketContain(enzObj.prodMinStruct))
        return;
    end
end

% rule 4: the product should not work on the children of specificed structure
if(~isempty(enzObj.prodMaxStruct))
    % prodMaxStruct does not have bracket, glycanObj can have
    % bracket
    if(~enzObj.prodMaxStruct.bracketContain(glycanObj))
        return;
    end
end

% rule 5: the enzyme should remove the residue from the target branch
% containing specified structure
if(~isempty(enzObj.targetbranchcontain))
    chartargetbranch = regexprep(enzObj.targetbranchcontain.name,'freeEnd','');
    charresidue      = glycanObj.getresiduetoroot(nreResidue.getParent);
    charresidue      = regexprep(charresidue,'freeEnd','');
    isValidRxnPos    = (~isempty(strfind(chartargetbranch,charresidue)));
    if(~isValidRxnPos)
        return;
    end
end

% rule 6: the enzyme should not act on the specified branch
if(~isempty(enzObj.targetNABranch))
    if(isa(enzObj.targetNABranch,'GlycanStruct'))
        chartargetbranch    = enzObj.targetNABranch.name;
        charresidue         = glycanObj.getresiduetoroot(...
            nreResidue.getParent);
        isValidRxnPos       = isValidRxnPos && (~strcmp(...
            chartargetbranch,charresidue));
    elseif(isa(enzObj.targetNABranch,'CellArrayList'))
        for ii = 1 : length(enzObj.targetNABranch)
            chartargetbranch = enzObj.targetNABranch.name;
            charresidue      = glycanObj.getresiduetoroot(...
                nreResidue.getParent);
            isValidRxnPos    = isValidRxnPos && ...
                (~strcmp(chartargetbranch,charresidue));
        end
    end
    
    if(~isValidRxnPos)
        return;
    end
end
nBondChoice = length(enzObj.linkFG.bond);
glycanSubstStructArray = CellArrayList;
isOneResiduedeleted = false;

if(nBondChoice>=1)
    for j = 1 : nBondChoice
        bond                 = enzObj.linkFG.bond(j,1);
        glycanSubstrObj      = glycanObj.clone;
        nreResidues          = glycanSubstrObj.getAllResidues;
        nreResidue           = nreResidues{1,i};
        isResidueAdded       = glycanSubstrObj.addResidue(nreResidue,...
            enzObj.resfuncgroup,enzObj.linkFG.anomer,...
            bond);
        if(isResidueAdded)
            isOneResiduedeleted = true;
            glycanSubstStructArray.add(glycanSubstrObj);
        end
    end
end

if(~isOneResiduedeleted)
    return;
end

for k = 1: glycanSubstStructArray.length
    glycanSubstrObj = glycanSubstStructArray.get(k);
    
    % rule 7 the substrate should not contain the specified branch
    if(~isempty(enzObj.substNABranch))
        if(isa(enzObj.substNABranch,'GlycanStruct'))
            isValidRxnPos = ~glycanSubstrObj.contains(enzObj.substNABranch);
        elseif(isa(enzObj.substNABranch,'CellArrayList'))
            isValidRxnPos = true;
            for ii = 1 : length(enzObj.substNABranch)
                isValidRxnPos = isValidRxnPos&&...
                    (~glycanSubstrObj.contains(enzObj.substNABranch.get(ii)));
            end
        end
        
        if(~isValidRxnPos)
            continue;
        end
    end
    
    % rule 8, the enzyme should not work on the parent of specified structure
    if(~isempty(enzObj.substMaxStruct))  % if the rule exists
        if(~enzObj.substMaxStruct.bracketContain(glycanSubstrObj))
            continue
        end
    end
    
    % rule 9, the enzyme should not work on the children of specified structure
    if(~isempty(enzObj.substMinStruct))  % if the rule exists
        if(~glycanSubstrObj.bracketContain(enzObj.substMinStruct))
            continue;
        end
    end
    
    substr    = GlycanSpecies(glycanSubstrObj);
    substrSpecies.add(substr);
end

end
