function substrSpecies = inferGTGlySubstr(prodObj,enzObj,usebracket)
% inferGTGlySubstr infer the substrate based on the GT enzyme and product.
%
%  substrSpecies = inferGTGlySubstr(prodObj,enz,usebracket)
%   returns the substrSpecies if the enzyme acts to form the product prodObj. If
%   no substrate is found, substrSpecies is set empty.
%
%      Example 1:
%             residueMap=load('residueTypes.mat');
%             b4GalI                   = GTEnz([2;4;1;38]);
%             b4GalI.resfuncgroup      = residueMap.allresidues('Gal');
%             b4GalI.resAtt2FG         = residueMap.allresidues('GlcNAc');
%             fucbond                  = GlycanBond('4','1');
%             b4GalI.linkFG            = struct('anomer','b','bond',fucbond);
%             gnbond                   = GlycanBond('?','1');
%             b4GalI.linkresAtt2FG     = struct('bond',gnbond,'anomer','b');
%             gn2m3gn                  = glycanMLread('1416.71.glycoct_xml');
%             b4GalI.substMinStruct    = gn2m3gn;
%             b4GalI.targetNABranch    = glycanMLread('b4galitargetNAbranch.glycoct_xml');
%             tetratest1_bracket       = GlycanSpecies(glycanMLread('tetratest1_bracket.glycoct_xml'));
%             substrSpecies = inferGTGlySubstr(tetratest1_bracket,b4GalI,1);
%             options    = displayset('showmass',true,'showLinkage',true,...
%                'showRedEnd',true);
%             glycanViewer(tetratest1_bracket,options)
%             for i = 1: length(substrSpecies)
%                glycanViewer(substrSpecies.get(i).glycanStruct,options);
%             end
%
%     Example 2:
%             residueMap=load('residueTypes.mat'); 
%             Futa23                     = GTEnz([2;4;1;152]);
%             Futa23.isTerminalTarget    = false;
%             Futa23.resfuncgroup        = residueMap.allresidues('Fuc');
%             fuctbond                   = GlycanBond('3','1');
%             Futa23.linkFG              = struct('anomer','a','bond',fuctbond);
%             glcnacResType              = residueMap.allresidues('GlcNAc');
%             glcnacBond                 = GlycanBond('?','1');
%             Futa23.resAtt2FG           = glcnacResType;
%             Futa23.linkresAtt2FG       = struct('bond', glcnacBond,'anomer','b');
%             Futa23.targetbranchcontain = glycanMLread('lactose.glycoct_xml');
%             Futa23.isTerminalTarget    = false; 
%             fuctestglycan              = GlycanSpecies(glycanMLread('fuctest_bracket.glycoct_xml'));
%             substrSpecies              = inferGTGlySubstr(fuctestglycan,Futa23,1);
%             options    = displayset('showmass',true,'showLinkage',true,...
%                'showRedEnd',true);
%             glycanViewer(fuctestglycan,options)
%             for i = 1: length(substrSpecies)
%                glycanViewer(substrSpecies.get(i).glycanStruct,options);
%             end
%
%
% See also inferGHGlyProd,inferGlyRevrPath.

% Author: Gang Liu,Yusen Zhou
% Date Last Updated: 11/18/14

substrSpecies = CellArrayList;
numSubstr     = 0;

glycanObj     = prodObj.glycanStruct;
isSubstrValid = true;
% rule 1: N-, O-, glycolipids
if(~isempty(enzObj.glycanTypeSpec))
    requirementglycantype     = enzObj.glycanTypeSpec;
    acceptortype    = glycanObj.glycanTypeSpec;
    isSubstrValid   = isequal(acceptortype,requirementglycantype);
end

if(~isempty(glycanObj.bracket))
    prodbracket   = 1;
    substrbracket = 1;
elseif(usebracket)
    prodbracket   = 0;
    substrbracket = 1; % inference
else
    prodbracket   = 0;
    substrbracket = 0; % multiple substrates if possible
end

nreResidues = glycanObj.getNonRedEndResidue;
numTermRes  =  length(nreResidues);

if(isSubstrValid)
    funcresrequirement    = enzObj.dispFuncResLink;
    for i = 1 : numTermRes
        nreResidue          = nreResidues{1,i};
        
        % rule 2: terminal residue must be function group defined in the
        % enzyme obj
        terminalInfo   = nreResidue.dispResidueInfo;
        isValidRxnPos = fuzzymatch(terminalInfo,funcresrequirement);
        if(~isValidRxnPos)
            continue;
        end
        
        if(prodbracket)
            bracketresidue    = getNonRedEndResidue(glycanObj,'bracket');
            inBracket = 0;
            for j = 1 : length(bracketresidue)
                jthnreResidue = bracketresidue{1,j};
                if(isequal(jthnreResidue.IDres,nreResidue.IDres))
                    inBracket = 1;
                end
            end
            if(inBracket)
               [numSubstr,substrSpecies] = ProceedByRule_withbracket(nreResidue,...
                   enzObj,glycanObj,i,substrSpecies,numSubstr,bracketresidue); 
            else
                for j = 1 : length(bracketresidue)
                    jthnreResidue  = bracketresidue{1,j};
                    isincorestruct = 1;
                    if(isempty(enzObj.isTerminalTarget))
                        enzObj.isTerminalTarget = 1;
                    end
                    if(isequal(jthnreResidue.residueType.name,enzObj.resfuncgroup.name))&&...
                            (enzObj.isTerminalTarget)
                        isincorestruct = 0;
                        break
                    end
                end
                if(isincorestruct)
                    [numSubstr,substrSpecies] = ProceedByRule(nreResidue,enzObj,glycanObj,i,substrSpecies,numSubstr);
                else
                    continue
                end
            end
        else
            [numSubstr,substrSpecies] = ProceedByRule(nreResidue,enzObj,glycanObj,i,substrSpecies,numSubstr);
        end
    end
    
    % Check if need to add bracket.
    if(length(substrSpecies)>1)&&(~prodbracket)&&(substrbracket)
        substrSpecies = isaddbracket(substrSpecies,numSubstr,enzObj);
    elseif(prodbracket)&&(length(substrSpecies)>1)
        substrSpecies = reorgnizestructure(substrSpecies,enzObj);
    end
end
end

function ismatch=fuzzymatch(str1,str2)
if(length(str1)~=length(str2))
    ismatch = 0;
    return
end

if(isempty(strfind(str1,'?')))&&(isempty(strfind(str2,'?')))
    ismatch  = strcmpi(str1,str2);
else
    str2(strfind(str1,'?')) = '?';
    str1(strfind(str2,'?')) = '?';
    ismatch  = strcmpi(str2,str1);
end

end

function [numSubstr,substrSpecies] = ProceedByRule_withbracket(nreResidue,enzObj,glycanObj,...
    i,substrSpecies,numSubstr,bracketresidue)
% rule 4: if the the substrate should target the terminal residue
if(~isempty(enzObj.isTerminalTarget))
    if(enzObj.isTerminalTarget)
        parent = nreResidue.getParent;
        if(isequal(parent.residueType.name,'#bracket'))
            isValidRxnPos = true;
        else
            if(length(nreResidue.getParent.linkageChildren)==1);
                isValidRxnPos = true;
            else
                isValidRxnPos = false;
            end
        end
    else
        isValidRxnPos = true;
    end
else
    isValidRxnPos = true;
end

if(~isValidRxnPos)
    return;
end

% rule 5:resiude next to terminal residue is the same as
% acceptorResidue
requirement        = enzObj.dispAttachResLink;
nexterminalResidue = nreResidue.getParent;
if(isequal(nexterminalResidue.residueType.name,'#bracket'))
    isValidRxnPos = 1;
else
    if(~isempty(nexterminalResidue))
        acceptorTerminal     = nexterminalResidue.dispResidueInfo;
    else
        acceptorTerminal    = [];
    end
    isValidRxnPos = fuzzymatch(acceptorTerminal,requirement);
end

if(~isValidRxnPos)
    return;
end

% rule 6: the product should not be greater than the specified structure
if(~isempty(enzObj.prodMinStruct))
    isValidRxnPos = isValidRxnPos && ...
        glycanObj.bracketContain(enzObj.prodMinStruct);
    if(~isValidRxnPos)
        return;
    end
end

% rule 7: the enzyme should not work on the children of specificed structure
if(~isempty(enzObj.prodMaxStruct))
    isValidRxnPos = isValidRxnPos && ...
        enzObj.prodMaxStruct.bracketContain(glycanObj);
    if(~isValidRxnPos)
        return;
    end
end

% rule 8: the enzyme should transfer the residue to the target branch containing
% specified structure
if(~isempty(enzObj.targetbranchcontain)) % if the rule exists
    chartargetbranch = enzObj.targetbranchcontain.name;
    chartargetbranch = regexprep(chartargetbranch,'freeEnd--',' ');
    chartargetbranch = regexprep(chartargetbranch,'?','');
    chartargetbranch = regexprep(chartargetbranch,'\d','');
    
    if(isempty(enzObj.isTerminalTarget))
        enzObj.isTerminalTarget = true;
    end
    
    if(~enzObj.isTerminalTarget)
        parentresidue = nreResidue.getParent;
        if(~isequal(parentresidue.residueType.name,'#bracket'))
            charresidue   = glycanObj.getresiduetoterminal(parentresidue,nreResidue);
            charresidue   = regexprep(charresidue,'freeEnd--','');
            charresidue   = regexprep(charresidue,'?','');
            isValidRxnPos = (~isempty(strfind(strtrim(charresidue),...
                strtrim(chartargetbranch))));
        else
            for j = 1 : length(bracketresidue)
                if(strcmp(bracketresidue{1,j}.residueType.name,nreResidue.residueType.name));
                    continue
                end
                Objstruct     = glycanObj.clone;
                Objresidues   = Objstruct.getNonRedEndResidue('bracket');
                isValidRxnPos = checkavailability(Objstruct,Objresidues{1,j},enzObj);
                if(isValidRxnPos)
                    break
                end
            end
        end
    else
        Objstruct     = glycanObj.clone;
        Objresidues   = Objstruct.getNonRedEndResidue;
        isValidRxnPos = checkavailability(Objstruct,Objresidues{1,i},enzObj);
    end
end

if(~isValidRxnPos)
    return;
end


glycanSubstrObj     = glycanObj.clone;
nreResidues         = glycanSubstrObj.getNonRedEndResidue;
nreResiduetoremove  = nreResidues{1,i};

isResidueDeleted    = glycanSubstrObj.removeNonRedEndResidue(...
    nreResiduetoremove);

%rule 11: It does not contain the specified residue
if(~isempty(enzObj.substNAResidue))
    NAResiduename = enzObj.substNAResidue.name;
    composition   = glycanSubstrObj.getComposition;
    isValidRxnPos  =  ~(isfield(composition,NAResiduename) ...
        && (composition.(NAResiduename)>=1));
    
    if(~isValidRxnPos)
        return;
    end
end

% rule 12, the enzyme should not work on the parent of specified structure
if(~isempty(enzObj.substMaxStruct)) % if the rule exists
    isValidRxnPos = isValidRxnPos && ...
        enzObj.substMaxStruct.bracketContain(glycanSubstrObj);
    
    if(~isValidRxnPos)
        return;
    end
end

% rule 13, the enzyme should not work on the children of specified structure
if(~isempty(enzObj.substMinStruct))  % if the rule exists
    isValidRxnPos = isValidRxnPos &&...
        glycanSubstrObj.bracketContain(enzObj.substMinStruct);
    if(~isValidRxnPos)
        return;
    end
end

numSubstr = numSubstr+1;
substr = GlycanSpecies(glycanSubstrObj);
substrSpecies.add(substr);
end

function [numSubstr,substrSpecies] = ProceedByRule(nreResidue,enzObj,glycanObj,i,substrSpecies,numSubstr)
% rule 4: if the the substrate should target the terminal residue
if(~isempty(enzObj.isTerminalTarget))
    if(enzObj.isTerminalTarget)
        if(length(nreResidue.getParent.linkageChildren)==1);
            isValidRxnPos = true;
        else
            isValidRxnPos = false;
        end
    else
        isValidRxnPos = true;
    end
else
    isValidRxnPos = true;
end

if(~isValidRxnPos)
    return;
end

% rule 5:resiude next to terminal residue is the same as
% acceptorResidue
requirement        = enzObj.dispAttachResLink;
nexterminalResidue = nreResidue.getParent;
if(~isempty(nexterminalResidue))
    acceptorTerminal     = nexterminalResidue.dispResidueInfo;
else
    acceptorTerminal    = '';
end

isValidRxnPos = fuzzymatch(acceptorTerminal,requirement);
if(~isValidRxnPos)
    return;
end

% rule 6: the product should not be greater than the specified structure
if(~isempty(enzObj.prodMinStruct))
    isValidRxnPos = isValidRxnPos && ...
        glycanObj.bracketContain(enzObj.prodMinStruct);
    if(~isValidRxnPos)
        return;
    end
end

% rule 7: the enzyme should not work on the children of specificed structure
if(~isempty(enzObj.prodMaxStruct))
    isValidRxnPos = isValidRxnPos && ...
        enzObj.prodMaxStruct.bracketContain(glycanObj);
    if(~isValidRxnPos)
        return;
    end
end

% rule 8: the enzyme should transfer the residue to the target branch
if(~isempty(enzObj.targetBranch))  % if the rule exists
    chartargetbranch   = enzObj.targetBranch.name;
    chartargetbranch   = regexprep(chartargetbranch,'freeEnd--','');
    chartargetbranch   = regexprep(chartargetbranch,'redEnd--','');
    charresidue        = glycanObj.getresiduetoroot(nreResidue.getParent);
    charresidue        = regexprep(charresidue,'freeEnd--','');
    charresidue        = regexprep(charresidue,'redEnd--','');
    
    isValidRxnPos      = fuzzymatch(chartargetbranch,charresidue);
    if(~isValidRxnPos)
        return;
    end
end
% rule 9: the enzyme should transfer the residue to the target branch containing
% specified structure
if(~isempty(enzObj.targetbranchcontain)) % if the rule exists
    chartargetbranch = enzObj.targetbranchcontain.name;
    chartargetbranch = regexprep(chartargetbranch,'freeEnd--',' ');
    chartargetbranch = regexprep(chartargetbranch,'?','');
    
    if(isempty(enzObj.isTerminalTarget))
        enzObj.isTerminalTarget = true;
    end
    
    if(~enzObj.isTerminalTarget)
        parentresidue = nreResidue.getParent;
        charresidue   = glycanObj.getresiduetoterminal(parentresidue,nreResidue);
        charresidue   = regexprep(charresidue,'freeEnd--','');
        charresidue   = regexprep(charresidue,'?','');
        
        isValidRxnPos = isValidRxnPos && (~isempty(strfind(strtrim(charresidue),...
            strtrim(chartargetbranch))));
    else
        parentresidue = nreResidue.getParent;
        charresidue   = glycanObj.getresiduetoroot(parentresidue);
        charresidue   = regexprep(charresidue,'freeEnd--','');
        charresidue   = regexprep(charresidue,'?','');
        
        isValidRxnPos = isValidRxnPos && (~isempty(strfind(strtrim(charresidue),...
            strtrim(chartargetbranch))));
    end
    
    if(~isValidRxnPos)
        return;
    end
end

% rule 10: the enzyme should not act on the specificed branch
if(~isempty(enzObj.targetNABranch))  % if the rule exists
    if(isa(enzObj.targetNABranch,'GlycanStruct'))
        chartargetbranch = enzObj.targetNABranch.name;
        charresidue           = glycanObj.getresiduetoroot(...
            nreResidue.getParent);
        isValidRxnPos         = isValidRxnPos && ...
            (~strcmp(chartargetbranch,charresidue));
    elseif(isa(enzObj.targetNABranch,'CellArrayList'))
        for ii = 1 : length(enzObj.targetNABranch)
            chartargetbranch  = enzObj.targetNABranch.name;
            charresidue       = glycanObj.getresiduetoroot(...
                nreResidue.getParent);
            isValidRxnPos     = isValidRxnPos && ...
                (~strcmp(chartargetbranch,charresidue));
        end
    end
    
    if(~isValidRxnPos)
        return;
    end
end

glycanSubstrObj     = glycanObj.clone;
nreResidues         = glycanSubstrObj.getNonRedEndResidue;
nreResiduetoremove  = nreResidues{1,i};

isResidueDeleted    = glycanSubstrObj.removeNonRedEndResidue(...
    nreResiduetoremove);

%rule 11: It does not contain the specified residue
if(~isempty(enzObj.substNAResidue))
    NAResiduename = enzObj.substNAResidue.name;
    composition   = glycanSubstrObj.getComposition;
    isValidRxnPos  =  ~(isfield(composition,NAResiduename) ...
        && (composition.(NAResiduename)>=1));
    
    if(~isValidRxnPos)
        return;
    end
end


% rule 12: the substrate should not contain the specificed branch
if(~isempty(enzObj.substNABranch))  % if the rule exists
    if(isa(enzObj.substNABranch,'GlycanStruct'))
        isValidRxnPos = isValidRxnPos &&...
            (~glycanSubstrObj.contains(enzObj.substNABranch));
    elseif(isa(enzObj.substNABranch,'CellArrayList'))
        for ii = 1 : length(enzObj.substNABranch)
            isValidRxnPos = isValidRxnPos && ...
                (~glycanSubstrObj.contains(enzObj.substNABranch.get(ii)));
        end
    end
    
    if(~isValidRxnPos)
        return;
    end
end

% rule 13, the enzyme should not work on the parent of specified structure
if(~isempty(enzObj.substMaxStruct)) % if the rule exists
    isValidRxnPos = isValidRxnPos && ...
        enzObj.substMaxStruct.bracketContain(glycanSubstrObj);
    
    if(~isValidRxnPos)
        return;
    end
end

% rule 14, the enzyme should not work on the children of specified structure
if(~isempty(enzObj.substMinStruct))  % if the rule exists
    isValidRxnPos = isValidRxnPos &&...
        glycanSubstrObj.bracketContain(enzObj.substMinStruct);
    if(~isValidRxnPos)
        return;
    end
end

numSubstr = numSubstr+1;
substr = GlycanSpecies(glycanSubstrObj);
substrSpecies.add(substr);
end

function substrSpecies = isaddbracket(substrSpecies,numSubstr,enzObj)
if(numSubstr>1)
    isbracket = 0;
    branches  = getBranch(substrSpecies.glycanStruct);
    isbranchdepthvalid = 1;
    for j = 1 : length(branches)
        if(branches(1,j).depth<6)
           isbranchdepthvalid = 0;
           break
        end
    end   
    if(numSubstr==length(branches))&&(isbranchdepthvalid)
        isbracket = 1;
    end
    
    if(isbracket)
        substr        = substrSpecies.get(1);
        substrSpecies = CellArrayList;
        substrSpecies.add(substr);
        substrnreResidues = substrSpecies.get(1).glycanStruct.getNonRedEndResidue;
        counter = 0;
        for i = 1 : length(substrnreResidues)
            ithnreresidue = substrnreResidues{1,i};
            if(isequal(ithnreresidue.residueType.name,'Fuc'))&&...
                    (isequal(ithnreresidue.linkageParent.bonds.posParent,'6'))
                continue
            end
            
            if(isequal(ithnreresidue.residueType.name,enzObj.resfuncgroup.name))
                residuestobracket(counter+1,1) = ithnreresidue;
                counter = counter+1;
            end
        end
        substrSpecies.get(1).glycanStruct.bracketResidue(residuestobracket);
    end
end
end

function glycanstruct = merge2exact(glycanstruct,addReidue)
objnreResidue = glycanstruct.getNonRedEndResidue('root');
glycanstruct.bracket = '';
for i = 1 : length(objnreResidue)
    if(~strcmp(objnreResidue{1,i}.residueType.name,'Fuc'))
        isResidueAdded = glycanstruct.addResidue(objnreResidue{1,i},...
            addReidue);
    end 
end
glycanstruct.resetjava;
end

function addReidue = findFirstReisude(GlycanResidue)
parentReisdue = GlycanResidue.getParent;
if(strcmp(parentReisdue.residueType.name,'#bracket'))
    addReidue = GlycanResidue;
    return
end
addReidue = findFirstReisude(parentReisdue);

end

function isValidRxnPos = checkavailability(Objstruct,Objresidue,enzObj)
addReidue   = findFirstReisude(Objresidue);
Objstruct   = merge2exact(Objstruct,addReidue);

allnreResidues = Objstruct.getNonRedEndResidue;
isValidRxnPos  = 1;
for j = 1 : length(allnreResidues)
    if(~strcmp(allnreResidues{1,j}.residueType.name,...
            Objresidue.residueType.name))
        continue
    end
    chartargetbranch = regexprep(enzObj.targetbranchcontain.name,'freeEnd','');
    chartargetbranch = regexprep(chartargetbranch,'\d','');
    chartargetbranch = regexprep(chartargetbranch,'?','');
    charresidue      = Objstruct.getresiduetoroot(allnreResidues{1,j});
    charresidue      = regexprep(charresidue,'freeEnd','');
    charresidue      = regexprep(charresidue,'\d','');
    charresidue      = regexprep(charresidue,'?','');
    isValid          = (~isempty(strfind(charresidue,chartargetbranch)));
    if(~isValid)
        isValidRxnPos = 0;
        break
    end
end
end





