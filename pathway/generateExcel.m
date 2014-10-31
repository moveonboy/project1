function generateExcel(varargin)
if(nargin>1)
    error('MATLAB:GNAT:ERRORNONCOMPLEX','WRONG NUMBER OF INPUT');
end

if(isempty(varargin))
    purpose = 'HL60WTN_glycandata';
end

if(varargin{1}==1)
    purpose = 'HL60WTN_glycandata';
elseif(varargin{1}==2)
    purpose = 'HL60WT_glycanspeceis';
end

glycanpath       = load('HL60WTN_glycanPathway1_25Red.mat');
HL60WTGlycanPath = glycanpath.nlinkedpath;
listofSpecies    = HL60WTGlycanPath.theSpecies;

% delete the Species that has terminal GlcNAc
newlistofSpecies = checkforTerGlc(listofSpecies);

% generate the composition
composition = strcomp(newlistofSpecies);

% generate monoisotopic mass
Monomass = strmass(newlistofSpecies);

% generate enzyme information
Enzinfo  = createEnzdata(newlistofSpecies);

% lump the information above
GlycanSpecies = struct('composition',[],'mass',[],'enzinfo',[]);
for i = 1 :length(composition)
    GlycanSpecies(i,1).composition = composition{i,1};
    GlycanSpecies(i,1).mass        = Monomass(i,1);
    GlycanSpecies(i,1).enzinfo     = Enzinfo(i,1);
end

% combine the glycanspecies with the same composition
newGlycanSpecies = combineglySpecies(GlycanSpecies,purpose);

% output the Excel file
outputExcel(newGlycanSpecies,purpose)

% buildlibrary
HL60grouplibrary = containers.Map;
for i = 1 : length(newGlycanSpecies)
    HL60grouplibrary(num2str(i)) = newGlycanSpecies(i).composition;
end

end

function composition = strcomp(listofSpecies)
composition = '';
counter     = 0;
for i = 1 : length(listofSpecies)
    ithglycanspecies  = listofSpecies.get(i);
    ithglySpeciesname = ithglycanspecies.glycanStruct.name;
    ithcompostion     = '';
    
    NeuAcposition = strfind(ithglySpeciesname,'NeuAc');
    if(~isempty(NeuAcposition))
        numofNeuAc    = length(NeuAcposition);
        ithcompostion = ['NeuAc' num2str(numofNeuAc)];
    end
    
    Fucposition = strfind(ithglySpeciesname,'Fuc');
    if(~isempty(Fucposition))
        numofFuc      = length(Fucposition);
        ithcompostion = [ithcompostion 'Fuc' num2str(numofFuc)];
    end
    
    Manposition = strfind(ithglySpeciesname,'Man');
    Galposition = strfind(ithglySpeciesname,'Gal');
    if(~isempty(Manposition))||(~isempty(Galposition))
        numofHex      = length(Manposition)+length(Galposition);
        ithcompostion = [ithcompostion 'Hex' num2str(numofHex)];
    end
    
    GlcNAcposition = strfind(ithglySpeciesname,'GlcNAc');
    if(~isempty(GlcNAcposition))
        numofHexNAc   = length(GlcNAcposition);
        ithcompostion = [ithcompostion 'HexNAc' num2str(numofHexNAc)];
    end
    
    composition{counter+1,1} = ithcompostion;
    counter                = counter+1;
end
end

function Monomass = strmass(listofSpecies)
Monomass = [];
counter  = 0;
for i = 1 : length(listofSpecies)
    ithglycanspecies  = listofSpecies.get(i);
    ithglySpeciesname = ithglycanspecies.glycanStruct.name;
    ithcompostion     = '';
    
    NeuAcposition = strfind(ithglySpeciesname,'NeuAc');
    if(~isempty(NeuAcposition))
        numofNeuAc = length(NeuAcposition);
        for j = 1 : numofNeuAc
            ithcompostion = [ithcompostion 's'];
        end
    end
    
    Fucposition = strfind(ithglySpeciesname,'Fuc');
    if(~isempty(Fucposition))
        numofFuc = length(Fucposition);
        for j = 1 : numofFuc
            ithcompostion = [ithcompostion 'f'];
        end
    end
    
    Manposition = strfind(ithglySpeciesname,'Man');
    Galposition = strfind(ithglySpeciesname,'Gal');
    if(~isempty(Manposition))||(~isempty(Galposition))
        numofHex = length(Manposition)+length(Galposition);
        for j = 1 : numofHex
            ithcompostion = [ithcompostion 'h'];
        end
    end
    
    GlcNAcposition = strfind(ithglySpeciesname,'GlcNAc');
    if(~isempty(GlcNAcposition))
        numofHexNAc = length(GlcNAcposition);
        for j = 1 : numofHexNAc
            ithcompostion = [ithcompostion 'n'];
        end
    end
    ithmonomass = glycanMolWt(ithcompostion);
    Monomass(counter+1,1) = ithmonomass;
    counter               = counter+1;
end
end

function Enzinfo  = createEnzdata(listofSpecies)
for i = 1 : length(listofSpecies)
    Enzinfo(i,1) = struct('mgat1',0,'mgat2',0,'mgat3',0,'mgat4',0,'mgat5',0,'b4GalI',0,...
        'iGnt',0,'STGalT',0,'Fut8',0,'Futa23',0);
    ithglycanSpecies = listofSpecies.get(i);
    allResidues      = getAllResidues(ithglycanSpecies.glycanStruct);
    for j = 1 : length(allResidues)
        jthresidue = allResidues{j};
        if(isequal(jthresidue.residueType.name,'freeEnd'))
            continue;
        end
        
        if(isequal(jthresidue.residueType.name,'GlcNAc'))
            parentresidue = jthresidue.linkageParent.parent;
            if(~isequal(parentresidue.residueType.name,'Man'))&&(~isequal(parentresidue.residueType.name,'Gal'))
                continue;
            end
            
            if(isequal(parentresidue.residueType.name,'Man'))
                linkage2parent = jthresidue.linkageParent.bonds.posParent;
                parentlinkage  = parentresidue.linkageParent.bonds.posParent;
                if(isequal(linkage2parent,'2'))&&(isequal(parentlinkage,'3'))
                    Enzinfo(i,1).mgat1 = 1;
                end
                
                if(isequal(linkage2parent,'2'))&&(isequal(parentlinkage,'6'))
                    Enzinfo(i,1).mgat2 = 1;
                end
                
                if(isequal(linkage2parent,'4'))&&(isequal(parentlinkage,'4'))
                    Enzinfo(i,1).mgat3 = 1;
                end
                
                if(isequal(linkage2parent,'4'))&&(isequal(parentlinkage,'3'))
                    Enzinfo(i,1).mgat4 = 1;
                end
                
                if(isequal(linkage2parent,'6'))&&(isequal(parentlinkage,'6'))
                    Enzinfo(i,1).mgat5 = 1;
                end
            elseif(isequal(parentresidue.residueType.name,'Gal'))
                Enzinfo(i,1).iGnt = 1;
            end
        end
        
        if(isequal(jthresidue.residueType.name,'Gal'))
            Enzinfo(i,1).b4GalI = 1;
        end
        
        if(isequal(jthresidue.residueType.name,'Fuc'))
            likage2parent = jthresidue.linkageParent.bonds.posParent;
            if(isequal(likage2parent,'6'))
                Enzinfo(i,1).Fut8 = 1;
            elseif(isequal(likage2parent,'3'))
                Enzinfo(i,1).Futa23 = 1;
            end
        end
        
        if(isequal(jthresidue.residueType.name,'NeuAc'))
            Enzinfo(i,1).STGalT = 1;
        end
    end
end
end

function newGlycanSpecies = combineglySpecies(GlycanSpecies,purpose)
newGlycanSpecies = struct('composition',[],'mass',[],'enzinfo',[]);
if(isequal(purpose,'HL60WTN_glycandata'))
    counter    = 0;
    for i = 1 : length(GlycanSpecies)
        ithglycanspeices = GlycanSpecies(i);
        isnewgroup = 1;
        for j = 1 : length(newGlycanSpecies)
            if(isequal(newGlycanSpecies(j,1).composition,ithglycanspeices.composition)) &&...
                  (isequal(newGlycanSpecies(j,1).enzinfo(1).mgat3,ithglycanspeices.enzinfo.mgat3))
                newGlycanSpecies(j,1).enzinfo = [newGlycanSpecies(j,1).enzinfo;...
                    ithglycanspeices.enzinfo];
                isnewgroup = 0;
                break
            end
        end
        if(isnewgroup)
            newGlycanSpecies(counter+1,1).composition = ithglycanspeices.composition;
            newGlycanSpecies(counter+1,1).mass        = ithglycanspeices.mass;
            newGlycanSpecies(counter+1,1).enzinfo     = ithglycanspeices.enzinfo;
            counter = counter+1;
        end
    end
elseif(isequal(purpose,'HL60WT_glycanspeceis'))
    counter    = 0;
    for i = 1 : length(GlycanSpecies)
        ithglycanspeices = GlycanSpecies(i);
        isnewgroup = 1;
        for j = 1 : length(newGlycanSpecies)
            if(isequal(newGlycanSpecies(j,1).composition,ithglycanspeices.composition))
                newGlycanSpecies(j,1).enzinfo = [newGlycanSpecies(j,1).enzinfo;...
                    ithglycanspeices.enzinfo];
                isnewgroup = 0;
                break
            end
        end
        
        if(isnewgroup)
            newGlycanSpecies(counter+1,1).composition = ithglycanspeices.composition;
            newGlycanSpecies(counter+1,1).mass        = ithglycanspeices.mass;
            newGlycanSpecies(counter+1,1).enzinfo     = ithglycanspeices.enzinfo;
            counter = counter+1;
        end
    end
end
end

function outputExcel(newGlycanSpecies,purpose)
if(isequal(purpose,'HL60WTN_glycandata'))
    composition = '';
    monomass    = [];
    mgat1       = [];
    mgat2       = [];
    mgat3       = [];
    mgat4       = [];
    mgat5       = [];
    b4GalI      = [];
    iGnt        = [];
    STGalT      = [];
    Fut8        = [];
    Futa23      = [];
    counter     = 0;
    for i = 1 : length(newGlycanSpecies)
        ithglycanspeceis         = newGlycanSpecies(i);
        composition{counter+1,1} = ithglycanspeceis.composition;
        monomass(counter+1,1)    = ithglycanspeceis.mass;
        totalmgat1   = 0;
        totalmgat2   = 0;
        totalmgat3   = 0;
        totalmgat4   = 0;
        totalmgat5   = 0;
        totalb4GalI  = 0;
        totaliGnt    = 0;
        totalSTGalT  = 0;
        totalFut8    = 0;
        totalFuta23  = 0;
        for j = 1 : length(ithglycanspeceis.enzinfo)
            totalmgat1   = totalmgat1 + ithglycanspeceis.enzinfo(j).mgat1;
            totalmgat2   = totalmgat2 + ithglycanspeceis.enzinfo(j).mgat2;
            totalmgat3   = totalmgat3 + ithglycanspeceis.enzinfo(j).mgat3;
            totalmgat4   = totalmgat4 + ithglycanspeceis.enzinfo(j).mgat4;
            totalmgat5   = totalmgat5 + ithglycanspeceis.enzinfo(j).mgat5;
            totalb4GalI  = totalb4GalI + ithglycanspeceis.enzinfo(j).b4GalI;
            totaliGnt    = totaliGnt + ithglycanspeceis.enzinfo(j).iGnt;
            totalSTGalT  = totalSTGalT + ithglycanspeceis.enzinfo(j).STGalT;
            totalFut8    = totalFut8 + ithglycanspeceis.enzinfo(j).Fut8;
            totalFuta23  = totalFuta23 + ithglycanspeceis.enzinfo(j).Futa23;
        end
        mgat1(counter+1,1)   = totalmgat1/length(ithglycanspeceis.enzinfo);
        mgat2(counter+1,1)   = totalmgat2/length(ithglycanspeceis.enzinfo);
        mgat3(counter+1,1)   = totalmgat3/length(ithglycanspeceis.enzinfo);
        mgat4(counter+1,1)   = totalmgat4/length(ithglycanspeceis.enzinfo);
        mgat5(counter+1,1)   = totalmgat5/length(ithglycanspeceis.enzinfo);
        b4GalI(counter+1,1)  = totalb4GalI/length(ithglycanspeceis.enzinfo);
        iGnt(counter+1,1)    = totaliGnt/length(ithglycanspeceis.enzinfo);
        STGalT(counter+1,1)  = totalSTGalT/length(ithglycanspeceis.enzinfo);
        Fut8(counter+1,1)    = totalFut8/length(ithglycanspeceis.enzinfo);
        Futa23(counter+1,1)  = totalFuta23/length(ithglycanspeceis.enzinfo);
        counter = counter+1;
    end
    
    outputfilename = 'HL60WTN_glycandata';
    A1 = cellstr('Composition');
    B1 = cellstr('Monoisotopic mass');
    C1 = cellstr('Mgat1');
    D1 = cellstr('Mgat2');
    E1 = cellstr('Mgat3');
    F1 = cellstr('Mgat4');
    G1 = cellstr('Mgat5');
    H1 = cellstr('b4GalI');
    I1 = cellstr('iGnt');
    J1 = cellstr('STGalT');
    K1 = cellstr('Fut8');
    L1 = cellstr('Futa2,3');
    
    xlswrite(outputfilename,A1,1,'A1');
    xlswrite(outputfilename,composition,1,'A2');
    xlswrite(outputfilename,B1,1,'B1');
    xlswrite(outputfilename,monomass,1,'B2');
    xlswrite(outputfilename,C1,1,'C1');
    xlswrite(outputfilename,mgat1,1,'C2');
    xlswrite(outputfilename,D1,1,'D1');
    xlswrite(outputfilename,mgat2,1,'D2');
    xlswrite(outputfilename,E1,1,'E1');
    xlswrite(outputfilename,mgat3,1,'E2');
    xlswrite(outputfilename,F1,1,'F1');
    xlswrite(outputfilename,mgat4,1,'F2');
    xlswrite(outputfilename,G1,1,'G1');
    xlswrite(outputfilename,mgat5,1,'G2');
    xlswrite(outputfilename,H1,1,'H1');
    xlswrite(outputfilename,b4GalI,1,'H2');
    xlswrite(outputfilename,I1,1,'I1');
    xlswrite(outputfilename,iGnt,1,'I2');
    xlswrite(outputfilename,J1,1,'J1');
    xlswrite(outputfilename,STGalT,1,'J2');
    xlswrite(outputfilename,K1,1,'K1');
    xlswrite(outputfilename,Fut8,1,'K2');
    xlswrite(outputfilename,L1,1,'L1');
    xlswrite(outputfilename,Futa23,1,'L2');
elseif(isequal(purpose,'HL60WT_glycanspeceis'))
    composition = '';
    monomass    = [];
    counter     = 0;
    for i = 1 : length(newGlycanSpecies)
        ithglycanspeceis         = newGlycanSpecies(i);
        composition{counter+1,1} = ithglycanspeceis.composition;
        monomass(counter+1,1)    = ithglycanspeceis.mass;
        counter = counter+1;
    end
    
    outputfilename = 'HL60WT_glycanspeceis';
    A1 = cellstr('Composition');
    B1 = cellstr('Monoisotopic mass');
    
    xlswrite(outputfilename,A1,1,'A1');
    xlswrite(outputfilename,composition,1,'A2');
    xlswrite(outputfilename,B1,1,'B1');
    xlswrite(outputfilename,monomass,1,'B2');
end
end

function newlistofSpecies = checkforTerGlc(listofSpecies)
newlistofSpecies = CellArrayList;
for i = 1 : length(listofSpecies)
    ithspecies        = listofSpecies.get(i);
    ithnonredresidues = getNonRedEndResidue(ithspecies.glycanStruct);
    ithresidues       = getAllResidues(ithspecies.glycanStruct);
    isiGnt      = 0;
    for j = 1 : length(ithresidues)
        jthresidue = ithresidues{j};
        if(isequal(jthresidue.residueType.name,'GlcNAc'))&&...
                (isequal(jthresidue.linkageParent.bonds.posParent,'3'))
            isiGnt = 1;
            break
        elseif(isequal(jthresidue.residueType.name,'GlcNAc'))&&...
                (isequal(jthresidue.getParent.residueType.name,'#bracket'))
            isiGnt = 1;
            break
        end
    end
    
    isdelete    = 0;
    for j = 1 : length(ithnonredresidues)
        jthnonredresidue  = ithnonredresidues{j};
        if(isequal(jthnonredresidue.residueType.name,'GlcNAc'))&&(isiGnt)
            %             if(~isempty(ithnonredresidues.getParent))
            jthnonredrePARENT = jthnonredresidue.getParent;
            if(~isequal(jthnonredrePARENT.residueType.name,'Man'))
                isdelete = 1;
                break
            else
                if(~isequal(jthnonredrePARENT.linkageParent.bonds.posParent,'4'))
                    isdelete = 1;
                    break
                end
            end
        end
    end
    
    if(~isdelete)
        newlistofSpecies.add(ithspecies)
    end
end
end
