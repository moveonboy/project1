function glycaninMS = createLocalDataBase(pathwayfilename,MSfilename)
%createLocalDataBase load the glycanpathway and MS rawdata to build a local glycan
% database, glycans in which have corresponding peaks in MS rawdata.
%
% glycaninMS = createLocalDataBase(pathwayfilename,MSfilename) load the glycanpathway
% named by pathwayfilename and MS rawdata named by MSfilename to create the local database. 
%
% Example:
%     glycaninMS=createLocalDataBase('MSglycandatabase.mat','MSdata.mat')
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 9/22/14
load(MSfilename);
glycanpath       = load(pathwayfilename);
HL60WTGlycanPath = glycanpath.nlinkedpath;
listofSpecies    = HL60WTGlycanPath.theSpecies;
peaklist         = MSdata.peaklist;
peakwidth        = MSdata.pfwhh;
isGlycaninMS     = 0;
glycaninMS       = CellArrayList;
for i = 1 : length(listofSpecies)
    ithspecies         = listofSpecies.get(i);
    ithcompostion      = strcomp(ithspecies);
    ithspeciesmonomass = glycanMolWt(ithcompostion);
    for j = 1 : length(peaklist)
        jthpeak = peaklist(j,1);
        pwfhh   = peakwidth(j,2)-peakwidth(j,1);
        if((abs(jthpeak-ithspeciesmonomass))<0.5*(pwfhh))
            isGlycaninMS = 1;
            break
        end
    end
    if(isGlycaninMS)
        glycaninMS.add(ithspecies)
    end
end
save('MSglycandatabase','glycaninMS')
end

function composition = strcomp(species)
composition    = '';
glySpeciesname = species.glycanStruct.name;

NeuAcposition  = strfind(glySpeciesname,'NeuAc');
if(~isempty(NeuAcposition))
    numofNeuAc = length(NeuAcposition);
    for j = 1 : numofNeuAc
        composition = [composition 's'];
    end
end

NeuGcposition  = strfind(glySpeciesname,'NeuGc');
if(~isempty(NeuGcposition))
    numofNeuGc = length(NeuGcposition);
    for j = 1 : numofNeuGc
        composition = [composition 'g'];
    end
end

Fucposition = strfind(glySpeciesname,'Fuc');
if(~isempty(Fucposition))
    numofFuc = length(Fucposition);
    for j = 1 : numofFuc
        composition = [composition 'f'];
    end
end

Manposition = strfind(glySpeciesname,'Man');
Galposition = strfind(glySpeciesname,'Gal');
if(~isempty(Manposition))||(~isempty(Galposition))
    numofHex = length(Manposition)+length(Galposition);
    for j = 1 : numofHex
        composition = [composition 'h'];
    end
end

GlcNAcposition = strfind(glySpeciesname,'GlcNAc');
if(~isempty(GlcNAcposition))
    numofHexNAc = length(GlcNAcposition);
    for j = 1 : numofHexNAc
        composition = [composition 'n'];
    end
end

Xylposition = strfind(glySpeciesname,'Xyl');
if(~isempty(Xylposition))
    numofXyl = length(Xylposition);
    for j = 1 : numofHexNAc
        composition = [composition 'x'];
    end
end

Kdnposition = strfind(glySpeciesname,'Kdn');
if(~isempty(Kdnposition))
    numofKdn = length(Kdnposition);
    for j = 1 : numofKdn
        composition = [composition 'k'];
    end
end

ManAposition = strfind(glySpeciesname,'ManA');
GalAposition = strfind(glySpeciesname,'GalA');
if(~isempty(ManAposition))||(~isempty(GalAposition))
    numofHexA = length(ManAposition)+length(GalAposition);
    for j = 1 : numofHexA
        composition = [composition 'u'];
    end
end
end