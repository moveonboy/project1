function [glycaninMS,expecGlycan] = createLocalDataBase(pathwayfilename,varargin)
%createLocalDataBase load the glycanpathway to build a local glycan
% database, including all the potential glycans.
%
% glycaninMS = createLocalDataBase(pathwayfilename,outputfilename) load the glycanpathway
% named by pathwayfilename to create the local database and write the composition and monoisotopic 
% mass into Excel file.. 
%
% Example:
%     glycaninMS=createLocalDataBase('MSglycandatabase.mat','HL60glycanList.xls')
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 11/06/14

if(isempty(varargin))
    outputfilename = '';
elseif(length(varargin)==1)
    outputfilename=varargin{1};
end
% load(MSfilename);
glycanpath       = load(pathwayfilename);
HL60WTGlycanPath = glycanpath.nlinkedpath;
listofSpecies    = HL60WTGlycanPath.theSpecies;
% peaklist         = MSdata.peaklist;
% peakwidth        = MSdata.pfwhh;
% isGlycaninMS     = 0;
glycaninMS       = CellArrayList;
expecGlycan      = cell(length(listofSpecies),1);
monoisomw        = cell(length(listofSpecies),1);
for i = 1 : length(listofSpecies)
    ithspecies           = listofSpecies.get(i);
    ithcompostion        = struct2strcomp(ithspecies);
    ithspeciesmonomass   = glycanMolWt(ithcompostion);
    expecGlycan{count+1} = ithcompostion;
    monoisomw{count+1}   = ithspeciesmonomass;
    glycaninMS.add(ithspecies)
    count = count+1;
%     for j = 1 : length(peaklist)
%         jthpeak = peaklist(j,1);
%         pwfhh   = peakwidth(j,2)-peakwidth(j,1);
%         if((abs(jthpeak-ithspeciesmonomass))<0.5*(pwfhh))
%             isGlycaninMS = 1;
%             break
%         end
%     end
%     if(isGlycaninMS)
%         glycaninMS.add(ithspecies)
%     end
end
save('MSglycandatabase','glycaninMS')

if(~isempty(outputfilename))
    A1=cellstr('Composition');
    B1=cellstr('Monoisotopic mass');
    xlswrite(outputfilename,A1,1,'A1');
    xlswrite(outputfilename,B1,1,'B1');
    xlswrite(outputfilename,expecGlycan,1,'A2');
    xlswrite(outputfilename, monoisomw,1,'B2');
end
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