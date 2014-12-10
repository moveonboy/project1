function composition = struct2strcomp(species)
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