residueMap=load('residueTypes.mat');

% GCNT1
GCNTI                    = GTEnz([2;4;1;102]);
GCNTI.resfuncgroup       = residueMap.allresidues('GlcNAc');
GCNTI.resAtt2FG          = residueMap.allresidues('GalNAc');
galnacBond               = GlycanBond('?','1');
GCNTI.linkresAtt2FG      = struct('bond',galnacBond,'anomer','a');
glcnacbond                  = GlycanBond('6','1');
GCNTI.linkFG             = struct('anomer','b','bond',glcnacbond);
GCNTI.substMaxStruct     = glycanMLread('518.26.glycoct_xml');

a1  = GlycanSpecies(glycanMLread('763.38.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a1,GCNTI);