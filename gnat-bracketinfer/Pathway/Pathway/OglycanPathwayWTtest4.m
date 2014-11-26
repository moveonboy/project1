residueMap=load('residueTypes.mat');

% ST6GalNAc I
ST6GalNAc                    = GTEnz([2;4;99;3]);
ST6GalNAc.resfuncgroup       = residueMap.allresidues('NeuAc');
ST6GalNAc.resAtt2FG          = residueMap.allresidues('GalNAc');
glcnacBond                   = GlycanBond('?','1');
ST6GalNAc.linkresAtt2FG      = struct('bond',glcnacBond,'anomer','a');
siabond                      = GlycanBond('6','2');
ST6GalNAc.linkFG             = struct('anomer','a','bond',siabond);
ST6GalNAc.substMinStruct     = glycanMLread('314.16.glycoct_xml');
ST6GalNAc.substMaxStruct     = glycanMLread('879.43.glycoct_xml');

a1  = GlycanSpecies(glycanMLread('1240.60.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a1,ST6GalNAc);