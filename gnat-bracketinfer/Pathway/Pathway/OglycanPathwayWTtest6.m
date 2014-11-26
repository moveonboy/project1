residueMap=load('residueTypes.mat');

ST3GalIII                     = GTEnz([2;4;99;6]);
ST3GalIII.isTerminalTarget    = true;
ST3GalIII.resfuncgroup        = residueMap.allresidues('NeuAc');
siaTbond                      = GlycanBond('3','2');
ST3GalIII.linkFG              = struct('anomer','a','bond',siaTbond);
galResType                    = residueMap.allresidues('Gal');
galBond                       = GlycanBond('4','1');
ST3GalIII.resAtt2FG           = galResType;
ST3GalIII.linkresAtt2FG       = struct('bond', galBond,'anomer','b');
ST3GalIII.targetBranch        = glycanMLread('lactosetagetBranch.glycoct_xml');

a26  = GlycanSpecies(glycanMLread('6820tetra_5lac.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a26,ST3GalIII);