residueMap=load('residueTypes.mat');

% ST3GalNAc I
ST3GalNAc                     = GTEnz([2;4;99;4]);
ST3GalNAc.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalNAc.resAtt2FG           = residueMap.allresidues('Gal');
galBond                       = GlycanBond('?','1'); 
ST3GalNAc. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                       = GlycanBond('3','2');
ST3GalNAc.linkFG              = struct('anomer','a','bond',siabond);
ST3GalNAc.isTerminalTarget    = true;

a2  = GlycanSpecies(glycanMLread('1124.56.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a2,ST3GalNAc);