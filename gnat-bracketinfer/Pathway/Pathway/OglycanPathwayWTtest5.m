residueMap=load('residueTypes.mat');

ST3GalIII                     = GTEnz([2;4;99;6]);
ST3GalIII.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalIII.resAtt2FG           = residueMap.allresidues('Gal');
galBond                       = GlycanBond('?','1');  
ST3GalIII. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                       = GlycanBond('3','2');
ST3GalIII.linkFG              = struct('anomer','a','bond',siabond);
ST3GalIII.targetBranch        = glycanMLread('763.38tb.glycoct_xml');
ST3GalIII.isTerminalTarget    = true;

a1  = GlycanSpecies(glycanMLread('1863.92.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a1,ST3GalIII);