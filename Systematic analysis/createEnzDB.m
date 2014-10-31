function enzDB = createEnzDB()
%createEnzDB create the enzyme information database by using the
% containers.Map
%
% enzDB = createEnzDB() using the containers.Map to build the database, 
% the key is the name of the enzyme and the mapObj is the enzyme
% information
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 9/23/14

residueMap=load('residueTypes.mat');

%mgat1
mgat1                    = GTEnz([2;4;1;101]);
mgat1.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('3','1');
mgat1.resAtt2FG          = manResType;
mgat1. linkresAtt2FG     = struct('bond',manBond,'anomer','a');
glcnacbond               = GlycanBond('2','1');
mgat1.linkFG             = struct('anomer','b','bond',glcnacbond);
a2                       = glycanMLread('1579.78.glycoct_xml');
mgat1.substMinStruct     = a2;
mgat1.substMaxStruct     = a2;
mgat1.targetBranch       = glycanMLread('mgat1targetbranch.glycoct_xml');

%mgat2
mgat2                    = GTEnz([2;4;1;143]);
mgat2.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('6','1');
mgat2.resAtt2FG          = manResType;
mgat2. linkresAtt2FG     = struct('bond',manBond,'anomer','a');
glcnacbond               = GlycanBond('2','1');
mgat2.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('gnm3gn.glycoct_xml');
mgat2.substMinStruct     = gn2m3gn;
mgat2.substNABranch      = glycanMLread('mgat2NAbranch2.glycoct_xml');
mgat2.substNAResidue     = residueMap.allresidues('Gal');

%mgat3
mgat3                    = GTEnz([2;4;1;144]);
mgat3.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('4','1');
mgat3.resAtt2FG          = manResType;
mgat3.linkresAtt2FG      = struct('bond', manBond,'anomer','b');
glcnacbond               = GlycanBond('4','1');
mgat3.linkFG             = struct('anomer','b','bond',glcnacbond);
m3gn                     = glycanMLread('m3gn.glycoct_xml');
mgat3.substMinStruct     = m3gn;
mgat3.targetBranch       = glycanMLread('mgat3targetbranch.glycoct_xml');
mgat3.substNAResidue     = residueMap.allresidues('Gal');

%mgat4
mgat4                    = GTEnz([2;4;1;145]);
mgat4.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('3','1');
mgat4.resAtt2FG          = manResType;
mgat4.linkresAtt2FG      = struct('bond', manBond,'anomer','a');
glcnacbond               = GlycanBond('4','1');
mgat4.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('1416.71.glycoct_xml');
mgat4.substMinStruct     = gn2m3gn;
mgat4.substNABranch      = glycanMLread('mgat4NAbranch1.glycoct_xml');
mgat4.substNAResidue     = residueMap.allresidues('Gal');

%mgat5
mgat5                    = GTEnz([2;4;1;155]);
mgat5.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('6','1');
mgat5.resAtt2FG          = manResType;
mgat5.linkresAtt2FG      = struct('bond', manBond,'anomer','a');
glcnacbond               = GlycanBond('6','1');
mgat5.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('mgat5minimal.glycoct_xml');
mgat5.substMinStruct     = gn2m3gn;
mgat5.substNABranch      = glycanMLread('mgat5NAbranch2.glycoct_xml');
mgat5.substNAResidue     = residueMap.allresidues('Gal');    

%manii
manResType               = residueMap.allresidues('Man');
manii                    = GHEnz([3;2;1;114]);
manii.resfuncgroup       = manResType;
manBond(1,1)             = GlycanBond('3','1');
manBond(2,1)             = GlycanBond('6','1');
manii.linkFG             = struct('bond',manBond ,'anomer','a');
manii.resAtt2FG          = manResType;
resbond                  = GlycanBond('6','1');
manii.linkresAtt2FG      = struct('bond', resbond,'anomer','a');
manii.substNABranch      = glycanMLread('1008.51.glycoct_xml');
manii.substMaxStruct     = glycanMLread('1824.91.glycoct_xml');

%mani
mani                     = GHEnz([3;2;1;113]);
mani.resfuncgroup        = residueMap.allresidues('Man');
mani.resAtt2FG           = residueMap.allresidues('Man');
manBond                  = GlycanBond('2','1');
mani.linkFG              = struct('anomer','a','bond',manBond);
manunknownbond           = GlycanBond('?','1');
mani.linkresAtt2FG       = struct('bond', manunknownbond,'anomer','a');
mani.substMinStruct      = glycanMLread('maniminstruct.glycoct_xml');
mani.substMaxStruct      = glycanMLread('manimaxstruct.glycoct_xml');

%ST3GalIV/VI/III, ST6GalT
STGalTs                     = GTEnz([2;4;99;4]);
STGalTs.isTerminalTarget    = true;
STGalTs.resfuncgroup        = residueMap.allresidues('NeuAc');
siaTbond                     = GlycanBond('?','2');
STGalTs.linkFG              = struct('anomer','a','bond',siaTbond);
galResType                   = residueMap.allresidues('Gal');
galBond                      = GlycanBond('4','1');
STGalTs.resAtt2FG           = galResType;
STGalTs.linkresAtt2FG       = struct('bond', galBond,'anomer','b');
STGalTs.targetbranchcontain = glycanMLread('lactose.glycoct_xml');

%Fut8
Fut8                     = GTEnz([2;4;1;68]);
Fut8.resfuncgroup        = residueMap.allresidues('Fuc');
Fut8.resAtt2FG           = residueMap.allresidues('GlcNAc');
fucbond                  = GlycanBond('6','1');
Fut8.linkFG              = struct('anomer','a','bond',fucbond);
gnbond                   = GlycanBond('?','?');
Fut8.linkresAtt2FG       = struct('bond',gnbond,'anomer','b');
gnm3                     = glycanMLread('1416.71.glycoct_xml');
Fut8.substMinStruct      = gnm3;
Fut8.substNABranch       = glycanMLread('1008.51.glycoct_xml');
Fut8.substNAResidue      = residueMap.allresidues('Fuc');
Fut8.isTerminalTarget    = false;

%b4GalI
b4GalI                   = GTEnz([2;4;1;38]);
b4GalI.resfuncgroup      = residueMap.allresidues('Gal');
b4GalI.resAtt2FG         = residueMap.allresidues('GlcNAc');
fucbond                  = GlycanBond('4','1');
b4GalI.linkFG            = struct('anomer','b','bond',fucbond);
gnbond                   = GlycanBond('?','1');
b4GalI.linkresAtt2FG     = struct('bond',gnbond,'anomer','b');
gn2m3gn                  = glycanMLread('1416.71.glycoct_xml');
b4GalI.substMinStruct    = gn2m3gn;
b4GalI.targetNABranch    = glycanMLread('b4galitargetNAbranch.glycoct_xml');

%iGnt
iGnt                     = GTEnz([2;4;1;149]);
iGnt.resfuncgroup        = residueMap.allresidues('GlcNAc');
galResType               = residueMap.allresidues('Gal');
iGnt.resAtt2FG           = galResType;
manBond                  = GlycanBond('4','1');
iGnt. linkresAtt2FG      = struct('bond',manBond,'anomer','b');
glcnacbond               = GlycanBond('3','1');
iGnt.linkFG              = struct('anomer','b','bond',glcnacbond);
iGnt.substMinStruct      = glycanMLread('igntminstruct.glycoct_xml');

%Fut7/4
Futa23                     = GTEnz([2;4;1;152]);
Futa23.isTerminalTarget    = false;
Futa23.resfuncgroup        = residueMap.allresidues('Fuc');
fuctbond                   = GlycanBond('3','1');
Futa23.linkFG              = struct('anomer','a','bond',fuctbond);
glcnacResType              = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','1');
Futa23.resAtt2FG           = glcnacResType;
Futa23.linkresAtt2FG       = struct('bond', glcnacBond,'anomer','b');
Futa23.targetbranchcontain = glycanMLread('lactose.glycoct_xml');
Futa23.isTerminalTarget    = false;



enzDB = containers.Map;
enzDB('mgat5')     = mgat5;
enzDB('mgat4')     = mgat4;
enzDB('mgat3')     = mgat3;
enzDB('mgat2')     = mgat2;
enzDB('mgat1')     = mgat1;
enzDB('manii')     = manii;
enzDB('mani')      = mani;
enzDB('Fut8')      = Fut8;
enzDB('iGnt')      = iGnt;
enzDB('b4GalI')    = b4GalI;
enzDB('STGalTs')   = STGalTs;
enzDB('Futa23')    = Futa23;

save('HL60enzDB.mat','enzDB');
end