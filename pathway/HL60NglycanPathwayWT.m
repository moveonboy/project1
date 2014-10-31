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


enzArray=CellArrayList;
enzArray.add(mgat1);
enzArray.add(mgat2);
enzArray.add(mgat3);
enzArray.add(mgat4);
enzArray.add(mgat5);
enzArray.add(manii);
enzArray.add(mani);
enzArray.add(iGnt);
enzArray.add(b4GalI);
enzArray.add(Fut8);
enzArray.add(STGalTs);
enzArray.add(Futa23);

%define the end prod
a1   = GlycanSpecies(glycanMLread('m3gn.glycoct_xml'));
% a2   = GlycanSpecies(glycanMLread('m5gn2.glycoct_xml'));
% a3   = GlycanSpecies(glycanMLread('2635hybrid.glycoct_xml'));
% a4   = GlycanSpecies(glycanMLread('2966bi_coreFuc_NeuAc.glycoct_xml'));
% a5   = GlycanSpecies(glycanMLread('3211bi_bisec_coreFuc_NeuAc.glycoct_xml'));
% a6   = GlycanSpecies(glycanMLread('3314bi_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a7   = GlycanSpecies(glycanMLread('3559bi_bisec_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a8   = GlycanSpecies(glycanMLread('3776tri_b14_coreFuc_NeuAc.glycoct_xml'));
% a9   = GlycanSpecies(glycanMLread('3776tri_b16_coreFuc_NeuAc.glycoct_xml'));
% a10  = GlycanSpecies(glycanMLread('4021tri_bisec_b14_coreFuc_NeuAc.glycoct_xml'));
% a11  = GlycanSpecies(glycanMLread('4021tri_bisec_b16_coreFuc_NeuAc.glycoct_xml'));
% a12  = GlycanSpecies(glycanMLread('4299tri_b14_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a13  = GlycanSpecies(glycanMLread('4299tri_b16_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a14  = GlycanSpecies(glycanMLread('4544tri_bisec_b14_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a15  = GlycanSpecies(glycanMLread('4544tri_bisec_b16_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a16  = GlycanSpecies(glycanMLread('4587tetra_coreFuc_NeuAc.glycoct_xml'));
% a17  = GlycanSpecies(glycanMLread('4574tetra_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a18  = GlycanSpecies(glycanMLread('4832tetra_bisec_coreFuc_NeuAc.glycoct_xml'));
% a19  = GlycanSpecies(glycanMLread('5472tetra_2lac.glycoct_xml'));
% a20  = GlycanSpecies(glycanMLread('5268tetra_bisec_1lac.glycoct_xml'));
% a21  = GlycanSpecies(glycanMLread('5384tetra_1lac.glycoct_xml'));
% a22  = GlycanSpecies(glycanMLread('4819tetra_bisec_coreFuc_terFuc_NeuAc.glycoct_xml'));
% a23  = GlycanSpecies(glycanMLread('5717tetra_bisec_2lac.glycoct_xml'));
a24  = GlycanSpecies(glycanMLread('5921tetra_3lac.glycoct_xml'));
a25  = GlycanSpecies(glycanMLread('6371tetra_4lac.glycoct_xml'));
% a26  = GlycanSpecies(glycanMLread('6820tetra_5lac.glycoct_xml'));
% a27  = GlycanSpecies(glycanMLread('5837tetra_6lac.glycoct_xml'));
% a27  = GlycanSpecies(glycanMLread('3142tetra_coreFuc.glycoct_xml'));

glycanArray = CellArrayList;
glycanArray.add(a1);
% glycanArray.add(a2);
% glycanArray.add(a3);
% glycanArray.add(a4);
% glycanArray.add(a5);
% glycanArray.add(a6);
% glycanArray.add(a7);
% glycanArray.add(a8);
% glycanArray.add(a9);
% glycanArray.add(a10);
% glycanArray.add(a11);
% glycanArray.add(a12);
% glycanArray.add(a13);
% glycanArray.add(a14);
% glycanArray.add(a15);
% glycanArray.add(a16);
% glycanArray.add(a17);
% glycanArray.add(a18);
% glycanArray.add(a19);
% glycanArray.add(a20);
% glycanArray.add(a21);
% glycanArray.add(a22);
% glycanArray.add(a23);
glycanArray.add(a24);
glycanArray.add(a25);
% glycanArray.add(a26);
% glycanArray.add(a27);

%Perform reaction
displayOptions = displayset('showMass',true,'showLinkage',true,'showRedEnd',true);
fprintf(1,'Input of glycan product structure is \n');
% glycanViewer(a1.glycanStruct,displayOptions);
% glycanViewer(a2.glycanStruct,displayOptions);
% glycanViewer(a3.glycanStruct,displayOptions);
% glycanViewer(a4.glycanStruct,displayOptions);
% glycanViewer(a5.glycanStruct,displayOptions);
% glycanViewer(a6.glycanStruct,displayOptions);
% glycanViewer(a7.glycanStruct,displayOptions);
% glycanViewer(a8.glycanStruct,displayOptions);
% glycanViewer(a9.glycanStruct,displayOptions);
% glycanViewer(a10.glycanStruct,displayOptions);
% glycanViewer(a11.glycanStruct,displayOptions);
% glycanViewer(a12.glycanStruct,displayOptions);
% glycanViewer(a13.glycanStruct,displayOptions);
% glycanViewer(a14.glycanStruct,displayOptions);
% glycanViewer(a15.glycanStruct,displayOptions);
% glycanViewer(a16.glycanStruct,displayOptions);
% glycanViewer(a17.glycanStruct,displayOptions);
% glycanViewer(a18.glycanStruct,displayOptions);
% glycanViewer(a19.glycanStruct,displayOptions);
% glycanViewer(a20.glycanStruct,displayOptions);
% glycanViewer(a21.glycanStruct,displayOptions);
% glycanViewer(a22.glycanStruct,displayOptions);
% glycanViewer(a23.glycanStruct,displayOptions);
% glycanViewer(a24.glycanStruct,displayOptions);
% glycanViewer(a25.glycanStruct,displayOptions);
% glycanViewer(a26.glycanStruct,displayOptions);
% glycanViewer(a27.glycanStruct,displayOptions);
[isPath,nlinkedpath]=inferGlyConnPath(glycanArray,enzArray);
fprintf(1,'Inferred network is shown below:\n');
glycanPathViewer(nlinkedpath);