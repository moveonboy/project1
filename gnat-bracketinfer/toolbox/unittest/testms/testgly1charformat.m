function tests = testgly1charformat
    tests = functiontests(localfunctions);
end

function testfunc1(testCase)
glycanresiduestring ='Hex2HexNAc5';
objtest1=gly1charformat(glycanresiduestring);
calcvalue=objtest1;
expvalue='nnnnnhh';
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc2(testCase)
glycanresiduestring ='Hex2HexNAc5Fuc10';
objtest2=gly1charformat(glycanresiduestring);
calcvalue=objtest2;
expvalue='nnnnnhhffffffffff';
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc3(testCase)
glycanresiduestring ='Hex1HexNAc4Fuc3';
objtest3=gly1charformat(glycanresiduestring);
calcvalue=objtest3;
expvalue='nnnnhfff';
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc4(testCase)
glycanresiduestring ='Hex10HexNAc4Fuc12';
objtest4=gly1charformat(glycanresiduestring);
calcvalue=objtest4;
expvalue='nnnnhhhhhhhhhhffffffffffff';
verifyEqual(testCase,calcvalue,expvalue);
end