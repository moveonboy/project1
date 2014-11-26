function tests = testglycanFormula
    tests = functiontests(localfunctions);
end

function testfunc1(testCase)
glycanstring='hhhhhnn';
calcvalue= glycanFormula(glycanstring);
expvalue = struct('C',69,'H',124,'O',36,'N',2,'S',0,'P',0,'Na',1);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc2(testCase)
options.methylation =false;
options.ion='none';
glycanstring='hhhhhnn';
calcvalue= glycanFormula(glycanstring,options);
expvalue = struct('C',46,'H',78,'O',36,'N',2,'S',0,'P',0);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc3(testCase)
options.methylation =false;
options.ion='none';
glycanstring='hhhhhnnfff';
calcvalue= glycanFormula(glycanstring,options);
expvalue = struct('C',64,'H',108,'O',48,'N',2,'S',0,'P',0);
verifyEqual(testCase,calcvalue,expvalue);
end


function testfunc4(testCase)
options.methylation =false;
options.ion='none';
glycanstring='hhhhhhhhhhhh';
calcvalue= glycanFormula(glycanstring,options);
expvalue = struct('C',72,'H',122,'O',61,'N',0,'S',0,'P',0);
verifyEqual(testCase,calcvalue,expvalue);
end


function testfunc5(testCase)
options.methylation =true;
options.ion='none';
glycanstring='hfsnn';
calcvalue= glycanFormula(glycanstring,options);
expvalue = struct('C',57,'H',101,'O',28,'N',3,'S',0,'P',0);
verifyEqual(testCase,calcvalue,expvalue);
end
