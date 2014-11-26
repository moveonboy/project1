function tests = testglycanMolWt
    tests = functiontests(localfunctions);
end

function testfunc1(testCase)
options.ion         ='none'; 
options.methylation =false;
objtest1            =num2str(glycanMolWt('shn',options)); 
calcvalue           =objtest1;
expvalue            =num2str(674.2382);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc2(testCase)
options.ion         ='none'; 
options.methylation =false;
objtest2            =num2str(glycanMolWt('sa2,3hb1,3{s2,6}n',options));
calcvalue           =objtest2;
expvalue            =num2str(965.3336);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc3(testCase)
options.ion         ='none'; 
options.methylation =false;
objtest3            =num2str(glycanMolWt('h(s)n',options));
calcvalue           =objtest3;
expvalue            =num2str(674.2382);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc4(testCase)
options.ion         ='none'; 
options.methylation =false;
slex='fshn';
objtest4            =num2str(glycanMolWt(strcat(slex,'hn'),options));
calcvalue           =objtest4;
expvalue            =num2str(1185.4283);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc5(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest5            =num2str(glycanMolWt('hhhhhhnn',options));
calcvalue           =objtest5;
expvalue            =num2str(1783.9);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc6(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest6            =num2str(glycanMolWt('hhhhhhhhhhhhhhhhhnnnnnnnnnnnnnnnnf',options));
calcvalue           =objtest6;
expvalue            =num2str(7632.8);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc7(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest7            =num2str(glycanMolWt('sshhhhhnnnn',options));
calcvalue           =objtest7;
expvalue            =num2str(2792.3);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc8(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest8            =num2str(glycanMolWt('hhhhhhhhhhnnnnnnnnnf',options));
calcvalue           =objtest8;
expvalue            =num2str(4489.9);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc9(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest9            =num2str(glycanMolWt('hhhhhhhhhhhhhnnnnnnnnnnnnf',options));
calcvalue           =objtest9;
expvalue            =num2str(5837.0);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc10(testCase)
options.ion         ='Na'; 
options.methylation =true;
objtest9            =num2str(glycanMolWt('hhhhhhhhhhhhhnnnnnnnnnnnnss',options));
calcvalue           =objtest9;
expvalue            =num2str(6386.2);
verifyEqual(testCase,calcvalue,expvalue);
end

