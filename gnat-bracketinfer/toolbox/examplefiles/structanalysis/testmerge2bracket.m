clc;clear;
m7{1}=GlycanSpecies(glycanMLread('m7a.glycoct_xml'));
m7{2}=GlycanSpecies(glycanMLread('m7b.glycoct_xml'));
m7{3}=GlycanSpecies(glycanMLread('m7c.glycoct_xml'));
m7{4}=GlycanSpecies(glycanMLread('m7d.glycoct_xml'));
m7bracket=merge2bracket(m7);
glycanViewer(m7bracket);
m6{1}=GlycanSpecies(glycanMLread('m6a.glycoct_xml'));
m6{2}=GlycanSpecies(glycanMLread('m6b.glycoct_xml'));
m6{3}=GlycanSpecies(glycanMLread('m6b.glycoct_xml'));
m6bracket = merge2bracket(m6);
glycanViewer(m6bracket);

% for i=1: length(m6)
%    glycanViewer(m6{i}) 
% end

tetra{1}=GlycanSpecies(glycanMLread('tetra_test1.glycoct_xml'));
tetra{2}=GlycanSpecies(glycanMLread('tetra_test2.glycoct_xml'));
tetra{3}=GlycanSpecies(glycanMLread('tetra_test3.glycoct_xml'));
tetra{4}=GlycanSpecies(glycanMLread('tetra_test4.glycoct_xml'));
for i=1: length(tetra)
   glycanViewer(tetra{i}) 
end

tetrabracket = merge2bracket(tetra);
glycanViewer(tetrabracket);

% for i=1: length(tetra)
%    glycanViewer(tetra{i}) 
% end

disp('end');
