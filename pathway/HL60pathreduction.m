% test bracket species
clc;
clear;
storedvarpath = load('HL60WTN_glycanPathwayPart1.mat');
binglycanpath = storedvarpath.nlinkedpath;
binglycanpath.resetjava;
glycanPathViewer(binglycanpath);
fprintf(1,'network size: %d Species *%d Reactions\n',...
    binglycanpath.getNSpecies,binglycanpath.getNReactions);
newpathway = pathreduction(binglycanpath);
glycanPathViewer(newpathway);
fprintf(1,'network size: %d Species *%d Reactions\n',...
    newpathway.getNSpecies,newpathway.getNReactions);