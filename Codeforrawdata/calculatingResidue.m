function ResidueHL60WT = calculatingResidue(filename)
% input: the '-mat'file that contains the groups of isotopic peaks existing overlapping with
% mass,experiment data and calculation

% output: the residues for each of the groupe

% grouping criteria: based on the isotopic peaks of different glycanspecies
%   that existing overlapping. For example: species1(3140),speceis2(3142).The
%   isotopic peaks for species1(3140) is '3140','3141','3142','3143','3144'
%   and theisotopic peaks for species2(3142) is'3142','3143','3144','3145','3146'
%   so combine them into one group that contains
%   '3140','3141','3142','3143','3144','3145','3146' intotal 7 peaks.

fullfilename = [filename '.mat'];
load(fullfilename);
ResidueHL60WT = struct('monomass',[],'relativetotalres',[]);
counter = 0;
for i = 1 : length(isopeakareagroup)
    ithgroup                          = isopeakareagroup{i};
    ResidueHL60WT(counter+1).monomass = ithgroup(1,1);
    totaldifference                   = 0;
    totalareaexp                      = 0;
    for j = 1 : length(ithgroup)
        differnce       = ithgroup(j,2)-ithgroup(j,3);
        totaldifference = totaldifference+differnce;
        totalareaexp    = totalareaexp+ithgroup(j,2);
    end
    ResidueHL60WT(counter+1).relativetotalres = totaldifference/totalareaexp;
    counter = counter+1;
end
end