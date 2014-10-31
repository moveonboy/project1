function glycangroupexptUB = createGlycanHL60Input()

% 2431
group1.fraction        = 0.107;
group1.glycanlist(1,1) = GlycanSpecies(glycanMLread('2431bracket.glycoct_xml'));

% 2489
group2.fraction        = 0.928;
group2.glycanlist(1,1) = GlycanSpecies(glycanMLread('2489.25.glycoct_xml'));

% 2605
group3.fraction        = 0.108;
group3.glycanlist(1,1) = GlycanSpecies(glycanMLread('2605bracket.glycoct_xml'));

% 2676
group4.fraction        = 0.023;
group4.glycanlist(1,1) = GlycanSpecies(glycanMLread('2676.33_a13bracket.glycoct_xml'));
group4.glycanlist(2,1) = GlycanSpecies(glycanMLread('2676.33_a16bracket.glycoct_xml'));

% 2779
group5.fraction        = 0.042;
group5.glycanlist(1,1) = GlycanSpecies(glycanMLread('2779.39bracket.glycoct_xml'));

% 2792
group6.fraction        = 0.030;
group6.glycanlist(1,1) = GlycanSpecies(glycanMLread('2792.38.glycoct_xml'));

% 2850
group7.fraction        = 0.097;
group7.glycanlist(1,1) = GlycanSpecies(glycanMLread('2850bracket.glycoct_xml'));

% 2966
group8.fraction        = 0.074;
group8.glycanlist(1,1) = GlycanSpecies(glycanMLread('2966.47.glycoct_xml'));

% 3054
group9.fraction        = 0.032;
group9.glycanlist(1,1) = GlycanSpecies(glycanMLread('3054.52_b14bracket.glycoct_xml'));
group9.glycanlist(2,1) = GlycanSpecies(glycanMLread('3054.52_b16bracket.glycoct_xml'));

% 3140
group10.fraction        = 0.021;
group10.glycanlist(1,1) = GlycanSpecies(glycanMLread('3140bracket.glycoct_xml'));

% 3211
group11.fraction        = 0.055;
group11.glycanlist(1,1) = GlycanSpecies(glycanMLread('3211.60.glycoct_xml'));

% 3241
group12.fraction        = 0.026;
group12.glycanlist(1,1) = GlycanSpecies(glycanMLread('3241b14bracket.glycoct_xml'));
group12.glycanlist(2,1) = GlycanSpecies(glycanMLread('3241b16bracket.glycoct_xml'));

% 3299
group13.fraction        = 0.033;
group13.glycanlist(1,1) = GlycanSpecies(glycanMLread('3299b14bracket.glycoct_xml'));
group13.glycanlist(2,1) = GlycanSpecies(glycanMLread('3299b16bracket.glycoct_xml'));

% 3415
group14.fraction        = 0.058;
group14.glycanlist(1,1) = GlycanSpecies(glycanMLread('3415b14bracket.glycoct_xml'));
group14.glycanlist(2,1) = GlycanSpecies(glycanMLread('3415b16bracket.glycoct_xml'));

% 3589
group15.fraction        = 0.034;
group15.glycanlist(1,1) = GlycanSpecies(glycanMLread('3589b14bracket.glycoct_xml'));
group15.glycanlist(2,1) = GlycanSpecies(glycanMLread('3589b16bracket.glycoct_xml'));

% 3660
group16.fraction        = 0.048;
group16.glycanlist(1,1) = GlycanSpecies(glycanMLread('3660b14bracket.glycoct_xml'));
group16.glycanlist(2,1) = GlycanSpecies(glycanMLread('3660b16bracket.glycoct_xml'));

% 3834
group17.fraction        = 0.029;
group17.glycanlist(1,1) = GlycanSpecies(glycanMLread('3834b14bracket.glycoct_xml'));
group17.glycanlist(2,1) = GlycanSpecies(glycanMLread('3834b16bracket.glycoct_xml'));

% 3864
group18.fraction        = 0.029;
group18.glycanlist(1,1) = GlycanSpecies(glycanMLread('3864bracket.glycoct_xml'));

% 3024
group19.fraction        = 0.021;
group19.glycanlist(1,1) = GlycanSpecies(glycanMLread('3024bracket.glycoct_xml'));

% 3228
group20.fraction        = 0.022;
group20.glycanlist(1,1) = GlycanSpecies(glycanMLread('3228b14bracket.glycoct_xml'));
group20.glycanlist(2,1) = GlycanSpecies(glycanMLread('3228b16bracket.glycoct_xml'));

% 3602
group21.fraction        = 0.019;
group21.glycanlist(1,1) = GlycanSpecies(glycanMLread('3602.78b14.glycoct_xml'));
group21.glycanlist(2,1) = GlycanSpecies(glycanMLread('3602.78b16.glycoct_xml'));

glycangroupexptUB = CellArrayList;
glycangroupexptUB.add(group1);
glycangroupexptUB.add(group2);
glycangroupexptUB.add(group3);
glycangroupexptUB.add(group4);
glycangroupexptUB.add(group5);
glycangroupexptUB.add(group6);
glycangroupexptUB.add(group7);
glycangroupexptUB.add(group8);
glycangroupexptUB.add(group9);
glycangroupexptUB.add(group10);
glycangroupexptUB.add(group11);
glycangroupexptUB.add(group12);
glycangroupexptUB.add(group13);
glycangroupexptUB.add(group14);
glycangroupexptUB.add(group15);
glycangroupexptUB.add(group16);
glycangroupexptUB.add(group17);
glycangroupexptUB.add(group18);
glycangroupexptUB.add(group19);
glycangroupexptUB.add(group20);
glycangroupexptUB.add(group21);
end
