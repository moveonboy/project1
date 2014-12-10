function glycanList = FindDiffGlycan(glycanlist1,glycanlist2)
%FindDiffGlycan is to compare two glycanlist and find the glycans that are
% unique in each glycanlist.
% 
% glycanList = FindDiffGlycan(glycanlist1,glycanlist2) compare the two
%  glycanlist2 to glycanlist1, and find the differences between these two
%  glycanlist.
%
% Input  : Two glycanlist.
% Output : glycanlist with unique glycans.
%
% Example:
%     glycanList = FindDiffGlycan(HL60WT,HL60KO);
%
%
%Author: Yusen Zhou
%Date Lastly Updated: 11/10/14

if(isempty(glycanlist1))||...
        (isempty(glycanlist2))
    error(message('IncorrectInputs'));
end

comparegroup1 = getcompos(glycanlist1);
comparegroup2 = getcompos(glycanlist2);

glycansincommon1 = [];
count           = 0;
for i = 1 : length(comparegroup1);
    ithcomp    = comparegroup1{i};
    isincommon = sum(strcmpi(comparegroup2,ithcomp));
    if(isincommon)
        glycansincommon1(count+1) = i;
        count = count+1;
    end  
end
glycansincommon2 = [];
count            = 0;
for i = 1 : length(comparegroup2);
    ithcomp    = comparegroup2{i};
    isincommon  = sum(strcmpi(comparegroup1,ithcomp));
    if(isincommon)
        glycansincommon2(count+1) = i;
        count = count+1;
    end  
end
comparegroup1(glycansincommon1)='';
comparegroup2(glycansincommon2)='';
glycanList.ListA = comparegroup1;
glycanList.ListB = comparegroup2;
end

function comparegroup = getcompos(glycanlist)
comparegroup = {};
count = 0;
for i = 1 : length(glycanlist)
    comparegroup(count+1)   = glycanlist(i).glycancompos;
    count = count+1;
end
end