function groupabundance = relativeabundance()
% To calculate the total relative abundance in one group based on three
% library


load('HL60WTlibrary.mat')
load('HL60Groups.mat')
load('HL60grouplibrary.mat')
groupabundance = zeros(length(HL60groups),1);
for i = 1 : length(groupabundance)
    ithgroupindex     = HL60groups{i};
    ithgroupabundance = 0;
    groupmembermass   = '';
    counter           = 0;
    for j = 1 : length(ithgroupindex)
        massindex    = num2str(ithgroupindex(j));
        jthindexmass = num2str(HL60grouplibrary(massindex));
        isnew        = 1;
        for jj = 1 : length(groupmembermass)
            jjthmass = groupmembermass{jj};
            if(isequal(jjthmass,jthindexmass))
                isnew = 0;
                break
            end
        end
        if(isnew)
        ithgroupabundance = ithgroupabundance+HL60library(jthindexmass);   
        groupmembermass{counter+1,1} = jthindexmass;
        counter = counter+1;
        end
    end
    groupabundance(i,1) = ithgroupabundance;
end
end