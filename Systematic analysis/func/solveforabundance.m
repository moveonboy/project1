function [peaktotalarea,everypeak] = solveforabundance(md,otherisomd,eqnindex,eqnarea,relativeabundance)
dist             = md(:,2);
if(length(otherisomd)==1)
    overlappingindex = [];
    counter          = 0;
    otherdist        = otherisomd.isomwarray(:,2);
    for i = 1 : length(eqnindex)
        ithmd         = eqnindex{i};
        for j = 1 : length(ithmd)
            jthothermd = ithmd{j};
            if(~size(jthothermd,1)==0)
                overlappingindex(counter+1) = i;
                counter                     = counter+1;
            end
        end
    end
    
    isoverlap = (overlappingindex(1)==1);
    if(~isoverlap)
        peakareaforfive = 0;
        peakdistforfive = 0;
        for i = 1 : overlappingindex(1)-1
            peakareaforfive = peakareaforfive+eqnarea(i);
            peakdistforfive = peakdistforfive+dist(i);
        end
        peaktotalarea   = peakareaforfive/peakdistforfive;
    else
        peakareaforfive     = 0;
        for i = 1 : length(eqnarea)
            peakareaforfive = peakareaforfive+eqnarea(i);
        end
        lastindex        = length(relativeabundance);
        try
            otherconc        = relativeabundance(lastindex);
        catch err
            debug
        end
        for i = 1 : length(eqnindex)
            otherindex    = eqnindex{i}{1,1}.index;
            peakareaforfive = peakareaforfive-otherconc*otherdist(otherindex);
        end
        if(peakareaforfive<0)
            peaktotalarea = 0;
        else
            peaktotalarea = peakareaforfive;
        end
    end
elseif(length(otherisomd)==2)
    otherisomd1 = otherisomd(1);
    otherisomd2 = otherisomd(2);
    overlappingindex1 = [];
    overlappingindex2 = [];
    counter1          = 0;
    counter2          = 0;
    otherdist1        = otherisomd1.isomwarray(:,2);
    otherdist2        = otherisomd2.isomwarray(:,2);
    for i = 1 : length(eqnindex)
        ithmd         = eqnindex{i}{1,1};
        if(~isempty(ithmd))
            overlappingindex1(counter1+1) = i;
            counter1                      = counter1+1;
        end
    end
    
    for i = 1 : length(eqnindex)
        ithmd         = eqnindex{i}{1,2};
        if(~isempty(ithmd))
            overlappingindex2(counter2+1) = i;
            counter2                      = counter2+1;
        end
    end
    
    if(overlappingindex1(1)~=1)&&(overlappingindex2(1)~=1)
        firstoverlappeak = min(overlappingindex1(1),overlappingindex2(1));
        peakareaforfive = 0;
        peakdistforfive = 0;
        for  i = 1 : firstoverlappeak-1
            peakareaforfive = peakareaforfive+eqnarea(i);
            peakdistforfive = peakdistforfive+dist(i);
        end
        peaktotalarea   = peakareaforfive/peakdistforfive;
    elseif((overlappingindex1(1)==1)&&(overlappingindex2(1)~=1))
        lastindex        = length(relativeabundance);
        otherconc        = relativeabundance(lastindex);
        peakareaforfive = 0;
        peakdistforfive = 0;
        for i = 1 : overlappingindex2(1)-1
            otherindex       = eqnindex{i}{1,1}.index;
            peakareaforother = otherconc*otherdist1(otherindex);
            peakareaforfive  = peakareaforfive+eqnarea(i)-peakareaforother;
            peakdistforfive  = peakdistforfive+dist(i);
        end
        if(peakareaforfive<0)
            peaktotalarea = 0;
        else
            peaktotalarea = peakareaforfive/peakdistforfive;
        end
    elseif((overlappingindex1(1)==1)&&(overlappingindex2(1)==1))
        peakareaforfive     = 0;
        for i = 1 : length(eqnarea)
            peakareaforfive = peakareaforfive+eqnarea(i);
        end
        lastindex1        = length(relativeabundance)-1;
        otherconc1        = relativeabundance(lastindex1);
        for i = 1 : length(eqnindex)
            otherindex1    = eqnindex{i}{1,1}.index;
            peakareaforfive = peakareaforfive-otherconc1*otherdist1(otherindex1);
        end
        lastindex2        = length(relativeabundance);
        otherconc2        = relativeabundance(lastindex2);
        for i = 1 : length(eqnindex)
            otherindex2    = eqnindex{i}{1,2}.index;
            peakareaforfive = peakareaforfive-otherconc2*otherdist2(otherindex2);
        end
        if(peakareaforfive<0)
            peaktotalarea = 0;
        else
            peaktotalarea = peakareaforfive;
        end
    end
end

if(length(eqnarea)==6)
    everypeak = zeros(6,3);
elseif(length(eqnarea)==5)
    everypeak = zeros(5,3);
end
for i = 1 : length(eqnarea)
    everypeak(i,1) = md(i,1);
    everypeak(i,2) = eqnarea(i);
    everypeak(i,3) = peaktotalarea*md(i,2);
end
end
    
