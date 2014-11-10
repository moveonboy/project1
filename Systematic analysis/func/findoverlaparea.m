function  glycantotalarea= findoverlaparea(md,peaklist,pfwhh)
numofglycans  = length(md);
glycantotalarea = '';
count         = 0;
for i = 1 : numofglycans
    if(i==1)
    ithglycanmz      = md{i}(:,1);
    ithglycandist    = md{i}(:,2);
    afterglycanmz    = md{i+1}(:,1);
    overlappeakindex = findoverlappeak(ithglycanmz,afterglycanmz);
    partialdist      = 0;
    partialpeakarea  = 0;
    for j = 1 : (overlappeakindex-1)
        dist          = ithglycandist(j,1);
        isotopicmz    = ithglycanmz(j,1);
        isotopicindex = findclosetpeak(isotopicmz,peaklist);
        if(isempty(isotopicindex))
            break
        end
        ithpeakarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
            *peaklist(isotopicindex,2);
        partialpeakarea = partialpeakarea + ithpeakarea;
        partialdist     = partialdist+dist;
    end
    glycantotalarea{count+1} = partialpeakarea/partialdist;
    count = count+1;
    elseif(i==numofglycans)
        glycan1mz     = md{1}(:,1);
        maxnumpeaks   = 6*length(md)-(length(md)-1);
        peak1         = glycan1mz(1,1);
        isotopicindex = findclosetpeak(peak1,peaklist);
        try
        peakarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
            *peaklist(isotopicindex,2);
        catch err
            debug
        end
        for j = 1 : maxnumpeaks-1
            jthpeakmz = peak1+j;
            isotopicindex = findclosetpeak(jthpeakmz,peaklist);
            if(~isempty(isotopicindex))
                jthpeakarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
            *peaklist(isotopicindex,2);
            else
                break
            end
            peakarea = peakarea+jthpeakarea;
        end
        ithglycanarea = peakarea;
        for j = 1 : i-1
            jthglycanarea = glycantotalarea{j};
            ithglycanarea = ithglycanarea - jthglycanarea;
        end
        if(ithglycanarea<0)
            glycantotalarea{count+1} = 0;
            count = count+1;
        else
            glycantotalarea{count+1} = ithglycanarea;
            count = count+1;
        end
    else
        ithglycanmz      = md{i}(:,1);
        ithglycandist    = md{i}(:,2);
        dist             = ithglycandist(1,1);
        isotopicmz       = ithglycanmz(1,1);
        isotopicindex    = findclosetpeak(isotopicmz,peaklist);
        ithpeakarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
            *peaklist(isotopicindex,2);
        for j = 1 : i-1
            jthglycanmz   = md{j}(:,1);
            jthglycandist = md{j}(:,1);
            prepeakindex  = findoverlappeak(jthglycanmz,ithglycanmz);
            peakareafromj = glycantotalarea{j}*jthglycandist(prepeakindex);
            ithpeakarea   = ithpeakarea-peakareafromj;
        end
        if(ithpeakarea<0)
            glycantotalarea{count+1} = 0;
            count = count+1;
        else
            glycantotalarea{count+1} = ithpeakarea/dist;
            count = count+1;
        end
    end
    
end
end

function overlappeakindex = findoverlappeak(glycan1mz,glycan2mz)
overlappeakindex = [];
for i = 1 : length(glycan1mz)
    ithisomz = round(glycan1mz(i,1));
    monoisopeak = round(glycan2mz(1,1));
    if(ithisomz==monoisopeak)
        overlappeakindex = i;
    end
end
end

function isotopicindex = findclosetpeak(isotopicmz,peaklist)
checkisotopicmz = abs(peaklist(:,1)-isotopicmz)<0.99;
isotopicindex   = find(checkisotopicmz);
% if multiple peaks are found close to isotopic mass, select the closet peak as the match
if(length(isotopicindex)>1)
    isotopicmzsdif         = abs(peaklist(checkisotopicmz)-isotopicmz);
    isotopicmultipleindex  = find(checkisotopicmz);
    isotopicnewindex       = find(~abs(isotopicmzsdif-min(isotopicmzsdif)));
    isotopicindex          = isotopicmultipleindex(isotopicnewindex);
end
end

