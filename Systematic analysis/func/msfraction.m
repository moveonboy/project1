function [glycanlist,HL60library]= msfraction(peaklist,pfwhh,expecGlycan,glycanmwarray,varargin)
%msfraction
% 
% [glycanlist]= msfraction(peaklist,pfwhh,expecGlycan,glycanmwarray) uses default options which
%   are 0.3 for oversegmentationfilter, 5 for heightfilter, true for
%   showplot. The peak list p is a matrix with two columns of peak
%   locations, and peak intensities. The matrix pfwhh has two columns
%   indicating the left and right location of the full width height.
%
% Example:
%     
%See also readMS,readMsprocess.

%Date Lastly Updated: 02/17/14
outputfilename = '';
if(length(varargin)==1)
    outputfilename=varargin{1};
    StoNratio = 3;
elseif(length(varargin)==2)
    outputfilename=varargin{1};
    StoNratio = varargin{2};
end

relativeabundance  = [];
allpeaktotalarea   =0;
monoisomw=[];
abundantmw=[];

index = 0;
poitnumber = 0;
deleteindex =[];
for i=1:length(expecGlycan)
    glycanindex = '';
    ithmddouble = '';
    if(poitnumber>=i)
        continue
    end
    isoverlapping  = 0;
    ithmddouble{1}    = glycanmwarray{i,1};
    isglycaninMS = checkglycaninMS(ithmddouble{1},peaklist,StoNratio);
    if(~isglycaninMS)
        deleteindex(end+1) = i;
        continue
    end
    minithglycanmw = ithmddouble{1}(1,1);
    maxithglycanmw = ithmddouble{1}(7,1);
    counter        = 1;
    glycanindex{1} = i;
    for j = i+1 : length(expecGlycan)
        jthmddouble    = glycanmwarray{j,1};
        minjthglycanmw = jthmddouble(1,1);
        maxjthglycanmw = jthmddouble(7,1);
        if((maxjthglycanmw-minithglycanmw)*(maxithglycanmw-minjthglycanmw)>=0)
            isoverlapping = 1;
            ithmddouble{counter+1} = jthmddouble;
            glycanindex{end+1} = j;
            counter                = counter+1;
        end
    end
    if(isoverlapping)
        peakarea = findoverlaparea(ithmddouble,peaklist,pfwhh);
    else
        peakarea = findpeakarea(ithmddouble{1},peaklist,pfwhh);
    end
    for j = 1 : length(peakarea)
        num = glycanindex{j};
        jthpeakarea = peakarea{j};
        relativeabundance(index+1,1)         = jthpeakarea;
        allpeaktotalarea                     = allpeaktotalarea+jthpeakarea;
        monoisomw{index+1,1}=glycanmwarray{num,1}(1,1);
        maxfraction=max(glycanmwarray{num,1}(:,2));
        numrow=find(~(glycanmwarray{num,1}(:,2)-maxfraction));
        abundantmw{index+1,1}=glycanmwarray{num,1}(numrow,1);
        glycanresiduestring=expecGlycan{num,1};
        glycan1letstring=gly1charformat(glycanresiduestring);
        options.mono=false;
        avaragemw{index+1,1}=glycanMolWt(glycan1letstring,options);
        index                                = index+1;
    end
    poitnumber = i+length(peakarea)-1;
end    

expecGlycan(deleteindex) = '';
for i = 1 : length(expecGlycan)
    relativeabundance(i,1)       = relativeabundance(i,1)/allpeaktotalarea;
end

for i=1:length(expecGlycan)
    glycanlist(i).relativeabundance=relativeabundance(i,1);
    glycanlist(i).glycancompos     = expecGlycan(i);
end

if(~isempty(outputfilename))
    A1=cellstr('Composition');
    B1=cellstr('Fraction');
    C1=cellstr('Monoisotopic mass');
    D1=cellstr('Abundant mass');
    E1=cellstr('Avarage mass');
    xlswrite(outputfilename,relativeabundance,1,'B2');
    xlswrite(outputfilename,expecGlycan,1,'A2');
    xlswrite(outputfilename, monoisomw,1,'C2');
    xlswrite(outputfilename, abundantmw,1,'D2');
    xlswrite(outputfilename, avaragemw,1,'E2');
    xlswrite(outputfilename,A1,1,'A1');
    xlswrite(outputfilename,B1,1,'B1');
    xlswrite(outputfilename,C1,1,'C1');
    xlswrite(outputfilename,D1,1,'D1');
    xlswrite(outputfilename,E1,1,'E1');
end

HL60library = buildlibrary(monoisomw,relativeabundance);

end

function HL60library = buildlibrary(monoisomw,relativeabundance)
HL60library = containers.Map;
for i = 1 : length(monoisomw)
    ithmass       = num2str(round(monoisomw{i}));
    ithrelativeabundance = relativeabundance(i);
    HL60library(ithmass) = ithrelativeabundance;
end
end

function isglycaninMS = checkglycaninMS(md,peaklist,StoNratio)
mz   = md(:,1);
dist = md(:,2);

% two criteria to define if the glycan exist in MS.

% First criteria: find the monoiso
isglycanexist   = 0;
isopeakexistnum = 0;
count           = 0;
for i = 1 : length(dist)
    ithdist = dist(i);
    if(ithdist>0.05)||(i==1)
        count = count+1;
        isotopicmz = mz(i,1);
        isotopicindex = findclosetpeak(isotopicmz,peaklist);
        if(~isempty(isotopicindex))
            isopeakexistnum =isopeakexistnum+1;
        end
    end
end
if(isopeakexistnum==count)
    isglycanexist = 1;
end
    

% Second criteria: the peak we find should to be the first peak in vicinity
isfirstpeak    = 0;
monoisotopicmz = mz(1,1);
preisotopicmz = monoisotopicmz - 1;
preisotopicindex = findclosetpeak(preisotopicmz,peaklist);
if(isempty(preisotopicindex))
    isfirstpeak = 1;
else
    twopeakratio = peaklist(preisotopicindex+1,2)/peaklist(preisotopicindex,2);
    if(twopeakratio>StoNratio)
        isfirstpeak = 1;
    end     
end

isglycaninMS = (isglycanexist)&&(isfirstpeak);
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



    
    



