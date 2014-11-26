function [retTime,retUnit] = retrieveMS1RetT(obj,scannum)
% retrieveActMethod read MS activation method from mzXML file 
%   
% Example 1: mzXMLobj  = mzXML('test1.mzXML',0,'memsave');
%            mod       = mzXMLobj.retrieveMS1RetT(2305);
%
% Example 2: mzXMLobj  = mzXML('test1.mzxml',0,'memsave');
%            msspectra = mzXMLobj.retrieveMS1RetT(3705);
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 06/17/14
    retTime   = -1;
    numscan   = obj.mzxmljava.getScanCount;
    ms1scannum = -1;
    for i = 1 : numscan
        if(obj.mzxmljava.rap(i).getNum==scannum)
           if(obj.mzxmljava.rap(i).getMsLevel==1)    
                retentionTime = char(obj.mzxmljava.rap(i).getRetentionTime);
                retention     = regexp(retentionTime,'(?<ret>\d+\.\d+)(?<unit>[a-z_A-Z]+)','names');
                if(~isempty(retention))
                     retTime = str2double(retention.ret);
                     retUnit = retention.unit;    
                else
                     retTime = -1;
                     retUnit = '';
                end
            return
           elseif(obj.mzxmljava.rap(i).getMsLevel==2)  
               ms1scannum = char(obj.mzxmljava.rap(i).getPrecursorScanNum);
               break
           end
        end
    end
    
    if(ms1scannum==-1)
        error('MATLAB:GPAT:ERRORSCANNUM','MS1 SCAN NOT FOUND');
    end
       
    for i = 1 : numscan
       if(obj.mzxmljava.rap(i).getNum==ms1scannum)
           retentionTime = char(obj.mzxmljava.rap(i).getRetentionTime);
           retention     = regexp(retentionTime,'(?<ret>\d+\.\d+)(?<unit>[a-z_A-Z]+)','names');
           if(~isempty(retention))
                 retTime = str2double(retention.ret);
                 retUnit = retention.unit;    
           else
                 retTime = -1;
                 retUnit = '';
           end
           return;
       end
    end
end