classdef mzXML < handle
    %MZXML: MS data stored in MZXML file format
    %
    % Example: mzXMLobj = mzXML('test1.mzxml');
    %          mzXMLobj = mzXML('test1.mzxml',1);
    %          mzXMLobj = mzXML('test1.mzxml',0,'memsave');
    %          mzXMLobj = mzXML('test1.mzxml',1,'fastload');
    %See Also getMSSpectra.
    
    % Author: Gang Liu
    % Date Lastly updated: 1/26/14
    properties
        mzxmljava;
        extmzxmljava;
        scandata;
        instrumentinfo;
        dataprocessinginfo;
    end
    
    methods
        function obj=mzXML(varargin)
            import org.systemsbiology.jrap.MSXMLParser;
            import java.io.File;
            import org.systemsbiology.jrap.extension.ExtMSXMLParser;
            
            if(nargin==1) || (nargin==2) || (nargin==3)
                fileexist = exist(varargin{1},'file')==2;
                if(~fileexist)
                    error('MATLAB:GPAT:FILENOTFOUND','mzXML FILE IS NOT FOUND');
                end
                
                mzxmlfilename = varargin{1};
                [pathstr,name,ext] = fileparts(mzxmlfilename);
                if(~strcmpi(ext,'.mzxml')&& ~strcmpi(ext,'.xml'))
                    error('MATLAB:GPAT:WRONGINPUTFILETYPE','WRONG INPUT FILE TYPE');
                end
                
                if(isempty(pathstr))
                    mzxmlfilename = which(mzxmlfilename);
                end
                
                if(nargin==2)||(nargin==3)
                    showpro = varargin{2};
                else
                    showpro = 0;
                end
            end
            
            if(nargin==3)
                if(strcmpi(varargin{3},'memsave'))
                    memObjsave = 1;
                    fastload  = 0;
                elseif(strcmpi(varargin{3},'fastload'))
                    memObjsave = 0;
                    fastload =  1;
                end
            else
                memObjsave = 0;
                fastload   = 0;
            end
            
            try
                obj.mzxmljava = MSXMLParser(mzxmlfilename);
                obj.extmzxmljava = ExtMSXMLParser(mzxmlfilename);
            catch err
                error('MATLAB:GlycoPAT:INCORRECTFILE','INCORRECT mzXML FILE');
            end
            
            if(~memObjsave)
                if(~fastload)
                    obj.setScanData(showpro);
                    obj.setInstrumentInfo();
                    obj.setDataProcessInfo();
                else
                    obj.setScanData(showpro,1);
                end
            end            
        end
    end
end

