function varargout = msrelativeabundance(MSdata,glycancomp,varargin)
%msrelativeabundance: Calculate the relative abundance for each potential
% glycans based on the MSdata and potential glycans, and write the result
% into a Excel file.
%
%glycanlist = msrelativeabundance(MSdata,nglycanname,abundanceExcelname)
% Output the relative abundance Excel file and glycanlist of MSdata
%
%[glycanlist abundancelibrary] = msrelativeabundance(MSdata,nglycanname,abundanceExcelname)
% Output the relative abundance Excel file, potential glycanlist of MSdata
% and containers.Map for glycan relative abundance.
%
% For example:
%         [glycanlist,abundancelibrary] = msrelativeabundance('MSdata.mat',..
%                                   'HL60WT_glycanspeceis.xls','HL60WTN_relativeabundance.xls');
%
%
%
% Author: Yusen Zhou 
% Date Lastly Updated: 10/27/14
nargoutchk(1,2);
if(isempty(varargin))
    abundanceExcelname = '';
    StoNratio          = '';
elseif(length(varargin)==1)
    if(isnumeric(varargin{1}))
        abundanceExcelname = '';
        StoNratio          = varargin{1};
    elseif(ischar(varargin{1}))
        abundanceExcelname = varargin{1};
        StoNratio          = '';
    else
        error('WRONG INPUT TYPE')
    end
elseif(length(varargin)==2)
    abundanceExcelname = varargin{1};
    StoNratio = varargin{2};
end

load(MSdata);
peaklist = MSdata.peaklist;
pfwhh    = MSdata.pfwhh;
if(ischar(glycancomp))
    [chowildtypecomps,glycanmwarray] = readCellNGlycanFromExcel(1,glycancomp);
elseif(iscellstr(glycancomp))
    [chowildtypecomps,glycanmwarray] = Calisomwarray(glycancomp);
end
if(~isempty(abundanceExcelname))&&(~isempty(StoNratio))
    [glycanlist,abundancelibrary]= msfraction(peaklist,pfwhh,chowildtypecomps,glycanmwarray,abundanceExcelname,StoNratio);
elseif(isempty(abundanceExcelname))&&(~isempty(StoNratio))
    [glycanlist,abundancelibrary]= msfraction(peaklist,pfwhh,chowildtypecomps,glycanmwarray,StoNratio);
elseif(~isempty(abundanceExcelname))&&(isempty(StoNratio))
    [glycanlist,abundancelibrary]= msfraction(peaklist,pfwhh,chowildtypecomps,glycanmwarray,abundanceExcelname);
else
    [glycanlist,abundancelibrary]= msfraction(peaklist,pfwhh,chowildtypecomps,glycanmwarray);
end

if(nargout==1)
    varargout{1} = glycanlist;
elseif(nargout==2)
    varargout{1} = glycanlist;
    varargout{2} = abundancelibrary;
    save('HL60WTlibrary.mat','HL60library');
end
end