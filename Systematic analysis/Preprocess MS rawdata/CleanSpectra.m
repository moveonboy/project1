function MSdata=CleanSpectra(msrawdatafilename,varargin)
% CleanSpectra: loads raw MS data and use msprocess code to convert raw peak 
% data to a peak list and their half height width column 
%
% MSdata = CleanSpectra(msrawdatafilename) use the defalut value of msrawdatapath
% and options to convert 'msrawdatafilename' file, a '.txt' file, into matlab form 
% and then preprocess it by using a four-step method to finally get the peak 
% list and their half height width column.
%
% MSdata = CleanSpectra(msrawdatafilename,'msrawdatapathstr',msrawdatapath) find the 'msrawdatafilename' 
% file in msrawdatapath and convert it into matlab form and then preprocess it 
% by using a four-step method to finally get the peak list and their half height 
% width column. 
%
% MSdata = CleanSpectra(msrawdatafilename,'optionsstr',options) set up your 
% own options and convert msrawdatafilename file into matlab form and then 
% preprocess it by using a four-step method to finally get the peak list and 
% their half height width column. 
% Example:
%     MSdata=CleanSpectra('HL60WTN_Glycans.txt')
%
% See also MSprocess.

% Author: Yusen Zhou & Gang Liu
% Date Lastly Updated: 9/22/14

narginchk(0,4)

% defalut value
msrawdatapath     = 'C:\Users\Yusen\Dropbox\HL60cells\N_glycan\Rawdata';
options.OverSegmentationFilter =0.9;
options.HEIGHTFILTER           =0.2;
options.showplot               =true;

nvarargin =  numel(varargin);
if rem(nvarargin,2)
    error(message('IncorrectNumberOfArguments'));
end
okargs    = {'msrawdatapathstr','optionsstr'};
for j=1:2:nvarargin
    pname = varargin{j};
    pval = varargin{j+1};
    k = find(strncmpi(pname, okargs,length(pname)));
    if isempty(k)
        error(message('UnknownParameterName'));
    elseif length(k)>1
        error(message('AmbiguousParameterName'));
    else
        switch(k)
            case 1 % msrawdatapathstr
                msrawdatapath = pval;
            case 2 % optionsstr
                options.OverSegmentationFilter = pval.OverSegmentationFilter;
                options.HEIGHTFILTER           = pval.HEIGHTFILTER;
                options.showplot               = pval.showplot;       
        end
    end
end

% read MS raw data
msrawdatafullfilename=fullfile(msrawdatapath,msrawdatafilename);
data = readMS(msrawdatafullfilename);
save('rawdata.mat','data');
% Preprocess of MS data
[peaklist,pfwhh] = msprocess(data,options);
MSdata = struct('peaklist',peaklist,'pfwhh',pfwhh);
save('MSdata.mat','MSdata');
end