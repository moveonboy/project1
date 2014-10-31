function MSdata=CleanSpectra(varargin)
%CleanSpectra load the rawdata and use msprocess code to convert raw peak 
% data to a peak list and their half height width column 
%
% MSdata = CleanSpectra(varargin) load the rawdata and the defalut
% options and use the funciton 'msprocess' to convert raw peak data to a 
% peak list and their half height width column using a four-step method. 
%
% Example:
%     MSdata=CleanSpectra('rawdata.mat')
%
%See also MSprocess.

% Author: Yusen Zhou & Gang Liu
% Date Lastly Updated: 9/22/14

narginchk(0,2)

options.OverSegmentationFilter =0.9;
options.HEIGHTFILTER           =0.5;
options.showplot               =true;

if(isempty(varargin))
    rawdata = 'rawdata.mat';
elseif(length(varargin)==1)
    rawdata = varargin{1};
elseif(length(varargin)==2)
    rawdata = varargin{1};
    options.OverSegmentationFilter = varargin{2}.OverSegmentationFilter;
    options.HEIGHTFILTER           = varargin{2}.HEIGHTFILTER;
    options.showplot               = varargin{2}.showplot;
end

load(rawdata);

[peaklist,pfwhh] = msprocess(data,options);
MSdata = struct('peaklist',peaklist,'pfwhh',pfwhh);
save('MSdata.mat','MSdata');
end