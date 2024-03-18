function [data,params] = read_cruna_data(fname,component,images)
%
% READ_CRUNA_DATA
%  reads image data files created by the cruna code
%  and creates a mono-image data array
%
%  (in) fname     : file name of the first image
%  (in) component : (optional) component to be read, default = all
%  (in) images    : (optional) list of images to be read, default = all
%
%  (out) data  : data      set of fname (all images)
%  (out) params: parameter set of fname (given image)

%% read params of fname
params = read_cruna_params( fname );

%% determine dnames (data files)
if ~exist('images','var')
    images = 1:(params.parallelism.i1*params.parallelism.i2*params.parallelism.i3);
end

fname_split_ending = strsplit(fname,'.')                 ;
fname_split_type   = strsplit(fname_split_ending{1},'__');
fname_split_image  = strsplit(fname_split_type{2}  ,'_' );

for i = 1:length(images)
    block_nbr = fname_split_image{1};
    image_nbr = num2str(images(i),'%0.5d')  ;
    
    file_idents = [];
    for m = 3:length(fname_split_image)
        file_idents = [file_idents '_' fname_split_image{m}]; %#ok<AGROW>
    end
    
    dnames{i} = [fname_split_type{1} '__' block_nbr '_' image_nbr file_idents '.h5'];   %#ok<SAGROW>
end

%% determine gnames (index files)
for i = 1:length(images)
    block_nbr = fname_split_image{1};
    image_nbr = num2str(images(i),'%0.5d')  ;
    
    gnames{i} = ['geometry_indices__' block_nbr '_' image_nbr '.h5'];                   %#ok<SAGROW>
end

%% create data field
% default layout 3d field
ns = [params.geom.n1,params.geom.n2,params.geom.n3,1,1];

% get dimension of dataset data
[~, DatasetDims, ~] = give_DatasetInfo(dnames{1},'/data');

% extent data layout corresponding to dataset dimension in file
for d = 4:size(DatasetDims,2)
    ns(d) = DatasetDims(d);
end

% reduce number of dimension if only a component is read
if exist('component','var')
    ns(4) = 1;
end

% init field
display(['Estimated variable size: ' num2str(prod(ns)*8/(1000^3),'%2.2f') 'Gb'])
data  = zeros(ns);

%% read data/params
for r = 1:length(dnames)
    disp(['read data file: ' dnames{r}])
    
    params_image = read_cruna_params( dnames{r} );
    
    [~, DatasetDims, DatasetName] = give_DatasetInfo(dnames{r},'/data');    
    start = ones(1,size(DatasetDims,2));
    count = DatasetDims;
    
    if size(DatasetDims,2) > 3
        if exist('component','var')
            start(4) = component;
            count(4) = 1;
        end
    end

    data_image = h5read(dnames{r},DatasetName,start,count);
    
    % assignment
    if isfield(params_image.geom,'xi10i')
        min_xi = params_image.geom.xi10i;
        max_xi = params_image.geom.xi11i;
        min_yi = params_image.geom.xi20i;
        max_yi = params_image.geom.xi21i;
        min_zi = params_image.geom.xi30i;
        max_zi = params_image.geom.xi31i;
        
    elseif exist(gnames{i},'file')
        idxs_image = hdf5read(gnames{r},'geometry_indices') ;
        
        min_xi = min(min(min((idxs_image(:,:,:,1)))));
        max_xi = max(max(max((idxs_image(:,:,:,1)))));
        min_yi = min(min(min((idxs_image(:,:,:,2)))));
        max_yi = max(max(max((idxs_image(:,:,:,2)))));
        min_zi = min(min(min((idxs_image(:,:,:,3)))));
        max_zi = max(max(max((idxs_image(:,:,:,3)))));
    else
        error('No global block assignment available: please code')
    end

    data(min_xi:max_xi,min_yi:max_yi,min_zi:max_zi,:,:) = squeeze(data_image);
end


end

function [DatasetIndex, DatasetDims, DatasetName] = give_DatasetInfo(fname,DatasetName)
    if ~exist('DatasetName','var')
        warning('use default DatasetName: /data')
        DatasetName = '/data';
    end

    finfo = hdf5info(fname);
    
    for i = 1,size(finfo.GroupHierarchy.Datasets,2);
        if strcmp(finfo.GroupHierarchy.Datasets(i).Name, DatasetName)
            DatasetIndex = i;
            DatasetName  = finfo.GroupHierarchy.Datasets(i).Name;
        end
    end
    
    DatasetDims = finfo.GroupHierarchy.Datasets(DatasetIndex).Dims; 
    
end

























