function [ params ] = read_cruna_params( fname )

[~,~,ext] = fileparts(fname);

comment_symbols = {'%','#'};


% for hdf5 file
if strcmp(ext,'.h5')
    try
        params_raw = hdf5read(fname,'params');
    catch
        params = [];
        return
    end
    
else % assume ascii file (preprocess)
    fid = fopen(fname);
    p_i = 0;
    while ~feof(fid)
        row=fgetl(fid);
        if isempty(row)
        elseif any(strcmp(row(1),comment_symbols)) % ignore comments
        else
            row=row(find(~isspace(row))); % rm whitespaces
            if any(strcmp(row(1),comment_symbols)) % ignore comments (after whitespaces)
            else
                p_i = p_i+1;
                %disp(strtok(row,comment_symbols))
                params_raw(p_i).Data=strtok(row,comment_symbols); % save part before comment
            end
        end
    end
    fclose(fid);
    
end

keys=cell(0,0);
for p_i = 1:length(params_raw)
    
    % split row
    key_value_pair = strsplit(params_raw(p_i).Data,'=');
    
    if isempty(key_value_pair{1})
        continue
    end
    
    % ignore comments
    if any(strcmp(key_value_pair{1}(1),comment_symbols))
        key_value_pair{1} = [];
        continue
    end
    
    % parameter name
    keys{p_i} = key_value_pair{1};
    
    % set value
    if length(key_value_pair) > 1
        value = strtrim(strsplit(key_value_pair{2},';','CollapseDelimiters',true));  
        double_value = str2double(value);
        
        for i = 1:length(value)
            if isnan(double_value(i)) % is char
                if ~isempty(value{i}) % no emty char
                    vals{p_i,i} = value{i};
                end
            else % save as number
                vals{p_i,i} = double_value(i);
            end
        end
    end
end

% transfer to struct
for p_i = 1:length(keys)
    if isempty(keys{p_i})
        continue
    end
    %p_i
    %if p_i==24; disp(keys{p_i}); keyboard; end
    key = strsplit(keys{p_i},'.');
    level_of_para_structure = length(key);
    %vals{p_i,1}
    if ischar(vals{p_i,1}) % keep char type
        char_or_array = vals(p_i,:);
        emptyCells                = cellfun(@isempty,char_or_array);
        char_or_array(emptyCells) = []; % remove empty cells
    else % convert to double array
        char_or_array = cell2mat(vals(p_i,:));
    end
    if isempty(key{1}) % empty row
        continue
    elseif strcmp(key{1}(1),comment_symbols) % comment row
        continue
    elseif level_of_para_structure == 1
        params.(key{1}) = char_or_array;      
    elseif level_of_para_structure == 2
        params.(key{1}).(key{2}) = char_or_array;
    elseif level_of_para_structure == 3
        params.(key{1}).(key{2}).(key{3}) = char_or_array;
    elseif level_of_para_structure > 3
        warning('parameter structure has more than 3 levels, not expected & ignored')
        keyboard
    end
end

end

