function [ channel_idx,names ] = determine_channels(nexFile,field_name,str_pattern)
%determine_channels
%   Determines the indices of the struct 'nexFile', within the desired field
% (given by 'field_name') that match the string 'str_pattern'

field2look = nexFile.(field_name);

channel_idx = [];
names = {};
name_iter = 1;

for i = 1:length(field2look)
    if ~isempty(strfind(field2look{i}.name,str_pattern))
        channel_idx = [channel_idx,i];
        names{name_iter} = field2look{i}.name;
        name_iter = name_iter + 1;
    end
end
       

end

