function [ master_dir,all_files ] = choose_files( fdir )
%CHOOSE_FILES Function that instigates loop for choosing files for
%analysis, compiling in a cell array of filenames (including their path
%names)
%   Detailed explanation goes here

if nargin < 1
    
    fprintf('Please choose master directory that all data files are found within\n');
    master_dir = uigetdir();
    
else
    
    master_dir = fdir;
    
end

terminate_flag = false;

all_files = {};

counter = 1;

while ~terminate_flag
    
    fprintf('Please navigate and select desired file or files\n');
    [fnam,sub_dir] = uigetfile(fullfile(master_dir,'*.nex'),sprintf('Choose one of more files'),'MultiSelect','On');
    sub_dir = strrep(sub_dir,master_dir,'');

    if iscell(fnam) % in case when multiple files have been selected
        for f_i = 1:length(fnam)
            all_files{counter} = fullfile(sub_dir,fnam{f_i});
            counter = counter + 1;
        end
    else  % when one file has been selected
        all_files{counter} = fullfile(sub_dir,fnam);
    end
     
    continue_string = input('Would you like to add more files? (Press n to end selection, press any key to continue adding files)\n','s');
    
    if strcmp(continue_string,'n')
        terminate_flag = true;
    end 

end

