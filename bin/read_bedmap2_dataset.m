function [ grid ] = read_bedmap2_dataset(dataset)
%READ_BEDMAP2_DATA Summary of this function goes here
%   Examples:
%   bed = read_bedmap2_dataset('bed');
%   surface = read_bedmap2_dataset('surface');
%   thickness = read_bedmap2_dataset('thickness');

fid = fopen([get_bedmap2_path(), '/bedmap2_', dataset, '.flt'],'r','l');
grid = fread(fid,[6667,6667],'float32');
fclose(fid);

end