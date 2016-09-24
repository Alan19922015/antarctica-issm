function bm2 = read_bedmap2()
    if ~exist('bm2.surface', 'var')
        bm2.surface = read_bedmap2_dataset('surface');
    end
    if ~exist('bm2.bed', 'var')
        bm2.bed = read_bedmap2_dataset('bed');
    end
    if ~exist('bm2.x', 'var')
        [bm2.x, bm2.y] = read_bedmap2_dataset('x,y');
    end
    if ~exist('bm2.X', 'var')
        [bm2.X, bm2.Y] = read_bedmap2_dataset('X,Y');
    end
    if ~exist('bm2.grounded', 'var')
        % original: -9999 == ocean; 0 == grounded; 1 == floating
        bm2.grounded = read_bedmap2_dataset('icemask_grounded_and_shelves');

        % new: 0 == ocean; 1 == grounded; -1 == floating
        bm2.grounded(bm2.grounded == 1) = -1;
        bm2.grounded(bm2.grounded == 0) = 1;
        bm2.grounded(bm2.grounded == -9999) = 0;
    end
    %bm2.velocity = read_bedmap2_dataset('velocity'); % Rignot velocity
end

function [ grid, grid2 ] = read_bedmap2_dataset(dataset)
%READ_BEDMAP2_DATA Summary of this function goes here
%   Examples:
%   bed = read_bedmap2_dataset('bed');
%   surface = read_bedmap2_dataset('surface');
%   thickness = read_bedmap2_dataset('thickness');

if strcmp(dataset, 'velocity')
    grid = read_bedmap2_grid('../data/rignot_velocity_bedmap2_grid.flt');
elseif strcmp(dataset, 'x,y')
    grid = -3333:3333;
    grid2 = -3333:3333;
    grid = grid .* 1000; % convert from km to m
    grid2 = grid2 .* 1000; % convert from km to m
elseif strcmp(dataset, 'X,Y')
    [grid, grid2] = meshgrid(-3333:3333, -3333:3333);
    grid = grid .* 1000; % convert from km to m
    grid2 = grid2 .* 1000; % convert from km to m
else
    grid = read_bedmap2_grid(['../data/bedmap2_bin/bedmap2_', ...
        dataset, '.flt']);
end

end

function grid = read_bedmap2_grid(path)
fid = fopen(path, 'r', 'l');
grid = fread(fid, [6667,6667], 'float32');
fclose(fid);
end