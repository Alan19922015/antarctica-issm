function [ rignot ] = read_rignot_velocity()
% Load Velocities
% http://nsidc.org/data/nsidc-0484.html
rignot.inputfile = '../data/antarctica_ice_velocity_900m.nc';
%rignot.inputfile = '../data/antarctica_ice_velocity_450m.nc';

% Get necessary data to build up the velocity grid
rignot.xmin    = ncreadatt(rignot.inputfile, '/', 'xmin');
rignot.ymax    = ncreadatt(rignot.inputfile, '/', 'ymax');
rignot.spacing = ncreadatt(rignot.inputfile, '/', 'spacing');

rignot.nx = double(ncreadatt(rignot.inputfile, '/', 'nx'));
rignot.ny = double(ncreadatt(rignot.inputfile, '/', 'ny'));
rignot.vx = double(ncread(rignot.inputfile, 'vx'));
rignot.vy = double(ncread(rignot.inputfile, 'vy'));

rignot.velocity_norm = sqrt(rignot.vx.^2 + rignot.vy.^2);

rignot.xmin = strtrim(rignot.xmin);
rignot.xmin = str2num(rignot.xmin(1:end-2));

rignot.ymax = strtrim(rignot.ymax);  
rignot.ymax = str2num(rignot.ymax(1:end-2));  

rignot.spacing = strtrim(rignot.spacing);
rignot.spacing = str2num(rignot.spacing(1:end-2));  

% Build the coordinates
rignot.x = rignot.xmin + (0:1:rignot.nx)' * rignot.spacing;
rignot.y = (rignot.ymax - rignot.ny*rignot.spacing) ...
    + (0:1:rignot.ny)' * rignot.spacing;

end
