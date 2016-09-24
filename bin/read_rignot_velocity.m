function [ rignot ] = read_rignot_velocity()
% Load Velocities
% http://nsidc.org/data/nsidc-0484.html
rignot.inputfile='../data/antarctica_ice_velocity_900m.nc'; 	

% Get necessary data to build up the velocity grid
rignot.xmin    = ncreadatt(rignot.inputfile,'/','xmin');
rignot.ymax    = ncreadatt(rignot.inputfile,'/','ymax');
rignot.spacing = ncreadatt(rignot.inputfile,'/','spacing');

rignot.nx = double(ncreadatt(rignot.inputfile,'/','nx'));
rignot.ny = double(ncreadatt(rignot.inputfile,'/','ny'));
rignot.vx = double(ncread(rignot.inputfile,'vx'));
rignot.vy = double(ncread(rignot.inputfile,'vy'));

rignot.velocity_norm = sqrt(rignot.vx.^2 + rignot.vy.^2);

rignot.xmin = strtrim(rignot.xmin);  % this is a string, and we need to recover the double value
rignot.xmin = str2num(rignot.xmin(1:end-2));  % get rid of the unit and convert to double

rignot.ymax = strtrim(rignot.ymax);  
rignot.ymax = str2num(rignot.ymax(1:end-2));  

rignot.spacing = strtrim(rignot.spacing);
rignot.spacing = str2num(rignot.spacing(1:end-2));  

% Build the coordinates
rignot.x = rignot.xmin+(0:1:rignot.nx)' * rignot.spacing;
rignot.y = (rignot.ymax - rignot.ny*rignot.spacing) ...
    + (0:1:rignot.ny)' * rignot.spacing;

%rignot.vx = flipud(rignot.vx);
%rignot.vy = flipud(rignot.vy);

%Limit the region to Pine Island
%posx  = find(x<=-12.0e5 & x>=-18.0e5);
%x_pig = x(posx);
%posy  = find(y<=1.0e5 & y>-4.0e5);
%y_pig = flipud(y(posy));

%vx_pig  = flipud(vx(posx,posy)');
%vy_pig  = flipud(vy(posx,posy)');
%vel_pig = sqrt(vx_pig.^2+vy_pig.^2);

%imagesc(x_pig,y_pig,log(vel_pig+1));
%axis xy equal tight;

end