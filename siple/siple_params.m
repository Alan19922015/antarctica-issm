% Parameters to change/Try
friction_coefficient = 10; % default [10]
Temp_change          =  0;  % default [0 K]

%Name and Coordinate system
md.miscellaneous.name = 'siple';
md.mesh.epsg = 3031;

%NetCdf Loading
disp('   Loading SeaRISE data from NetCDF');
ncdata = '../data/Antarctica_5km_withshelves_v0.75.nc';
x1    = ncread(ncdata, 'x1');
y1    = ncread(ncdata, 'y1');
usrf  = ncread(ncdata, 'usrf')';
topg  = ncread(ncdata, 'topg')';
temp  = ncread(ncdata, 'presartm')';
smb   = ncread(ncdata, 'presprcp')';
gflux = ncread(ncdata, 'bheatflx_fox')';

disp('   Loading velocities data from NetCDF');
nsidc_vel='../data/antarctica_ice_velocity_900m.nc';

xmin = ncreadatt(nsidc_vel, '/', 'xmin');
xmin = strtrim(xmin);  %this is a string, and we need to recover the double value
xmin = xmin(1:end-2);  %get rid of the unit
xmin = str2num(xmin);  %convert to double

ymax = ncreadatt(nsidc_vel, '/', 'ymax'); 
ymax = strtrim(ymax);  
ymax = ymax(1:end-2);  
ymax = str2num(ymax); 

nx = ncreadatt(nsidc_vel, '/', 'nx');
ny = ncreadatt(nsidc_vel, '/', 'ny');

spacing = ncreadatt(nsidc_vel, '/', 'spacing'); 
spacing = strtrim(spacing);  
spacing = spacing(1:end-2);  
spacing = str2num(spacing); 

velx = double(ncread(nsidc_vel, 'vx'));
vely = double(ncread(nsidc_vel, 'vy'));

x2 = xmin + (0:1:nx)' * spacing; 
x2 = double(x2);

y2 = (ymax - ny * spacing) + (0:1:ny)' * spacing; 
y2 = double(y2);

% Geometry
disp('   Interpolating surface and ice base');
md.geometry.base    = InterpFromGridToMesh( ...
    x1, y1, topg, ...
    md.mesh.x, md.mesh.y, 0);
md.geometry.surface = InterpFromGridToMesh( ...
    x1, y1, usrf, ...
    md.mesh.x, md.mesh.y, 0);
clear usrf topg;

disp('   Constructing thickness');
md.geometry.thickness = md.geometry.surface - md.geometry.base;

% ensure hydrostatic equilibrium on ice shelf: 
di = md.materials.rho_ice / md.materials.rho_water;

%Get the node numbers of floating nodes
pos = find(md.mask.groundedice_levelset < 0); 

%apply a flotation criterion on the precedingly defined nodes and
%redefine base and thickness accordingly
md.geometry.thickness(pos) = 1 / (1 - di) * md.geometry.surface(pos);
md.geometry.base(pos) = md.geometry.surface(pos) ...
    - md.geometry.thickness(pos);
md.geometry.hydrostatic_ratio = ones(md.mesh.numberofvertices, 1);

%Set min thickness to 1 meter
pos0=find(md.geometry.thickness<=0);
md.geometry.thickness(pos0) = 1;
md.geometry.surface = md.geometry.thickness+md.geometry.base;

%Initialization parameters
disp('   Interpolating temperatures');
md.initialization.temperature = InterpFromGridToMesh( ...
    x1, y1, ...
    temp, ...
    md.mesh.x, md.mesh.y, 0) + 273.15 + Temp_change;
clear temp;

disp('   Set observed velocities')
vx_obs=InterpFromGridToMesh( ...
    x2, y2, flipud(velx'), ...
    md.mesh.x, md.mesh.y, 0);
vy_obs=InterpFromGridToMesh( ...
    x2, y2, flipud(vely'), ...
    md.mesh.x, md.mesh.y, 0);
clear velx vely;

vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);
md.initialization.vx = vx_obs;
md.initialization.vy = vy_obs;
md.initialization.vz = zeros(md.mesh.numberofvertices, 1);
md.initialization.vel = vel_obs;

disp('   Set Pressure');
md.initialization.pressure = ...
    md.materials.rho_ice * md.constants.g * md.geometry.thickness;

disp('   Construct ice rheological properties');
md.materials.rheology_n = 3 * ones(md.mesh.numberofelements, 1);
md.materials.rheology_B = paterson(md.initialization.temperature);

%Forcings
disp('   Interpolating surface mass balance');
md.smb.mass_balance = InterpFromGridToMesh( ...
    x1, y1, smb, ...
    md.mesh.x, md.mesh.y, 0);
md.smb.mass_balance = ...
    md.smb.mass_balance * md.materials.rho_water / md.materials.rho_ice;
clear smb;

disp('   Set geothermal heat flux');
md.basalforcings.geothermalflux = InterpFromGridToMesh( ...
    x1, y1, gflux, ...
    md.mesh.x, md.mesh.y, 0);
clear gflux;

%Friction and inversion set up
disp('   Construct basal friction parameters');
md.friction.coefficient = ...
    friction_coefficient * ones(md.mesh.numberofvertices, 1);
md.friction.p = ones(md.mesh.numberofelements, 1);
md.friction.q = ones(md.mesh.numberofelements, 1);

%no friction applied on floating ice
pos = find(md.mask.groundedice_levelset < 0);
md.friction.coefficient(pos) = 0;

md.inversion = m1qn3inversion();
md.inversion.vx_obs = vx_obs;
md.inversion.vy_obs = vy_obs;
md.inversion.vel_obs = vel_obs;

disp('   Set boundary conditions');
md = SetMarineIceSheetBC(md);
md.basalforcings.floatingice_melting_rate = ...
    zeros(md.mesh.numberofvertices, 1);
md.basalforcings.groundedice_melting_rate = ...
    zeros(md.mesh.numberofvertices, 1);

%impose observed temperature on surface
md.thermal.spctemperature = [md.initialization.temperature; 1]; 
md.masstransport.spcthickness = NaN * ones(md.mesh.numberofvertices, 1);
