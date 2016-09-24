%% general script settings
plot_data_sets = false;
plot_meshes = true;
plot_geometry = true;
plot_grounding_line = true;
plot_friction_coefficient = true;

%% add to path
addpath('../bin/');    % my scripts
addpath('../../bin');  % issm/trunk/bin
addpath('../../lib');  % issm/trunk/lib

% create output folders
if ~exist('models', 'file')
    mkdir models
end

if ~exist('figures', 'file')
    mkdir figures
end

%% Read input data

% read Bedmap2 data if not loaded
if ~exist('bm2', 'var')
    bm2 = read_bedmap2();
end

% read Rignot velocity if not loaded
if ~exist('rignot', 'var')
    rignot = read_rignot_velocity();
end

% visualize input data
if plot_data_sets
    close all
    plot_bed(bm2)
    plot_surface(bm2)
    plot_grounded(bm2)
    plot_velocity(rignot, 1)
end


%% Define model limits
% 1. Press 'add a contour (closed)'
% 2. Click around area of interest (no need to close polygon)
% 3. Press <Enter>
% 4. Close exptool dialog box by pressing 'Quit' button
% 5. Close figure
domain = 'DomainOutline.exp';
if ~exist(domain, 'file')
    if ~plot_data_sets
        plot_velocity(rignot, 1)
    end
    exptool(domain)
end


%% Mesh generation

% mesh parameters
hinit = 10000;   % element size for the initial mesh
hmax = 40000;    % maximum element size of the final mesh
hmin = 5000;     % minimum element size of the final mesh
gradation = 1.7; % maximum size ratio between two neighboring elements
err = 8;         % maximum error between interpolated and control field

% generate an initial uniform mesh (resolution = hinit meters)
md = bamg(model, 'domain', domain, 'hmax', hinit);
%plotmodel(md,'data','mesh')

% TODO: check if the velocity data below is oriented correctly
% interpolate velocities onto coarse mesh
vx_obs = InterpFromGridToMesh( ...
    rignot.x, rignot.y, ...
    rignot.vx, ...
    md.mesh.x, md.mesh.y, ...
    0);

vy_obs = InterpFromGridToMesh( ...
    rignot.x, rignot.y, ...
    rignot.vy, ...
    md.mesh.x, md.mesh.y, ...
    0);

vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);

% adapt the mesh to minimize error in velocity interpolation
md = bamg(md, 'hmax', hmax, 'hmin', hmin, ...
    'gradation', gradation, 'field', vel_obs, 'err', err);

%clear vx_obs vy_obs vel_obs;

if plot_meshes
    plotmodel(md, 'data', 'mesh')
    saveas(gcf, 'figures/siple_mesh')
    saveas(gcf, 'figures/siple_mesh.pdf')
end

% save model
save models/siple_mesh_generation md;


%% Apply masks for grounded/floating ice

md = loadmodel('models/siple_mesh_generation');

% interpolate onto our mesh vertices
groundedice = double(InterpFromGridToMesh(...
    bm2.x', bm2.y', bm2.grounded, ...
    md.mesh.x, md.mesh.y, 0));

% fill in the md.mask structure
% ice is grounded for mask equal one
md.mask.groundedice_levelset = groundedice;
clear groundedice

% interpolate onto our mesh vertices
ice = double(InterpFromGridToMesh(...
    bm2.x', bm2.y', bm2.grounded, ...
    md.mesh.x, md.mesh.y, 0));
ice(ice > 1.0e-5) = -1;   % make ice shelves count as ice

% ice is present when negative
%md.mask.ice_levelset = -1 * ones(md.mesh.numberofvertices, 1);% all is ice
md.mask.ice_levelset = ice;
clear ice

if plot_grounding_line
    plotmodel(md, ...
        'data', md.mask.groundedice_levelset, ...
        'title', 'grounded/floating', ...
        'data', md.mask.ice_levelset, ...
        'title', 'ice/no-ice');
    saveas(gcf, 'figures/siple_grounding_line')
    saveas(gcf, 'figures/siple_grounding_line.pdf')
end

% Save model
save models/siple_set_mask md;


%% Parameterization
md = loadmodel('models/siple_set_mask');
md = parameterize(md, 'siple_params.m');

% Use a MacAyeal flow model
md = setflowequation(md, 'SSA', 'all');

% Save model
save models/siple_parameterization md;


%% Find stress balance (control method)

md = loadmodel('models/siple_parameterization');

% Control general
md.inversion.iscontrol = 1;
md.inversion.maxsteps = 20;
md.inversion.maxiter = 40;
md.inversion.dxmin = 0.1;
md.inversion.gttol = 1.0e-4;
md.verbose=verbose('solution', true, 'control', true);

% Cost functions
md.inversion.cost_functions = [101 103 501];
md.inversion.cost_functions_coefficients = ...
    ones(md.mesh.numberofvertices, 3);
md.inversion.cost_functions_coefficients(:,1) = 1;
md.inversion.cost_functions_coefficients(:,2) = 1;
md.inversion.cost_functions_coefficients(:,3) = 8e-15;

% Controls
md.inversion.control_parameters = {'FrictionCoefficient'};
md.inversion.min_parameters = 1 * ones(md.mesh.numberofvertices, 1);
md.inversion.max_parameters = 200 * ones(md.mesh.numberofvertices, 1);

% Additional parameters
md.stressbalance.restol = 0.01;
md.stressbalance.reltol = 0.1;
md.stressbalance.abstol = NaN;

% Solve
md.toolkits = toolkits;
md.cluster = generic('name', oshostname, 'np', 2);
md = solve(md, StressbalanceSolutionEnum);

% Update model friction fields accordingly
md.friction.coefficient = ...
    md.results.StressbalanceSolution.FrictionCoefficient;

if plot_friction_coefficient
    plotmodel(md, 'data', md.friction.coefficient)
    saveas(gcf, 'figures/siple_friction')
    saveas(gcf, 'figures/siple_friction.pdf')
end

% Save model
save models/siple_control_drag md;

%% Post visualization
plotmodel(md, 'nlines', 2, 'ncols', 2, ...
    'unit#all', 'km', 'axis#all', 'equal', ...
    'xlim#all', [min(md.mesh.x) max(md.mesh.x)] / 10^3, ...
    'ylim#all', [min(md.mesh.y) max(md.mesh.y)] / 10^3, ...
    'FontSize#all', 12, ...
    'data', md.initialization.vel, ...
    'title', 'Observed velocity', ...
    'data', md.results.StressbalanceSolution.Vel, ...
    'title', 'Modeled Velocity', ...
    'data', md.geometry.base, ...
    'title', 'Bed elevation', ...
    'data', md.results.StressbalanceSolution.FrictionCoefficient, ...
    'title', 'Friction Coefficient', ...
    'colorbar#all', 'on', 'colorbartitle#1-2', '[m/yr]', ...
    'caxis#1-2', ([1.5,4000]), ...
    'colorbartitle#3', '[m]', 'log#1-2', 10);
saveas(gcf, 'figures/siple_combined')
saveas(gcf, 'figures/siple_combined.pdf')



