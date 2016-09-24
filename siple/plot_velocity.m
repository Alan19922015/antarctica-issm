function plot_velocity(dataset, log)

figure

if ~exist('log', 'var')
    imagesc(dataset.x, dataset.y, dataset.velocity_norm);
else
    imagesc(dataset.x, dataset.y, log10(dataset.velocity_norm + 1));
end

axis xy equal tight

title('Surface velocity')
xlabel('x [m]')
ylabel('y [m]')

cb = colorbar();
ylabel(cb,'log_1_0(v) [m/s]')

end