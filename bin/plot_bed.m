function plot_bed(dataset)

figure

imagesc(dataset.x, dataset.y, dataset.bed);
axis xy equal tight

title('Bed elevation')
xlabel('x [m]')
ylabel('y [m]')

cb = colorbar();
ylabel(cb,'z [m]')

if ~exist('figures', 'file')
    mkdir figures
end
saveas(gcf, 'figures/bed')
saveas(gcf, 'figures/bed.pdf')

end