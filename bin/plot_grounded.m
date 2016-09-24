function plot_grounded(dataset)

figure

imagesc(dataset.x, dataset.y, dataset.grounded);
axis xy equal tight

title('Grounded ice')
xlabel('x [m]')
ylabel('y [m]')

cb = colorbar();
ylabel(cb,'Flag')

if ~exist('figures', 'file')
    mkdir figures
end
saveas(gcf, 'figures/grounded')
saveas(gcf, 'figures/grounded.pdf')

end