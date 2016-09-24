function plot_grounded(dataset)

figure

imagesc(dataset.x, dataset.y, dataset.grounded);
axis xy equal tight

title('Grounded ice')
xlabel('x [m]')
ylabel('y [m]')

cb = colorbar();
ylabel(cb,'Flag')

end