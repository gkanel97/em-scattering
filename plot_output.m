function plot_output (filePath1, filePath2)

    max_err = importdata(filePath1);
    rcs = importdata(filePath2);

    radius = linspace(0.8, 1.0, length(rcs));
    figure();
    plot(radius, rcs, "LineWidth", 2);
    grid();
    title("RCS against cylinder radius", 'FontSize', 14);
    xlabel("Cylinder radius", "FontSize", 12);
    ylabel("RCS", "FontSize", 12);

    radius = linspace(0.01, 0.99, length(max_err));
    figure();
    plot(radius, max_err, "LineWidth", 2);
    grid();
    title("Maximum absolute error against surface radius", "FontSize", 13);
    xlabel("Auxiliary source radius", "FontSize", 12);
    ylabel("Maximum error", "FontSize", 12);
end
