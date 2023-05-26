function h_fig = get_figure_h(position)

h_fig = findobj('Position',position);
if (isempty(h_fig))
    h_fig = figure;
    set(h_fig,'Position',position)
end