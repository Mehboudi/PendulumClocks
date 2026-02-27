% Script to overlay Hist_filtered_M2.fig and Hist_unfiltered_M2.fig
% The filtered histogram is placed on top with alpha = 0.6

% Load the unfiltered histogram figure
fig1 = openfig('Figs/Hist_unfiltered_M2.fig', 'invisible');
ax1 = gca;
h_unfiltered = findobj(ax1, 'Type', 'histogram');
if isempty(h_unfiltered)
    h_unfiltered = findobj(ax1, 'Type', 'bar');
end
if isempty(h_unfiltered)
    h_unfiltered = findobj(ax1, 'Type', 'patch');
end

% Load the filtered histogram figure
fig2 = openfig('Figs/Hist_filtered_M2.fig', 'invisible');
ax2 = gca;
h_filtered = findobj(ax2, 'Type', 'histogram');
if isempty(h_filtered)
    h_filtered = findobj(ax2, 'Type', 'bar');
end
if isempty(h_filtered)
    h_filtered = findobj(ax2, 'Type', 'patch');
end

% Create new figure for overlay
fig3 = figure;
ax3 = axes(fig3);
hold on;

% Copy unfiltered histogram first (bottom layer)
if ~isempty(h_unfiltered)
    h_unfilt_copy = copyobj(h_unfiltered, ax3);
    % Set properties for unfiltered (opaque)
    if isprop(h_unfilt_copy, 'FaceAlpha')
        set(h_unfilt_copy, 'FaceAlpha', 1.0);
    end
    if isprop(h_unfilt_copy, 'DisplayName')
        set(h_unfilt_copy, 'DisplayName', 'Unfiltered');
    end
end

% Copy filtered histogram on top with transparency
if ~isempty(h_filtered)
    h_filt_copy = copyobj(h_filtered, ax3);
    % Set transparency for filtered
    if isprop(h_filt_copy, 'FaceAlpha')
        set(h_filt_copy, 'FaceAlpha', 0.6);
    end
    if isprop(h_filt_copy, 'DisplayName')
        set(h_filt_copy, 'DisplayName', 'Filtered');
    end
end

% Add manual xlines at x=0, x=1, x=2 with colors matching the image
xline(ax3, 0, 'Color', [0.5, 0.6, 0.7], 'LineWidth', 1.5);
xline(ax3, 1, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5);
xline(ax3, 2, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5);

% Copy axis labels with all properties including interpreter
xlabel(ax3, ax1.XLabel.String, 'Interpreter', ax1.XLabel.Interpreter, ...
    'FontSize', ax1.XLabel.FontSize, 'FontName', ax1.XLabel.FontName, ...
    'FontWeight', ax1.XLabel.FontWeight, 'FontAngle', ax1.XLabel.FontAngle, ...
    'Color', ax1.XLabel.Color);
ylabel(ax3, ax1.YLabel.String, 'Interpreter', ax1.YLabel.Interpreter, ...
    'FontSize', ax1.YLabel.FontSize, 'FontName', ax1.YLabel.FontName, ...
    'FontWeight', ax1.YLabel.FontWeight, 'FontAngle', ax1.YLabel.FontAngle, ...
    'Color', ax1.YLabel.Color);

% Copy axis appearance properties
ax3.FontSize = ax1.FontSize;
ax3.FontName = ax1.FontName;
ax3.FontWeight = ax1.FontWeight;
ax3.FontAngle = ax1.FontAngle;
ax3.XLim = ax1.XLim;
ax3.YLim = ax1.YLim;
ax3.XScale = ax1.XScale;
ax3.YScale = ax1.YScale;
ax3.XTick = ax1.XTick;
ax3.YTick = ax1.YTick;
ax3.XTickLabel = ax1.XTickLabel;
ax3.YTickLabel = ax1.YTickLabel;
ax3.XGrid = ax1.XGrid;
ax3.YGrid = ax1.YGrid;
ax3.GridLineStyle = ax1.GridLineStyle;
ax3.GridAlpha = ax1.GridAlpha;
ax3.LineWidth = ax1.LineWidth;
ax3.Box = ax1.Box;
ax3.XColor = ax1.XColor;
ax3.YColor = ax1.YColor;
ax3.Color = ax1.Color;

% Update title to indicate overlay
if ~isempty(ax1.Title.String)
    title(ax3, 'Overlay: Filtered and Unfiltered Histograms', ...
        'Interpreter', ax1.Title.Interpreter, ...
        'FontSize', ax1.Title.FontSize, 'FontName', ax1.Title.FontName, ...
        'FontWeight', ax1.Title.FontWeight, 'FontAngle', ax1.Title.FontAngle, ...
        'Color', ax1.Title.Color);
end

% Add legend
legend(ax3, 'show');

% Save the overlay figure
savefig(fig3, 'Figs/Hist_overlay_M2.fig');
saveas(fig3, 'Figs/Hist_overlay_M2.png');

% Close invisible figures
close(fig1);
close(fig2);

disp('Overlay histogram created and saved to Figs/Hist_overlay_M2.fig');
