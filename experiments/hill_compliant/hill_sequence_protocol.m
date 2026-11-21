function hill_sequence_protocol(ax, t, labels)
    % HILL_SEQUENCE_PROTOCOL Shade the original waits and label load stages.
    % Axis padding prevents protocol labels from covering force peaks.
    limits = ylim(ax);
    if labels
        traces = findobj(ax, 'Type', 'line');
        peak = max(arrayfun(@(h) max(h.YData(:)), traces));
        limits(2) = max(limits(2), peak + .15 * diff(limits));
    end
    waiting = logical(t.waiting);
    edges = diff([false; waiting; false]);
    starts = find(edges == 1);
    ends = find(edges == -1) - 1;
    for k = 1:numel(starts)
        left = t.t(starts(k));
        right = t.t(ends(k));
        h = patch(ax, [left right right left], limits([1 1 2 2]), [.965 .970 .975], ...
                  'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(h, 'bottom');
    end
    xlim(ax, [t.t(1), t.t(end)]);
    ylim(ax, limits);
    if labels
        for stage = 1:3
            s = t.segments(stage);
            text(ax, mean([s.startTime, s.endTime]), limits(2), sprintf('%g kg | speed code %g', s.loadKg, s.speedCode), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 12);
        end
    end
end
