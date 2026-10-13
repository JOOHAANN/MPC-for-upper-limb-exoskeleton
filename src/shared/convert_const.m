clear all;
load ('regression_const');
data = zeros(15, 24);

for mode = 1:2
    if mode == 1
        varMode = 'elevation';
    else
        varMode = 'speed';
    end

    for comp = 1:4
        varComp = strcat('Comp', int2str(comp));
        for musc = 1:15
            varMusc = strcat('muscle_', int2str(musc));
            start = (comp - 1) * 3 + 1 + (mode - 1) * 12;
            data(musc, start:start + 2) = regression_const.data.(varComp).table.(varMusc).(varMode);
        end
    end
end

csvwrite('regression.txt', data);

load('weights_generic');
csvwrite('weightings.txt', weights_generic);

load('gaussXP_new', 'XP_FITN'); % // with the gaussians
XP_FIT = XP_FITN;
XP_FIT(3, :) = ((XP_FIT(3, :) - min(XP_FIT(3, :))) / (max(XP_FIT(3, :) - min(XP_FIT(3, :)))));
XP_data(1, :) = XP_FIT(1, :) / max(XP_FIT(1, :));
XP_data(2, :) = XP_FIT(2, :) / max(XP_FIT(2, :));
XP_data(3, :) = XP_FIT(3, :) / max(XP_FIT(3, :));
XP_data(4, :) = XP_FIT(4, :) / max(XP_FIT(4, :));

csvwrite('gaussXP.txt', XP_data);
