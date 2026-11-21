function diagnostics = inspect_bse_model(subjectNames, showFigures, modelVersion)
%INSPECT_BSE_MODEL Expose muscle recruitment, torque and elastic soft limits.
%   inspect_bse_model                  % All four subjects
%   inspect_bse_model("725")
%   diagnostics = inspect_bse_model([],false);
%
% Uses the same selected model version as run_bse_demo. No fitting, reference smoothing,
% trajectory clipping or changes to the original data occur here.
if nargin < 1, subjectNames = []; end
if nargin < 2, showFigures = true; end
base = fileparts(mfilename('fullpath'));
if nargin < 3
    modelVersion = 'current';
end
folder = fullfile(base,'diagnostics');
if ~isfolder(folder), mkdir(folder); end
useNative = build_bse_solver;
subjects = bse_subjects(subjectNames,modelVersion);
rows = struct([]);
diagnostics.subjects = cell(size(subjects));
diagnostics.modelVersion = string(modelVersion);

for k = 1:numel(subjects)
    subject = subjects{k};
    cases = calculate_cases(subject,useNative);
    for j = 1:numel(cases)
        item = cases{j};
        for source = 1:2
            out = item.outputs(:,:,source);
            plant = bse_model(subject.models{source},subject.balanceController, ...
                subject.sequences{1}.standingAngle,useNative);
            [lower,upper] = soft_limits(plant.physics,out(:,1));
            % The plotted split must reproduce the torque used by the plant.
            assert(max(abs(lower+upper-out(:,26))) < 1e-9,'Boundary accounting failed.');
            assert(max(abs(out(:,3)-out(:,4)+out(:,5))) < 1e-9,'Muscle accounting failed.');
            valid = isfinite(item.referenceTorque);
            inputCorrelation = corrcoef(item.inputs(valid,1,source),item.referenceTorque(valid));
            torqueCorrelation = corrcoef(out(valid,3),item.referenceTorque(valid));
            rows = [rows; struct('Subject',string(subject.id),'Source',source, ...
                'LoadKg',item.loadKg,'SpeedCode',item.speedCode, ...
                'BackPeakNm',max(out(:,4)),'AbdomenPeakNm',max(out(:,5)), ...
                'MeanAbdomenToMeanBackRatio',mean(out(:,5))/mean(out(:,4)), ...
                'LowerBoundaryPeakNm',max(lower),'UpperBoundaryPeakMagnitudeNm',max(-upper), ...
                'MinimumThetaRad',min(out(:,1)),'MaximumThetaRad',max(out(:,1)), ...
                'ReferenceSupportFraction',mean(valid), ...
                'InputReferenceCorrelation',inputCorrelation(1,2), ...
                'NetTorqueReferenceCorrelation',torqueCorrelation(1,2))]; %#ok<AGROW>
        end
    end
    plot_muscles(subject,cases,folder,showFigures);
    plot_boundaries(subject,cases,folder,showFigures);
    plot_boundary_traces(subject,cases,folder,showFigures);
    plot_signal_path(subject,cases,folder,showFigures);
    plot_torque_budget(subject,cases,folder,showFigures);
    diagnostics.subjects{k} = struct('id',subject.id,'models',{subject.models},'cases',{cases});
    fprintf('Inspected %s: nine cases, two fixed model vectors.\n',subject.id);
end
plot_structure(folder,showFigures);
diagnostics.metrics = struct2table(rows);
writetable(diagnostics.metrics,fullfile(folder,'muscle_boundary_metrics.csv'));
save(fullfile(folder,'bse_model_diagnostics.mat'),'diagnostics','-v7');
fprintf('BSE diagnostic figures and traces: %s\n',folder);
end

function cases = calculate_cases(subject,useNative)
cases = {};
for group = 1:numel(subject.sequences)
    sequence = subject.sequences{group};
    traces = cell(1,2); inputs = cell(1,2);
    for source = 1:2
        model = subject.models{source};
        inputs{source} = bse_inputs(sequence,model,source);
        plant = bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
        traces{source} = plant.simulate(inputs{source},sequence.mass,sequence.waiting, ...
            sequence.initialState,sequence.dt,10);
    end
    for stage = 1:numel(sequence.segments)
        part = sequence.segments{stage}; index = part.startIndex:part.endIndex;
        item = struct('loadKg',part.loadKg,'speedCode',part.speedCode, ...
            'time',sequence.time(index)-sequence.time(index(1)), ...
            'referenceTheta',sequence.referenceTheta(index), ...
            'referenceTorque',sequence.referenceTorque(index), ...
            'inputs',cat(3,inputs{1}(index,:),inputs{2}(index,:)), ...
            'outputs',cat(3,traces{1}(index,:),traces{2}(index,:)));
        if isfield(part,'cadenceBPM'), item.cadenceBPM = part.cadenceBPM; end
        item.referenceTorque(~sequence.torqueValid(index)) = NaN;
        cases{end+1} = item; %#ok<AGROW>
    end
end
end

function plot_muscles(subject,cases,folder,showFigures)
[fig,layout] = new_figure(subject.id+" | Muscle torques",3,3,showFigures,[1500 960]);
backColor = [0 .447 .741]; abdomenColor = [.85 .35 .12];
for k = 1:numel(cases)
    item = cases{k}; ax = nexttile(layout); hold(ax,'on');
    plot(ax,item.time,item.outputs(:,4,1),'Color',backColor,'LineWidth',1.5);
    plot(ax,item.time,item.outputs(:,4,2),'--','Color',backColor,'LineWidth',1.5);
    plot(ax,item.time,item.outputs(:,5,1),'Color',abdomenColor,'LineWidth',1.5);
    plot(ax,item.time,item.outputs(:,5,2),'--','Color',abdomenColor,'LineWidth',1.5);
    title(ax,condition_title(item,k));
    ylabel(ax,'Muscle torque magnitude (N m)'); xlabel(ax,'Action time (s)');
    grid(ax,'on'); box(ax,'on'); xlim(ax,item.time([1 end]));
end
key = legend(ax,{'Back: experimental input','Back: MAP input', ...
    'Abdomen: experimental input','Abdomen: MAP input'},'Orientation','horizontal');
key.Layout.Tile = 'south';
title(layout,subject.id+" | Back and abdominal torque in all nine conditions",'Interpreter','none');
subtitle(layout,'Positive magnitudes are shown; abdominal torque enters the joint balance with a minus sign.');
save_figure(fig,folder,subject.shortName+"_muscle_torques",showFigures);
end

function plot_boundaries(subject,cases,folder,showFigures)
[fig,layout] = new_figure(subject.id+" | Soft boundaries",2,2,showFigures,[1400 900]);
sourceNames = ["Experimental input","MAP input"];
example = find_case(cases,5,45);
for source = 1:2
    plant = bse_model(subject.models{source},subject.balanceController,1.5);
    h = plant.physics;
    q = linspace(h.boundaryLower-.45,h.boundaryUpper+.45,600)';
    [lower,upper] = soft_limits(h,q);
    allTheta = cellfun(@(c) c.outputs(:,1,source),cases,'UniformOutput',false);
    allTheta = vertcat(allTheta{:});
    ax = nexttile(layout,source); hold(ax,'on');
    patch(ax,[min(allTheta) max(allTheta) max(allTheta) min(allTheta)], ...
        [-h.boundaryCapacity -h.boundaryCapacity h.boundaryCapacity h.boundaryCapacity], ...
        [.85 .93 .85],'EdgeColor','none','DisplayName','Angles visited in nine cases');
    plot(ax,q,lower,'Color',[0 .447 .741],'LineWidth',1.8,'DisplayName','Lower boundary');
    plot(ax,q,upper,'Color',[.85 .35 .12],'LineWidth',1.8,'DisplayName','Upper boundary');
    xline(ax,h.boundaryLower,':','HandleVisibility','off');
    xline(ax,h.boundaryUpper,':','HandleVisibility','off');
    title(ax,sourceNames(source)+" | elastic characteristic");
    ylabel(ax,'Signed boundary torque (N m)'); xlabel(ax,'Theta (rad)');
    grid(ax,'on'); xlim(ax,q([1 end])); legend(ax,'Location','southwest');
    subtitle(ax,sprintf('Dead zone %.2f-%.2f rad; stiffness %.0f N m/rad; capacity %.0f N m', ...
        h.boundaryLower,h.boundaryUpper,h.boundaryStiffness,h.boundaryCapacity));
    ax = nexttile(layout,source+2); hold(ax,'on');
    [lower,upper] = soft_limits(h,example.outputs(:,1,source));
    plot(ax,example.time,lower,'Color',[0 .447 .741],'LineWidth',1.8);
    plot(ax,example.time,upper,'--','Color',[.85 .35 .12],'LineWidth',1.8);
    actual = cellfun(@(c) max(abs(c.outputs(:,26,source))),cases);
    title(ax,sprintf('Actual %g kg trace | peak over all nine cases = %.6g N m',example.loadKg,max(actual)));
    subtitle(ax,sprintf('All-nine visited theta: %.3f to %.3f rad. No boundary damping.',min(allTheta),max(allTheta)));
    ylabel(ax,'Signed boundary torque (N m)'); xlabel(ax,'Action time (s)');
    extent = max(1,1.15*max(actual));
    grid(ax,'on'); ylim(ax,[-extent extent]); xlim(ax,example.time([1 end]));
    legend(ax,{'Lower boundary','Upper boundary'},'Location','best');
end
title(layout,subject.id+" | Capacity is not the torque actually applied",'Interpreter','none');
save_figure(fig,folder,subject.shortName+"_soft_boundaries",showFigures);
end

function plot_boundary_traces(subject,cases,folder,showFigures)
[fig,layout] = new_figure(subject.id+" | Actual boundary torques",3,3,showFigures,[1500 960]);
styles = {'-','--'}; colors = [0 .447 .741; .85 .35 .12];
for k = 1:numel(cases)
    item = cases{k}; ax = nexttile(layout); hold(ax,'on');
    for source = 1:2
        plant = bse_model(subject.models{source},subject.balanceController,1.5);
        [lower,upper] = soft_limits(plant.physics,item.outputs(:,1,source));
        plot(ax,item.time,lower,styles{source},'Color',colors(1,:),'LineWidth',1.5);
        plot(ax,item.time,upper,styles{source},'Color',colors(2,:),'LineWidth',1.5);
    end
    yline(ax,0,':','Color',[.5 .5 .5],'HandleVisibility','off');
    title(ax,condition_title(item,k));
    xlabel(ax,'Action time (s)'); ylabel(ax,'Signed soft-limit torque (N m)');
    grid(ax,'on'); box(ax,'on'); xlim(ax,item.time([1 end]));
end
key = legend(ax,{'Lower: experimental input','Upper: experimental input', ...
    'Lower: MAP input','Upper: MAP input'},'Orientation','horizontal');
key.Layout.Tile = 'south';
title(layout,subject.id+" | Actual upper/lower torques in all nine cases",'Interpreter','none');
save_figure(fig,folder,subject.shortName+"_boundary_traces",showFigures);
end

function plot_torque_budget(subject,cases,folder,showFigures)
% Show both muscles and the actual soft-limit contribution on the same figure.
item = find_case(cases,15,45);
[fig,layout] = new_figure(subject.id+" | Torque budget",2,2,showFigures,[1450 850]);
sourceNames = ["Experimental input","MAP input"];
for source = 1:2
    out = item.outputs(:,:,source);
    ax = nexttile(layout,source); hold(ax,'on');
    plot(ax,item.time,item.referenceTorque,'Color',[.75 .15 .15],'LineWidth',1.5);
    plot(ax,item.time,out(:,4),'Color',[0 .447 .741],'LineWidth',1.6);
    plot(ax,item.time,out(:,5),'Color',[.85 .35 .12],'LineWidth',1.6);
    plot(ax,item.time,out(:,3),'Color',[.49 .18 .56],'LineWidth',1.6);
    title(ax,sourceNames(source)); ylabel(ax,'Muscle torque (N m)');
    xlabel(ax,'Action time (s)'); grid(ax,'on'); xlim(ax,item.time([1 end]));
    legend(ax,{'Available CEINMS net estimate','Back magnitude','Abdomen magnitude','Net = back - abdomen'},'Location','best');
    plant = bse_model(subject.models{source},subject.balanceController,1.5);
    [lower,upper] = soft_limits(plant.physics,out(:,1));
    ax = nexttile(layout,source+2); hold(ax,'on');
    plot(ax,item.time,lower,'Color',[0 .447 .741],'LineWidth',1.6);
    plot(ax,item.time,upper,'--','Color',[.85 .35 .12],'LineWidth',1.6);
    title(ax,sprintf('Actual soft-limit torque | peak %.3f N m',max(abs(lower+upper))));
    ylabel(ax,'Signed limit torque (N m)'); xlabel(ax,'Action time (s)');
    grid(ax,'on'); xlim(ax,item.time([1 end]));
    extent = max(1,1.15*max(abs([lower;upper]))); ylim(ax,[-extent extent]);
    legend(ax,{'Lower boundary','Upper boundary'},'Location','best');
end
title(layout,subject.id+" | "+subject.modelVersion+" | "+condition_title(item,1),'Interpreter','none');
subtitle(layout,'Abdominal torque is a positive magnitude here and enters joint balance with a minus sign.');
save_figure(fig,folder,subject.shortName+"_torque_budget",showFigures);
end

function plot_signal_path(subject,cases,folder,showFigures)
item = find_case(cases,5,45);
[fig,layout] = new_figure(subject.id+" | Input to torque",4,2,showFigures,[1450 1160]);
sourceNames = ["Experimental input","Trained MAP input"];
colors = [.12 .53 .39; .49 .18 .56; .89 .43 .1];
for source = 1:2
    out = item.outputs(:,:,source); input = item.inputs(:,:,source);
    model = subject.models{source}; r = model.neuralParameters; p = model.parameters;
    ax = nexttile(layout,source); hold(ax,'on');
    plot(ax,item.time,input(:,1),'Color',colors(1,:),'LineWidth',1.5);
    plot(ax,item.time,out(:,18),'Color',colors(2,:),'LineWidth',1.4);
    plot(ax,item.time,out(:,14),'Color',colors(3,:),'LineWidth',1.4);
    title(ax,sourceNames(source)+" | back recruitment");
    subtitle(ax,sprintf('Feedback blend %.3f; causal lead %.3f s; exponent %.3f',r(7),r(9),p(6)));
    ylabel(ax,'Back activation / command');
    legend(ax,{'Prepared input','Neural command','Effective Hill activation'},'Location','best');
    ax = nexttile(layout,source+2); hold(ax,'on');
    plot(ax,item.time,input(:,2),'Color',colors(1,:),'LineWidth',1.5);
    plot(ax,item.time,out(:,19),'Color',colors(2,:),'LineWidth',1.4);
    plot(ax,item.time,out(:,15),'Color',colors(3,:),'LineWidth',1.4);
    title(ax,sprintf('Abdominal recruitment | exponent %.3f',p(10)));
    ylabel(ax,'Abdominal activation / command');
    legend(ax,{'Recruited input','Neural command','Effective Hill activation'},'Location','best');
    ax = nexttile(layout,source+4); hold(ax,'on');
    plot(ax,item.time,item.referenceTheta,'--','Color',[.5 .5 .5],'LineWidth',1.4);
    plot(ax,item.time,out(:,24),'Color',colors(2,:),'LineWidth',1.4);
    plot(ax,item.time,out(:,1),'Color',[0 .447 .741],'LineWidth',1.5);
    title(ax,'Angle from the same Hill / joint solve'); ylabel(ax,'Theta (rad)');
    legend(ax,{'Original IK reference','Decoded motion intent','BSE theta'},'Location','best');
    ax = nexttile(layout,source+6); hold(ax,'on');
    plot(ax,item.time,item.referenceTorque,'Color',[.85 .2 .15],'LineWidth',1.4);
    plot(ax,item.time,out(:,4),'Color',[0 .447 .741],'LineWidth',1.4);
    plot(ax,item.time,out(:,5),'Color',colors(3,:),'LineWidth',1.4);
    plot(ax,item.time,out(:,3),'Color',[.15 .15 .15],'LineWidth',1.5);
    title(ax,'CEINMS estimate versus actual Hill torques'); ylabel(ax,'Torque (N m)');
    referenceLabel = 'CEINMS estimate';
    if endsWith(subject.id,'A01'), referenceLabel = 'CEINMS estimate, A=0.1'; end
    legend(ax,{referenceLabel,'Back magnitude','Abdominal magnitude','Net = back - abdomen'},'Location','best');
end
axesList = findall(fig,'Type','axes');
for ax = axesList'
    grid(ax,'on'); box(ax,'on'); xlim(ax,item.time([1 end])); xlabel(ax,'Action time (s)');
end
title(layout,subject.id+" | 5 kg / 45 | Prepared activation is not the final Hill activation",'Interpreter','none');
save_figure(fig,folder,subject.shortName+"_signal_path",showFigures);
end

function plot_structure(folder,showFigures)
[fig,layout] = new_figure("BSE model structure",1,2,showFigures,[1500 820]);
ax = nexttile(layout); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
xlim(ax,[-.75 1.35]); ylim(ax,[-.30 1.65]);
q = 1; tip = [cos(q),sin(q)];
plot(ax,[-.65 1.2],[0 0],'Color',[.5 .5 .5],'LineWidth',1.2);
plot(ax,[0 tip(1)],[0 tip(2)],'Color',[.24 .28 .35],'LineWidth',12);
plot(ax,0,0,'ko','MarkerFaceColor','k','MarkerSize',7);
plot(ax,tip(1),tip(2),'o','Color',[.24 .28 .35],'MarkerSize',15,'LineWidth',2);
arc = linspace(0,q,70); plot(ax,.32*cos(arc),.32*sin(arc),'k-','LineWidth',1.1);
text(ax,.38,.12,'\theta','FontSize',16);
text(ax,-.15,-.12,'L5/S1','FontSize',12);
quiver(ax,.5*tip(1),.5*tip(2),0,-.35,0,'Color',[.6 .15 .15],'LineWidth',1.8,'MaxHeadSize',.35);
text(ax,.5*tip(1)+.04,.5*tip(2)-.24,'Body weight','FontSize',11);
quiver(ax,tip(1)+.13,tip(2),0,-.38,0,'Color',[.6 .15 .15],'LineWidth',1.8,'MaxHeadSize',.35);
text(ax,tip(1)+.17,tip(2)-.17,'Load: 5 / 10 / 15 kg','FontSize',11);
torque_arrow(ax,.48,[0 .447 .741],1);
torque_arrow(ax,.58,[.85 .35 .12],-1);
text(ax,-.70,.59,'Back: +\tau_b','Color',[0 .447 .741],'FontSize',13);
text(ax,.64,.59,'Abdomen: -\tau_a','Color',[.85 .35 .12],'FontSize',13);
text(ax,-.67,1.45,{'Equivalent sagittal-plane model', ...
    'Two antagonistic Hill muscle groups', ...
    'CE || PE, in series with elastic tendon'},'FontSize',12,'VerticalAlignment','top');
hill_symbol(ax,-.65,1.02,.72,'Back Hill',[0 .447 .741]);
hill_symbol(ax,.55,1.02,.72,'Abdominal Hill',[.85 .35 .12]);
text(ax,-.67,-.20,{'Schematic only; not anatomical scale.', ...
    'Constant effective moment arms; one rotational joint.'},'FontSize',10);
title(ax,'Mechanical structure','FontSize',15);
ax = nexttile(layout); hold(ax,'on'); axis(ax,[0 1 0 1]); axis(ax,'off');
boxes = {
    'Prepared experimental activation OR trained MAP activation';
    'Causal intent decoder + muscle proprioceptive feedback';
    'Neural command blend -> power mapping -> activation dynamics';
    'Back and abdominal Hill fibers + compliant tendons';
    'Coupled joint balance -> theta and muscle torques'};
for k = 1:5
    y = .90-(k-1)*.145;
    rectangle(ax,'Position',[.035 y-.075 .91 .105],'Curvature',.08, ...
        'FaceColor',[.94 .96 .98],'EdgeColor',[.4 .5 .6]);
    text(ax,.49,y-.023,boxes{k},'HorizontalAlignment','center','FontSize',10.5,'Interpreter','none');
    if k < 5
        plot(ax,[.49 .49],[y-.077 y-.103],'k-','LineWidth',1.2);
        patch(ax,[.48 .50 .49],[y-.103 y-.103 y-.115],'k','EdgeColor','none');
    end
end
plot(ax,[.96 .99 .99 .96],[.30 .30 .73 .73],'Color',[.49 .18 .56],'LineWidth',1.4);
text(ax,.98,.53,'theta / velocity feedback','Rotation',90,'HorizontalAlignment','center','FontSize',9);
text(ax,.04,.13,'Net torque = back - abdomen','FontSize',10.5);
text(ax,.04,.085,'$J\ddot{\theta}=\tau_b-\tau_a+\tau_g+\tau_s+\tau_{exo}-b\dot{\theta}$', ...
    'Interpreter','latex','FontSize',12);
text(ax,.04,.03,'Reference theta and CEINMS torque are comparison signals only.','FontSize',10.5);
title(ax,'Shared signal path for both activation sources','FontSize',15);
title(layout,'Current BSE | same Hill system produces theta and torque');
save_figure(fig,folder,"bse_model_structure",showFigures);
end

function torque_arrow(ax,radius,color,direction)
angle = linspace(.12,.80,45);
if direction < 0, angle = fliplr(angle); end
x = radius*cos(angle); y = radius*sin(angle);
plot(ax,x,y,'Color',color,'LineWidth',2);
tip = [x(end) y(end)];
direction = tip-[x(end-1) y(end-1)]; direction = direction/norm(direction);
base = tip-.07*direction; side = .025*[-direction(2) direction(1)];
patch(ax,[tip(1) base(1)+side(1) base(1)-side(1)], ...
    [tip(2) base(2)+side(2) base(2)-side(2)],color,'EdgeColor','none');
end

function hill_symbol(ax,x,y,width,label,color)
% A small CE/PE parallel pair followed by a series elastic tendon.
text(ax,x+width/2,y+.22,label,'Color',color,'FontSize',10.5,'HorizontalAlignment','center');
plot(ax,x+width*[0 .10 .10 .64 .64],y+[.06 .06 .13 .13 .06],'Color',color);
plot(ax,x+width*[.10 .10 .64 .64],y+[.06 -.02 -.02 .06],'Color',color);
rectangle(ax,'Position',[x+.22*width y+.08 .29*width .10],'FaceColor','w','EdgeColor',color);
text(ax,x+.365*width,y+.13,'CE','FontSize',8.5,'HorizontalAlignment','center','Color',color);
text(ax,x+.365*width,y-.07,'PE','FontSize',8.5,'HorizontalAlignment','center','Color',color);
springX = x+width*[.64 .70 .75 .80 .85 .90 .96 1];
springY = y+[.06 .06 .10 .02 .10 .02 .06 .06];
plot(ax,springX,springY,'Color',color,'LineWidth',1.3);
text(ax,x+.83*width,y+.15,'SE','FontSize',8.5,'HorizontalAlignment','center','Color',color);
end

function [lower,upper] = soft_limits(h,q)
low = max(h.boundaryLower-q,0); high = max(q-h.boundaryUpper,0);
low = low.^2./(low+.05); high = high.^2./(high+.05);
lower = h.boundaryCapacity*tanh(h.boundaryStiffness*low/h.boundaryCapacity);
upper = -h.boundaryCapacity*tanh(h.boundaryStiffness*high/h.boundaryCapacity);
end

function item = find_case(cases,load,speed)
index = find(cellfun(@(c) c.loadKg == load && c.speedCode == speed,cases),1);
if isempty(index), index = find(cellfun(@(c) c.loadKg == load,cases),1); end
if isempty(index), index = 1; end
item = cases{index};
end

function text = condition_title(item,index)
text = sprintf('%g kg | speed code %g',item.loadKg,item.speedCode);
if item.speedCode < 0 && isfield(item,'cadenceBPM')
    text = sprintf('%g kg | %.2f BPM | bout %g',item.loadKg,item.cadenceBPM,index);
end
end

function [fig,layout] = new_figure(name,rows,columns,showFigures,sizePixels)
visible = 'off'; if showFigures, visible = 'on'; end
fig = figure('Name',name,'Color','w','Position',[50 50 sizePixels],'Visible',visible);
layout = tiledlayout(fig,rows,columns,'TileSpacing','compact','Padding','compact');
end

function save_figure(fig,folder,name,showFigures)
exportgraphics(fig,fullfile(folder,name+".png"),'Resolution',160);
savefig(fig,fullfile(folder,name+".fig"));
if ~showFigures, close(fig); end
end
