function metrics=audit_hill_reference(root)
%AUDIT_HILL_REFERENCE Diagnostic only; never used to generate model predictions.
% Test whether the supplied reference pair satisfies the assumed joint mechanics.
    out=fullfile(root,'results','hill_forward'); d=load(fullfile(out,'hill_forward.mat')); r=d.result;
    records=struct([]);
    fig=figure('Visible','off','Color','w','Position',[30 30 1400 1100]); tiledlayout(3,3,'TileSpacing','compact');
    for k=1:9
        t=r.trials(k); q=sgolayfilt(t.q,3,21); v=gradient(q,.01); acc=gradient(v,.01);
        I=t.m*.760606^2+4*22.6006*.345^2/3;
        required=I.*acc+9.81*22.6006*.345*cos(q)+9.81*t.m*.760606.*cos(q-.936847);
        ix=11:numel(q)-10;
        % Report both ends of the permitted viscosity range, without fitting
        % a new torque model or assuming a particular activation source.
        rec=struct('loadKg',t.loadKg,'speed',t.speed, ...
            'referenceConsistencyRMSE_B0',sqrt(mean((required(ix)-t.torque(ix)).^2)), ...
            'referenceConsistencyRMSE_B8',sqrt(mean((required(ix)+8*v(ix)-t.torque(ix)).^2)), ...
            'negativeRequiredTorquePercent_B0',100*mean(required(ix)<0));
        if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
        nexttile; plot(t.t,[t.torque required required+8*v],'LineWidth',1.1); grid on;
        title(sprintf('%g kg / speed %g',t.loadKg,t.speed)); ylabel('Torque (Nm)'); xlabel('Time (s)');
        if k==1, legend('Recorded reference','Required by theta: B=0','Required by theta: B=8','Location','best'); end
    end
    sgtitle('Reference consistency audit ONLY | No curve here is a Hill prediction');
    exportgraphics(fig,fullfile(out,'reference_consistency.png'),'Resolution',160);
    exportgraphics(fig,fullfile(out,'reference_consistency.pdf'),'ContentType','vector'); close(fig);
    metrics=struct2table(records); writetable(metrics,fullfile(out,'reference_consistency.csv')); disp(metrics);
end
