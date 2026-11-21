% One experimental tuning across all nine recorded conditions.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=1;
if ~exist('Release_rate','var'), Release_rate=28; end
if ~exist('Release_band','var'), Release_band=5; end
if ~exist('Angle_lead','var'), Angle_lead=0; end
if ~exist('Early_recovery','var'), Early_recovery=false; end
if ~exist('Guard_floor','var'), Guard_floor=1; end
for group=1:3
    sequence=subject.sequences{group}; model=subject.models{source};
    plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
    inputs=bse_inputs(sequence,model,source); planned=bse_inputs(sequence,model,2);
    cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528', ...
        sprintf('nmpc_group%d_source%d.mat',group,source)));
    options=cache.options; options.blocks=8;
    options.predictionErrorGuardScale=.3; options.predictionErrorGuardFloorNm=Guard_floor;
    options.backTorqueMarginNm=.5; options.releaseRateNmPerSec=Release_rate;
    options.releaseTrackingBandNm=Release_band;
    options.maximumAngleLeadRad=Angle_lead;
    if Angle_lead>0, options.angleGovernorUpperRad=1.7; end
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=Early_recovery;
    fprintf('\nEXPERIMENTAL TRACKING group %d | release rate %g\n',group,Release_rate);
    run=simulate_threshold_nmpc(plant,sequence,inputs,planned,cache.threshold,options);
    rows=[];
    for stage=1:3
        seg=sequence.segments{stage}; ix=seg.startIndex:seg.endIndex;
        demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
        gap=cache.threshold(demand)-run.actual(demand,4);
        rows=[rows;seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
            mean(max(gap,0)) max(run.actual(ix,25))]; %#ok<AGROW>
    end
    fprintf('EXPERIMENTAL RESULT group %d [load speed maximum-excess mean-demand-gap assist-peak]\n',group);
    disp(rows);
    save(fullfile(folder,sprintf('experimental_release%g_band%g_lead%g_early%d_floor%g_group%d.mat',Release_rate,Release_band,Angle_lead,Early_recovery,Guard_floor,group)), ...
        'run','options','rows','source','group','-v7');
end
