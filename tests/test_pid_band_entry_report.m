function test_pid_band_entry_report
%TEST_PID_BAND_ENTRY_REPORT Check exported events against saved physical traces.
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'Demo','NMPC','output','current_BSE_528');
d=load(fullfile(folder,'current_BSE_NMPC_528.mat'));r=d.results;
t=readtable(fullfile(folder,'PID_NMPC_band_entry_events.csv'),'TextType','string');
summary=readtable(fullfile(folder,'PID_NMPC_band_entry_summary.csv'),'TextType','string');
expected=0;
for group=1:3
    q=r.groups{group}.sequence;
    for source=1:2
        pair=r.groups{group}.runs{source};
        for stage=1:3
            seg=q.segments{stage};ix=(seg.startIndex:seg.endIndex)';
            active=pair.baseline(ix,4)>pair.threshold(ix);
            edge=diff([false;active;false]);begins=find(edge==1);ends=find(edge==-1)-1;
            expected=expected+2*sum((ends-begins+1)*q.dt>=.30-1e-10);
        end
    end
end
assert(height(t)==expected,'An eligible demand event was added or omitted.');
for i=1:height(t)
    q=r.groups{t.Group(i)}.sequence;pair=r.groups{t.Group(i)}.runs{t.SourceIndex(i)};
    seg=q.segments{t.Stage(i)};first=t.DemandStartIndex(i);
    last=first+round(t.DemandDurationSeconds(i)/q.dt)-1;
    assert(first>=seg.startIndex && last<=seg.endIndex);
    assert(all(pair.baseline(first:last,4)>pair.threshold(first:last)));
    assert(abs(t.ThresholdNm(i)-pair.threshold(first))<1e-9);
    assert(abs(t.BandHalfWidthNm(i)-max(5,.05*t.ThresholdNm(i)))<1e-9);
    run=pair.(lower(char(t.Controller(i))));error=run.actual(first:last,4)-t.ThresholdNm(i);
    [~,peak]=max(error);inside=abs(error)<=t.BandHalfWidthNm(i);holdSamples=ceil(.15/q.dt-1e-10)+1;
    smoothed=movmean(error,max(1,round(.05/q.dt)));
    assert(t.OvershootExcursions(i)==sum(diff([false;smoothed>t.BandHalfWidthNm(i);false])==1));
    assert(t.UndershootExcursions(i)==sum(diff([false;smoothed<-t.BandHalfWidthNm(i);false])==1));
    assert(abs(t.BandFractionAfterPeak(i)-mean(inside(peak:end)))<1e-12);
    lastOutside=find(~inside(peak:end),1,'last');
    if isempty(lastOutside),finalIndex=peak;else,finalIndex=peak+lastOutside;end
    settled=numel(inside)-finalIndex+1>=holdSamples;
    assert(t.SettledByDemandEnd(i)==settled);
    if settled
        assert(abs(t.FinalSettlingTimeSeconds(i)-(finalIndex-1)*q.dt)<1e-10);
    else
        assert(isnan(t.FinalSettlingTimeSeconds(i)));
    end
    if t.SustainedEntry(i)
        entry=round(t.FirstEntrySeconds(i)/q.dt)+1;
        assert(entry>=peak && entry+holdSamples-1<=numel(inside));
        assert(all(inside(entry:entry+holdSamples-1)));
        assert(t.StayedUntilDemandEnd(i)==all(inside(entry:end)));
        assert(abs(t.ConfirmationTimeSeconds(i)-t.FirstEntrySeconds(i)-(holdSamples-1)*q.dt)<1e-10);
        for earlier=peak:entry-1
            if earlier+holdSamples-1<=numel(inside)
                assert(~all(inside(earlier:earlier+holdSamples-1)),'Entry is not the first qualifying hold.');
            end
        end
    else
        assert(isnan(t.FirstEntrySeconds(i)) && ~t.StayedUntilDemandEnd(i));
        for trial=peak:numel(inside)-holdSamples+1
            assert(~all(inside(trial:trial+holdSamples-1)),'A qualifying hold was omitted.');
        end
    end
end
for i=1:height(summary)
    x=t(t.Source==summary.Source(i) & t.Controller==summary.Controller(i),:);
    assert(summary.EligibleEvents(i)==height(x));
    assert(summary.SustainedEntries(i)==sum(x.SustainedEntry));
    assert(summary.StayedUntilDemandEnd(i)==sum(x.StayedUntilDemandEnd));
    assert(summary.SettledByDemandEnd(i)==sum(x.SettledByDemandEnd));
    assert(summary.AtMostTwoOvershootsAndSettled(i)==sum(x.OvershootExcursions<=2 & x.SettledByDemandEnd));
end
fprintf('PASS: post-peak band entries, confirmation times, unsettled events and source summaries match the physical traces.\n');
end
