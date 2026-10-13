function filename = exo_data_file(root, trial, name)
%EXO_DATA_FILE Resolve both archived and supplied recording folder layouts.
%   Original scripts refer to a CalibrationData level that is absent in the
%   supplied 240528Pilot dataset. Prefer the exact path, then try that layout.
    parts = regexp(char(trial), '[\\/]+', 'split');
    parts = parts(~cellfun('isempty', parts));
    filename = fullfile(root, parts{:}, name);
    if ~isfile(filename)
        parts(strcmp(parts, 'CalibrationData')) = [];
        candidate = fullfile(root, parts{:}, name);
        if isfile(candidate)
            filename = candidate;
        end
    end
    assert(isfile(filename), 'Exo:MissingRecording', ...
        'Recording not found: %s', filename);
end
