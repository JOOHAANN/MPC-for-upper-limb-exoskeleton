function folder = exo_data_root()
%EXO_DATA_ROOT Locate raw recordings, optionally outside the repository.
%   Set EXO_DATA_ROOT to override the default data/raw directory.
    folder = getenv('EXO_DATA_ROOT');
    if isempty(folder)
        root = fileparts(fileparts(mfilename('fullpath')));
        folder = fullfile(root, 'data', 'raw');
    end
    assert(isfolder(folder), 'Exo:MissingData', ...
        'Raw data directory not found: %s. Set EXO_DATA_ROOT.', folder);
end
