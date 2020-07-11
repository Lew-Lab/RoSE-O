folders = split(pwd, filesep);
homepath = join(folders(1:end - 1), filesep);
addpath(homepath{:});

folders = {'utils', 'phasemask'};

for ii = 1:length(folders)

    addpath(genpath(fullfile(homepath{:}, folders{ii})));

end