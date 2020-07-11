% Tianben Ding 190226

% This code returns subfolder names in a folder. No '.' or '..' from dir
% function

function [subFolderName] = subdir(folderName)
allFolderCont = dir(folderName);
isdirInd = [allFolderCont(1:end).isdir];

allFolderCont(isdirInd == 0) = [];

inds = zeros(length(allFolderCont), 1);
n = 0;
k = 1;

while n < 2 && k <= length(allFolderCont)
    if any(strcmp(allFolderCont(k).name, {'.', '..'}))
        inds(k) = 1;
        n = n + 1;
    end
    k = k + 1;
end

allFolderCont(logical(inds)) = [];
subFolderName = cell(length(allFolderCont), 1);
for i = 1:length(allFolderCont)
    subFolderName{i} = allFolderCont(i).name;
end

end