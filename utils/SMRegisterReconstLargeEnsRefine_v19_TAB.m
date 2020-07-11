% 180611 - v19 TD - changed integration area for intensity calculation from
% sinple square to rectangle

% 180530 - v18 TD - Changed start and end frame indexes to start index and
% total frame number.

% 180403 - v17 TD - added BinLimits to h3 figure for the all possible
% distance distribution.

% 171025 - v16 TD - Save a csv file contains localization number in each
% frame.

% 171024 - v15 TD - Added a parameter for tuning bias toralance of two
% channel registration.

% 171021 - v14 TD - Changed ROI definition. Manual selection with roipoly
% -> automatic selection.

% 171021 - v13 TD - Shifted bin positions of on-time histograms

% 171017 - v12 TD - Fixed the way to calculate locaNumROIBinShift and
% burstNumROIBinShift.

% 171013 - v11 TD - Added mean and std information of some critical
% analysis resutls. Debugged localization correlation part by removing
% over-counting of unpaired localizations.

% 171005 - v10 TD - Added ROI using roipoly function. Tuned visualization
% and output saving.

% 171004 - v9 TD - Changed the way to get paired precision from mean to
% max. Also changed the final precision histogram to only contain
% localizations after thresholding.

% 170927 - v8 TD - Added startFrame. Changed totFrame -> endFrame.

% 170921 - v7 TD - Correlate bursts across multiple frames and calculate
% total photon and on-time per burst.

% 170919 - v6 TD - Unified the burst number analysis and turned off Figure
% h15. Now, the code accepts frame number not only multiples of any
% specific number. Also, added an ensemble factor for adjusting
% registration tolerance. Debugged angle filtering. Changed window -> bin
% for reconstruction. Changed colormap input.

% 170903 - v5 TD - Changed angular filtering (for re-generating transform
% function) from the combination of sin and cos to an angle (radian) based
% method. Tuned histogram bin sizes.

% 170830 - v4 TD - Tuned the algorithm a little. Changed restriction for
% sin and cos fitting part in case we need to re-generate the transforma
% function. Also changed the restriciton for distance population
% (windowSize -> disSigma). Also, fixed the final pairing and unpairing
% algorithm. Now, the code not only working on the frames containing
% localizations in both channels. Adjusted the way to deal with
% channelRatio.

% 170827 - v3 TD - If there is any obvious bias on the two channel
% registration, compensate it and re-generate a transform function based on
% the single-molecule data. Otherwise, go forward with the original
% calibrated transform function.

% 170824 Tianben Ding

% This code matches SM standard PSF in two channels on a detector, creates
% a reconstruction, and calculates the photon numbers and the burst rate.
% You should have a localization list of the standard PSFs e.g., from
% ThunderSTORM, and a rough transform function for registrating two
% channels.

% This code assumes your images are from the polarization sensitive 4f
% system with pyramidal mirror.

% This code is based on the first and the visualization part of
% linearDichroism_r11.m

% The raw data should be a single tif sequence.

% The calculated localization precision is based on the LS part of the
% reference shown below
% 1. Rieger, B. & Stallinga, S. The lateral and axial localization
% uncertainty in super-resolution light microscopy. ChemPhysChem 15,
% 664–670 (2014).

clear;
close all;
clc; % clear everything

% Get output figure size (for large figures)
P0 = get(0, 'ScreenSize');
P1(3:4) = P0(3:4) - 150;

% Default figure setting
set(0, 'DefaultFigureColor', 'White', ...
    'DefaultTextFontName', 'Arial', 'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontSize', 24, 'DefaultTextFontSize', 24, ...
    'DefaultAxesLabelFontSizeMultiplier', 1.2, ...
    'DefaultAxesTickLength', [0.01, 0.01], 'DefaultAxesTickDir', 'out', ...
    'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input parameters

fFormat = '.tif'; % file format
conversionF = 0.49; % conversion factor (electron/count)
pixelSize = 58.5; % pixel size on objective plane (nm)
% channelRatio = 0.9275;% FF458-Di02,490/60, BP-3
% channelRatio = 0.8873;% Di03-R488/561,523/610, intensity ration between two channels (intensity of left channel/intensity of right channel)
% channelRatio = 0.8272;% Di03-R488/561,525/50, BP-2
channelRatio = 0.8663; % Di02-R514,550/49, BP-4
% channelRatio = 0.8978;% Di02-R514,582/64, BP-5
% channelRatio = 1.1613;%Di03-R488/561,523/610, BP-1 with 561 excitation
expTime = 20; % ms, exposure time on a detector

% geometric transform information for channel registration (calibrated by bright nanostructures)
tformData = '180530_controlPointsMatching_r5\180530 YG beads,Di03-R488-561,523_610,vS100 geometric transformation data.mat';

% SM parameters (partial chip imaging)
boundaryXPosSM = 1024; % boudary pixel of two channels in x axis for single molecule imaging (pixel)
imageSizeVSM = 100; % vertical image size(pixel)
imageSizeHSM = 2048; % horizontal image size(pixel)
squareSV = 7; % square size for intensity integration (pixel, vertical, odd)
squareSH = 7; % square size for intensity integration (pixel, horizontal, odd)
startFrame = 1; % First frame for analysis
totFrame = 5000; % total frame number for analysis

sigmaLowThreshold = 50; % threshold for removing localization results with very small sigma (nm)
sigmaHighThreshold = 150; % threshold for removing localization results with very large sigma (nm)

photonLowThreshold = 100; % threshold for removing localization results with very dim intensity (photon number)
unpairedPhotonLowThreshold = photonLowThreshold; % photon threshold for unpaired bursts
photonLowThresholdReg = 0; % photon threshold for two channel registration
zoomRegion = [27.8, 29.3, 2, 4.1] * 1e+3; % Zoomed in region with target for final output figures (nm) [xStart xEnd yStart yEnd]
workingRegion = [] * 1e+3; % Working region with target for analysis (nm) [xStart xEnd yStart yEnd], this matrix can be empty.

% Raw data and localization data address for single molecule
dataPlaceSM = 'C:\Users_NotBackedUp\tding\RawData\5. Transient amyloid binding imaging\180726_2158BBHP TAB performance with different buffers\offsetSubtracted\Data12';
% save name
saveName = '180726 Buffer-3 488 100mW Data12';
saveFormat = '.png';

binSize = 20; % nm, bin size for reconstruction
bmin = [0]; % nm
bmax = [500]; %nm
scaleBar = 300; % nm, scale bar length
histBinSize = 100; % frame, histogram bin size for burst ratio plots
maxPhotonTick = 2500; % photon number, maximum photon number for visualizations of photons detected

% precision thresholding
preThre = binSize;

% bias tolerance
biasThre = 15; %preThre;

% factor of most common ensemble
ensFactor = 1;

% colormap information
colorMap = parula_hdr;
close all;

% Saturate population
satPop = 0.98; % 0-1

% threshold of ROI selection
threROI = 3;

% Input part end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get information for saving

t = datetime('today');

dirName = datestr(t, 'yymmdd');
dirName = [dirName, '_', mfilename];

if exist(dirName, 'dir') ~= 7
    mkdir(dirName);
end

save([dirName, '\', saveName, ' initial setting', '.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load geometric transform information

load(tformData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start analysis

% Read localized data
num = csvread([dataPlaceSM, '.csv'], 1, 0);

% Clear the localized emitters which are too close to the edge
index = (ceil(num(:, 2) ./ pixelSize) - (squareSH + 1) / 2) < 1;
num(index, :) = [];
index = (ceil(num(:, 3) ./ pixelSize) - (squareSV + 1) / 2) < 1;
num(index, :) = [];
index = (ceil(num(:, 2) ./ pixelSize) + (squareSH + 1) / 2) > imageSizeHSM;
num(index, :) = [];
index = (ceil(num(:, 3) ./ pixelSize) + (squareSV + 1) / 2) > imageSizeVSM;
num(index, :) = [];

% Remove very small sigma localization
index = num(:, 4) < sigmaLowThreshold;
num(index, :) = [];

% Remove very large sigma localization
index = num(:, 4) > sigmaHighThreshold;
num(index, :) = [];

% Select working and zoomed FOV regions
xPosRightTemp = num(num(:, 2) > pixelSize*boundaryXPosSM, 2);
xPosRightTemp = xPosRightTemp - pixelSize * boundaryXPosSM;
yPosRightTemp = num(num(:, 2) > pixelSize*boundaryXPosSM, 3);

% Figure coordinate
imgAxisX = (1:imageSizeHSM - boundaryXPosSM) - 0.5;
imgAxisY = (1:imageSizeVSM) - 0.5;
imgAxisX = imgAxisX .* pixelSize;
imgAxisY = imgAxisY .* pixelSize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot in right channel
figure('Position', P1);
hold on
scatter(xPosRightTemp, yPosRightTemp, 3, 'w', 'filled');
axis image ij
grid on;
grid minor;
set(gca, 'GridColor', 'white')
set(gca, 'MinorGridColor', 'white')
% axis([imgAxisX(1) imgAxisX(end) imgAxisY(1) imgAxisY(end)])
axis([(imgAxisX(end) - imgAxisX(1)) / 2 - 1.5 * (imgAxisY(end) - imgAxisY(1)), (imgAxisX(end) - imgAxisX(1)) / 2 + 1.5 * (imgAxisY(end) - imgAxisY(1)), ...
    imgAxisY(1), imgAxisY(end)])
box off
xlabel('x position [nm]')
ylabel('y position [nm]')
title('Localized right channel emitters')
set(gca, 'Color', 'black')
hold off;
h0 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot in right channel (zoom)
figure('Position', P1);
hold on
scatter(xPosRightTemp, yPosRightTemp, 3, 'w', 'filled');
axis image ij
grid on;
grid minor;
set(gca, 'GridColor', 'white')
set(gca, 'MinorGridColor', 'white')
axis(zoomRegion)
box off
xlabel('x position [nm]')
ylabel('y position [nm]')
title({['Zoom region confirmation']; ['Localized right channel emitters']})
set(gca, 'Color', 'black')
hold off;
h0_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = 'Is the zoom region reasonable? \n Yes-> Enter / No-> xStart xEnd yStart yEnd [um]: \n';
zoomCon = input(prompt, 's');
while ~isempty(zoomCon)
    zoomRegion = str2num(zoomCon) * 1e+3; % Zoomed in region with target for final output figures (nm) [xStart xEnd yStart yEnd]
    close(h0_1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scatter plot in right channel (zoom, double-check)
    figure('Position', P1);
    hold on
    scatter(xPosRightTemp, yPosRightTemp, 3, 'w', 'filled');
    axis image ij
    grid on;
    grid minor;
    set(gca, 'GridColor', 'white')
    set(gca, 'MinorGridColor', 'white')
    axis(zoomRegion)
    xlabel('x position [nm]')
    ylabel('y position [nm]')
    title({['Zoom region confirmation']; ['Localized right channel emitters']})
    set(gca, 'Color', 'black')
    hold off;
    h0_1 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    prompt = 'Are you sure? \n Yes-> Enter / No-> xStart xEnd yStart yEnd [um]: \n';
    zoomCon = input(prompt, 's');
end

% Check of working FOV
if isempty(workingRegion)
    %     workingRegion = [zoomRegion(1)-(zoomRegion(2)-zoomRegion(1))*2/3, zoomRegion(2)+(zoomRegion(2)-zoomRegion(1))*2/3,...
    %         zoomRegion(3)-(zoomRegion(4)-zoomRegion(3))*2/3, zoomRegion(4)+(zoomRegion(4)-zoomRegion(3))*2/3];
    workingRegion = [zoomRegion(1) - 0.4 * 1e+3, zoomRegion(2) + 0.4 * 1e+3, zoomRegion(3) - 0.4 * 1e+3, zoomRegion(4) + 0.4 * 1e+3];
end
% working FOV with target for the following analysis (nm) [xStart xEnd yStart yEnd]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot (working FOV)
figure('Position', P1);
hold on
scatter(xPosRightTemp, yPosRightTemp, 3, 'w', 'filled');
axis image ij
grid on;
grid minor;
set(gca, 'GridColor', 'white')
set(gca, 'MinorGridColor', 'white')
axis(workingRegion)
box off
xlabel('x position [nm]')
ylabel('y position [nm]')
title({['Working region confirmation']; ['Localized right channel emitters']})
set(gca, 'Color', 'black')
hold off;
h0_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = 'Is the working region reasonable? \n Yes-> Enter / No-> xStart xEnd yStart yEnd [um]: \n';
workingCon = input(prompt, 's');
while ~isempty(workingCon)
    workingRegion = str2num(workingCon) * 1e+3; % Zoomed in region with target for final output figures (nm) [xStart xEnd yStart yEnd]
    close(h0_2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scatter plot (zoom, double-check)
    figure('Position', P1);
    hold on
    scatter(xPosRightTemp, yPosRightTemp, 3, 'w', 'filled');
    axis image ij
    grid on;
    grid minor;
    set(gca, 'GridColor', 'white')
    set(gca, 'MinorGridColor', 'white')
    axis(workingRegion)
    box off
    xlabel('x position [nm]')
    ylabel('y position [nm]')
    title({['Working region confirmation']; ['Localized right channel emitters']})
    set(gca, 'Color', 'black')
    hold off;
    h0_2 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    prompt = 'Are you sure? \n Yes-> Enter / No-> xStart xEnd yStart yEnd [um]: \n';
    workingCon = input(prompt, 's');
end
export_fig(h0, [dirName, '\', saveName, ' h0, scatter, localizations in r-cha, large FOV', saveFormat], '-transparent')
export_fig(h0_1, [dirName, '\', saveName, ' h0_1, scatter, localizations in r-cha, zoom', saveFormat], '-transparent')
export_fig(h0_2, [dirName, '\', saveName, ' h0_2, scatter, localizations in r-cha, working FOV', saveFormat], '-transparent')
close(h0, h0_1, h0_2)

workingRegionRight = [workingRegion(1) + pixelSize * boundaryXPosSM, workingRegion(2) + pixelSize * boundaryXPosSM, ...
    workingRegion(3), workingRegion(4)];
workingRegionCoor = [workingRegion(1), workingRegion(3); workingRegion(2), workingRegion(4)];
workingRegionCoor = transformPointsInverse(tformRight2Left, workingRegionCoor);

workingRegionCoor(:, 1) = -workingRegionCoor(:, 1) + pixelSize * boundaryXPosSM;
workingRegionLeft = [workingRegionCoor(2, 1), workingRegionCoor(1, 1), workingRegionCoor(1, 2), workingRegionCoor(2, 2)];
% Be careful about the mirrored coordinate

% % Checking workingRegionLeft in the left channel
% xPosLeftTemp = num(num(:,2) < pixelSize*boundaryXPosSM,2);
% yPosLeftTemp = num(num(:,2) < pixelSize*boundaryXPosSM,3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scatter plot in left channel
% figure('Position',P0);
% hold on
% scatter(xPosLeftTemp,yPosLeftTemp,3,'w','filled');
% axis image ij
% grid on;
% grid minor;
% set(gca,'GridColor','white')
% set(gca,'MinorGridColor','white')
% axis(workingRegionLeft)
% box off
% xlabel('x position [nm]','FontSize',24)
% ylabel('y position [nm]','FontSize',24)
% title('Working region confirmation in left channel','FontSize',24)
% set(gca,'FontSize',20)
% set(gca,'Color','black')
% hold off;
% set(gcf, 'Color', 'w');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the localized emitters which are not in the working FOV (right channel)
indexXL = num(:, 2) > workingRegionRight(1);
indexXU = num(:, 2) < workingRegionRight(2);
indexYL = num(:, 3) > workingRegionRight(3);
indexYU = num(:, 3) < workingRegionRight(4);
indexAll = indexXL + indexXU + indexYL + indexYU;

numRight = num(indexAll > 3, :);

% Clear the localized emitters which are not in the working FOV (left channel)
indexXL = num(:, 2) > workingRegionLeft(1);
indexXU = num(:, 2) < workingRegionLeft(2);
indexYL = num(:, 3) > workingRegionLeft(3);
indexYU = num(:, 3) < workingRegionLeft(4);
indexAll = indexXL + indexXU + indexYL + indexYU;

numLeft = num(indexAll > 3, :);

% Only keep the emitters in the working FOV
num = [numLeft; numRight];

% Read raw data folder information
dataFolderInfo = dir([dataPlaceSM, '\*', fFormat]);

% Initialize memory for the next process
xPosAll = [];
yPosAll = [];
sigmaAll = [];
idAll = [];
backgroundAllVectorAll = [];
intensAllVectorAll = [];

xPosRightPreRegAll = [];
yPosRightPreRegAll = [];
xPosLeftPreRegAll = [];
yPosLeftPreRegAll = [];

precRightAll = [];
precLeftAll = [];

allVectorDisAll = [];
allVectorRadAll = [];
% allVectorCosAll = [];
% allVectorSinAll = [];
allVectorPosAll = [];
allVectorPrecAll = [];

for h = startFrame:(startFrame + totFrame - 1)
    id = transpose(num(num(:, 1) == h, 1));
    xPos = transpose(num(num(:, 1) == h, 2));
    yPos = transpose(num(num(:, 1) == h, 3));
    sigma = transpose(num(num(:, 1) == h, 4));

    if ~isempty(id)
        % Read the images
        capturedImage = tiffread2([dataPlaceSM, '\', dataFolderInfo(h).name], 1, 1);
        capturedImage = double(capturedImage.data);

        xPosCopy = xPos;
        yPosCopy = yPos;
        for i = 1:size(xPos, 1)
            distMat = sqrt((xPos(:)-xPos(i)).^2+(yPos(:) - yPos(i)).^2);
            % ind = distMat < sqrt(squareSV^2+squareSH^2)*pixelSize;
            ind = distMat < max(squareSV, squareSH) * pixelSize;
            if sum(ind) > 1 % if there is too close emitters
                % remove too close emitters
                id(ind) = nan;
                xPosCopy(ind) = nan;
                yPosCopy(ind) = nan;
                sigma(ind) = nan;
            end
        end
        xPos = xPosCopy;
        yPos = yPosCopy;

        xPos(isnan(xPos)) = [];
        yPos(isnan(yPos)) = [];
        sigma(isnan(sigma)) = [];
        id(isnan(id)) = [];

        % Subtraction of background
        % insideMatrix: contains all emitter candidates' PSF information
        insideMatrix = ...
            capturedImage(reshape(repmat(ceil(yPos ./ pixelSize), [squareSV, 1]), [1, (squareSV) * length(yPos)])+repmat(-(squareSV - 1) / 2:(squareSV - 1) / 2, [1, length(yPos)]), ...
            reshape(repmat(ceil(xPos ./ pixelSize), [squareSH, 1]), [1, (squareSH) * length(xPos)])+repmat(-(squareSH - 1) / 2:(squareSH - 1) / 2, [1, length(xPos)]) ...
            );
        insideIndex = kron(eye(length(xPos)), ones(squareSV, squareSH));
        insideMatrix(~logical(insideIndex)) = 0;
        insideVector = sum(insideMatrix, 1);
        insideVector = reshape(insideVector, [squareSH, length(xPos)]);
        insideVector = sum(insideVector, 1);

        % largeMatrix: contains all emitter candidates' PSF information (larger region)
        largeMatrix = ...
            capturedImage(reshape(repmat(ceil(yPos ./ pixelSize), [squareSV + 2, 1]), [1, (squareSV + 2) * length(yPos)])+repmat(-(squareSV + 2 - 1) / 2:(squareSV + 2 - 1) / 2, [1, length(yPos)]), ...
            reshape(repmat(ceil(xPos ./ pixelSize), [squareSH + 2, 1]), [1, (squareSH + 2) * length(xPos)])+repmat(-(squareSH + 2 - 1) / 2:(squareSH + 2 - 1) / 2, [1, length(xPos)]) ...
            );
        largeIndex = kron(eye(length(xPos)), ones(squareSV + 2, squareSH + 2));
        largeMatrix(~logical(largeIndex)) = 0;
        largeVector = sum(largeMatrix, 1);
        largeVector = reshape(largeVector, [squareSH + 2, length(xPos)]);
        largeVector = sum(largeVector, 1);

        backgroundVector = largeVector - insideVector;
        backgroundVector = backgroundVector ./ (squareSV * 2 + (squareSH + 2) * 2);
        % backgroundVector = round(backgroundVector);

        % backgroundMatrix: contains background information at all emitters
        % candidates' positions
        backgroundMatrix = reshape(repmat(backgroundVector, [length(xPos) * squareSV * squareSH, 1]), [length(xPos) * squareSV, length(xPos) * squareSH]);
        backgroundMatrix(~logical(insideIndex)) = 0;

        % Subtract background from inside box
        % intensMatrix: contains all emitter candidates' PSF information after
        % background subtraction
        intensMatrix = insideMatrix - backgroundMatrix;
        % intensMatrix(intensMatrix < 0) = 0;

        % figure('Position',P1);
        % imagesc(intensMatrix)
        % caxis([0 100])
        % colormap hot;
        % colorbar;

        intensAllVector = sum(intensMatrix, 1);
        intensAllVector = reshape(intensAllVector, [squareSH, length(xPos)]); % not rounded yet
        %intensAllVector = round(sum(intensAllVector,1).*conversionF);
        intensAllVector = sum(intensAllVector, 1) .* conversionF;
        backgroundAllVector = backgroundVector .* conversionF;

        % Remove localizations with negative intensities
        xPos(intensAllVector <= 0) = [];
        yPos(intensAllVector <= 0) = [];
        sigma(intensAllVector <= 0) = [];
        id(intensAllVector <= 0) = [];
        backgroundAllVector(intensAllVector <= 0) = [];
        intensAllVector(intensAllVector <= 0) = [];

        % Record for the next process
        xPosAll = [xPosAll, xPos];
        yPosAll = [yPosAll, yPos];
        sigmaAll = [sigmaAll, sigma];
        idAll = [idAll, id];
        backgroundAllVectorAll = [backgroundAllVectorAll, backgroundAllVector];
        intensAllVectorAll = [intensAllVectorAll, intensAllVector];

        % Distinguish two channels
        rightInd = xPos > boundaryXPosSM * pixelSize;
        leftInd = xPos <= boundaryXPosSM * pixelSize;

        % Store the right channel information
        xPosRight = xPos(rightInd);
        yPosRight = yPos(rightInd);
        sigmaRight = sigma(rightInd);
        intensRight = intensAllVector(rightInd);
        backgroundRight = backgroundAllVector(rightInd);
        if channelRatio > 1
            intensRight = intensRight .* channelRatio;
            backgroundRight = backgroundRight .* channelRatio;
        end
        photonRight = round(intensRight);
        backgroundRight = round(backgroundRight);
        idRight = id(rightInd);

        % Remove localization if its photon number is smaller than
        % threshold (registration)
        ind = photonRight <= photonLowThresholdReg;
        xPosRight(ind) = [];
        yPosRight(ind) = [];
        sigmaRight(ind) = [];
        photonRight(ind) = [];
        backgroundRight(ind) = [];
        idRight(ind) = [];

        % Store the left channel information
        xPosLeft = xPos(leftInd);
        yPosLeft = yPos(leftInd);
        sigmaLeft = sigma(leftInd);
        intensLeft = intensAllVector(leftInd);
        backgroundLeft = backgroundAllVector(leftInd);
        if channelRatio <= 1
            intensLeft = intensLeft ./ channelRatio;
            backgroundLeft = backgroundLeft ./ channelRatio;
        end
        photonLeft = round(intensLeft);
        backgroundLeft = round(backgroundLeft);
        idLeft = id(leftInd);

        % Remove localization if its photon number is smaller than
        % threshold (registration)
        ind = photonLeft <= photonLowThresholdReg;
        xPosLeft(ind) = [];
        yPosLeft(ind) = [];
        sigmaLeft(ind) = [];
        photonLeft(ind) = [];
        backgroundLeft(ind) = [];
        idLeft(ind) = [];

        % Calculate localization precision based on captured photon number
        % (Least square)
        tauRight = 2 * pi * backgroundRight .* (sigmaRight.^2 + pixelSize^2 / 12) ./ photonRight ./ (pixelSize^2);
        precRight = (sigmaRight.^2 + pixelSize^2 / 12) ./ photonRight .* (16 / 9 + 4 * tauRight);
        % precRight = (sigmaRight.^2+pixelSize^2/12)./photonRight.*(1+4*tauRight+sqrt(2*tauRight./(1+4*tauRight)));
        precRight = sqrt(precRight);

        tauLeft = 2 * pi * backgroundLeft .* (sigmaLeft.^2 + pixelSize^2 / 12) ./ photonLeft ./ (pixelSize^2);
        precLeft = (sigmaLeft.^2 + pixelSize^2 / 12) ./ photonLeft .* (16 / 9 + 4 * tauLeft);
        % precLeft = (sigmaLeft.^2+pixelSize^2/12)./photonLeft.*(1+4*tauLeft+sqrt(2*tauLeft./(1+4*tauLeft)));
        precLeft = sqrt(precLeft);

        % thresholding using the calculated precision
        ind = precRight > preThre;
        xPosRight(ind) = [];
        yPosRight(ind) = [];
        sigmaRight(ind) = [];
        photonRight(ind) = [];
        backgroundRight(ind) = [];
        idRight(ind) = [];
        precRight(ind) = [];
        precRightAll = [precRightAll, precRight];

        ind = precLeft > preThre;
        xPosLeft(ind) = [];
        yPosLeft(ind) = [];
        sigmaLeft(ind) = [];
        photonLeft(ind) = [];
        backgroundLeft(ind) = [];
        idLeft(ind) = [];
        precLeft(ind) = [];
        precLeftAll = [precLeftAll, precLeft];

        if ~isempty(xPosRight) && ~isempty(xPosLeft)
            xPosRight = xPosRight - pixelSize * boundaryXPosSM;
            xPosLeft = -xPosLeft + pixelSize * boundaryXPosSM;

            xPosRightPreRegAll = [xPosRightPreRegAll, xPosRight];
            yPosRightPreRegAll = [yPosRightPreRegAll, yPosRight];
            xPosLeftPreRegAll = [xPosLeftPreRegAll, xPosLeft];
            yPosLeftPreRegAll = [yPosLeftPreRegAll, yPosLeft];

            [xPosLeftTrans, yPosLeftTrans] = transformPointsInverse(tform, xPosLeft.', yPosLeft.');
            xPosLeftTrans = xPosLeftTrans.';
            yPosLeftTrans = yPosLeftTrans.';

            rightMeshX = repmat(xPosRight, length(xPosLeftTrans), 1);
            rightMeshY = repmat(yPosRight, length(yPosLeftTrans), 1);
            leftMeshX = repmat(xPosLeft.', 1, length(xPosRight));
            leftMeshY = repmat(yPosLeft.', 1, length(yPosRight));
            leftTransMeshX = repmat(xPosLeftTrans.', 1, length(xPosRight));
            leftTransMeshY = repmat(yPosLeftTrans.', 1, length(yPosRight));
            precRightMesh = repmat(precRight, length(precLeft), 1);
            precLeftMesh = repmat(precLeft.', 1, length(precRight));

            % calculate Euclidean distance between all SMs in two channels
            allVectorDis = hypot(rightMeshX-leftTransMeshX, rightMeshY-leftTransMeshY);
            allVectorCos = (rightMeshX - leftTransMeshX) ./ allVectorDis;
            allVectorSin = (rightMeshY - leftTransMeshY) ./ allVectorDis;
            allVectorRad = angle(allVectorCos+1j*allVectorSin);
            allVectorDisAll = [allVectorDisAll, reshape(allVectorDis, [1, size(allVectorDis, 1) * size(allVectorDis, 2)])];
            allVectorRadAll = [allVectorRadAll, reshape(allVectorRad, [1, size(allVectorRad, 1) * size(allVectorRad, 2)])];
            %             allVectorCosAll = [allVectorCosAll reshape(allVectorCos,[1,size(allVectorCos,1)*size(allVectorCos,2)])];
            %             allVectorSinAll = [allVectorSinAll reshape(allVectorSin,[1,size(allVectorSin,1)*size(allVectorSin,2)])];
            allVectorPosAll = [allVectorPosAll, ...
                [reshape(rightMeshX, [1, size(rightMeshX, 1) * size(rightMeshX, 2)]); ...
                reshape(rightMeshY, [1, size(rightMeshY, 1) * size(rightMeshY, 2)]); ...
                reshape(leftMeshX, [1, size(leftMeshX, 1) * size(leftMeshX, 2)]); ...
                reshape(leftMeshY, [1, size(leftMeshY, 1) * size(leftMeshY, 2)])]];
            allVectorPrecAll = [allVectorPrecAll, ...
                [reshape(precRightMesh, [1, size(precRightMesh, 1) * size(precRightMesh, 2)]); ...
                reshape(precLeftMesh, [1, size(precLeftMesh, 1) * size(precLeftMesh, 2)])]];
        end
    end
end

[xPosLeftPreRegAllTrans, yPosLeftPreRegAllTrans] = transformPointsInverse(tform, xPosLeftPreRegAll.', yPosLeftPreRegAll.');
xPosLeftPreRegAllTrans = xPosLeftPreRegAllTrans.';
yPosLeftPreRegAllTrans = yPosLeftPreRegAllTrans.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all localization used for registration check (before pairing, with tform, only ones with high precision)
figure('Position', P1);
hold on;
scatter(xPosRightPreRegAll, yPosRightPreRegAll, 'LineWidth', 1)
scatter(xPosLeftPreRegAllTrans, yPosLeftPreRegAllTrans, 'LineWidth', 1)
axis image ij
axis(zoomRegion)
% axis(workingRegion)
grid on
grid minor
legend('rightPoints', 'leftPointsTrans')
xlabel('x position [nm]')
ylabel('y position [nm]')
title({['Geometric transform of SM images']; ...
    ['before pairing, localizations with precision < ', num2str(preThre), ' nm']})
hold off
h1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localization precision of the all localization (histogram, before pairing, only ones with high precision)
figure('Position', P1);
hold on;
hPR = histogram(precRightAll);
hPL = histogram(precLeftAll);
hPL.BinWidth = hPR.BinWidth;
legend(['right cha, median = ', num2str(median(precRightAll)), ' nm'], ['left cha, median = ', num2str(median(precLeftAll)), ' nm'])
xlabel('localization precision [nm]')
% ylabel('Frequency')
title({['Localization precision distribution of emitters']; ...
    ['after ', num2str(photonLowThresholdReg), ' photon and ', num2str(preThre), ' nm precision thresholding']})
hold off
h2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all possible pairing distance (histogram, only ones with high precision)
figure('Position', P1);
hold on;
if ~isempty(bmin) && ~isempty(bmax)
    histDist = histogram(allVectorDisAll, 'BinLimits', [bmin, bmax]);
else
    histDist = histogram(allVectorDisAll);
end
[~, histDistMaxInd] = max(histDist.Values);
histDist.BinEdges = 0:binSize / 4:histDist.BinEdges(histDistMaxInd+1) * 2;
% Gaussian equation for one dimentional fitting
gaussEqn = 'a*exp(-(x-b)^2/(2*c^2))+d';
x = histDist.BinEdges + histDist.BinWidth / 2;
x = x(1:end-1);
[maxi, maxiInd] = max(histDist.Values);
startPoints = [maxi, x(maxiInd), binSize, 0];
f = fit(x.', (histDist.Values).', gaussEqn, 'Start', startPoints);
plot(f, x, (histDist.Values));
xlabel('Pairing distance [nm]')
% ylabel('Frequency')
title({['Distribution of all possible pairing distance']; ...
    ['after ', num2str(photonLowThresholdReg), ' photon and ', num2str(preThre), ' nm precision thresholding']})
hold off
h3 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disPeak = f.b;
disSigma = (f.c) / ensFactor;

export_fig(h1, [dirName, '\', saveName, ' h1, scatter, localizations with high prec, without pair and corr', saveFormat], '-transparent')
export_fig(h2, [dirName, '\', saveName, ' h2, hist, prec distribution after thresholding', saveFormat], '-transparent')
export_fig(h3, [dirName, '\', saveName, ' h3, hist, all possible dist after prec thresholding before pairing', saveFormat], '-transparent')
close(h1, h2, h3)

if disPeak < hypot(biasThre, biasThre)
    % If there is no bias on the registration, go forward with the
    % calibratyed transform function
    tformSM = tform;
    tformSMRight2Left = tformRight2Left;
    matchErrorMeanSM = matchErrorMean;
    matchErrorStdSM = matchErrorStd;
    matchErrorMaxSM = matchErrorMax;
else
    % If there is any obvious linear bias, compensate it and re-generate a
    % new transform function
    indL = allVectorDisAll > (disPeak - disSigma);
    indU = allVectorDisAll < (disPeak + disSigma);
    indLU = (indL + indU) > 1;

    % Radian filtering
    workingRad = allVectorRadAll(indLU);
    resVec = sum(exp(1j * workingRad));
    meanRad = angle(resVec);
    resVecL = abs(resVec) / length(workingRad);
    stdRad = sqrt(-2*log(resVecL));
    xRadMin = meanRad - pi / 2;
    if xRadMin < -pi
        workingRad((pi-abs(-pi - xRadMin)) < workingRad) = workingRad((pi-abs(-pi - xRadMin)) < workingRad) - 2 * pi;
    end
    xRadMax = meanRad + pi / 2;
    if xRadMax > pi
        workingRad(workingRad < (-pi + abs(xRadMax - pi))) = workingRad(workingRad < (-pi + abs(xRadMax - pi))) + 2 * pi;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Aangles between all possible pairing localization, after distance thresholding (histogram)
    figure('Position', P1);
    hold on;
    histRad = histogram(workingRad);
    histRad.BinEdges = xRadMin:0.05:xRadMax;
    % Gaussian equation for one dimentional fitting
    gaussEqn = 'a*exp(-(x-b)^2/(2*c^2))+d';
    x = histRad.BinEdges + histRad.BinWidth / 2;
    x = x(1:end-1);
    [maxRad, maxRadInd] = max(histRad.Values);
    startPointsRad = [maxRad, x(maxRadInd), pi / 4, 0];
    fRad = fit(x.', (histRad.Values).', gaussEqn, 'Start', startPointsRad);
    plot(fRad, x, (histRad.Values));
    xlabel('Radian')
    ylabel('Frequency')
    title({['Angles of all possible pairs, after ', num2str(photonLowThresholdReg), ' photon,']; ...
        [num2str(preThre), ' nm precision, and distance thresholding']})
    hold off
    h4 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    radPeak = fRad.b;
    radSigma = (fRad.c) / ensFactor;
    indRadU = allVectorRadAll < (radPeak + radSigma);
    indRadL = allVectorRadAll > (radPeak - radSigma);
    indLUAll = (indL + indU + indRadL + indRadU) > 3;

    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % cos and sin of angles between all possible pairing localization, after distance thresholding (histogram)
    %     figure('Position',P1);
    %     hold on;
    %     histCos = histogram(allVectorCosAll(indLU));
    %     histSin = histogram(allVectorSinAll(indLU));
    %     histCos.BinEdges = -1:0.01:1;
    %     histSin.BinEdges = -1:0.01:1;
    %     % Gaussian equation for one dimentional fitting
    %     gaussEqn = 'a*exp(-(x-b)^2/(2*c^2))+d';
    %     x = (-1:0.01:1)+0.01/2;
    %     x = x(1:end-1);
    %     [maxCos,maxCosInd] = max(histCos.Values);
    %     [maxSin,maxSinInd] = max(histSin.Values);
    %     startPointsCos = [maxCos x(maxCosInd) (1-(-1))/4 0];
    %     startPointsSin = [maxSin x(maxSinInd) (1-(-1))/4 0];
    %     fCos = fit(x.',(histCos.Values).',gaussEqn,'Start',startPointsCos);
    %     fSin = fit(x.',(histSin.Values).',gaussEqn,'Start',startPointsSin);
    %     plot(fCos,x,(histCos.Values));
    %     plot(fSin,x,(histSin.Values));
    %     legend('cos','sin')
    %     xlabel('sin/cos')
    %     ylabel('Frequency')
    %     title({['Cos and sin of all possible pairs, after ' num2str(photonLowThresholdReg) ' photon,'];...
    %         [num2str(preThre) ' nm precision, and distance thresholding']})
    %     hold off
    %     h4 = gcf;
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     cosPeak = fCos.b;
    %     cosSigma = fCos.c;
    %     sinPeak = fSin.b;
    %     sinSigma = fSin.c;
    %     if sinPeak < -1
    %         sinPeak = -1;
    %     elseif sinPeak > 1
    %         sinPeak = 1;
    %     end
    %     if cosPeak < -1
    %         cosPeak = -1;
    %     elseif cosPeak > 1
    %         cosPeak = 1;
    %     end
    %     if (cosPeak - cosSigma) < -1
    %         indCosL = allVectorCosAll >= -1;
    %         indCosU = allVectorCosAll < (cosPeak + cosSigma + (-1 - (cosPeak - cosSigma)));
    %     elseif (cosPeak + cosSigma) > 1
    %         indCosL = allVectorCosAll > (cosPeak - cosSigma - ((cosPeak + cosSigma) - 1));
    %         indCosU = allVectorCosAll <= 1;
    %     else
    %         indCosL = allVectorCosAll > (cosPeak - cosSigma);
    %         indCosU = allVectorCosAll < (cosPeak + cosSigma);
    %     end
    %     if (sinPeak - sinSigma) < -1
    %         indSinL = allVectorSinAll >= -1;
    %         indSinU = allVectorSinAll < (sinPeak + sinSigma + (-1 - (sinPeak - sinSigma)));
    %     elseif (sinPeak + sinSigma) > 1
    %         indSinL = allVectorSinAll > (sinPeak - sinSigma - ((sinPeak + sinSigma) - 1));
    %         indSinU = allVectorSinAll <= 1;
    %     else
    %         indSinL = allVectorSinAll > (sinPeak - sinSigma);
    %         indSinU = allVectorSinAll < (sinPeak + sinSigma);
    %     end
    %     indLUAll = (indL+indU+indCosL+indCosU+indSinL+indSinU) > 5;

    pairedPos = allVectorPosAll(:, indLUAll);
    pairedPrec = allVectorPrecAll(:, indLUAll);

    % Fit the old geometric transformation to the control point pairs
    movingPoints = pairedPos(1:2, :).';
    fixedPoints = pairedPos(3:4, :).';
    fixedPointsTransUncorr = transformPointsInverse(tform, fixedPoints);
    medXDir = median(movingPoints(:, 1)-fixedPointsTransUncorr(:, 1));
    medYDir = median(movingPoints(:, 2)-fixedPointsTransUncorr(:, 2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform without correction (paired localizations)
    figure('Position', P1);
    hold on;
    scatter(movingPoints(:, 1), movingPoints(:, 2), 'LineWidth', 1)
    scatter(fixedPointsTransUncorr(:, 1), fixedPointsTransUncorr(:, 2), 'LineWidth', 1)
    for k = 1:size(movingPoints, 1)
        plot([movingPoints(k, 1), fixedPointsTransUncorr(k, 1)], [movingPoints(k, 2), fixedPointsTransUncorr(k, 2)], '-g');
    end
    axis image ij
    axis(zoomRegion)
    % axis(workingRegion)
    grid on
    grid minor
    legend('rightPoints', 'leftPointsTrans')
    xlabel('x position [nm]')
    ylabel('y position [nm]')
    title({['Geometric transform of SM images']; ...
        ['without corrected transform, Paired # = ', num2str(size(movingPoints, 1))]})
    hold off
    h5 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % localization precision used for pairing (histogram)
    figure('Position', P1);
    hold on;
    histPrec1 = histogram(pairedPrec(1, :));
    histPrec2 = histogram(pairedPrec(2, :));
    histPrec2.BinWidth = histPrec1.BinWidth;
    legend(['right cha, #', num2str(length(pairedPrec(1, :))), ' median = ', num2str(median(pairedPrec(1, :))), ' nm'], ...
        ['left cha, #', num2str(length(pairedPrec(2, :))), ' median = ', num2str(median(pairedPrec(2, :))), ' nm'])
    xlabel('localization precision [nm]')
    % ylabel('Frequency')
    title('Localization precision used for pairing and re-calibration')
    hold off
    h6 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform without correction (histogram)
    figure('Position', P1);
    hold on;
    histX = histogram(movingPoints(:, 1)-fixedPointsTransUncorr(:, 1));
    histY = histogram(movingPoints(:, 2)-fixedPointsTransUncorr(:, 2));
    histY.BinWidth = histX.BinWidth;
    legend(['xDir, median:', num2str(medXDir), 'nm'], ...
        ['yDir, median:', num2str(medYDir), 'nm'])
    xlabel('Error [nm]')
    % ylabel('Frequency')
    title({['Error distribution of geometric transform']; ...
        ['without corrected transform']})
    hold off
    h7 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lengthReg = hypot((movingPoints(:, 1)-fixedPointsTransUncorr(:, 1)), (movingPoints(:, 2) - fixedPointsTransUncorr(:, 2)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform without correction (histogram of distance)
    figure('Position', P1);
    hold on;
    histTestDis = histogram(lengthReg);
    legend(['median= ', num2str(median(lengthReg)), ' nm']);
    xlabel('Error [nm]')
    % ylabel('Frequency')
    title({['Error distribution of geometric transform']; ...
        ['without corrected transform (distance)']})
    hold off
    h8 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tformSM = fitgeotrans(movingPoints, fixedPoints, 'polynomial', 2);
    tformSMRight2Left = fitgeotrans(fixedPoints, movingPoints, 'polynomial', 2);

    % Check the quality of fitgeotrans
    fixedPointsTrans = transformPointsInverse(tformSM, fixedPoints);
    movingPointsTrans = transformPointsInverse(tformSMRight2Left, movingPoints);

    disTrans = movingPoints - fixedPointsTrans;
    disTrans = hypot(disTrans(:, 1), disTrans(:, 2));

    matchErrorMeanSM = mean(disTrans);
    matchErrorStdSM = std(disTrans);

    matchErrorMaxSM = max(disTrans);

    disTransRight2Left = movingPointsTrans - fixedPoints;
    disTransRight2Left = hypot(disTransRight2Left(:, 1), disTransRight2Left(:, 2));

    matchErrorRight2LeftMeanSM = mean(disTransRight2Left);
    matchErrorRight2LeftStdSM = std(disTransRight2Left);

    matchErrorRight2LeftMaxSM = max(disTransRight2Left);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform with correction (paired localizations)
    figure('Position', P1);
    hold on;
    scatter(movingPoints(:, 1), movingPoints(:, 2), 'LineWidth', 1)
    scatter(fixedPointsTrans(:, 1), fixedPointsTrans(:, 2), 'LineWidth', 1)
    for k = 1:size(movingPoints, 1)
        plot([movingPoints(k, 1), fixedPointsTrans(k, 1)], [movingPoints(k, 2), fixedPointsTrans(k, 2)], '-g');
    end
    axis image ij
    axis(zoomRegion)
    % axis(workingRegion)
    grid on
    grid minor
    legend('rightPoints', 'leftPointsTrans')
    xlabel('x position [nm]')
    ylabel('y position [nm]')
    title({['Geometric transform of SM images']; ...
        ['with corrected transform, Paired # = ', num2str(size(movingPoints, 1))]})
    hold off
    h9 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform with correction (histogram)
    figure('Position', P1);
    hold on;
    histTest1 = histogram(movingPoints(:, 1)-fixedPointsTrans(:, 1));
    histTest2 = histogram(movingPoints(:, 2)-fixedPointsTrans(:, 2));
    histTest2.BinWidth = histTest1.BinWidth;
    legend(['xDir, median:', num2str(median(movingPoints(:, 1) - fixedPointsTrans(:, 1))), 'nm'], ...
        ['yDir, median:', num2str(median(movingPoints(:, 2) - fixedPointsTrans(:, 2))), 'nm'])
    xlabel('Error [nm]')
    % ylabel('Frequency')
    title({['Error distribution of geometric transform']; ...
        ['with corrected transform (FRE)']})
    hold off
    h10 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lengthRegCorr = hypot((movingPoints(:, 1)-fixedPointsTrans(:, 1)), (movingPoints(:, 2) - fixedPointsTrans(:, 2)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test figure transform with correction (histogram of distance)
    figure('Position', P1);
    hold on;
    histTestDis = histogram(lengthRegCorr);
    legend(['median= ', num2str(median(lengthRegCorr)), ' nm']);
    xlabel('Error [nm]')
    % ylabel('Frequency')
    title({['Error distribution of geometric transform']; ...
        ['with corrected transform (distance, FRE)']})
    hold off
    h11 = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    export_fig(h4, [dirName, '\', saveName, ' h4, hist, all possible CosSin after prec, dist thresholding', saveFormat], '-transparent')
    export_fig(h5, [dirName, '\', saveName, ' h5, scatter, localization after pairing, without corr', saveFormat], '-transparent')
    export_fig(h6, [dirName, '\', saveName, ' h6, hist, prec of localization for pairing', saveFormat], '-transparent')
    export_fig(h7, [dirName, '\', saveName, ' h7, hist, x,y bias for paired localizations, without corr', saveFormat], '-transparent')
    export_fig(h8, [dirName, '\', saveName, ' h8, hist, dist between paired localizations, without corr', saveFormat], '-transparent')

    export_fig(h9, [dirName, '\', saveName, ' h9, scatter, localization after pairing, with corr', saveFormat], '-transparent')
    export_fig(h10, [dirName, '\', saveName, ' h10, hist, x,y bias for paired localizations, with corr', saveFormat], '-transparent')
    export_fig(h11, [dirName, '\', saveName, ' h11, hist, dist between paired localizations, with corr', saveFormat], '-transparent')
    close(h4, h5, h6, h7, h8, h9, h10, h11)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use the fine tuned geometric transformation to pair and analyze bursts

% again
% Distinguish two channels
rightInd = xPosAll > boundaryXPosSM * pixelSize;
leftInd = xPosAll <= boundaryXPosSM * pixelSize;

% Store the right channel information
xPosAllRight = xPosAll(rightInd);
yPosAllRight = yPosAll(rightInd);
sigmaAllRight = sigmaAll(rightInd);
intensAllRight = intensAllVectorAll(rightInd);
backgroundAllRight = backgroundAllVectorAll(rightInd);
if channelRatio > 1
    intensAllRight = intensAllRight .* channelRatio;
    backgroundAllRight = backgroundAllRight .* channelRatio;
end
photonAllRight = round(intensAllRight);
backgroundAllRight = round(backgroundAllRight);
idAllRight = idAll(rightInd);

% Remove localization if its photon number is smaller than
% threshold
ind = photonAllRight < photonLowThreshold;
xPosAllRight(ind) = [];
yPosAllRight(ind) = [];
sigmaAllRight(ind) = [];
photonAllRight(ind) = [];
backgroundAllRight(ind) = [];
idAllRight(ind) = [];

% Store the left channel information
xPosAllLeft = xPosAll(leftInd);
yPosAllLeft = yPosAll(leftInd);
sigmaAllLeft = sigmaAll(leftInd);
intensAllLeft = intensAllVectorAll(leftInd);
backgroundAllLeft = backgroundAllVectorAll(leftInd);
if channelRatio <= 1
    intensAllLeft = intensAllLeft ./ channelRatio;
    backgroundAllLeft = backgroundAllLeft ./ channelRatio;
end
photonAllLeft = round(intensAllLeft); % round! photon number
backgroundAllLeft = round(backgroundAllLeft);
idAllLeft = idAll(leftInd);

% Remove localization if its photon number is smaller than
% threshold
ind = photonAllLeft < photonLowThreshold;
xPosAllLeft(ind) = [];
yPosAllLeft(ind) = [];
sigmaAllLeft(ind) = [];
photonAllLeft(ind) = [];
backgroundAllLeft(ind) = [];
idAllLeft(ind) = [];

% Initialize the memory
rightPointPairSMAll = [];
leftPointPairSMAll = [];
posPairSMAll = [];
sigmaRightPairSMAll = [];
sigmaLeftPairSMAll = [];
photonNumPairSMAll = [];
backgroundPairSMAll = [];
idPairSMAll = [];
precPairSMAll = [];

posUnpairedSMRightAll = [];
posUnpairedSMLeftAll = [];
sigmaUnpairedSMRightAll = [];
sigmaUnpairedSMLeftAll = [];
photonNumUnpairedSMRightAll = [];
photonNumUnpairedSMLeftAll = [];
backgroundUnpairedSMRightAll = [];
backgroundUnpairedSMLeftAll = [];
idUnpairedSMRightAll = [];
idUnpairedSMLeftAll = [];
precUnpairedSMRightAll = [];
precUnpairedSMLeftAll = [];

precRightAll = [];
precLeftAll = [];

posPairUnpairOld = [];
photonNumPairUnpairOld = [];
precPairUnpairOld = [];
onTimeOld = [];

posPairUnpairCorrAll = [];
photonNumPairUnpairCorrAll = [];
precPairUnpairCorrAll = [];
idPairUnpairCorrAll = [];
onTimeCorrAll = [];

corrPos = []; % store all correlated localizations

for h = startFrame:(startFrame + totFrame - 1)
    indIdRight = idAllRight == h;
    indIdLeft = idAllLeft == h;

    xPosRight = xPosAllRight(indIdRight);
    yPosRight = yPosAllRight(indIdRight);
    sigmaRight = sigmaAllRight(indIdRight);
    photonRight = photonAllRight(indIdRight);
    backgroundRight = backgroundAllRight(indIdRight);

    xPosLeft = xPosAllLeft(indIdLeft);
    yPosLeft = yPosAllLeft(indIdLeft);
    sigmaLeft = sigmaAllLeft(indIdLeft);
    photonLeft = photonAllLeft(indIdLeft);
    backgroundLeft = backgroundAllLeft(indIdLeft);

    % Calculate localization precision based on captured photon number
    % (Least square)
    tauRight = 2 * pi * backgroundRight .* (sigmaRight.^2 + pixelSize^2 / 12) ./ photonRight ./ (pixelSize^2);
    precRight = (sigmaRight.^2 + pixelSize^2 / 12) ./ photonRight .* (16 / 9 + 4 * tauRight);
    % precRight = (sigmaRight.^2+pixelSize^2/12)./photonRight.*(1+4*tauRight+sqrt(2*tauRight./(1+4*tauRight)));
    precRight = sqrt(precRight);
    precRightAll = [precRightAll, precRight];

    tauLeft = 2 * pi * backgroundLeft .* (sigmaLeft.^2 + pixelSize^2 / 12) ./ photonLeft ./ (pixelSize^2);
    precLeft = (sigmaLeft.^2 + pixelSize^2 / 12) ./ photonLeft .* (16 / 9 + 4 * tauLeft);
    % precLeft = (sigmaLeft.^2+pixelSize^2/12)./photonLeft.*(1+4*tauLeft+sqrt(2*tauLeft./(1+4*tauLeft)));
    precLeft = sqrt(precLeft);
    precLeftAll = [precLeftAll, precLeft];

    xPosRight = xPosRight - pixelSize * boundaryXPosSM;
    xPosLeft = -xPosLeft + pixelSize * boundaryXPosSM;

    [xPosLeftTrans, yPosLeftTrans] = transformPointsInverse(tformSM, xPosLeft.', yPosLeft.');
    xPosLeftTrans = xPosLeftTrans.';
    yPosLeftTrans = yPosLeftTrans.';

    unpairedIndSMRight = ones(1, length(xPosRight));
    unpairedIndSMLeft = ones(1, length(xPosLeft));

    if ~isempty(xPosRight) && ~isempty(xPosLeft)
        %     imgAxisX = (1:boundaryXPosSM)-0.5;
        %     imgAxisY = (1:imageSizeVSM)-0.5;
        %     imgAxisX = imgAxisX.*pixelSize;
        %     imgAxisY = imgAxisY.*pixelSize;
        %     figure('Position',P1);
        %     hold on;
        %     scatter(xPosRight,yPosRight,'LineWidth',1)
        %     scatter(xPosLeftTrans,yPosLeftTrans,'LineWidth',1)
        %     axis image ij
        %     axis([imgAxisX(1) imgAxisX(end) imgAxisY(1) imgAxisY(end)])
        %     legend('posRight','transformed posLeft')
        %     xlabel('x position [nm]','FontSize',24)
        %     ylabel('y position [nm]','FontSize',24)
        %     title('Geometric transformed images','FontSize',24)
        %     set(gca,'FontSize',20)
        %     hold off

        rightMeshX = repmat(xPosRight, length(xPosLeftTrans), 1);
        rightMeshY = repmat(yPosRight, length(yPosLeftTrans), 1);
        leftTransMeshX = repmat(xPosLeftTrans.', 1, length(xPosRight));
        leftTransMeshY = repmat(yPosLeftTrans.', 1, length(yPosRight));

        % calculate Euclidean distance between all SMs in two channels
        allVectorDis = hypot(rightMeshX-leftTransMeshX, rightMeshY-leftTransMeshY);

        % Changed precision margin from simple sum to element-wise max
        precMatrix = max(repmat(precRight, length(precLeft), 1), repmat(precLeft.', 1, length(precRight)));

        % Tight condition (sigma)
        %candInd = allVectorDis < (precMatrix+matchErrorMean+matchErrorStd);
        % Loose condition (3*sigma)
        candInd = allVectorDis < (3 * precMatrix + matchErrorMeanSM + 3 * matchErrorStdSM);

        while sum(sum(candInd)) > 0
            candidate = nan(size(allVectorDis));
            candidate(candInd) = allVectorDis(candInd);
            [miniR, minIndR] = nanmin(candidate, [], 2);
            [mini, minIndC] = nanmin(miniR);
            pairIndSM = [minIndR(minIndC); minIndC];
            rightPointPairSM = [xPosRight(pairIndSM(1)); yPosRight(pairIndSM(1))];
            rightPointPairSMAll = [rightPointPairSMAll, rightPointPairSM];
            leftPointPairSM = [xPosLeftTrans(pairIndSM(2)); yPosLeftTrans(pairIndSM(2))];
            leftPointPairSMAll = [leftPointPairSMAll, leftPointPairSM];
            posPairSM = mean([rightPointPairSM, leftPointPairSM], 2);
            posPairSMAll = [posPairSMAll, posPairSM];
            sigmaRightPairSMAll = [sigmaRightPairSMAll, sigmaRight(pairIndSM(1))];
            sigmaLeftPairSMAll = [sigmaLeftPairSMAll, sigmaLeft(pairIndSM(2))];
            photonNumPairSM = photonRight(pairIndSM(1)) + photonLeft(pairIndSM(2));
            photonNumPairSMAll = [photonNumPairSMAll, photonNumPairSM];
            backgroundPairSM = backgroundRight(pairIndSM(1)) + backgroundLeft(pairIndSM(2));
            backgroundPairSMAll = [backgroundPairSMAll, backgroundPairSM];
            idPairSMAll = [idPairSMAll, h];
            % precPairSM = mean([precRight(pairIndSM(1)) precLeft(pairIndSM(2))]);
            precPairSM = [precRight(pairIndSM(1)); precLeft(pairIndSM(2))];
            precPairSMAll = [precPairSMAll, precPairSM];
            unpairedIndSMRight(pairIndSM(1)) = 0;
            unpairedIndSMLeft(pairIndSM(2)) = 0;
            candInd(minIndC, :) = 0;
            candInd(:, minIndR(minIndC)) = 0;
        end
    end

    if sum([unpairedIndSMLeft, unpairedIndSMRight]) ~= 0
        % Read the images for photon calculation of unpaired localizations
        capturedImage = tiffread2([dataPlaceSM, '\', dataFolderInfo(h).name], 1, 1);
        capturedImage = double(capturedImage.data);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-process for unpaired localization in right channel
    posUnpairedSMRight = [];
    photonNumUnpairedSMRight = [];
    precUnpairedSMRight = [];
    if sum(unpairedIndSMRight) ~= 0
        unpairedIndSMRight = logical(unpairedIndSMRight);

        [xPosUnpairedRight2Left, yPosUnpairedRight2Left] = ...
            transformPointsInverse(tformSMRight2Left, xPosRight(unpairedIndSMRight).', yPosRight(unpairedIndSMRight).');
        xPosUnpairedRight2Left = xPosUnpairedRight2Left.';
        yPosUnpiredRight2Left = yPosUnpairedRight2Left.';
        xPosUnpairedRight2Left = -xPosUnpairedRight2Left + pixelSize * boundaryXPosSM;
        % Subtraction of background
        % insideMatrix: contains all emitter candidates' PSF information
        insideMatrix = ...
            capturedImage(reshape(repmat(ceil(yPosUnpairedRight2Left ./ pixelSize), [squareSV, 1]), [1, (squareSV) * length(yPosUnpairedRight2Left)])+repmat(-(squareSV - 1) / 2:(squareSV - 1) / 2, [1, length(yPosUnpairedRight2Left)]), ...
            reshape(repmat(ceil(xPosUnpairedRight2Left ./ pixelSize), [squareSH, 1]), [1, (squareSH) * length(xPosUnpairedRight2Left)])+repmat(-(squareSH - 1) / 2:(squareSH - 1) / 2, [1, length(xPosUnpairedRight2Left)]) ...
            );
        insideIndex = kron(eye(length(xPosUnpairedRight2Left)), ones(squareSV, squareSH));
        insideMatrix(~logical(insideIndex)) = 0;
        insideVector = sum(insideMatrix, 1);
        insideVector = reshape(insideVector, [squareSH, length(xPosUnpairedRight2Left)]);
        insideVector = sum(insideVector, 1);

        % largeMatrix: contains all emitter candidates' PSF information (larger region)
        largeMatrix = ...
            capturedImage(reshape(repmat(ceil(yPosUnpairedRight2Left ./ pixelSize), [squareSV + 2, 1]), [1, (squareSV + 2) * length(yPosUnpairedRight2Left)])+repmat(-(squareSV + 2 - 1) / 2:(squareSV + 2 - 1) / 2, [1, length(yPosUnpairedRight2Left)]), ...
            reshape(repmat(ceil(xPosUnpairedRight2Left ./ pixelSize), [squareSH + 2, 1]), [1, (squareSH + 2) * length(xPosUnpairedRight2Left)])+repmat(-(squareSH + 2 - 1) / 2:(squareSH + 2 - 1) / 2, [1, length(xPosUnpairedRight2Left)]) ...
            );
        largeIndex = kron(eye(length(xPosUnpairedRight2Left)), ones(squareSV + 2, squareSH + 2));
        largeMatrix(~logical(largeIndex)) = 0;
        largeVector = sum(largeMatrix, 1);
        largeVector = reshape(largeVector, [squareSH + 2, length(xPosUnpairedRight2Left)]);
        largeVector = sum(largeVector, 1);

        backgroundVector = largeVector - insideVector;
        backgroundVector = backgroundVector ./ (squareSV * 2 + (squareSH + 2) * 2);
        % backgroundVector = round(backgroundVector);

        % backgroundMatrix: contains background information at all emitters
        % candidates' positions
        backgroundMatrix = reshape(repmat(backgroundVector, [length(xPosUnpairedRight2Left) * squareSV * squareSH, 1]), [length(xPosUnpairedRight2Left) * squareSV, length(xPosUnpairedRight2Left) * squareSH]);
        backgroundMatrix(~logical(insideIndex)) = 0;

        % Subtract background from inside box
        % intensMatrix: contains all emitter candidates' PSF information after
        % background subtraction
        intensMatrix = insideMatrix - backgroundMatrix;
        % intensMatrix(intensMatrix < 0) = 0;

        intensAllVector = sum(intensMatrix, 1);
        intensAllVector = reshape(intensAllVector, [squareSH, length(xPosUnpairedRight2Left)]); % not rounded yet
        %intensAllVector = round(sum(intensAllVector,1).*conversionF);
        intensAllVector = sum(intensAllVector, 1) .* conversionF;
        backgroundAllVector = backgroundVector .* conversionF;

        % Reset negative intensities as 0
        intensAllVector(intensAllVector < 0) = 0;

        intensUnpairedRight2Left = intensAllVector;
        if channelRatio <= 1
            intensUnpairedRight2Left = intensUnpairedRight2Left ./ channelRatio;
            backgroundAllVector = backgroundAllVector ./ channelRatio;
        end

        photonUnpairedRight2Left = round(intensUnpairedRight2Left);
        backgroundUnpairedRight2Left = round(backgroundAllVector);

        posUnpairedSMRight = [xPosRight(unpairedIndSMRight); yPosRight(unpairedIndSMRight)];
        sigmaUnpairedSMRight = sigmaRight(unpairedIndSMRight);
        photonNumUnpairedSMRight = photonRight(unpairedIndSMRight) + photonUnpairedRight2Left;
        backgroundUnpairedSMRight = backgroundRight(unpairedIndSMRight) + backgroundUnpairedRight2Left;
        idUnpairedSMRight = h * ones(1, sum(unpairedIndSMRight));
        precUnpairedSMRight = precRight(unpairedIndSMRight);

        ind = photonNumUnpairedSMRight < unpairedPhotonLowThreshold;
        posUnpairedSMRight(:, ind) = [];
        sigmaUnpairedSMRight(ind) = [];
        photonNumUnpairedSMRight(ind) = [];
        backgroundUnpairedSMRight(ind) = [];
        idUnpairedSMRight(ind) = [];
        precUnpairedSMRight(ind) = [];

        posUnpairedSMRightAll = [posUnpairedSMRightAll, posUnpairedSMRight];
        sigmaUnpairedSMRightAll = [sigmaUnpairedSMRightAll, sigmaUnpairedSMRight];
        photonNumUnpairedSMRightAll = [photonNumUnpairedSMRightAll, photonNumUnpairedSMRight];
        backgroundUnpairedSMRightAll = [backgroundUnpairedSMRightAll, backgroundUnpairedSMRight];
        idUnpairedSMRightAll = [idUnpairedSMRightAll, idUnpairedSMRight];
        precUnpairedSMRightAll = [precUnpairedSMRightAll, precUnpairedSMRight];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-process for unpaired localization in left channel
    posUnpairedSMLeft = [];
    photonNumUnpairedSMLeft = [];
    precUnpairedSMLeft = [];
    if sum(unpairedIndSMLeft) ~= 0
        unpairedIndSMLeft = logical(unpairedIndSMLeft);

        [xPosUnpairedLeft2Right, yPosUnpairedLeft2Right] = ...
            transformPointsInverse(tformSM, xPosLeft(unpairedIndSMLeft).', yPosLeft(unpairedIndSMLeft).');
        xPosUnpairedLeft2Right = xPosUnpairedLeft2Right.';
        yPosUnpairedLeft2Right = yPosUnpairedLeft2Right.';
        xPosUnpairedLeft2Right = xPosUnpairedLeft2Right + pixelSize * boundaryXPosSM;

        % Subtraction of background
        % insideMatrix: contains all emitter candidates' PSF information
        insideMatrix = ...
            capturedImage(reshape(repmat(ceil(yPosUnpairedLeft2Right ./ pixelSize), [squareSV, 1]), [1, (squareSV) * length(yPosUnpairedLeft2Right)])+repmat(-(squareSV - 1) / 2:(squareSV - 1) / 2, [1, length(yPosUnpairedLeft2Right)]), ...
            reshape(repmat(ceil(xPosUnpairedLeft2Right ./ pixelSize), [squareSH, 1]), [1, (squareSH) * length(xPosUnpairedLeft2Right)])+repmat(-(squareSH - 1) / 2:(squareSH - 1) / 2, [1, length(xPosUnpairedLeft2Right)]) ...
            );
        insideIndex = kron(eye(length(xPosUnpairedLeft2Right)), ones(squareSV, squareSH));
        insideMatrix(~logical(insideIndex)) = 0;
        insideVector = sum(insideMatrix, 1);
        insideVector = reshape(insideVector, [squareSH, length(xPosUnpairedLeft2Right)]);
        insideVector = sum(insideVector, 1);

        % largeMatrix: contains all emitter candidates' PSF information (larger region)
        largeMatrix = ...
            capturedImage(reshape(repmat(ceil(yPosUnpairedLeft2Right ./ pixelSize), [squareSV + 2, 1]), [1, (squareSV + 2) * length(yPosUnpairedLeft2Right)])+repmat(-(squareSV + 2 - 1) / 2:(squareSV + 2 - 1) / 2, [1, length(yPosUnpairedLeft2Right)]), ...
            reshape(repmat(ceil(xPosUnpairedLeft2Right ./ pixelSize), [squareSH + 2, 1]), [1, (squareSH + 2) * length(xPosUnpairedLeft2Right)])+repmat(-(squareSH + 2 - 1) / 2:(squareSH + 2 - 1) / 2, [1, length(xPosUnpairedLeft2Right)]) ...
            );
        largeIndex = kron(eye(length(xPosUnpairedLeft2Right)), ones(squareSV + 2, squareSH + 2));
        largeMatrix(~logical(largeIndex)) = 0;
        largeVector = sum(largeMatrix, 1);
        largeVector = reshape(largeVector, [squareSH + 2, length(xPosUnpairedLeft2Right)]);
        largeVector = sum(largeVector, 1);

        backgroundVector = largeVector - insideVector;
        backgroundVector = backgroundVector ./ (squareSV * 2 + (squareSH + 2) * 2);
        % backgroundVector = round(backgroundVector);

        % backgroundMatrix: contains background information at all emitters
        % candidates' positions
        backgroundMatrix = reshape(repmat(backgroundVector, [length(xPosUnpairedLeft2Right) * squareSV * squareSH, 1]), [length(xPosUnpairedLeft2Right) * squareSV, length(xPosUnpairedLeft2Right) * squareSH]);
        backgroundMatrix(~logical(insideIndex)) = 0;

        % Subtract background from inside box
        % intensMatrix: contains all emitter candidates' PSF information after
        % background subtraction
        intensMatrix = insideMatrix - backgroundMatrix;
        % intensMatrix(intensMatrix < 0) = 0;

        intensAllVector = sum(intensMatrix, 1);
        intensAllVector = reshape(intensAllVector, [squareSH, length(xPosUnpairedLeft2Right)]); % not rounded yet
        %intensAllVector = round(sum(intensAllVector,1).*conversionF);
        intensAllVector = sum(intensAllVector, 1) .* conversionF;
        backgroundAllVector = backgroundVector .* conversionF;

        % Reset negative intensities as 0
        intensAllVector(intensAllVector < 0) = 0;

        intensUnpairedLeft2Right = intensAllVector;
        if channelRatio > 1
            intensUnpairedLeft2Right = intensUnpairedLeft2Right .* channelRatio;
            backgroundAllVector = backgroundAllVector .* channelRatio;
        end

        photonUnpairedLeft2Right = round(intensUnpairedLeft2Right);
        backgroundUnpairedLeft2Right = round(backgroundAllVector);

        posUnpairedSMLeft = [xPosLeftTrans(unpairedIndSMLeft); yPosLeftTrans(unpairedIndSMLeft)];
        sigmaUnpairedSMLeft = sigmaLeft(unpairedIndSMLeft);
        photonNumUnpairedSMLeft = photonUnpairedLeft2Right + photonLeft(unpairedIndSMLeft);
        backgroundUnpairedSMLeft = backgroundUnpairedLeft2Right + backgroundLeft(unpairedIndSMLeft);
        idUnpairedSMLeft = h * ones(1, sum(unpairedIndSMLeft));
        precUnpairedSMLeft = precLeft(unpairedIndSMLeft);

        ind = photonNumUnpairedSMLeft < unpairedPhotonLowThreshold;
        posUnpairedSMLeft(:, ind) = [];
        sigmaUnpairedSMLeft(ind) = [];
        photonNumUnpairedSMLeft(ind) = [];
        backgroundUnpairedSMLeft(ind) = [];
        idUnpairedSMLeft(ind) = [];
        precUnpairedSMLeft(ind) = [];

        posUnpairedSMLeftAll = [posUnpairedSMLeftAll, posUnpairedSMLeft];
        sigmaUnpairedSMLeftAll = [sigmaUnpairedSMLeftAll, sigmaUnpairedSMLeft];
        photonNumUnpairedSMLeftAll = [photonNumUnpairedSMLeftAll, photonNumUnpairedSMLeft];
        backgroundUnpairedSMLeftAll = [backgroundUnpairedSMLeftAll, backgroundUnpairedSMLeft];
        idUnpairedSMLeftAll = [idUnpairedSMLeftAll, idUnpairedSMLeft];
        precUnpairedSMLeftAll = [precUnpairedSMLeftAll, precUnpairedSMLeft];
    end

    %%

    % Correlate localization across multiple frames
    posPairUnpairNew = posPairSMAll(:, idPairSMAll == h);
    photonNumPairUnpairNew = photonNumPairSMAll(idPairSMAll == h);
    precPairUnpairNew = max(precPairSMAll(:, idPairSMAll == h), [], 1);
    if exist('posUnpairedSMRight', 'var') == 1
        posPairUnpairNew = [posPairUnpairNew, posUnpairedSMRight];
        photonNumPairUnpairNew = [photonNumPairUnpairNew, photonNumUnpairedSMRight];
        precPairUnpairNew = [precPairUnpairNew, precUnpairedSMRight];
    end
    if exist('posUnpairedSMLeft', 'var') == 1
        posPairUnpairNew = [posPairUnpairNew, posUnpairedSMLeft];
        photonNumPairUnpairNew = [photonNumPairUnpairNew, photonNumUnpairedSMLeft];
        precPairUnpairNew = [precPairUnpairNew, precUnpairedSMLeft];
    end
    %     backgroundPairUnpairNew = [backgroundPairSMAll(idPairSMAll==h) backgroundUnpairedSMRight backgroundUnpairedSMLeft];

    unpairedIndSMOld = ones(1, length(photonNumPairUnpairOld));
    unpairedIndSMNew = ones(1, length(photonNumPairUnpairNew));

    posPairUnpairOldCorr = [];
    photonNumPairUnpairOldCorr = [];
    precPairUnpairOldCorr = [];
    onTimeOldCorr = [];
    if ~isempty(photonNumPairUnpairOld) && ~isempty(photonNumPairUnpairNew)
        posOldMeshX = repmat(posPairUnpairOld(1, :), length(posPairUnpairNew(1, :)), 1);
        posOldMeshY = repmat(posPairUnpairOld(2, :), length(posPairUnpairNew(2, :)), 1);
        posNewMeshX = repmat(posPairUnpairNew(1, :).', 1, length(posPairUnpairOld(1, :)));
        posNewMeshY = repmat(posPairUnpairNew(2, :).', 1, length(posPairUnpairOld(2, :)));
        % calculate Euclidean distance between all SMs in two successive
        % frames
        allVectorDis = hypot(posOldMeshX-posNewMeshX, posOldMeshY-posNewMeshY);

        precMatrix = max(repmat(precPairUnpairOld, length(precPairUnpairNew), 1), repmat(precPairUnpairNew.', 1, length(precPairUnpairOld)));

        candInd = allVectorDis < (3 * precMatrix);

        while sum(sum(candInd)) > 0
            candidate = nan(size(allVectorDis));
            candidate(candInd) = allVectorDis(candInd);
            [miniR, minIndR] = nanmin(candidate, [], 2);
            [mini, minIndC] = nanmin(miniR);
            pairIndSM = [minIndR(minIndC); minIndC];
            corrPos = [corrPos, [posPairUnpairOld(:, pairIndSM(1)); posPairUnpairNew(:, pairIndSM(2))]];
            posPairUnpairOld(:, pairIndSM(1)) = posPairUnpairNew(:, pairIndSM(2));
            posPairUnpairOldCorr = [posPairUnpairOldCorr, posPairUnpairOld(:, pairIndSM(1))];
            photonNumPairUnpairOld(pairIndSM(1)) = photonNumPairUnpairOld(pairIndSM(1)) + photonNumPairUnpairNew(pairIndSM(2));
            photonNumPairUnpairOldCorr = [photonNumPairUnpairOldCorr, photonNumPairUnpairOld(pairIndSM(1))];
            precPairUnpairOld(pairIndSM(1)) = precPairUnpairNew(pairIndSM(2));
            precPairUnpairOldCorr = [precPairUnpairOldCorr, precPairUnpairOld(pairIndSM(1))];
            onTimeOld(pairIndSM(1)) = onTimeOld(pairIndSM(1)) + 1;
            onTimeOldCorr = [onTimeOldCorr, onTimeOld(pairIndSM(1))];
            unpairedIndSMOld(pairIndSM(1)) = 0;
            unpairedIndSMNew(pairIndSM(2)) = 0;
            candInd(minIndC, :) = 0;
            candInd(:, minIndR(minIndC)) = 0;
        end
    end

    % store information of localizations that are not correlated with localizations in the new frame
    posPairUnpairCorrAll = [posPairUnpairCorrAll, posPairUnpairOld(:, logical(unpairedIndSMOld))];
    photonNumPairUnpairCorrAll = [photonNumPairUnpairCorrAll, photonNumPairUnpairOld(logical(unpairedIndSMOld))];
    precPairUnpairCorrAll = [precPairUnpairCorrAll, precPairUnpairOld(logical(unpairedIndSMOld))];
    idPairUnpairCorrAll = [idPairUnpairCorrAll, (h - 1) * ones(1, sum(unpairedIndSMOld))];
    onTimeCorrAll = [onTimeCorrAll, onTimeOld(logical(unpairedIndSMOld))];

    % replace "old" frame information with correlated localizations and the new frame information
    posPairUnpairOld = [posPairUnpairOldCorr, posPairUnpairNew(:, logical(unpairedIndSMNew))];
    photonNumPairUnpairOld = [photonNumPairUnpairOldCorr, photonNumPairUnpairNew(logical(unpairedIndSMNew))];
    precPairUnpairOld = [precPairUnpairOldCorr, precPairUnpairNew(logical(unpairedIndSMNew))];
    onTimeOld = [onTimeOldCorr, ones(1, sum(unpairedIndSMNew))];

end
% store information of localizations in the last frame
posPairUnpairCorrAll = [posPairUnpairCorrAll, posPairUnpairOld];
photonNumPairUnpairCorrAll = [photonNumPairUnpairCorrAll, photonNumPairUnpairOld];
precPairUnpairCorrAll = [precPairUnpairCorrAll, precPairUnpairOld];
idPairUnpairCorrAll = [idPairUnpairCorrAll, h * ones(1, size(posPairUnpairOld, 2))];
onTimeCorrAll = [onTimeCorrAll, onTimeOld];
onTimeCorrAll = onTimeCorrAll .* expTime;

xPosRightAll = [posPairSMAll(1, :), posUnpairedSMRightAll(1, :)];
yPosRightAll = [posPairSMAll(2, :), posUnpairedSMRightAll(2, :)];
sigmaRightAll = [sigmaRightPairSMAll, sigmaUnpairedSMRightAll];
idRightAll = [idPairSMAll, idUnpairedSMRightAll];

xPosLeftAll = [posPairSMAll(1, :), posUnpairedSMLeftAll(1, :)];
yPosLeftAll = [posPairSMAll(2, :), posUnpairedSMLeftAll(2, :)];
sigmaLeftAll = [sigmaLeftPairSMAll, sigmaUnpairedSMLeftAll];
idLeftAll = [idPairSMAll, idUnpairedSMLeftAll];

xPosPairUnpair = [posPairSMAll(1, :), posUnpairedSMRightAll(1, :), posUnpairedSMLeftAll(1, :)];
yPosPairUnpair = [posPairSMAll(2, :), posUnpairedSMRightAll(2, :), posUnpairedSMLeftAll(2, :)];
photonNumPairUnpair = [photonNumPairSMAll, photonNumUnpairedSMRightAll, photonNumUnpairedSMLeftAll];
backgroundPairUnapir = [backgroundPairSMAll, backgroundUnpairedSMRightAll, backgroundUnpairedSMLeftAll];
precPairUnpair = [precPairSMAll, [precUnpairedSMRightAll; nan(1, length(precUnpairedSMRightAll))], ...
    [nan(1, length(precUnpairedSMLeftAll)); precUnpairedSMLeftAll]];
idPairUnpair = [idPairSMAll, idUnpairedSMRightAll, idUnpairedSMLeftAll];

% Get zoomed region information
% Matched up emitters in zoomed region
indZoom = zeros(size(zoomRegion, 2), size(posPairSMAll, 2));
indZoom(1, :) = posPairSMAll(1, :) > zoomRegion(1);
indZoom(2, :) = posPairSMAll(1, :) < zoomRegion(2);
indZoom(3, :) = posPairSMAll(2, :) > zoomRegion(3);
indZoom(4, :) = posPairSMAll(2, :) < zoomRegion(4);
indZoomWhole = sum(indZoom, 1) > 3;
xPosPairZoom = posPairSMAll(1, indZoomWhole);
yPosPairZoom = posPairSMAll(2, indZoomWhole);
photonNumPairZoom = photonNumPairSMAll(indZoomWhole);
backgroundPairZoom = backgroundPairSMAll(indZoomWhole);
precPairZoom = precPairSMAll(:, indZoomWhole);
idPairZoom = idPairSMAll(indZoomWhole);

% Unpaired emitters in right
indZoom = zeros(size(zoomRegion, 2), size(posUnpairedSMRightAll, 2));
indZoom(1, :) = posUnpairedSMRightAll(1, :) > zoomRegion(1);
indZoom(2, :) = posUnpairedSMRightAll(1, :) < zoomRegion(2);
indZoom(3, :) = posUnpairedSMRightAll(2, :) > zoomRegion(3);
indZoom(4, :) = posUnpairedSMRightAll(2, :) < zoomRegion(4);
indZoomWhole = sum(indZoom, 1) > 3;
xPosRightUnpairedZoom = posUnpairedSMRightAll(1, indZoomWhole);
yPosRightUnpairedZoom = posUnpairedSMRightAll(2, indZoomWhole);
photonNumUnpairedRightZoom = photonNumUnpairedSMRightAll(indZoomWhole);
backgroundUnpairedRightZoom = backgroundUnpairedSMRightAll(indZoomWhole);
precUnpairedRightZoom = precUnpairedSMRightAll(indZoomWhole);
idUnpairedRightZoom = idUnpairedSMRightAll(indZoomWhole);

% Unpaired emitters in left
indZoom = zeros(size(zoomRegion, 2), size(posUnpairedSMLeftAll, 2));
indZoom(1, :) = posUnpairedSMLeftAll(1, :) > zoomRegion(1);
indZoom(2, :) = posUnpairedSMLeftAll(1, :) < zoomRegion(2);
indZoom(3, :) = posUnpairedSMLeftAll(2, :) > zoomRegion(3);
indZoom(4, :) = posUnpairedSMLeftAll(2, :) < zoomRegion(4);
indZoomWhole = sum(indZoom, 1) > 3;
xPosLeftUnpairedZoom = posUnpairedSMLeftAll(1, indZoomWhole);
yPosLeftUnpairedZoom = posUnpairedSMLeftAll(2, indZoomWhole);
photonNumUnpairedLeftZoom = photonNumUnpairedSMLeftAll(indZoomWhole);
backgroundUnpairedLeftZoom = backgroundUnpairedSMLeftAll(indZoomWhole);
precUnpairedLeftZoom = precUnpairedSMLeftAll(indZoomWhole);
idUnpairedLeftZoom = idUnpairedSMLeftAll(indZoomWhole);

% Localized all emitters in right
xPosRightAllZoom = [xPosPairZoom, xPosRightUnpairedZoom];
yPosRightAllZoom = [yPosPairZoom, yPosRightUnpairedZoom];
idRightAllZoom = [idPairZoom, idUnpairedRightZoom];

% Localized all emitters in left
xPosLeftAllZoom = [xPosPairZoom, xPosLeftUnpairedZoom];
yPosLeftAllZoom = [yPosPairZoom, yPosLeftUnpairedZoom];
idLeftAllZoom = [idPairZoom, idUnpairedLeftZoom];

% Correlated bursts
indZoom = zeros(size(zoomRegion, 2), size(posPairUnpairCorrAll, 2));
indZoom(1, :) = posPairUnpairCorrAll(1, :) > zoomRegion(1);
indZoom(2, :) = posPairUnpairCorrAll(1, :) < zoomRegion(2);
indZoom(3, :) = posPairUnpairCorrAll(2, :) > zoomRegion(3);
indZoom(4, :) = posPairUnpairCorrAll(2, :) < zoomRegion(4);
indZoomWhole = sum(indZoom, 1) > 3;
posPairUnpairCorrAllZoom = posPairUnpairCorrAll(:, indZoomWhole);
photonNumPairUnpairCorrAllZoom = photonNumPairUnpairCorrAll(indZoomWhole);
precPairUnpairCorrAllZoom = precPairUnpairCorrAll(indZoomWhole);
idPairUnpairCorrAllZoom = idPairUnpairCorrAll(indZoomWhole);
onTimeCorrAllZoom = onTimeCorrAll(indZoomWhole);

% Get background region information
% Matched up emitters in background region
indBackground = zeros(size(zoomRegion, 2), size(posPairSMAll, 2));
indBackground(1, :) = posPairSMAll(1, :) < zoomRegion(1);
indBackground(2, :) = posPairSMAll(1, :) > zoomRegion(2);
indBackground(3, :) = posPairSMAll(2, :) < zoomRegion(3);
indBackground(4, :) = posPairSMAll(2, :) > zoomRegion(4);
indBackgroundWhole = sum(indBackground, 1) > 0;
xPosPairBackground = posPairSMAll(1, indBackgroundWhole);
yPosPairBackground = posPairSMAll(2, indBackgroundWhole);
idPairBackground = idPairSMAll(indBackgroundWhole);

% Unpaired emitters in right
indBackground = zeros(size(zoomRegion, 2), size(posUnpairedSMRightAll, 2));
indBackground(1, :) = posUnpairedSMRightAll(1, :) < zoomRegion(1);
indBackground(2, :) = posUnpairedSMRightAll(1, :) > zoomRegion(2);
indBackground(3, :) = posUnpairedSMRightAll(2, :) < zoomRegion(3);
indBackground(4, :) = posUnpairedSMRightAll(2, :) > zoomRegion(4);
indBackgroundWhole = sum(indBackground, 1) > 0;
xPosRightUnpairedBackground = posUnpairedSMRightAll(1, indBackgroundWhole);
yPosRightUnpairedBackground = posUnpairedSMRightAll(2, indBackgroundWhole);
photonNumUnpairedRightbackground = photonNumUnpairedSMRightAll(indBackgroundWhole);
backgroundUnpairedRightbackground = backgroundUnpairedSMRightAll(indBackgroundWhole);
idUnpairedRightBackground = idUnpairedSMRightAll(indBackgroundWhole);

% Unpaired emitters in left
indBackground = zeros(size(zoomRegion, 2), size(posUnpairedSMLeftAll, 2));
indBackground(1, :) = posUnpairedSMLeftAll(1, :) < zoomRegion(1);
indBackground(2, :) = posUnpairedSMLeftAll(1, :) > zoomRegion(2);
indBackground(3, :) = posUnpairedSMLeftAll(2, :) < zoomRegion(3);
indBackground(4, :) = posUnpairedSMLeftAll(2, :) > zoomRegion(4);
indBackgroundWhole = sum(indBackground, 1) > 0;
xPosLeftUnpairedBackground = posUnpairedSMLeftAll(1, indBackgroundWhole);
yPosLeftUnpairedBackground = posUnpairedSMLeftAll(2, indBackgroundWhole);
photonNumUnpairedLeftBackground = photonNumUnpairedSMLeftAll(indBackgroundWhole);
backgroundUnpairedLeftBackground = backgroundUnpairedSMLeftAll(indBackgroundWhole);
idUnpairedLeftBackground = idUnpairedSMLeftAll(indBackgroundWhole);

% Localized all emitters in right
xPosRightAllBackground = [xPosPairBackground, xPosRightUnpairedBackground];
yPosRightAllBackground = [yPosPairBackground, yPosRightUnpairedBackground];
idRightAllBackground = [idPairBackground, idUnpairedRightBackground];

% Localized all emitters in left
xPosLeftAllBackground = [xPosPairBackground, xPosLeftUnpairedBackground];
yPosLeftAllBackground = [yPosPairBackground, yPosLeftUnpairedBackground];
idLeftAllBackground = [idPairBackground, idUnpairedLeftBackground];

PhotonArray = [photonNumPairZoom, photonNumUnpairedRightZoom, photonNumUnpairedLeftZoom];
BackgroundArray = [backgroundPairZoom, backgroundUnpairedRightZoom, backgroundUnpairedLeftZoom];
PrecArray = [reshape(precPairZoom, [1, size(precPairZoom, 1) * size(precPairZoom, 2)]), precUnpairedRightZoom, precUnpairedLeftZoom];

% Check correlation
% test figure transform with correction (paired localizations)
figure('Position', P1);
hold on;
for k = 1:size(corrPos, 2)
    plot([corrPos(1, k), corrPos(3, k)], [corrPos(2, k), corrPos(4, k)], '-r');
end
axis image ij
axis(zoomRegion)
% axis(workingRegion)
grid on
grid minor
xlabel('x position [nm]')
ylabel('y position [nm]')
title('Correlation check')
hold off
hCorr1 = gcf;

corrPosDist = hypot((corrPos(1, :)-corrPos(3, :)), (corrPos(2, :) - corrPos(4, :)));
figure('Position', P1);
hold on;
hCorrDist = histogram(corrPosDist);
hCorrDist.BinWidth = binSize;
grid on
grid minor
xlabel('Distance [nm]')
title(['Correlated distance, median = ', num2str(median(corrPosDist))])
hold off
hCorr2 = gcf;

export_fig(hCorr1, [dirName, '\', saveName, ' hCorr1, lines, correlated localizations', saveFormat], '-transparent')
export_fig(hCorr2, [dirName, '\', saveName, ' hCorr2, hist, distance between correlated localizations', saveFormat], '-transparent')

% Get ROI
Xedges = (zoomRegion(1):binSize:zoomRegion(2));
Yedges = (zoomRegion(3):binSize:zoomRegion(4));
[reconstructed, ~, ~] = histcounts2(xPosPairUnpair, yPosPairUnpair, Xedges, Yedges);
reconstructed = reconstructed.';
figure('Position', P1);
hold on
imagesc(reconstructed)
N3 = reconstructed;
N3 = reshape(N3, [size(N3, 1) * size(N3, 2), 1]);
N3(N3 == 0) = [];
sNum3 = 0;
indN3 = 0;
while sNum3 < length(N3) * satPop
    indN3 = indN3 + 1;
    sNum3 = sNum3 + sum(N3 == indN3);
end
tickMinBoth = 0;
tickMaxBoth = indN3;
colormap(colorMap)
caxis([tickMinBoth, tickMaxBoth]);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
xLim = xlim;
yLim = ylim;
plotX = [xLim(2) - scaleBar / binSize - (xLim(2) - xLim(1)) / 20, xLim(2) - (xLim(2) - xLim(1)) / 20];
plotY = [yLim(2) - (yLim(2) - yLim(1)) / 100 - (yLim(2) - yLim(1)) / 20, yLim(2) - (yLim(2) - yLim(1)) / 100 - (yLim(2) - yLim(1)) / 20];
plot(plotX, plotY, '-w', 'LineWidth', 3)
title({['Localizations after two channel registration']; ['scale bar = ', num2str(scaleBar), 'nm']}, 'FontSize', 24)
hold off
hTest = gcf;

binRec = reconstructed > threROI;
binRecFill = imfill(binRec, 'holes');
CC = bwconncomp(binRecFill);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
fibROI = zeros(size(reconstructed));
fibROI(CC.PixelIdxList{idx}) = 1;
fibROIAlpha = fibROI .* 0.3;
fibROI = logical(fibROI);
fibROIBW = ~fibROI;
pixIdx = 1;
while fibROIBW(pixIdx) == 1
    pixIdx = pixIdx + 1;
end
pixIdx = pixIdx - 1;
if mod(pixIdx, size(fibROIBW, 1)) == 0
    pixR = size(fibROIBW, 1);
else
    pixR = mod(pixIdx, size(fibROIBW, 1));
end
pixC = ceil(pixIdx/size(fibROIBW, 1));
fibROIBound = bwtraceboundary(fibROIBW, [pixR, pixC], 'E', 4, Inf, 'clockwise');
xROI = fibROIBound(:, 2);
yROI = fibROIBound(:, 1);
recSat = reconstructed;
recSat(reconstructed > indN3) = indN3;
recSat = round(recSat./indN3.*length(colorMap));
I = ind2rgb(recSat, colorMap);
green = cat(3, zeros(size(recSat)), ones(size(recSat)), zeros(size(recSat)));
figure('Position', P1);
hold on
imagesc(I)
h = imagesc(green);
set(h, 'AlphaData', fibROIAlpha);
plot(xROI, yROI, 'w', 'LineWidth', 1)
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
axis image ij
title({['Selected ROI']; ['transparent green: ROI, white line: boundary']}, 'FontSize', 24)
hold off
hTest1 = gcf;

in = inpolygon(xPosPairUnpair, yPosPairUnpair, (xROI - 0.5).*binSize+zoomRegion(1), (yROI - 0.5).*binSize+zoomRegion(3));
inCorr = inpolygon(posPairUnpairCorrAll(1, :), posPairUnpairCorrAll(2, :), (xROI - 0.5).*binSize+zoomRegion(1), (yROI - 0.5).*binSize+zoomRegion(3));
pixelROI = nan(size(reconstructed));
pixelROI(fibROI) = reconstructed(fibROI);
pixelROI = reshape(pixelROI, [1, size(pixelROI, 1) * size(pixelROI, 2)]);
pixelROI(isnan(pixelROI)) = [];

figure('Position', P1);
hold on
plot((xROI-0.5).*binSize+zoomRegion(1), (yROI - 0.5).*binSize+zoomRegion(3))
plot(xPosPairUnpair(in), yPosPairUnpair(in), 'r+')
plot(xPosPairUnpair(~in), yPosPairUnpair(~in), 'bo')
grid on
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'on');
axis image ij
axis(zoomRegion)
title({['Selected ROI vs localizations']; ['red cross: localizations inside, blue circle: localizations outside']}, 'FontSize', 24)
hold off
hTest2 = gcf;

% figure('Position',P1);
% hold on
% plot((xROI-0.5).*binSize+zoomRegion(1),(yROI-0.5).*binSize+zoomRegion(3))
% plot(posPairUnpairCorrAll(1,inCorr),posPairUnpairCorrAll(2,inCorr),'r+')
% plot(posPairUnpairCorrAll(1,~inCorr),posPairUnpairCorrAll(2,~inCorr),'bo')
% axis(zoomRegion)
% grid on
% set(gca, 'XTickLabel', []);
% set(gca, 'YTickLabel', []);
% set(gca, 'Box', 'on');
% axis image ij
% hold off
% hTest3 = gcf;

% figure('Position',P1);
% hold on
% imagesc(reconstructed)
% plot(xROI,yROI)
% colormap(colorMap)
% caxis([tickMinBoth,tickMaxBoth]);
% set(gca, 'XTickLabel', []);
% set(gca, 'YTickLabel', []);
% set(gca, 'Box', 'off');
% axis image ij
% %title({['Localized emission in both channels']; ['scale bar = ' num2str(100) 'nm']},'FontSize',24)
% hold off
% hTest4 = gcf;

export_fig(hTest, [dirName, '\', saveName, ' hTest, 2Dhist, reconstruction', saveFormat], '-transparent')
export_fig(hTest1, [dirName, '\', saveName, ' hTest1, 2Dhist, check ROI vs localizations, image', saveFormat], '-transparent')
export_fig(hTest2, [dirName, '\', saveName, ' hTest2, scatter, check ROI vs localizations', saveFormat], '-transparent')

close(hCorr1, hCorr2, hTest, hTest1)

photonNumPairUnpairROI = photonNumPairUnpair(in);
backgroundPairUnpairROI = backgroundPairUnapir(in);
precPairUnpairROI = reshape(precPairUnpair(:, in), [1, size(precPairUnpair(:, in), 1) * size(precPairUnpair(:, in), 2)]);
precPairUnpairROI(isnan(precPairUnpairROI)) = [];
idPairUnpairROI = idPairUnpair(in);

photonNumPairUnpairCorrROI = photonNumPairUnpairCorrAll(inCorr);
% precPairUnpairCorrROI = precPairUnpairCorrAll(inCorr);
onTimeCorrROI = onTimeCorrAll(inCorr);
idPairUnpairCorrROI = idPairUnpairCorrAll(inCorr);

% Count localization and burst number per frame
locaNumZoom = nan(totFrame, 1);
locaNumROI = nan(totFrame, 1);
burstNumZoom = nan(totFrame, 1);
burstNumROI = nan(totFrame, 1);
hhId = 1;
for hh = startFrame:(startFrame + totFrame - 1)
    locaNumPairZoom = idPairZoom == hh;
    locaNumUnpairedRightZoom = idUnpairedRightZoom == hh;
    locaNumUnpairedLeftZoom = idUnpairedLeftZoom == hh;
    locaNumZoom(hhId) = sum(double(locaNumPairZoom)) + sum(double(locaNumUnpairedRightZoom)) + sum(double(locaNumUnpairedLeftZoom));
    locaNumROI(hhId) = sum(idPairUnpairROI == hh);
    burstNumZoom(hhId) = sum(idPairUnpairCorrAllZoom == hh);
    burstNumROI(hhId) = sum(idPairUnpairCorrROI == hh);
    hhId = hhId + 1;
end

% bin number based on specified window size
binNum = [(zoomRegion(2) - zoomRegion(1)) / binSize, (zoomRegion(4) - zoomRegion(3)) / binSize];
binNumAll = binNum(1) * binNum(2);
% ROI area
areaROI = polyarea((xROI-0.5).*binSize, (yROI - 0.5).*binSize);

% Localization number in a shifting window
locaNumZoomBinShift = nan(1, totFrame-histBinSize+1);
locaNumROIBinShift = nan(1, totFrame-histBinSize+1);
burstNumZoomBinShift = nan(1, totFrame-histBinSize+1);
burstNumROIBinShift = nan(1, totFrame-histBinSize+1);
for hh = 1:totFrame - histBinSize + 1
    locaNumZoomBinShift(hh) = sum(locaNumZoom(hh:hh + histBinSize - 1)) / binNumAll;
    locaNumROIBinShift(hh) = sum(locaNumROI(hh:hh + histBinSize - 1)) / areaROI;
    burstNumZoomBinShift(hh) = sum(burstNumZoom(hh:hh + histBinSize - 1)) / binNumAll;
    burstNumROIBinShift(hh) = sum(burstNumROI(hh:hh + histBinSize - 1)) / areaROI;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final check of figure transform with correction (histogram)
figure('Position', P1);
hold on;
histTest1 = histogram(rightPointPairSMAll(1, :)-leftPointPairSMAll(1, :));
histTest2 = histogram(rightPointPairSMAll(2, :)-leftPointPairSMAll(2, :));
histTest2.BinWidth = histTest1.BinWidth;
legend(['xDir, median:', num2str(median(rightPointPairSMAll(1, :) - leftPointPairSMAll(1, :))), 'nm'], ...
    ['yDir, median:', num2str(median(rightPointPairSMAll(2, :) - leftPointPairSMAll(2, :))), 'nm'])
xlabel('Error [nm]')
% ylabel('Frequency')
title({['Error distribution of geometric transform']; ...
    ['with corrected transform (all paired localizations for reconstruction, TRE)']})
hold off
h12 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengthRegCorr = hypot((rightPointPairSMAll(1, :)-leftPointPairSMAll(1, :)), (rightPointPairSMAll(2, :) - leftPointPairSMAll(2, :)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final check of figure transform with correction (histogram of distance)
figure('Position', P1);
hold on;
histTestDis = histogram(lengthRegCorr);
legend(['paired #: ', num2str(length(lengthRegCorr)), ', median= ', num2str(median(lengthRegCorr)), ' nm']);
xlabel('Error [nm]')
% ylabel('Frequency')
title({['Error distribution of geometric transform']; ...
    ['with corrected transform (distance, all paired localizations for reconstruction, TRE)']})
hold off
h13 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of all localized emitter number in both channel vs paired emitters
figure('Position', P1);
hold on
hRight = histogram(idRightAll);
hLeft = histogram(idLeftAll);
hId = histogram(idPairSMAll);
hRight.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hLeft.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hId.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
legend('Right channel (x-pol)', 'Left channel (y-pol)', 'matched up emissions')
xlabel('Frame number')
ylabel('Accumulated frequency per bin')
title(['Localized SM number per bin (before and after two channel matchup, bin size =', num2str(histBinSize), ')'])
hold off
h14 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of emitter number in zoomed region (both channel vs paired emitters)
figure('Position', P1);
hold on
hRightZ = histogram(idRightAllZoom);
hLeftZ = histogram(idLeftAllZoom);
hIdZ = histogram(idPairZoom);
hRightZ.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hLeftZ.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hIdZ.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
legend('Right channel (x-pol)', 'Left channel (y-pol)', 'matched up emissions')
xlabel('Frame number')
ylabel('Accumulated frequency per bin')
title(['Localized SM number per bin (Zoom, before and after two channel matchup, bin size =', num2str(histBinSize), ')'])
hold off
h14_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of emitter number in background region (both channel vs paired emitters)
figure('Position', P1);
hold on
hRightB = histogram(idRightAllBackground);
hLeftB = histogram(idLeftAllBackground);
hIdB = histogram(idPairBackground);
hRightB.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hLeftB.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
hIdB.BinEdges = startFrame:histBinSize:(startFrame + totFrame - 1);
legend('Right channel (x-pol)', 'Left channel (y-pol)', 'matched up emissions')
xlabel('Frame number')
ylabel('Accumulated frequency per bin')
title(['Localized SM number per bin (Background, before and after two channel matchup, bin size =', num2str(histBinSize), ')'])
hold off
h14_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

export_fig(h12, [dirName, '\', saveName, ' h12, hist, x,y bias of paired localizations, for reconstruction', saveFormat], '-transparent')
export_fig(h13, [dirName, '\', saveName, ' h13, hist, dist between paired localizations, for reconstruction', saveFormat], '-transparent')
export_fig(h14, [dirName, '\', saveName, ' h14, hist, loca per bin', saveFormat], '-transparent')
export_fig(h14_1, [dirName, '\', saveName, ' h14_1, hist, loca per bin, zoom', saveFormat], '-transparent')
export_fig(h14_2, [dirName, '\', saveName, ' h14_2, hist, loca per bin, background', saveFormat], '-transparent')
close(h12, h13, h14, h14_1, h14_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Position',P1);
% frameTick = histBinSize/2:histBinSize:totFrames-histBinSize/2;
% plot(frameTick,burstNumZoomBinPerArea)
% xlabel('Frame')
% ylabel('Burst num [/unit area]')
% title(['Accumulated localization in every ' num2str(histBinSize) ' frames per unit area (unit area = ' num2str(binSize) '\times' num2str(binSize) 'nm)'])
% hold off
% h15=gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', P1);
frameTick = startFrame:(startFrame + totFrame - 1) - histBinSize + 1;
plot(frameTick, locaNumZoomBinShift)
xlabel('Frame')
ylabel('Localization num [/unit area]')
title(['Accumulated localization in every ', num2str(histBinSize), ' frames per unit area (unit area = ', num2str(binSize), '\times', num2str(binSize), 'nm, shift window, zoom)'])
hold off
h15_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', P1);
frameTick = startFrame:(startFrame + totFrame - 1) - histBinSize + 1;
plot(frameTick, locaNumROIBinShift)
xlabel('Frame')
ylabel('Localization num [/unit area]')
title(['Accumulated localization in every ', num2str(histBinSize), ' frames per unit area (unit area = ', num2str(binSize), '\times', num2str(binSize), 'nm, shift window, ROI)'])
hold off
h15_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', P1);
frameTick = startFrame:(startFrame + totFrame - 1) - histBinSize + 1;
plot(frameTick, burstNumZoomBinShift)
xlabel('Frame')
ylabel('Burst num [/unit area]')
title(['Accumulated bursts in every ', num2str(histBinSize), ' frames per unit area (unit area = ', num2str(binSize), '\times', num2str(binSize), 'nm, shift window, zoom)'])
hold off
h15_3 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', P1);
frameTick = startFrame:(startFrame + totFrame - 1) - histBinSize + 1;
plot(frameTick, burstNumROIBinShift)
xlabel('Frame')
ylabel('Burst num [/unit area]')
title(['Accumulated bursts in every ', num2str(histBinSize), ' frames per unit area (unit area = ', num2str(binSize), '\times', num2str(binSize), 'nm, shift window, ROI)'])
hold off
h15_4 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of 2d histograms and LD map
[N1, ~, ~] = histcounts2(xPosLeftAll, yPosLeftAll, Xedges, Yedges);
N1 = reshape(N1, [size(N1, 1) * size(N1, 2), 1]);
N1(N1 == 0) = [];
sNum1 = 0;
indN1 = 0;
while sNum1 < length(N1) * satPop
    indN1 = indN1 + 1;
    sNum1 = sNum1 + sum(N1 == indN1);
end
[N2, ~, ~] = histcounts2(xPosRightAll, yPosRightAll, Xedges, Yedges);
N2 = reshape(N2, [size(N2, 1) * size(N2, 2), 1]);
N2(N2 == 0) = [];
sNum2 = 0;
indN2 = 0;
while sNum2 < length(N2) * satPop
    indN2 = indN2 + 1;
    sNum2 = sNum2 + sum(N2 == indN2);
end
tickMin = 0;
tickMax = max(indN1, indN2);
subplotTight = @(m, n, p)subtightplot(m, n, p, [0, 0.001], [0, 0], [0, 0]);
figure('Position', P1);
ax1 = subplotTight(2, 3, 1);
hold on
hist1 = histogram2(xPosLeftAll, yPosLeftAll, Xedges, Yedges, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor', 'none');
colormap(colorMap)
caxis([tickMin, tickMax]);
axis image ij
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
%title({['Localized emission in y channel']; ['scale bar = ' num2str(100) 'nm']},'FontSize',24)
hold off
ax2 = subplotTight(2, 3, 2);
hold on
hist2 = histogram2(xPosRightAll, yPosRightAll, Xedges, Yedges, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor', 'none');
colormap(colorMap)
caxis([tickMin, tickMax]);
axis image ij
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
%title({['Localized emission in x channel']; ['scale bar = ' num2str(100) 'nm']},'FontSize',24)
hold off
ax3 = subplotTight(2, 3, 3);
hold on
hist3 = histogram2(xPosPairUnpair, yPosPairUnpair, Xedges, Yedges, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on', 'EdgeColor', 'none');
N3 = hist3.Values;
N3 = reshape(N3, [size(N3, 1) * size(N3, 2), 1]);
N3(N3 == 0) = [];
sNum3 = 0;
indN3 = 0;
while sNum3 < length(N3) * satPop
    indN3 = indN3 + 1;
    sNum3 = sNum3 + sum(N3 == indN3);
end
tickMinBoth = 0;
tickMaxBoth = indN3;
colormap(colorMap)
caxis([tickMinBoth, tickMaxBoth]);
axis image ij
ax = gca;
rec = rectangle('Position', [ax.XLim(2) - scaleBar - (ax.XLim(2) - ax.XLim(1)) / 20, ...
    ax.YLim(2) - (ax.YLim(2) - ax.YLim(1)) / 100 - (ax.YLim(2) - ax.YLim(1)) / 20, ...
    scaleBar, (ax.YLim(2) - ax.YLim(1)) / 100], 'FaceColor', 'w', 'EdgeColor', 'w');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'Box', 'off');
%title({['Localized emission in both channels']; ['scale bar = ' num2str(100) 'nm']},'FontSize',24)
hold off
ax4 = subplotTight(2, 3, 4);
set(gca, 'visible', 'off');
ax5 = subplotTight(2, 3, 5);
set(gca, 'visible', 'off');
colormap(ax5, colorMap)
caxis([tickMin, tickMax]);
c12 = colorbar('north', 'Ticks', [tickMin, tickMax], ...
    'TickLabels', {num2str(tickMin), num2str(tickMax)});
c12.Position(1) = c12.Position(1) + c12.Position(3) / 4;
c12.Position(3) = c12.Position(3) * 1 / 2;
c12.Box = 'off';
% c12.Color = 'w';
ax6 = subplotTight(2, 3, 6);
set(gca, 'visible', 'off');
colormap(ax6, colorMap)
caxis([tickMinBoth, tickMaxBoth]);
c3 = colorbar('north', 'Ticks', [tickMinBoth, tickMaxBoth], ...
    'TickLabels', {num2str(tickMinBoth), num2str(tickMaxBoth)});
c3.Position(1) = c3.Position(1) + c3.Position(3) / 4;
c3.Position(3) = c3.Position(3) * 1 / 2;
c3.Box = 'off';
% c3.Color = 'w';
h16 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of photon num
figure('Position', P1);
hold on
hPhotonNumAll = histogram(photonNumPairUnpair);
hPhotonNum = histogram(photonNumPairSMAll);
hPhotonNumAll.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
hPhotonNum.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All localizations: ', num2str(length(photonNumPairUnpair)), ...
    ' , median = ', num2str(median(photonNumPairUnpair))], ...
    ['Paired localizations: ', num2str(length(photonNumPairSMAll)), ' , median = ', num2str(median(photonNumPairSMAll))])
xlabel('Photon number/localization')
% ylabel('Frequency')
title('Photon number per localization')
hold off
h17 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of correlated photon num
figure('Position', P1);
hold on
hCorrPhotonNum = histogram(photonNumPairUnpairCorrAll);
hCorrPhotonNum.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrAll)), ...
    ' , median = ', num2str(median(photonNumPairUnpairCorrAll))])
xlabel('Photon number/burst')
% ylabel('Frequency')
title('Total photon number per burst')
hold off
h17_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% export_fig(h15,[dirName '\' saveName ' h15, line, loca per bin' saveFormat],'-transparent')
export_fig(h15_1, [dirName, '\', saveName, ' h15_1, line, loca per bin, slide window, zoom', saveFormat], '-transparent')
export_fig(h15_2, [dirName, '\', saveName, ' h15_2, line, loca per bin, slide window, ROI', saveFormat], '-transparent')
export_fig(h15_3, [dirName, '\', saveName, ' h15_3, line, burst per bin, slide window, zoom', saveFormat], '-transparent')
export_fig(h15_4, [dirName, '\', saveName, ' h15_4, line, burst per bin, slide window, ROI', saveFormat], '-transparent')
export_fig(h16, [dirName, '\', saveName, ' h16, 2Dhist, reconstructed structure after all process', saveFormat], '-transparent')
export_fig(h17, [dirName, '\', saveName, ' h17, hist, photon number used for reconstruction', saveFormat], '-transparent')
export_fig(h17_1, [dirName, '\', saveName, ' h17_1, hist, photon number used for reconstruction per burst', saveFormat], '-transparent')
close(h15_1, h15_2, h15_3, h17, h17_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of photon num (zoom)
figure('Position', P1);
hold on
hPhotonNumAllZoom = histogram(PhotonArray);
hPhotonNumZoom = histogram(photonNumPairZoom);
hPhotonNumAllZoom.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
hPhotonNumZoom.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All localizations: ', num2str(length(PhotonArray)), ' , median = ', num2str(median(PhotonArray))], ...
    ['Paired localizations: ', num2str(length(photonNumPairZoom)), ' , median = ', num2str(median(photonNumPairZoom))])
xlabel('Photon number/frame')
ylabel('Frequency')
title('Photon number per frame (zoom)')
hold off
h18 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of correlated photon num (zoom)
figure('Position', P1);
hold on
hCorrPhotonNumZoom = histogram(photonNumPairUnpairCorrAllZoom);
hCorrPhotonNumZoom.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrAllZoom)), ' , median = ', num2str(median(photonNumPairUnpairCorrAllZoom))])
xlabel('Photon number/burst')
ylabel('Frequency')
title('Total photon number per burst (zoom)')
hold off
h18_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of photon num (ROI)
figure('Position', P1);
hold on
hPhotonNumAllROI = histogram(photonNumPairUnpairROI);
hPhotonNumAllROI.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All localizations: ', num2str(length(photonNumPairUnpairROI)), ' , median = ', num2str(median(photonNumPairUnpairROI))])
xlabel('Photon number/frame')
ylabel('Frequency')
title(['Photon number per frame (ROI), mean = ', num2str(mean(photonNumPairUnpairROI)), ', std = ', num2str(std(photonNumPairUnpairROI))])
hold off
h18_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of correlated photon num (ROI)
figure('Position', P1);
hold on
hCorrPhotonNumROI = histogram(photonNumPairUnpairCorrROI);
hCorrPhotonNumROI.BinEdges = 0:maxPhotonTick / 100:maxPhotonTick;
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrROI)), ' , median = ', num2str(median(photonNumPairUnpairCorrROI))])
xlabel('Photon number/burst')
ylabel('Frequency')
title(['Total photon number per burst (ROI), mean = ', num2str(mean(photonNumPairUnpairCorrROI)), ', std = ', num2str(std(photonNumPairUnpairCorrROI))])
hold off
h18_3 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of sigma
figure('Position', P1);
hold on;
hSigmaRight = histogram(sigmaRightAll);
hSigmaLeft = histogram(sigmaLeftAll);
hSigmaLeft.BinWidth = hSigmaRight.BinWidth;
% hSigmaRight.BinEdges = 0:(sigmaHighThreshold+sigmaHighThreshold/5)/100:sigmaHighThreshold+sigmaHighThreshold/5;
% hSigmaLeft.BinEdges = 0:(sigmaHighThreshold+sigmaHighThreshold/5)/100:sigmaHighThreshold+sigmaHighThreshold/5;
legend('Right channel (x-pol)', 'Left channel (y-pol)')
xlabel('\sigma [nm]')
% ylabel('Frequency')
title('Estimated \sigma in both channels')
hold off
h19 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of precision
figure('Position', P1);
hold on;
hprecRight = histogram(precRightAll);
hprecLeft = histogram(precLeftAll);
hprecRight.BinEdges = 0:binSize / 2:binSize * 20;
hprecLeft.BinEdges = 0:binSize / 2:binSize * 20;
legend(['x-pol prec, med:', num2str(median(precRightAll)), 'nm'], ['y-pol prec, med:', num2str(median(precLeftAll)), 'nm'])
xlabel('Precision [nm]')
% ylabel('Frequency')
title('Calcualted precision distribution in both channels')
hold off
h20 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of precision (zoom)
figure('Position', P1);
hold on
hPrecAllZoom = histogram(PrecArray);
hPrecAllZoom.BinEdges = 0:binSize / 2:binSize * 20;
legend(['All localizations: ', num2str(length(PrecArray)), ' , median = ', num2str(median(PrecArray))])
xlabel('Precision [nm]')
% ylabel('Frequency')
title('Calcualted precision distribution in both channels after all filtering (zoom)')
hold off
h20_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of precision (ROI)
figure('Position', P1);
hold on
hPrecAllROI = histogram(precPairUnpairROI);
hPrecAllROI.BinEdges = 0:binSize / 2:binSize * 20;
legend(['All localizations: ', num2str(length(precPairUnpairROI)), ' , median = ', num2str(nanmedian(precPairUnpairROI))])
xlabel('Precision [nm]')
% ylabel('Frequency')
title({['Calcualted precision distribution in both channels after all filtering (ROI)']; ...
    ['mean = ', num2str(mean(precPairUnpairROI)), ', std = ', num2str(std(precPairUnpairROI))]})
hold off
h20_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of on-time
figure('Position', P1);
hold on
hOnFrame = histogram(onTimeCorrAll-expTime/2);
hOnFrame.BinEdges = 0:expTime:max(onTimeCorrAll) + expTime;
xlim([0, max(onTimeCorrAll) + expTime])
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrAll)), ...
    ' , median = ', num2str(median(onTimeCorrAll))])
xlabel('On-time [ms]')
% ylabel('Frequency')
title('Total on frames per burst')
hold off
h21 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of on-time (zoom)
figure('Position', P1);
hold on
hOnFrame = histogram(onTimeCorrAllZoom-expTime/2);
hOnFrame.BinEdges = 0:expTime:max(onTimeCorrAllZoom) + expTime;
% Exponential distribution for fitting
expEqn = 'a*exp(-b*x)';
x = hOnFrame.BinEdges;
x = x(2:end-1);
startPoints = [2 * hOnFrame.Values(1), 1 / median(onTimeCorrAllZoom)];
f = fit(x.', (hOnFrame.Values(1:end - 1)).', expEqn, 'Start', startPoints);
plot(f, x, (hOnFrame.Values(1:end - 1)));
xlim([0, max(onTimeCorrAllZoom) + expTime])
onTimeCorrAllZoomTau = 1 / f.b;
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrAllZoom)), ...
    ' , fit tau = ', num2str(round(onTimeCorrAllZoomTau, 3, 'significant'))])
xlabel('On-time [ms]')
% ylabel('Frequency')
title('Total on frames per burst (zoom)')
hold off
h21_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of on-time (ROI)
figure('Position', P1);
hold on
hOnFrame = histogram(onTimeCorrROI-expTime/2);
hOnFrame.BinEdges = 0:expTime:max(onTimeCorrROI) + expTime;
% Exponential distribution for fitting
expEqn = 'a*exp(-b*x)';
x = hOnFrame.BinEdges;
x = x(2:end-1);
startPoints = [2 * hOnFrame.Values(1), 1 / median(onTimeCorrROI)];
f = fit(x.', (hOnFrame.Values(1:end - 1)).', expEqn, 'Start', startPoints);
plot(f, x, (hOnFrame.Values(1:end - 1)));
xlim([0, max(onTimeCorrROI) + expTime])
onTimeCorrROITau = 1 / f.b;
legend(['All bursts: ', num2str(length(photonNumPairUnpairCorrROI)), ...
    ' , fit tau = ', num2str(round(onTimeCorrROITau, 3, 'significant'))])
xlabel('On-time [ms]')
% ylabel('Frequency')
title(['Total on frames per burst (ROI), mean = ', num2str(mean(onTimeCorrROI)), ', std = ', num2str(std(onTimeCorrROI))])
hold off
h21_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of background
figure('Position', P1);
hold on
hBackgroundAll = histogram(backgroundPairUnapir);
% hBackgroundAll.BinEdges = 0:maxPhotonTick/100:maxPhotonTick;
legend(['All localizations: ', num2str(length(backgroundPairUnapir)), ' , median = ', num2str(median(backgroundPairUnapir))])
xlabel('Photon number/pixel/localization')
% ylabel('Frequency')
title('Background per localization')
hold off
h22 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of background (zoom)
figure('Position', P1);
hold on
hBackgroundZoom = histogram(BackgroundArray);
hBackgroundZoom.BinEdges = hBackgroundAll.BinEdges;
legend(['All localizations: ', num2str(length(BackgroundArray)), ' , median = ', num2str(median(BackgroundArray))])
xlabel('Photon number/pixel/localization')
% ylabel('Frequency')
title('Background per localization (zoom)')
hold off
h22_1 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of background (ROI)
figure('Position', P1);
hold on
hBackgroundROI = histogram(backgroundPairUnpairROI);
hBackgroundROI.BinEdges = hBackgroundAll.BinEdges;
legend(['All bursts: ', num2str(length(backgroundPairUnpairROI)), ' , median = ', num2str(median(backgroundPairUnpairROI))])
xlabel('Photon number/pixel/localization')
% ylabel('Frequency')
title(['Background per localization (ROI), mean = ', num2str(mean(backgroundPairUnpairROI)), ', std = ', num2str(std(backgroundPairUnpairROI))])
hold off
h22_2 = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

export_fig(h18, [dirName, '\', saveName, ' h18, hist, photon number used for reconstruction, zoom', saveFormat], '-transparent')
export_fig(h18_1, [dirName, '\', saveName, ' h18_1, hist, photon number used for reconstruction per burst, zoom', saveFormat], '-transparent')
export_fig(h18_2, [dirName, '\', saveName, ' h18_2, hist, photon number used for reconstruction, ROI', saveFormat], '-transparent')
export_fig(h18_3, [dirName, '\', saveName, ' h18_3, hist, photon number used for reconstruction per burst, ROI', saveFormat], '-transparent')
export_fig(h19, [dirName, '\', saveName, ' h19, hist, sigma estimation for both channels', saveFormat], '-transparent')
export_fig(h20, [dirName, '\', saveName, ' h20, hist, calculated precision', saveFormat], '-transparent')
export_fig(h20_1, [dirName, '\', saveName, ' h20_1, hist, calculated precision (after filter, zoom)', saveFormat], '-transparent')
export_fig(h20_2, [dirName, '\', saveName, ' h20_2, hist, calculated precision (after filter, ROI)', saveFormat], '-transparent')
export_fig(h21, [dirName, '\', saveName, ' h21, hist, on time per burst', saveFormat], '-transparent')
export_fig(h21_1, [dirName, '\', saveName, ' h21_1, hist, on time per burst, zoom', saveFormat], '-transparent')
export_fig(h21_2, [dirName, '\', saveName, ' h21_2, hist, on time per burst, ROI', saveFormat], '-transparent')
export_fig(h22, [dirName, '\', saveName, ' h22, hist, background', saveFormat], '-transparent')
export_fig(h22_1, [dirName, '\', saveName, ' h22_1, hist, background, zoom', saveFormat], '-transparent')
export_fig(h22_2, [dirName, '\', saveName, ' h22_2, hist, background, ROI', saveFormat], '-transparent')
close(h18, h18_1, h20, h20_1, h21, h21_1, h22, h22_1)

reconstructed = uint16(reconstructed);
imwrite(reconstructed, [dirName, '\', saveName, ' reconstructed zoom', fFormat], 'Compression', 'none')

% save localization number per frame as a csv file
dataLocNum = [(startFrame:(startFrame + totFrame - 1)).', locaNumZoom, locaNumROI];
cHeader = {'Frame', 'zoom', 'ROI'};
commaHeader = [cHeader; repmat({','}, 1, numel(cHeader))];
commaHeader = commaHeader(:).';
textHeader = cell2mat(commaHeader);
fcsv = fopen([dirName, '\', saveName, ' localization number per frame.csv'], 'w');
fprintf(fcsv, '%s\n', textHeader);
fclose(fcsv);
dlmwrite([dirName, '\', saveName, ' localization number per frame.csv'], dataLocNum, '-append');

% save([dirName '\' saveName ' analysis output.mat'],'zoomRegion','xPosPairUnpair','yPosPairUnpair','PhotonArray','BackgroundArray','burstNumZoom','burstNumZoomBinPerArea')
save([dirName, '\', saveName, ' analysis output.mat'], ...
    'zoomRegion', 'workingRegion', 'xROI', 'yROI', 'pixelROI', 'fibROI', ...
    'xPosPairUnpair', 'yPosPairUnpair', 'photonNumPairUnpair', 'backgroundPairUnapir', 'precPairUnpair', 'idPairUnpair', ...
    'PhotonArray', 'BackgroundArray', 'PrecArray', ...
    'photonNumPairUnpairROI', 'backgroundPairUnpairROI', 'precPairUnpairROI', 'idPairUnpairROI', ...
    'posPairUnpairCorrAll', 'photonNumPairUnpairCorrAll', 'idPairUnpairCorrAll', 'onTimeCorrAll', ...
    'posPairUnpairCorrAllZoom', 'photonNumPairUnpairCorrAllZoom', 'idPairUnpairCorrAllZoom', 'onTimeCorrAllZoom', ...
    'photonNumPairUnpairCorrROI', 'idPairUnpairCorrROI', 'onTimeCorrROI', ...
    'locaNumZoom', 'locaNumROI', 'burstNumZoom', 'burstNumROI', ...
    'locaNumZoomBinShift', 'locaNumROIBinShift', 'burstNumZoomBinShift', 'burstNumROIBinShift', ...
    'expTime', 'onTimeCorrAllZoomTau', 'onTimeCorrROITau', 'binSize', 'pixelSize')