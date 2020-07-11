%% this script checks the correctness of the forward model used in RoSE-O

% set up a nanoscope object with default settings
imageSize = 81;
if ~exist('n1')
    n1 = Nanoscope('imageSize', imageSize);
end

% create PSF
[FPSFx, FPSFy] = n1.createPSFstruct(n1);

% clear the job queue
if exist('F', 'var')
    clear F
end

%% set up parameters

distance2grid = [-10]; %nm
number_cases = 10; %
mean_brightness = 5e3;
mean_backg = 5;


%orientation parameters
theta = [90] * pi / 180;
phi = [90] * pi / 180;
rotMobility = [1];
pos_x = distance2grid;
pos_y = distance2grid;

Emitter = struct('position_para', {}, 'theta', {}, 'phi', {}, 'rotMobility', {});

Emitter(1).position_para.x = pos_x;
Emitter(1).position_para.y = pos_y;
Emitter(1).position_para.z = 0;

Emitter(1).theta = theta(1);
Emitter(1).phi = phi(1);
Emitter(1).rotMobility = rotMobility(1);

num_char = 0;

%create handle

formImg_h = @(E)n1.formImage(E, 'channel_mismatch', [0, 0]);

p = gcp('nocreate');
if isempty(p)
    p = parpool('local'); %current pool
end

for k = 1:number_cases
    F(k) = parfeval(p, formImg_h, 2, Emitter);
end
%collect the results and monitor the progress
cnt = 1;
for k = 1:number_cases
    % fetchNext blocks until next results are available.
    [completedIndx, imgx, imgy] = fetchNext(F);
    img_t = [imgx, imgy];
    % apply Poisson noise

    SMLM_img(:, :, (completedIndx - 1)*1+1:(completedIndx)*1) = ...
        poissrnd(bsxfun(@times, img_t, ... .
        mean_brightness)+mean_backg);
    num_char = progress_bar(cnt/number_cases, num_char, 20);
    cnt = cnt + 1;
end

%% analyze

num_frames = size(SMLM_img, 3);
imgSize = size(SMLM_img, 1);
backg = [mean_backg, mean_backg];
%handle for recovery
RoSEO_h = @(img)RoSEO(n1, img, backg, FPSFx, FPSFy);

num_char = 0;

for k = 1:number_cases
    %creat emitter structure
    F(k) = parfeval(p, RoSEO_h, 2, SMLM_img(:, :, k));
end

%collect the results and monitor the progress
loc_data = cell(1, num_frames);
cnt = 1;

for k = 1:number_cases

    % fetchNext blocks until next results are available.
    [completedIndx, gammaf, recovS] = fetchNext(F);
    if ~isempty(gammaf)
        loc_data{completedIndx} = get_loc_data2(gammaf, recovS);
    end
    num_char = progress_bar(cnt/num_frames, num_char, 20);
    cnt = cnt + 1;
end

%% plot results

x = [];
y = [];
for i = 1:number_cases
    x = [x; loc_data{i}(1, 2)];
    y = [y; loc_data{i}(1, 3)];
end

figure;
h = histogram(x);
hold on
h2 = histogram(y);
h.BinWidth = 1;
h2.BinWidth = 1;

legend({'x', 'y'})
axis square
