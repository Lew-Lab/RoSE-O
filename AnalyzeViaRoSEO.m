%% Necessary information

% place following files in a directory called source_dir
% 1- transmissionRatio_name: a mat file that contains transmit_ratio_L2R as
% a variable
% transmit_ratio_L2R is a scalar indicating the light transmission ratio
% between left and right channels
%
% 2- zerothOrder_name: a mat file containing zerothOrder as a variable
% zerothOrder_RL is a 1*2 array specifying the zeroth order factor
% in the pupil plane according to [x_channel y_channel] or [right channel left
% channel] formt
%
% 3- image stack and background stack in binary format
% You may want to write your image and background stack into binary format
% with names img_stack_name.bin and backg_stack_name.bin accordingly
% use writeSMLMbackg2bin function

%% add folders to search path

%----------------------------------------------
folders = {'utils', 'phasemask'};

for ii = 1:length(folders)

    addpath(genpath(folders{ii}));

end

%% data information

source_dir = 'test_data'; % directory containing required data

% transmission ratio
transmissionRatio_name = 'PcCDt3_0815c1r06_NR_transratio_L2R'; % name of the transmission ratio file

% zeroth-order leakage of PSF
zerothOrder_name = 'PcCDt3_0815c1r06_zerothOrder_RL'; % name of the zeroth order file

img_stack_name = 'SMLM_img.bin'; %name of the image stack to be analyzed (add .bin)

backg_stack_name = 'backg.bin'; %name of the background stack (add .bin)

%% load necessary data

%----------------------------------------------
% transmission ratio

try
    load(fullfile(source_dir, transmissionRatio_name))
catch
    msg = [source_dir, ' or', transmissionRatio_name, ' dose not exist'];
    causeException = MException('MATLAB:AnylzeViaRoSEO:TransmissionRatio', msg);
    ME = addCause(ME, causeException);
    rethrow(ME)
end

% zeroth-order leakage of PSF
try
    load(fullfile(source_dir, zerothOrder_name))
catch
    msg = [source_dir, ' or', zerothOrder_name, ' dose not exist'];
    causeException = MException('MATLAB:AnalyzeViaRoSEO:zerothOrder', msg);
    ME = addCause(ME, causeException);
    rethrow(ME)
end

%% format data and computational resources

%----------------------------------------------
imgSize = 121; % side length of the image to be analyzed (in camera pixels)
num_frames = 450; % total number of frames in your stack
num_frames_backg = 450; % number of background frames; it can be equal to ONE frame or the same as num_frames
num_frames_p = 50; % number of frames to be analyzed in parallel by all workers
number_of_workers = 2; % number of workers


% initialize a parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool('local', number_of_workers)
end
num_full_stack = floor(num_frames/num_frames_p);

% map the address of the data to m_img & m_backg
try
    m_img = memmapfile(fullfile(source_dir, img_stack_name), ...
        'Format', {'single', [imgSize * imgSize * 2, num_frames_p], 'SMLM_img'}, 'Repeat', num_full_stack);
catch ME
    msg = [img_stack_name, ' dose not exist or', ' the dimension of ', img_stack_name, ' does not match ', num2str(imgSize)];
    causeException = MException('MATLAB:AnylzeViaRoSEO:formatData', msg);
    ME = addCause(ME, causeException);
    rethrow(ME)
end


if num_frames_backg == 1
    num_full_stack_backg = 1;
    num_frames_backg_p = 1;
else
    num_full_stack_backg = num_full_stack;
    num_frames_backg_p = num_frames_p;
end
try
    m_backg = memmapfile(fullfile(source_dir, backg_stack_name), ...
        'Format', {'single', [imgSize * imgSize * 2, num_frames_backg_p], 'backg'}, 'Repeat', num_full_stack_backg);
catch ME
    msg = ['backg.mat dose not exist or', ' the dimension of backg does not match ', num2str(imgSize)];
    causeException = MException('MATLAB:AnylzeViaRoSEO:formatData', msg);
    ME = addCause(ME, causeException);
    rethrow(ME)
end

%% analyze via RoSEO

%----------------------------------------------
% create nanoscope object analysis
maskName = fullfile('phasemask', 'tri-spot');
emitter_wavelength = 610; %nm
left_to_right_trans_ratio = transmit_ratio_L2R;
zerothOrder = zerothOrder_RL;
%construct phasemaskpara
phasemaskpara.zeroorder = zerothOrder;
phasemaskpara.maskname = maskName;

n1 = Nanoscope('imageSize', imgSize, 'ADcount', 1, 'emissWavelength', emitter_wavelength, ...
    'phasemaskpara', phasemaskpara);
% create PSF matrix accounting for channel transmission ratio
[FPSFx, FPSFy] = n1.createPSFstruct(n1, 'ytoxchanneltransratio', left_to_right_trans_ratio);

%%

loc_data_tot = cell(1, num_full_stack*num_frames_p);


for ll = 1:num_full_stack

    num_char = 0; %initialize the progress bar
    cnt = 1;

    if exist('F', 'var')
        cancel(F);
        clear F
    end

    % load data on demand
    SMLM_img = reshape(m_img.Data(ll).SMLM_img, imgSize, imgSize*2, num_frames_p);
    backg = reshape(m_backg.Data(ll).backg, imgSize, imgSize*2, num_frames_backg_p);

    %handle for recovery

    if num_frames_backg_p == 1

        RoSEO_h = @(img)RoSEO(n1, img, backg, FPSFx, FPSFy, 'regval', .22);
        for k = 1:num_frames_p
            %creat emitter structure

            F(k) = parfeval(p, RoSEO_h, 2, SMLM_img(:, :, k));

        end
    else
        RoSEO_h = @(img, backg)RoSEO(n1, img, backg, FPSFx, FPSFy, 'regval', .22);
        for k = 1:num_frames_p
            %creat emitter structure

            F(k) = parfeval(p, RoSEO_h, 2, SMLM_img(:, :, k), backg(:, :, k));

        end
    end


    %collect the results and monitor the progress
    loc_data = cell(1, num_frames_p);
    for k = 1:num_frames_p
        % fetchNext blocks until next results are available.
        [completedIndx, gammaf, recovS] = fetchNext(F);
        if ~isempty(gammaf)
            loc_data{completedIndx} = get_loc_data2(gammaf, recovS);
        end
        %free the memory
        clear gammaf;

        %display progress
        num_char = progress_bar(cnt/num_frames_p, num_char, 20);
        cnt = cnt + 1;
    end

    loc_data_tot((ll-1)*num_frames_p+1:ll*num_frames_p) = loc_data;
end

%%

% filter localization far from center; and too close
%----------------------------------------------
minimum_distance = 300; % nm
pixelSize = 58.5;
distance_from_center = (imgSize - 50) / 2; %pixels
br_thres = .8e3;
num_frames_p = numel(loc_data_tot);


x_est = cell(1, num_frames_p);
y_est = cell(1, num_frames_p);
secM = cell(1, num_frames_p);
br = cell(1, num_frames_p);
for i = 1:num_frames_p

    if ~isempty(loc_data_tot{i})

        x_est_t = loc_data_tot{i}(:, 2);
        y_est_t = loc_data_tot{i}(:, 3);
        secM_t = loc_data_tot{i}(:, 5:end);
        br_t = loc_data_tot{i}(:, 4);
        %filter based on distance from center
        indx_valid = (x_est_t > -distance_from_center * pixelSize) .* (x_est_t < distance_from_center * pixelSize) .* ...
            (y_est_t > -distance_from_center * pixelSize) .* (y_est_t < distance_from_center * pixelSize);

        x_est_t = x_est_t(indx_valid > 0);
        y_est_t = y_est_t(indx_valid > 0);
        secM_t = secM_t(indx_valid > 0, :);
        br_t = br_t(indx_valid > 0);

        indx = zeros(1, numel(x_est_t));

        %         filter based on intera-distance
        if ~isempty(x_est_t)
            for j = 1:numel(x_est_t)


                dist_t = sqrt(((x_est_t(j) - x_est_t).^2+(y_est_t(j) - y_est_t).^2));
                indx_t = dist_t > 0;
                if ~any(indx_t) || min(dist_t(indx_t > 0)) > minimum_distance

                    if br_t(j) > br_thres
                        indx(j) = 1;
                    end
                end
            end

            x_est{i} = x_est_t(indx > 0);
            y_est{i} = y_est_t(indx > 0);
            secM{i} = secM_t(indx > 0, :);
            br{i} = br_t(indx > 0);
        end
    end
end

%% save results

% % issue an email when done
% %----------------------------------------------
%     props = java.lang.System.getProperties;
%     props.setProperty('mail.smtp.auth','true');
%     props.setProperty('mail.smtp.socketFactory.class', ...
%         'javax.net.ssl.SSLSocketFactory');
%     props.setProperty('mail.smtp.socketFactory.port','465');
%     setpref('Internet','SMTP_Username','hesam.creative@gmail.com');
%     setpref('Internet','SMTP_Password','');
%     setpref('Internet','E_mail','hesam.creative@gmail.com');
%     setpref('Internet','SMTP_Server','smtp.gmail.com');
%     sendmail('hesam.creative@gmail.com','Simulation is done!');
%
