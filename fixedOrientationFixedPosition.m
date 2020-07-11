%% INFO

% generates images of a molecule with fixed positions and a fixed
% molecular orientation

%% set up a nanoscope object with default settings

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
%

%% set up parameters

distance2grid = [30]; %nm
number_cases = 200; %
mean_brightness = 2e3;
mean_backg = 5;


%orientation parameters
theta = [0:10:90] * pi / 180;
phi = [0, pi / 4, pi / 2];
rotMobility = [.1:.1:1];

for ii = 1:numel(theta)
    for jj = 1:numel(phi)
        for kkk = 1:numel(rotMobility)

            SMLM_img = zeros(n1.imageSize, 2*n1.imageSize, number_cases);
            pos_x = distance2grid * rand(1);
            pos_y = distance2grid * rand(1);
            %set position of emitters
            Emitter = struct('position_para', {}, 'theta', {}, 'phi', {}, 'rotMobility', {});
            for l = 1:number_cases
                Emitter(l).position_para.x = pos_x;
                Emitter(l).position_para.y = pos_y;
                Emitter(l).position_para.z = 0;

                Emitter(l).theta = theta(ii);
                Emitter(l).phi = phi(jj);
                Emitter(l).rotMobility = rotMobility(kkk);
            end

            %% creat images on x-, y-channels

            % channel_mismatch set to 0

            num_char = 0;

            %create handle

            formImg_h = @(E)n1.formImage(E, 'channel_mismatch', [0, 0]);
            p = gcp('nocreate');
            if isempty(p)
                p = parpool('local', 14); %current pool
            end

            for k = 1:number_cases
                F(k) = parfeval(p, formImg_h, 2, Emitter(k));
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

            tic;
            num_char = 0;

            for k = 1:num_frames
                %creat emitter structure
                F(k) = parfeval(p, RoSEO_h, 2, SMLM_img(:, :, k));
            end

            %collect the results and monitor the progress
            loc_data = cell(1, num_frames);
            cnt = 1;

            for k = 1:num_frames

                % fetchNext blocks until next results are available.
                [completedIndx, gammaf, recovS] = fetchNext(F);
                if ~isempty(gammaf)
                    loc_data{completedIndx} = get_loc_data2(gammaf, recovS);
                end
                num_char = progress_bar(cnt/num_frames, num_char, 20);
                cnt = cnt + 1;
            end

            elapsed_time = toc;

            %% save results

            folder_path = fullfile('tests', ['fixedOrientationFixedPos', '_theta', num2str(theta(ii)), '_phi', ...
                num2str(phi(jj)), '_rotMobility', num2str(rotMobility(kkk))]);
            mkdir(folder_path)
            save(fullfile(folder_path, 'Emitter'), 'Emitter');

            save(fullfile(folder_path, 'loc_data'), 'loc_data');

            save(fullfile(folder_path, 'elapsed_time'), 'elapsed_time');

        end
    end
end

%%

% x_bias=zeros(numel(theta),numel(phi),numel(rotMobility));
% x_prec=zeros(numel(theta),numel(phi),numel(rotMobility));
%
% y_bias=zeros(numel(theta),numel(phi),numel(rotMobility));
% y_prec=zeros(numel(theta),numel(phi),numel(rotMobility));
%
% CRLB_x=zeros(numel(theta),numel(phi),numel(rotMobility));
% CRLB_y=zeros(numel(theta),numel(phi),numel(rotMobility));
%
% %%
% for ii=1:numel(theta)
% for jj=1:numel(phi)
% for kkk=1:numel(rotMobility)
%        file_path=fullfile('tests',['fixedOrientationFixedPos','_theta',num2str(theta(ii)),'_phi',...
% num2str(phi(jj)),'_rotMobility',num2str(rotMobility(kkk))],'loc_data');
%     load(file_path)
%  file_path=fullfile('tests',['fixedOrientationFixedPos','_theta',num2str(theta(ii)),'_phi',...
% num2str(phi(jj)),'_rotMobility',num2str(rotMobility(kkk))],'Emitter');
%     load(file_path)
%
%     x_est=[];
%     y_est=[];
%     x_gt=Emitter(1).position_para.x;
%     y_gt=Emitter(1).position_para.y;
%
%     for iter_loc=1:numel(loc_data)
%
%         if ~isempty(loc_data{iter_loc})
%
%             x_est=[x_est;loc_data{iter_loc}(2)];
%             y_est=[y_est;loc_data{iter_loc}(3)];
%
%         end
%     end
%
%     [CRLB_vector]=n1.CRB_orinet(Emitter(1),'paraType','angular');
%     x_bias(ii,jj,kkk)=-mean(x_est)+x_gt;
%     x_prec(ii,jj,kkk)=std(x_est);
%
%     y_bias(ii,jj,kkk)=-mean(y_est)+y_gt;
%     y_prec(ii,jj,kkk)=std(y_est);
%
%     CRLB_x(ii,jj,kkk)=CRLB_vector(1);
%     CRLB_y(ii,jj,kkk)=CRLB_vector(2);
%
%
% end
% end
% end
%
% disp('Analysis is finished!')
% %% visialize
