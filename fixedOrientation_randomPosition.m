%% INFO

% generates images of a molecule with random positions and a fixed
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

distance2grid = 58.5 / 2; %nm
number_cases = 200; %
mean_brightness = 2e3;
mean_backg = 5;


%orientation parameters
theta = (0:10:90) * pi / 180;
phi = [0, pi / 4, pi / 2];
rotMobility = [.2, .5, .8];

for ii = 1:5
    for jj = 1:5
        for kkk = 1:3

            SMLM_img = zeros(n1.imageSize, 2*n1.imageSize, number_cases);

            %set position of emitters
            Emitter = struct('position_para', {}, 'theta', {}, 'phi', {}, 'rotMobility', {});
            for l = 1:number_cases
                Emitter(l).position_para.x = distance2grid * rand(1);
                Emitter(l).position_para.y = distance2grid * rand(1);
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
                p = parpool('local'); %current pool
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

            folder_path = fullfile('tests', ['fixedOrientationRandomPos', '_theta', num2str(theta(ii)), '_phi', ...
                num2str(phi(jj)), '_rotMobility', num2str(rotMobility(kkk))]);
            mkdir(folder_path)
            save(fullfile(folder_path, 'Emitter'), 'Emitter');

            save(fullfile(folder_path, 'loc_data'), 'loc_data');

            save(fullfile(folder_path, 'elapsed_time'), 'elapsed_time');

        end
    end
end

%%

% phi_bias=zeros(5,2,3);
% phi_prec=zeros(5,2,3);
% theta_bias=zeros(5,2,3);
% theta_prec=zeros(5,2,3);
% rotMobil_bias=zeros(5,2,3);
% rotMobil_prec=zeros(5,2,3);
% MAE=zeros(5,2,3);
% CRLB_MAE=zeros(5,2,3);
% CRLB_rotMob=zeros(5,2,3);
%
% %%
% for ii=1:numel(theta)
% for jj=1:2
% for kkk=1
%     file_path=fullfile('tests',['fixedOrientationRandomPos','_theta',num2str(theta(ii)),'_phi',...
% num2str(phi(jj)),'_rotMobility',num2str(rotMobility(kkk))],'loc_data');
%     load(file_path)
%
%
%     mux=nan(numel(loc_data),1);
%     muz=nan(numel(loc_data),1);
%     muy=nan(numel(loc_data),1);
%     rotMobil=nan(numel(loc_data),1);
%
%
%
%     % map to orientation & position space
%     for iter_loc=1:numel(loc_data)
%
%         if ~isempty(loc_data{iter_loc})
%
%             try
%     [mux(iter_loc,1),...
%         muy(iter_loc,1),...
%         muz(iter_loc,1),...
%         rotMobil(iter_loc)]=secondM2SymmCone( loc_data{iter_loc}(5:end));
%             catch ME
%                 warning([ME.message,'Problem with secondM2SymmCone.']);
%             end
%
%         end
%     end
%
%      indxNan=isnan(mux);
%      mux(indxNan)=[];
%      muy(indxNan)=[];
%      muz(indxNan)=[];
%      rotMobil(indxNan)=[];
%
%
%     muz_gt=cos(theta(ii));
%     mux_gt=sin(theta(ii)).*cos(phi(jj));
%     muy_gt=sin(theta(ii)).*sin(phi(jj));
%
%
%     mu_gt=[mux_gt muy_gt muz_gt];
%
%     indx=muz<0;
%     muz(indx>0)=-muz(indx>0);
%     mux(indx>0)=-mux(indx>0);
%     muy(indx>0)=-muy(indx>0);
%     mu_est=[mux muy muz];
%
%     theta_est=acosd(muz);
%     phi_est_t=atan2d(muy,mux)+180;
%
%     phi_error_t=[phi_est_t-phi(jj)*180/pi,phi_est_t-(phi(jj)*180/pi+180),phi_est_t-(phi(jj)*180/pi-180)];
%     min_phi_error=min(abs(phi_error_t),[],2);
%     indx=abs(phi_error_t)==min_phi_error;
%     phi_error=phi_error_t(indx);
%
%     phi_bias(ii,jj,kkk)=mean(phi_error);
%     phi_prec(ii,jj,kkk)=std(phi_error);
%
%     theta_bias(ii,jj,kkk)=mean(theta_est)-theta(ii)*180/pi;
%     theta_prec(ii,jj,kkk)=std(theta_est);
%
%     %mean-square angular error
%     delta_mu=min([sqrt(sum((mu_est-mu_gt).^2,2)),sqrt(sum((mu_est+mu_gt).^2,2))],[],2);
%     MAE(ii,jj,kkk)=mean((2*asin(delta_mu/2)).^2);
%
%     rotMobil_prec(ii,jj,kkk)=std(rotMobil);
%     rotMobil_bias(ii,jj,kkk)=mean(rotMobil)-rotMobility(kkk);
%
%     % fundamental bound for mean-square angular error
%     Emitter=struct();
%     Emitter.position_para.x=0;
%     Emitter.position_para.y=0;
%      Emitter.theta=theta(ii);
%     Emitter.phi=phi(jj);
%     Emitter.rotMobility=rotMobility(kkk);
%
%     CRLB_vector=CRB_orinet(n1,Emitter,'brightness',mean_brightness,...
%         'background',mean_backg);
%
%
%     CRLB_MAE(ii,jj,kkk) =cos(pi/2-theta(ii))^2*CRLB_vector(4)^2+CRLB_vector(3)^2; % double-check
%
%     CRLB_rotMob(ii,jj,kkk)=CRLB_vector(end);
% end
% end
% end
% disp('Analysis is finished!')
%
% %% visualize
%
% %Mean angular error
%
%
% hf=figure;
%
% plot(theta*180/pi,sqrt(CRLB_MAE(:,1,1))*180/pi,'-k')
% hold on
% plot(theta*180/pi,sqrt(MAE(:,1,1))*180/pi,'--^k')
%
% xlabel('theta(degree)')
% ylabel('mean angular error(degree)')
%
% plot(theta*180/pi,sqrt(CRLB_MAE(:,2,1))*180/pi,'-m')
% hold on
% plot(theta*180/pi,sqrt(MAE(:,2,1))*180/pi,'--^m')
%
% xlabel('\theta(degree)')
% ylabel('mean angular error(degree)')
%
% legend({'Square root of MSAE_B (\phi=0)',' Estimated square root of MSAE_B (\phi=0)',...
%     'Square root of MSAE_B (\phi=90)',' Estimated square root of MSAE_B (\phi=90)'})
% axis square
% xlim([0 90])
%
% %%
% hf2=figure;
%
% errorbar(theta*180/pi,rotMobil_bias(:,1,3),rotMobil_prec(:,1,3),'LineWidth',1.6);
% hold on
%
% errorbar(theta*180/pi,0*rotMobil_bias(:,1,3),CRLB_rotMob(:,1,3),'LineWidth',1.6);
% xlabel('\theta(degree)')
% legend({'Estimated errors(\gamma=0.8,\phi=0)',' CRLB  (\gamma=0.8,\phi=0)'})
% axis square
% axis square
