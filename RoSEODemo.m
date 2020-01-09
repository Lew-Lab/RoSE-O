% 200109 Tianben Ding
% Load synthesized orientation data with standard polarized PSF and
% estimate phi, theta, Omega using RoSE-O for a quick demo.

%% Analysis configuration
clear; close all;clc;

% Add useful tools into path of MATLAB
folderPath = pwd;
subFolders = {'utils';'phasemask'};
for f = 1:length(subFolders)
    addpath(fullfile(folderPath,subFolders{f}))
    subSubFolders = subdir(fullfile(folderPath,subFolders{f}));
    if ~isempty(subSubFolders)
        for ff = 1:length(subSubFolders)
            addpath(fullfile(folderPath,subFolders{f},subSubFolders{ff}))
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parameters
synthDataPath = 'demodata';
synthDataName = 'standard,Omega0sr,theta90deg,phi0-180deg,sig1000,back2,iter100 synthesized data';

load([synthDataPath filesep synthDataName])
% pixelSize: pixel size in object space
% refractiveIndxSam: refractive index of sample medium
% sizeROI: image pixel size in a channel
% emitWaveL = emission wave length
% alpha_mean: 0th order photon leakage, set as 0 if you use standard PSF
% frameN: frame number of synthesized data at each orientation
% phasemaskpara: phase mask information for RoSEO
% theta_A: true theta angle in radian, first entry of SMLM_imgS, radian
% phi_A: true phi angle in radian, sencond entry of SMLM_imgS, radian
% rotMobil_A: true gamma (rotational constraint) calculated from Omega, third entry of SMLM_imgS 0-1
% sig_t: true photons detected
% backg_t: true background photons per pixel

%--
anaLoop = frameN; % partitioned frame numbers for preventing RAM overflow, should be smaller than frameN
%--

locN = round(frameN/10^round(log10(frameN)))*10^round(log10(frameN)); % only analyze and compare the first "locN" localizations from RoSEO

% penalty of localization
reg = 0.25;

photonLowThreshold = 200; % threshold for removing localization results with dimmer intensities than this value (photon number)

%% Start analysis
%% Preparation of bases for later estimation
phasemaskpara.zeroorder=alpha_mean; % add zero-order leakage into phasemaskpara
n1=Nanoscope('imageSize',sizeROI,'ADcount',1,'emissWavelength',emitWaveL,...
    'refractiveIndxSam',refractiveIndxSam,'phasemaskpara',phasemaskpara);
n1.visualizeBases;

% create PSF matrix accounting for channel transmission ratio
[FPSFx,FPSFy]=n1.createPSFstruct(n1,'ytoxchanneltransratio',1);

Emitter_t.position_para.x = 0.0; % dummy
Emitter_t.position_para.y = 0.0;
Emitter_t.position_para.z = 0.0;

[bx.XX,bx.YY,bx.ZZ, bx.XY,bx.XZ,bx.YZ,...
    by.XX,by.YY,by.ZZ,by.XY,by.XZ,by.YZ]=...
    Nanoscope.computeBases(n1,Emitter_t);

%handle for croping region of interest
up_sample=n1.pixelUpsample;
img_size=n1.imageSize;
N_pupil=size(n1.phaseMask,1);

roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,....
    -up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,:);

bx.XX = roi(bx.XX);
bx.YY = roi(bx.YY);
bx.ZZ = roi(bx.ZZ);
bx.XY = roi(bx.XY);
bx.XZ = roi(bx.XZ);
bx.YZ = roi(bx.YZ);

by.XX = roi(by.XX);
by.YY = roi(by.YY);
by.ZZ = roi(by.ZZ);
by.XY = roi(by.XY);
by.XZ = roi(by.XZ);
by.YZ = roi(by.YZ);

% normalization factor
Emitter_t.polar_para.phiD=pi/2;
Emitter_t.polar_para.thetaD=pi/2;
Emitter_t.position_para.x=0;
Emitter_t.position_para.y=0;
Emitter_t.position_para.z=0;

[brightness_scalingX,brightness_scalingY]=n1.simDipole_novotny(n1,Emitter_t);
brightness_scaling = brightness_scalingX + brightness_scalingY;
sumNorm = sum(sum(roi(brightness_scaling)));

% calculate Hadamard products of basis for FIM calculation in LS projection
XXx = n1.XXxBasis;
YYx = n1.YYxBasis;
ZZx = n1.ZZxBasis;
XYx = n1.XYxBasis;
XZx = n1.XZxBasis;
YZx = n1.YZxBasis;

XXy = n1.XXyBasis;
YYy = n1.YYyBasis;
ZZy = n1.ZZyBasis;
XYy = n1.XYyBasis;
XZy = n1.XZyBasis;
YZy = n1.YZyBasis;

% Hadamard products, for x channel
Bx.aa = (XXx).*(XXx);
Bx.ab = (XXx).*(YYx);
Bx.ac = (XXx).*(ZZx);
Bx.ad = (XXx).*(XYx);
Bx.ae = (XXx).*(XZx);
Bx.af = (XXx).*(YZx);

Bx.bb = (YYx).*(YYx);
Bx.bc = (YYx).*(ZZx);
Bx.bd = (YYx).*(XYx);
Bx.be = (YYx).*(XZx);
Bx.bf = (YYx).*(YZx);

Bx.cc = (ZZx).*(ZZx);
Bx.cd = (ZZx).*(XYx);
Bx.ce = (ZZx).*(XZx);
Bx.cf = (ZZx).*(YZx);

Bx.dd = (XYx).*(XYx);
Bx.de = (XYx).*(XZx);
Bx.df = (XYx).*(YZx);

Bx.ee = (XZx).*(XZx);
Bx.ef = (XZx).*(YZx);

Bx.ff = (YZx).*(YZx);

% for y channel
By.aa = (XXy).*(XXy);
By.ab = (XXy).*(YYy);
By.ac = (XXy).*(ZZy);
By.ad = (XXy).*(XYy);
By.ae = (XXy).*(XZy);
By.af = (XXy).*(YZy);

By.bb = (YYy).*(YYy);
By.bc = (YYy).*(ZZy);
By.bd = (YYy).*(XYy);
By.be = (YYy).*(XZy);
By.bf = (YYy).*(YZy);

By.cc = (ZZy).*(ZZy);
By.cd = (ZZy).*(XYy);
By.ce = (ZZy).*(XZy);
By.cf = (ZZy).*(YZy);

By.dd = (XYy).*(XYy);
By.de = (XYy).*(XZy);
By.df = (XYy).*(YZy);

By.ee = (XZy).*(XZy);
By.ef = (XZy).*(YZy);

By.ff = (YZy).*(YZy);

%% Estimation of PSFs
% initialize record variables

theta_std = nan(length(theta_A),length(phi_A),length(rotMobil_A)); % standard deviation of estimation parameters
phi_std = theta_std;
omega_std = theta_std;
x_std = theta_std;
y_std = theta_std;

theta_err = theta_std; % accumulated error of estimation parameters
phi_err = theta_std;
omega_err = theta_std;
x_err = theta_std;
y_err = theta_std;

detNum = theta_std; % detected molecule number
detFrameNum = theta_std; % frame numbers needed to reach "locN"

sigMed = theta_std; % median of detected signal

MAll = cell(length(theta_A),length(phi_A),length(rotMobil_A));
sigAll = MAll;
frmAll = MAll;
locPosAll = MAll;

muxAll = MAll;
muyAll = MAll;
muzAll = MAll;
thetaAll = MAll;
phiAll = MAll;
rotMobilAll = MAll;
omegaAll = MAll;
locPosXAll = MAll;
locPosYAll = MAll;

% generatin of background frames, assume you estimate background correctly
backg = repmat(backg_t,[sizeROI,2*sizeROI,frameN]);

% RoSEO estimation
for rInd = 1:length(rotMobil_A)
    rotMobil_t = rotMobil_A(rInd);
    omega_t = rotCon2omega(rotMobil_t);
    for phiInd = 1:length(phi_A)
        for thetaInd = 1:length(theta_A)
            phi_t = phi_A(phiInd);
            theta_t = theta_A(thetaInd);
            mux_t = sin(theta_t)*cos(phi_t);
            muy_t = sin(theta_t)*sin(phi_t);
            muz_t = cos(theta_t);
            
            disp('theta,phi,gamma=')
            disp([rad2deg(theta_t),rad2deg(phi_t),rotMobil_t])
            
            emitter.theta = theta_t;
            emitter.phi = phi_t;
            emitter.rotMobility = rotMobil_t;
            emitter.position_para.x = 0;
            emitter.position_para.y = 0;
            
            SMLM_img = SMLM_imgS{thetaInd,phiInd,rInd};
            
            %% RoSEO analysis
            % devide image stack into sub-stack for parfor analysis
            subLoopNum = frameN/anaLoop;
            
            % handle for recovery
            RoSEO_h=@(img,background)RoSEO(n1,img,background,FPSFx,FPSFy,'regVal',reg);
            
            loc_data=cell(1,frameN);
            p=gcp(); %current pool
            
            for l = 1:subLoopNum
                SMLM_imgTemp = SMLM_img(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
                backgTemp = backg(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
                
                for i=1:anaLoop
                    F(i)=parfeval(p,RoSEO_h,3,SMLM_imgTemp(:,:,i),backgTemp(:,:,i));
                end
                
                for i=1:anaLoop
                    % fetchNext blocks until next results are available.
                    [completedIndx,~,~,locD]=fetchNext(F);
                    
                    loc_data{(l-1)*anaLoop+completedIndx}=locD;
                end
                clear F
            end
            
            %% Projection of setimated second moments to physical space
            % accumulating estimated M vector
            M = [];
            sig = [];
            frm = [];
            locPos = [];
            
            i = 0;
            while size(M,1) < locN && i < frameN
                i = i + 1;
                if ~isempty(loc_data{i})
                    ind = loc_data{i}(:,4) < photonLowThreshold;
                    loc_data{i}(ind,:) = [];
                    M = [M; loc_data{i}(:,5:10)];
                    sig = [sig; loc_data{i}(:,4)];
                    frm = [frm; i*ones(size(loc_data{i},1),1)];
                    locPos = [locPos; loc_data{i}(:,2:3)];
                end
            end
            
            MAll{thetaInd,phiInd,rInd} = M;
            sigAll{thetaInd,phiInd,rInd} = sig;
            frmAll{thetaInd,phiInd,rInd} = frm;
            locPosAll{thetaInd,phiInd,rInd} = locPos;
            
            % projection to physical angles
            mux = nan(size(M,1),1);
            muy = nan(size(M,1),1);
            rotMobil = nan(size(M,1),1);
            muz = nan(size(M,1),1);
            
            for l = 1:subLoopNum
                backgTemp = backg(:,:, ((l-1)*anaLoop+1) : (l*anaLoop) );
                firstAnaLoopMInd = zeros(size(frm));
                lastAnaLoopMInd = zeros(size(frm));
                firstShift = 0;
                lastShift = 0;
                while sum(firstAnaLoopMInd) == 0
                    firstAnaLoopMInd = frm == ( (l-1)*anaLoop+1 + firstShift );
                    firstAnaLoopMInd = find(firstAnaLoopMInd);
                    firstShift = firstShift + 1;
                end
                firstAnaLoopMInd = firstAnaLoopMInd(1);
                while sum(lastAnaLoopMInd) == 0
                    lastAnaLoopMInd = frm == ( l*anaLoop - lastShift );
                    lastAnaLoopMInd = find(lastAnaLoopMInd);
                    lastShift = lastShift + 1;
                end
                lastAnaLoopMInd = lastAnaLoopMInd(end);
                
                parfor ii = firstAnaLoopMInd:lastAnaLoopMInd
                    [mux(ii),muy(ii),muz(ii),rotMobil(ii),~]= ...
                        secondM2SymmConeWeighted(bx,by,Bx,By,sumNorm,M(ii,:),sig(ii),backgTemp(:,:,frm(ii)));
                end
            end
            
            muR = sqrt(mux.^2+muy.^2+muz.^2);
            mux = mux./muR;
            muy = muy./muR;
            muz = muz./muR;
            
            % flip the estimated each molecule orientation based on its deviation from
            % the true orientation
            % this is flipping on dipole-like emitter
            locPosRefMu = [mux_t muy_t muz_t];
            locPosRefMu = repmat(locPosRefMu,[size(M,1),1]);
            dotMu2Ref = mux.*locPosRefMu(:,1) + muy.*locPosRefMu(:,2) + muz.*locPosRefMu(:,3);
            mux(dotMu2Ref < 0) = -mux(dotMu2Ref < 0);
            muy(dotMu2Ref < 0) = -muy(dotMu2Ref < 0);
            muz(dotMu2Ref < 0) = -muz(dotMu2Ref < 0);
            
            theta = acos(muz);
            
            phi = atan2(muy,mux);
            phi((emitter.phi - phi) > pi) = phi((emitter.phi - phi) > pi) + 2*pi;
            phi((emitter.phi - phi) < -pi) = phi((emitter.phi - phi) < -pi) - 2*pi;
            
            omega = rotCon2omega(rotMobil);
            
            muxAll{thetaInd,phiInd,rInd} = mux;
            muyAll{thetaInd,phiInd,rInd} = muy;
            muzAll{thetaInd,phiInd,rInd} = muz;
            thetaAll{thetaInd,phiInd,rInd} = theta;
            phiAll{thetaInd,phiInd,rInd} = phi;
            rotMobilAll{thetaInd,phiInd,rInd} = rotMobil;
            omegaAll{thetaInd,phiInd,rInd} = omega;
            locPosXAll{thetaInd,phiInd,rInd} = locPos(:,1);
            locPosYAll{thetaInd,phiInd,rInd} = locPos(:,2);
            
            theta_std(thetaInd,phiInd,rInd) = std(theta);
            phi_std(thetaInd,phiInd,rInd) = std(phi);
            omega_std(thetaInd,phiInd,rInd) = std(omega);
            x_std(thetaInd,phiInd,rInd) = std(locPos(:,1));
            y_std(thetaInd,phiInd,rInd) = std(locPos(:,2));
                       
            theta_err(thetaInd,phiInd,rInd) = sum(theta-emitter.theta)/length(sig);
            phi_err(thetaInd,phiInd,rInd) = sum(phi-emitter.phi)/length(sig);
            omega_err(thetaInd,phiInd,rInd) = sum(omega-omega_t)/length(sig);
            x_err(thetaInd,phiInd,rInd) = sum(locPos(:,1)-emitter.position_para.x)/length(sig);
            y_err(thetaInd,phiInd,rInd) = sum(locPos(:,2)-emitter.position_para.y)/length(sig);
            
            detNum(thetaInd,phiInd,rInd) = length(sig);
            detFrameNum(thetaInd,phiInd,rInd) = frm(end);
                        
            sigMed(thetaInd,phiInd,rInd) = median(sig);            
        end
    end
end
%% Figure visualization
figure;
subplot(2,3,1)
plot(rad2deg(phi_A),rad2deg(theta_std))
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\sigma_{\theta} [deg]')
title(['\theta precision'])
subplot(2,3,2)
plot(rad2deg(phi_A),rad2deg(phi_std))
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\sigma_{\phi} [deg]')
title(['\phi precision'])
subplot(2,3,3)
plot(rad2deg(phi_A),omega_std)
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\sigma_{\Omega} [sr]')
title(['\Omega precision'])
subplot(2,3,4)
plot(rad2deg(phi_A),rad2deg(theta_err))
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\theta-\theta_0 [deg]')
title(['\theta bias'])
subplot(2,3,5)
plot(rad2deg(phi_A),rad2deg(phi_err))
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\phi-\phi_0 [deg]')
title(['\phi bias'])
subplot(2,3,6)
plot(rad2deg(phi_A),omega_err)
xlim([rad2deg(phi_A(1)) rad2deg(phi_A(end))])
xlabel('\phi_0 [deg]')
ylabel('\Omega-\Omega_0 [sr]')
title(['\Omega bias'])
