function [FXXx, FYYx, FZZx, FXYx, FXZx, FYZx, lateral_grid_p] = FourierTriSpotPSF(imgPara, pmask)
% FourierTriSpotPSF produces the Fourier transforms of basis matrices. It uses the
% vectorial model and the theory of this function is based on  Adam S. Backer and
% W. E. Moerner, Opt. Express 23, 4255-4276 (2015)
%->---
%input
%->---
%img_par:     	 structure
%img_para.lambda--   	 emission wavelength
%img_para.n--       	 refractive index
%img_para.Mag--          imaging system mgnification
%img_para.NA--           numerical aperture
%img_para.image_size--   side length of the image (image must be an square
%region)
%
%pmask--  	          	 phase mask
%---->-
%output
%---->-
%FXXy--                  Foureir transform of the XXy basis
%FYYy--                  Foureir transform of the YYy basis
%FZZy--                  Foureir transform of the ZZy basis
%FXYy--                  Foureir transform of the XYy basis
%FXZy--                  Foureir transform of the XZy basis
%FYZy--                  Foureir transform of the YZy basis

position_para.z = 0;
position_para.x = 0;
position_para.y = 0;


% basis images of the imaging system
%-----------------------------------

% define a handle for convenience
simDipole_novotny_v2_h = @(polar_par)simDipole_novotny_v3(polar_par, position_para, imgPara, pmask);

polar_par.phiD = zeros(1, 1);
polar_par.thetaD = pi / 2 * ones(1, 1);

[XXy, XXx] = simDipole_novotny_v2_h(polar_par);

polar_par.phiD = pi / 2 * ones(1, 1);
polar_par.thetaD = pi / 2 * ones(1, 1);

[YYy, YYx] = simDipole_novotny_v2_h(polar_par);

polar_par.phiD = pi / 2 * ones(1, 1);
polar_par.thetaD = 0 * ones(1, 1);

[ZZy, ZZx] = simDipole_novotny_v2_h(polar_par);

polar_par.phiD = pi / 4 * ones(1, 1);
polar_par.thetaD = pi / 2 * ones(1, 1);

[XYytmp, XYxtmp] = simDipole_novotny_v2_h(polar_par);

polar_par.phiD = 0 * ones(1, 1);
polar_par.thetaD = pi / 4 * ones(1, 1);

[XZytmp, XZxtmp] = simDipole_novotny_v2_h(polar_par);

polar_par.phiD = pi / 2 * ones(1, 1);
polar_par.thetaD = pi / 4 * ones(1, 1);

[YZytmp, YZxtmp] = simDipole_novotny_v2_h(polar_par);

XYx = 2 * XYxtmp - XXx - YYx;
XZx = 2 * XZxtmp - XXx - ZZx;
YZx = 2 * YZxtmp - YYx - ZZx;

XYy = 2 * XYytmp - XXy - YYy;
XZy = 2 * XZytmp - XXy - ZZy;
YZy = 2 * YZytmp - YYy - ZZy;


% corp the basis images to match the desired image size
%------------------------------------------------------

%accounting for photon loss

polar_par.phiD = pi / 2;
polar_par.thetaD = pi / 2;

brightness_scaling = simDipole_novotny_v2(polar_par, position_para, imgPara, pmask);


img_size = imgPara.image_size;
up_sample = 1;
N_pupil = size(XXy, 1);
roi = @(img)img(-up_sample*(img_size - 1)/2+N_pupil/2+1:1:up_sample*(img_size - 1)/2+N_pupil/2+1, ... .
    -up_sample*(img_size - 1)/2+N_pupil/2+1:1:up_sample*(img_size - 1)/2+N_pupil/2+1, :);

sumnorm = sum(sum(roi(brightness_scaling)));

XXy_corp = roi(XXy) / sumnorm;
YYy_corp = roi(YYy) / sumnorm;
ZZy_corp = roi(ZZy) / sumnorm;
XYy_corp = roi(XYy) / sumnorm;
XZy_corp = roi(XZy) / sumnorm;
YZy_corp = roi(YZy) / sumnorm;

XXx_corp = roi(XXx) / sumnorm;
YYx_corp = roi(YYx) / sumnorm;
ZZx_corp = roi(ZZx) / sumnorm;
XYx_corp = roi(XYx) / sumnorm;
XZx_corp = roi(XZx) / sumnorm;
YZx_corp = roi(YZx) / sumnorm;
% compute Fourier transforms
% --------------------------

FXXy = single(fft2(ifftshift(up_sample^2 * XXy_corp)));
FYYy = single(fft2(ifftshift(up_sample^2 * YYy_corp)));
FZZy = single(fft2(ifftshift(up_sample^2 * ZZy_corp)));
FXYy = single(fft2(ifftshift(up_sample^2 * XYy_corp)));
FXZy = single(fft2(ifftshift(up_sample^2 * XZy_corp)));
FYZy = single(fft2(ifftshift(up_sample^2 * YZy_corp)));


FXXx = single(fft2(ifftshift(up_sample^2 * XXx_corp)));
FYYx = single(fft2(ifftshift(up_sample^2 * YYx_corp)));
FZZx = single(fft2(ifftshift(up_sample^2 * ZZx_corp)));
FXYx = single(fft2(ifftshift(up_sample^2 * XYx_corp)));
FXZx = single(fft2(ifftshift(up_sample^2 * XZx_corp)));
FYZx = single(fft2(ifftshift(up_sample^2 * YZx_corp)));
% computing  lateral grid points
%-------------------------------
pixel_size = imgPara.pixel_size;
Mag = imgPara.Mag;
pixel_size_obj = pixel_size / Mag;
lateral_grid_p = (-(img_size - 1) / 2:1 / up_sample:(img_size - 1) / 2) * pixel_size_obj;

end
