function img = simDipole_novotny_v2(polar_para, position_para, img_para, pmask)
%Simulation of a dipole at an interface--By Adam Backer, Jan 1 2013
%Adapted from Axelrod, Journal of Microscopy article, 2012.
%Modified by Hesam Mazidi and Oumeng Zhang @ Lew Lab
%->---
%input
%->---
%polar_par:      structure
%polar_par.phiD-- 		azimuthal orientation, size (1,num_molecule)
%polar_par.thetaD--		polar orientation, size (1,num_molecule)
%
%position_par:   structure
%position_par.z-- 		defocus    size (1,num_molecule) (z position)
%position_par.x-- 		x position size (1,num_molecule)
%position_par.y-- 		y position size (1,num_molecule)
%
%img_par:     	 structure
%img_para.lambda--       emission wavelength
%img_para.n--       	 refractive index
%img_para.Mag--          imaging system mgnification
%img_para.NA--           numerical aperture
%
%pmask--  		 phase mask
%
%---->--
%outputs
%---->--
%img-- 			image of dipole emitter.


% get position parameters
%------------------------
fileds = {'z', 'x', 'y'};

checkFileds(position_para, fileds); %validate fields
z = position_para.z;
deltax = position_para.x;
deltay = position_para.y;

% get molecular orientation parameters
%-------------------------------------
fileds = {'phiD', 'thetaD'};
checkFileds(polar_para, fileds); %validate fields
phiD = polar_para.phiD;
thetaD = polar_para.thetaD;

% set image formation parameters
%-------------------------------
zh = 0; % ?
z2 = 0; % ?
n1 = img_para.n; %refractive index
n2 = n1; %
nh = n1; % ?
lambda = img_para.lambda; %wavelength
NA = img_para.NA; %numerical aperture
M = img_para.Mag; %magnification
K = numel(z); %number of image stacks
pixel_size = img_para.pixel_size; %camera pixel size (nm)
upsampling = 1; %upsampling factor of image space

%calculate both pupil and image plane sampling,
%one will affect the other, so make sure not to introduce aliasing

N = size(pmask, 1);
% dx_true = (pixel_size*10^9/(upsampling*M))*1e-9; %image plane sampling
dx = n1 * (pixel_size * 10^9 / (upsampling * M)) * 1e-9; %due to Abbe sine...
% condition, scale by imaging medium r.i.
dv = 1 / (N * dx); %pupil sampling, related to image plane by FFT


% define pupil coordinates
%-------------------------
[eta, xi] = meshgrid(((-1 / (2 * dx)) + (1 / (2 * N * dx))):dv:(-(1 / (2 * N * dx)) ...
    +(1 / (2 * dx))), ((-1 / (2 * dx)) + (1 / (2 * N * dx))):dv:(-(1 / (N * 2 * dx)) + (1 / (2 * dx))));
x = lambda * (eta);
y = lambda * (xi);
[phi, rho] = cart2pol(x, y);
rho_max = NA / n1; %pupil region of support determined by NA and imaging medium r.i.
k1 = n1 * (2 * pi / lambda);
kh = nh * (2 * pi / lambda);
k2 = n2 * (2 * pi / lambda);
theta1 = asin(rho); %theta in matched medium
thetah = asin((n1/nh)*sin(theta1)); %theta in thin film
theta2 = asin((n1/n2)*sin(theta1)); %theta in mismatched medium

% Fresnel coefficients
%---------------------
tp_2h = 2 * n2 * cos(theta2) ./ (n2 * cos(thetah) + nh * cos(theta2));
ts_2h = 2 * n2 * cos(theta2) ./ (nh * cos(thetah) + n2 * cos(theta2));
tp_h1 = 2 * nh * cos(thetah) ./ (nh * cos(theta1) + n1 * cos(thetah));
ts_h1 = 2 * nh * cos(thetah) ./ (n1 * cos(theta1) + nh * cos(thetah));

rp_2h = (n2 * cos(theta2) - nh * cos(thetah)) ./ (n2 * cos(theta2) + nh * cos(thetah));
rs_2h = (nh * cos(theta2) - n2 * cos(thetah)) ./ (nh * cos(theta2) + n2 * cos(thetah));
rp_h1 = (nh * cos(thetah) - n1 * cos(theta1)) ./ (nh * cos(thetah) + n1 * cos(theta1));
rs_h1 = (n1 * cos(thetah) - nh * cos(theta1)) ./ (n1 * cos(thetah) + nh * cos(theta1));

% Axelrod's equations for E-fields at pupil plane
%------------------------------------------------
mux = reshape(sin(thetaD).*cos(phiD), 1, 1, K);
muy = reshape(sin(thetaD).*sin(phiD), 1, 1, K);
muz = reshape(cos(thetaD), 1, 1, K);
tp = tp_2h .* tp_h1 .* exp(1i*kh*cos(thetah)*zh) ./ (1 + rp_2h .* rp_h1 .* exp(2i * kh * zh * cos(thetah)));
ts = ts_2h .* ts_h1 .* exp(1i*kh*cos(thetah)*zh) ./ (1 + rs_2h .* rs_h1 .* exp(2i * kh * zh * cos(thetah)));

Es = bsxfun(@times, ts.*(cos(theta1) ./ cos(theta2)).*(n1 / n2), ...
    (bsxfun(@times, muy, cos(phi)) - bsxfun(@times, mux, sin(phi))));

Ep = bsxfun(@times, tp, bsxfun(@times, (n1 / n2) .* cos(theta1), (bsxfun(@times, mux, cos(phi)) + ...
    bsxfun(@times, muy, sin(phi))))- ...
    bsxfun(@times, bsxfun(@times, muz, sin(theta1)), (n1 / n2)^2 .* (cos(theta1) ./ cos(theta2))));

PupilFilt = (rho < rho_max);

Ey_common = bsxfun(@times, PupilFilt.*(1 ./ sqrt(cos(theta1))).*exp(1i * kh * zh * cos(thetah)).* ...
    exp(1i * k2 * z2 * cos(theta2)), ...
    (bsxfun(@times, cos(phi), Es) + bsxfun(@times, sin(phi), Ep)));

Ey_1 = exp(1i*k1*(bsxfun(@times, reshape(z, 1, 1, K), cos(theta1)))); %defocus term at pupil plane
Ey_2 = exp(1i*k1*(bsxfun(@times, reshape(deltax, 1, 1, K), sin(theta1) .* cos(phi)))); %phase shift
Ey_3 = exp(1i*k1*(bsxfun(@times, reshape(deltay, 1, 1, K), sin(theta1) .* sin(phi)))); %phase shift
Ey_t = Ey_1 .* Ey_2 .* Ey_3;
Ey = bsxfun(@times, Ey_t, Ey_common);

% for propagation from pupil plane E-field to image plane via tube-lens, paraxial
% approximation is in force.
%--------------------------------------------------------------------------------
imgEy = (fftshift(fft2(Ey .* repmat((pmask), 1, 1, K))));

% image on the camera is the amplitude squared of the electric field
%-------------------------------------------------------------------
img = abs(imgEy).^2;

end