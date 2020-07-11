% 190719 Tianben Ding
% Maps the second moments onto parameters of a symmetric cone model (theta and phi (mu_x, mu_y), and gamma (Omega))
% Weight LS estimation using FIM
function [mux, muy, muz, rotMobil, x0] = secondM2SymmConeWeighted(bx, by, Bx, By, sumNorm, secM, signal, backg)

%% define a least square problem structure

Ix = secM(1) .* (bx.XX) + secM(2) .* (bx.YY) + secM(3) .* (bx.ZZ) + ...
    secM(4) .* (bx.XY) + secM(5) .* (bx.XZ) + secM(6) .* (bx.YZ);

Iy = secM(1) .* (by.XX) + secM(2) .* (by.YY) + secM(3) .* (by.ZZ) + ...
    secM(4) .* (by.XY) + secM(5) .* (by.XZ) + secM(6) .* (by.YZ);

Ix = Ix / sumNorm;
Iy = Iy / sumNorm;

Ix = signal .* Ix;
Iy = signal .* Iy;

FIM = calFIMSecondM(Bx, By, Ix, Iy, signal, backg);

%objective function
%------------------------------------------------------------
objFcn = @(z)(symmCone2secM(z).' - secM) * FIM * (symmCone2secM(z) - secM.');

%  constraints in full 3D coordinate

lb = [-ones(3, 1); zeros(1, 1)];
ub = ones(4, 1);

options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'interior-point');

% Initial value based the largest eigen value
%------------------------------------------------------------
% construct the M matrix
M = [secM(1), secM(4), secM(5); ... .
    secM(4), secM(2), secM(6); ...
    secM(5), secM(6), secM(3)];
[V, D] = eig(M);

% initialization via SVD
x0SVD(1) = real(V(1, 3));
x0SVD(2) = real(V(2, 3));
x0SVD(3) = real(V(3, 3));
x0SVD(4) = 1.5 * real(D(3, 3)) - .5;

[solSVD, ~] = fmincon(objFcn, x0SVD, [], [], [], [], lb, ub, @mycon, options);

mux = solSVD(1);
muy = solSVD(2);
muz = solSVD(3);
rotMobil = solSVD(4);
x0 = x0SVD;

    function out = symmCone2secM(z)
        z = reshape(z, 1, 4);
        muz_t = z(:, 3);
        muxx = z(:, 4) .* z(:, 1).^2 + (1 - z(:, 4)) / 3;
        muyy = z(:, 4) .* z(:, 2).^2 + (1 - z(:, 4)) / 3;
        muzz = z(:, 4) .* muz_t.^2 + (1 - z(:, 4)) / 3;
        muxy = z(:, 4) .* z(:, 1) .* z(:, 2);
        muxz = z(:, 4) .* z(:, 1) .* muz_t;
        muyz = z(:, 4) .* z(:, 2) .* muz_t;
        out = [muxx, muyy, muzz, muxy, muxz, muyz]';
    end

    function [c, ceq] = mycon(z)
        c = [];
        ceq = sum(z(:, 1:3).^2, 2) - 1;
    end
end