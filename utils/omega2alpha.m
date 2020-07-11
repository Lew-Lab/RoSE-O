% 190918 Tianben Ding
% Convert dipole wobbling area (Omega, 0-4pi) into wobbling angle (alpha,
% 0-pi/2)

function alpha = omega2alpha(omega)
alpha = 2 * asin(sqrt(omega / 8 / pi));
end