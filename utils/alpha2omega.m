% 190918 Tianben Ding
% Convert wobbling angle (alpha, 0-pi/2) into dipole wobbling area (Omega,
% 0-4pi)

function omega = alpha2omega(alpha)
omega = 2 .* 4 .* pi .* sin(alpha./2).^2;
end