% 190918 Tianben Ding
% Convert dipole wobbling area (Omega, 0-4pi) into rotational constraint
% (gamma, 0-1)

function rotCon = omega2rotCon(omega)
alpha = omega2alpha(omega);
rotCon = alpha2rotCon(alpha);
end