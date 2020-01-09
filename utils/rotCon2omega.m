% 190918 Tianben Ding
% Convert rotational constraint (gamma, 0-1) into dipole wobbling area
% (Omega, 0-4pi)

function omega = rotCon2omega(rotCon)
alpha = acos( (  -1 + sqrt(1+8*rotCon)  )/2);
omega = alpha2omega(alpha);
end