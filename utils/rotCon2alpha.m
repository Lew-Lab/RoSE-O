% 190918 Tianben Ding
% Convert rotational constraint (gamma, 0-1) into wobbling angle (alpha,
% 0-pi/2)

function alpha = rotCon2alpha(rotCon)
alpha = acos((-1 + sqrt(1 + 8 * rotCon))/2);
end