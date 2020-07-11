% 190918 Tianben Ding
% Convert wobbling angle (alpha, 0-pi/2) into rotational constraint (gamma,
% 0-1)

function rotCon = alpha2rotCon(alpha)
rotCon = (cos(alpha).^2 + cos(alpha)) / 2;
end