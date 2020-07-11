% 190404 Tianben Ding
% Calculate Fisher information matric of second moments
% Using fully pixelated images

function FIM = calFIMSecondM(Bx, By, Ix, Iy, s, backg)

backgx = backg(:, 1:(size(backg, 2) / 2));
backgy = backg(:, (size(backg, 2) / 2)+1:size(backg, 2));

FIM = zeros(6, 6);

FIM(1, 1) = 0.5 .* (s^2) .* (sum(sum((Bx.aa) ./ (Ix + backgx) + (By.aa) ./ (Iy + backgy))));
FIM(2, 2) = 0.5 .* (s^2) .* (sum(sum((Bx.bb) ./ (Ix + backgx) + (By.bb) ./ (Iy + backgy))));
FIM(3, 3) = 0.5 .* (s^2) .* (sum(sum((Bx.cc) ./ (Ix + backgx) + (By.cc) ./ (Iy + backgy))));
FIM(4, 4) = 0.5 .* (s^2) .* (sum(sum((Bx.dd) ./ (Ix + backgx) + (By.dd) ./ (Iy + backgy))));
FIM(5, 5) = 0.5 .* (s^2) .* (sum(sum((Bx.ee) ./ (Ix + backgx) + (By.ee) ./ (Iy + backgy))));
FIM(6, 6) = 0.5 .* (s^2) .* (sum(sum((Bx.ff) ./ (Ix + backgx) + (By.ff) ./ (Iy + backgy))));

FIM(1, 2) = (s^2) .* (sum(sum((Bx.ab) ./ (Ix + backgx) + (By.ab) ./ (Iy + backgy))));
FIM(1, 3) = (s^2) .* (sum(sum((Bx.ac) ./ (Ix + backgx) + (By.ac) ./ (Iy + backgy))));
FIM(1, 4) = (s^2) .* (sum(sum((Bx.ad) ./ (Ix + backgx) + (By.ad) ./ (Iy + backgy))));
FIM(1, 5) = (s^2) .* (sum(sum((Bx.ae) ./ (Ix + backgx) + (By.ae) ./ (Iy + backgy))));
FIM(1, 6) = (s^2) .* (sum(sum((Bx.af) ./ (Ix + backgx) + (By.af) ./ (Iy + backgy))));

FIM(2, 3) = (s^2) .* (sum(sum((Bx.bc) ./ (Ix + backgx) + (By.bc) ./ (Iy + backgy))));
FIM(2, 4) = (s^2) .* (sum(sum((Bx.bd) ./ (Ix + backgx) + (By.bd) ./ (Iy + backgy))));
FIM(2, 5) = (s^2) .* (sum(sum((Bx.be) ./ (Ix + backgx) + (By.be) ./ (Iy + backgy))));
FIM(2, 6) = (s^2) .* (sum(sum((Bx.bf) ./ (Ix + backgx) + (By.bf) ./ (Iy + backgy))));

FIM(3, 4) = (s^2) .* (sum(sum((Bx.cd) ./ (Ix + backgx) + (By.cd) ./ (Iy + backgy))));
FIM(3, 5) = (s^2) .* (sum(sum((Bx.ce) ./ (Ix + backgx) + (By.ce) ./ (Iy + backgy))));
FIM(3, 6) = (s^2) .* (sum(sum((Bx.cf) ./ (Ix + backgx) + (By.cf) ./ (Iy + backgy))));

FIM(4, 5) = (s^2) .* (sum(sum((Bx.de) ./ (Ix + backgx) + (By.de) ./ (Iy + backgy))));
FIM(4, 6) = (s^2) .* (sum(sum((Bx.df) ./ (Ix + backgx) + (By.df) ./ (Iy + backgy))));

FIM(5, 6) = (s^2) .* (sum(sum((Bx.ef) ./ (Ix + backgx) + (By.ef) ./ (Iy + backgy))));

FIM = FIM + FIM.';
end
