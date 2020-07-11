function FIM = computeFIMSecM(nanosopeObj, backg, loc_data)

%% 1- create image

secM.muxx = loc_data(:, 5);
secM.muyy = loc_data(:, 6);
secM.muzz = loc_data(:, 7);
secM.muxy = loc_data(:, 8);
secM.muxz = loc_data(:, 9);
secM.muyz = loc_data(:, 10);

molecule_num = size(loc_data, 1);
s = reshape(loc_data(:, 4), 1, 1, molecule_num);


%function handle for computing images formed on the camera
img = @(bases, moments) bsxfun(@times, bases.XX, reshape(moments.muxx, 1, 1, molecule_num)) ...
    +bsxfun(@times, bases.YY, reshape(moments.muyy, 1, 1, molecule_num)) + ...
    +bsxfun(@times, bases.ZZ, reshape(moments.muzz, 1, 1, molecule_num)) + ...
    bsxfun(@times, bases.XY, reshape(moments.muxy, 1, 1, molecule_num)) + ...
    bsxfun(@times, bases.XZ, reshape(moments.muxz, 1, 1, molecule_num)) + ...
    bsxfun(@times, bases.YZ, reshape(moments.muyz, 1, 1, molecule_num));


bx.XX = nanosopeObj.XXxBasis;
bx.YY = nanosopeObj.YYxBasis;
bx.ZZ = nanosopeObj.ZZxBasis;
bx.XY = nanosopeObj.XYxBasis;
bx.XZ = nanosopeObj.XZxBasis;
bx.YZ = nanosopeObj.YZxBasis;


by.XX = nanosopeObj.XXyBasis;
by.YY = nanosopeObj.YYyBasis;
by.ZZ = nanosopeObj.ZZyBasis;
by.XY = nanosopeObj.XYyBasis;
by.XZ = nanosopeObj.XZyBasis;
by.YZ = nanosopeObj.YZyBasis;

Ix = s .* img(bx, secM);
Iy = s .* img(by, secM);


% cache hadamard product of basis images

B11x = bx.XX .* (bx.XX);
B12x = (bx.XX) .* (bx.YY);
B13x = (bx.XX) .* (bx.ZZ);
B14x = (bx.XX) .* (bx.XY);
B15x = (bx.XX) .* (bx.XZ);
B16x = (bx.XX) .* (bx.YZ);

B22x = (bx.YY) .* (bx.YY);
B23x = (bx.YY) .* (bx.ZZ);
B24x = (bx.YY) .* (bx.XY);
B25x = (bx.YY) .* (bx.XZ);
B26x = (bx.YY) .* (bx.YZ);


B33x = (bx.ZZ) .* (bx.ZZ);
B34x = (bx.ZZ) .* (bx.XY);
B35x = (bx.ZZ) .* (bx.XZ);
B36x = (bx.ZZ) .* (bx.YZ);


B44x = (bx.XY) .* (bx.XY);
B45x = (bx.XY) .* (bx.XZ);
B46x = (bx.XY) .* (bx.YZ);

B55x = (bx.XZ) .* (bx.XZ);
B56x = (bx.XZ) .* (bx.YZ);

B66x = (bx.YZ) .* (bx.YZ);


B11y = (by.XX) .* (by.XX);
B12y = (by.XX) .* (by.YY);
B13y = (by.XX) .* (by.ZZ);
B14y = (by.XX) .* (by.XY);
B15y = (by.XX) .* (by.XZ);
B16y = (by.XX) .* (by.YZ);

B22y = (by.YY) .* (by.YY);
B23y = (by.YY) .* (by.ZZ);
B24y = (by.YY) .* (by.XY);
B25y = (by.YY) .* (by.XZ);
B26y = (by.YY) .* (by.YZ);


B33y = (by.ZZ) .* (by.ZZ);
B34y = (by.ZZ) .* (by.XY);
B35y = (by.ZZ) .* (by.XZ);
B36y = (by.ZZ) .* (by.YZ);


B44y = (by.XY) .* (by.XY);
B45y = (by.XY) .* (by.XZ);
B46y = (by.XY) .* (by.YZ);

B55y = (by.XZ) .* (by.XZ);
B56y = (by.XZ) .* (by.YZ);

B66y = (by.YZ) .* (by.YZ);

%% 2- get the background estimate

bx = backg;
by = backg;

%% 3- compute FIM

FIM = zeros(6, 6, molecule_num);


FIM(1, 1, :) = .5 * (s.^2) .* (sum(sum((B11x./(Ix + bx)) + (B11y ./ (Iy + by)), 1), 2));
FIM(1, 2, :) = (s.^2) .* (sum(sum((B12x./(Ix + bx)) + (B12y ./ (Iy + by)), 1), 2));
FIM(1, 3, :) = (s.^2) .* (sum(sum((B13x./(Ix + bx)) + (B13y ./ (Iy + by)), 1), 2));
FIM(1, 4, :) = (s.^2) .* (sum(sum((B14x./(Ix + bx)) + (B14y ./ (Iy + by)), 1), 2));
FIM(1, 5, :) = (s.^2) .* (sum(sum((B15x./(Ix + bx)) + (B15y ./ (Iy + by)), 1), 2));
FIM(1, 6, :) = (s.^2) .* (sum(sum((B16x./(Ix + bx)) + (B16y ./ (Iy + by)), 1), 2));
FIM(2, 2, :) = .5 * (s.^2) .* (sum(sum((B22x./(Ix + bx)) + (B22y ./ (Iy + by)), 1), 2));
FIM(2, 3, :) = (s.^2) .* (sum(sum((B23x./(Ix + bx)) + (B23y ./ (Iy + by)), 1), 2));
FIM(2, 4, :) = (s.^2) .* (sum(sum((B24x./(Ix + bx)) + (B24y ./ (Iy + by)), 1), 2));
FIM(2, 5, :) = (s.^2) .* (sum(sum((B25x./(Ix + bx)) + (B25y ./ (Iy + by)), 1), 2));
FIM(2, 6, :) = (s.^2) .* (sum(sum((B26x./(Ix + bx)) + (B26y ./ (Iy + by)), 1), 2));
FIM(3, 3, :) = .5 * (s.^2) .* (sum(sum((B33x./(Ix + bx)) + (B33y ./ (Iy + by)), 1), 2));
FIM(3, 4, :) = (s.^2) .* (sum(sum((B34x./(Ix + bx)) + (B34y ./ (Iy + by)), 1), 2));
FIM(3, 5, :) = (s.^2) .* (sum(sum((B35x./(Ix + bx)) + (B35y ./ (Iy + by)), 1), 2));
FIM(3, 6, :) = (s.^2) .* (sum(sum((B36x./(Ix + bx)) + (B36y ./ (Iy + by)), 1), 2));
FIM(4, 4, :) = (s.^2) .* (sum(sum((B44x./(Ix + bx)) + (B44y ./ (Iy + by)), 1), 2));
FIM(4, 5, :) = (s.^2) .* (sum(sum((B45x./(Ix + bx)) + (B45y ./ (Iy + by)), 1), 2));
FIM(4, 6, :) = (s.^2) .* (sum(sum((B46x./(Ix + bx)) + (B46y ./ (Iy + by)), 1), 2));
FIM(5, 5, :) = .5 * (s.^2) .* (sum(sum((B55x./(Ix + bx)) + (B55y ./ (Iy + by)), 1), 2));
FIM(5, 6, :) = (s.^2) .* (sum(sum((B56x./(Ix + bx)) + (B56y ./ (Iy + by)), 1), 2));
FIM(6, 6, :) = .5 .* (s.^2) .* (sum(sum((B66x./(Ix + bx)) + (B66y ./ (Iy + by)), 1), 2));
FIM = permute(FIM, [2, 1, 3]) + FIM;


end
