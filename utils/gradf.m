function g = gradf(gamma, SMLM_img, b, FPSFx, FPSFy, recovStruct)
%gradf computes the gradient of the negative Poisson log-likelihood
%w.r.t. molecular parameters at the current estimate gamma
%->-----
% input
%->-----
% gamma:        array (N*3,n_f) -current estimate of the molecular parameters
% SMLM_img:     array (M,n_f)   -input camera images
% b:            array (M,n_f)   -background

% N:                            -number of grid points
% M:                            -image size in number of camera pixels
% n_f:                          -number of frames
% upsamplef:                    -object space upsample factor
%----->-
% output
%----->-
% g:                            -array (N*3,n_f) gradient of the negative Poisson
%                                log-likelihood w.r.t. gamma

%% 1- global parameters

N = recovStruct.n_grid_p;
M = recovStruct.img_size;
n_f = recovStruct.subframe_l;
upsamplef = recovStruct.upsample_factor;

%number of boundry pixels to guard against artifacts
n_boundry_p = 5;

%% 2- extract x-y channel parameters

%------------------------------------------------------------
%x_channel Fourier transfor of PSFs
FXXx = FPSFx.FXXx;
FYYx = FPSFx.FYYx;
FZZx = FPSFx.FZZx;
FXYx = FPSFx.FXYx;
FXZx = FPSFx.FXZx;
FYZx = FPSFx.FYZx;
%y_channel Fourier transfor of PSFs
FXXy = FPSFy.FXXy;
FYYy = FPSFy.FYYy;
FZZy = FPSFy.FZZy;
FXYy = FPSFy.FXYy;
FXZy = FPSFy.FXZy;
FYZy = FPSFy.FYZy;

% x-y channel images
SMLM_img_x = SMLM_img(1:M^2, :);
SMLM_img_y = SMLM_img(M^2+1:end, :);
% x-y channel background
b_x = b(1:M^2, :);
b_y = b(M^2+1:2*M^2, :);

FXX(:, :, 1) = FXXx;
FXX(:, :, 2) = FXXy;
FYY(:, :, 1) = FYYx;
FYY(:, :, 2) = FYYy;
FZZ(:, :, 1) = FZZx;
FZZ(:, :, 2) = FZZy;

FXY(:, :, 1) = FXYx;
FXY(:, :, 2) = FXYy;
FXZ(:, :, 1) = FXZx;
FXZ(:, :, 2) = FXZy;
FYZ(:, :, 1) = FYZx;
FYZ(:, :, 2) = FYZy;

%% 3- computations are performed in Fourier domain

%------------------------------------------------------------

% c_x_temp=(down_sample(fast_mul_fft_x(reshapeMat(gamma))));
%
% c_y_temp=(down_sample(fast_mul_fft_y(reshapeMat(gamma))));

c_temp = (down_sample(fast_mul_fft(reshapeMat(gamma))));
c_x_temp = c_temp(:, :, 1);
c_y_temp = c_temp(:, :, 2);

c_x = bsxfun(@plus, reshape(c_x_temp, M^2, n_f), b_x);
c_y = bsxfun(@plus, reshape(c_y_temp, M^2, n_f), b_y);

g_x = fast_transp_mul_fft_x((1-SMLM_img_x./c_x));

g_y = fast_transp_mul_fft_y((1-SMLM_img_y./c_y));

g = g_x + g_y;

%% utility functions

%------------------------------------------------------------
    function out = reshapeMat(A)
        out = reshape(A, sqrt(N), sqrt(N), 6, n_f);
    end

    function down_sampled = down_sample(x)
        down_sampled = x(1:upsamplef:end, 1:upsamplef:end, :);
    end


    function out_N = fast_mul_fft(x)

        out_N = real(ifft2(bsxfun(@times, FXX, fft2(xxgrid(x))) + bsxfun(@times, FYY, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZ, fft2(zzgrid(x))) + bsxfun(@times, FXY, fft2(xygrid(x))) + bsxfun(@times, FXZ, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZ, fft2(yzgrid(x)))));


    end

    function out_N2 = fast_mul_fft_x(x)

        out_N2 = real(ifft2(bsxfun(@times, FXXx, fft2(xxgrid(x))) + bsxfun(@times, FYYx, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZx, fft2(zzgrid(x))) + bsxfun(@times, FXYx, fft2(xygrid(x))) + bsxfun(@times, FXZx, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZx, fft2(yzgrid(x)))));
    end
    function out_N2 = fast_mul_fft_y(x)
        out_N2 = real(ifft2(bsxfun(@times, FXXy, fft2(xxgrid(x))) + bsxfun(@times, FYYy, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZy, fft2(zzgrid(x))) + bsxfun(@times, FXYy, fft2(xygrid(x))) + bsxfun(@times, FXZy, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZy, fft2(yzgrid(x)))));
    end


    function out_N1_inN2 = xxgrid(x)
        out_N1_inN2_t = (reshape(x(:, :, 1, :), sqrt(N), sqrt(N), n_f));
        out_N1_inN2 = (padarray(out_N1_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));

    end

    function out_N2_inN2 = yygrid(x)
        out_N2_inN2_t = (reshape(x(:, :, 2, :), sqrt(N), sqrt(N), n_f));
        out_N2_inN2 = (padarray(out_N2_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
    end

    function out_N3_inN2 = zzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 3, :), sqrt(N), sqrt(N), n_f));
        out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
    end

    function out_N3_inN2 = xygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 4, :), sqrt(N), sqrt(N), n_f));
        out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
    end

    function out_N3_inN2 = xzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 5, :), sqrt(N), sqrt(N), n_f));
        out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
    end

    function out_N3_inN2 = yzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 6, :), sqrt(N), sqrt(N), n_f));
        out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
    end


    function out_N3 = fast_transp_mul_fft_x(z)
        z_up = up_samp(z);
        fz = fft2(z_up);
        xx_temp = ifft2(bsxfun(@times, conj(FXXx), fz));
        yy_temp = ifft2(bsxfun(@times, conj(FYYx), fz));
        zz_temp = ifft2(bsxfun(@times, conj(FZZx), fz));
        xy_temp = ifft2(bsxfun(@times, conj(FXYx), fz));
        xz_temp = ifft2(bsxfun(@times, conj(FXZx), fz));
        yz_temp = ifft2(bsxfun(@times, conj(FYZx), fz));

        indx = [1:n_boundry_p, size(xx_temp, 1) - n_boundry_p + 1:size(xx_temp, 1)];
        xx_temp(indx, :) = 0;
        xx_temp(:, indx) = 0;
        yy_temp(indx, :) = 0;
        yy_temp(:, indx) = 0;
        zz_temp(indx, :) = 0;
        zz_temp(:, indx) = 0;
        xy_temp(indx, :) = 0;
        xy_temp(:, indx) = 0;
        xz_temp(indx, :) = 0;
        xz_temp(:, indx) = 0;
        yz_temp(indx, :) = 0;
        yz_temp(:, indx) = 0;

        out_N3 = real([reshape(xx_temp, N, n_f); ...
            reshape(yy_temp, N, n_f); ...
            reshape(zz_temp, N, n_f); ...
            reshape(xy_temp, N, n_f); ...
            reshape(xz_temp, N, n_f); ...
            reshape(yz_temp, N, n_f)]);

        %         out_N3= real([reshape(padarray(xx_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(yy_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(zz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(xy_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(xz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(yz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f)]);


        function out_N1_N3 = up_samp(x)
            out_N1_N3 = zeros(sqrt(N), sqrt(N), n_f);
            out_N1_N3(1:upsamplef:end, 1:upsamplef:end, :) = reshape(x, M, M, n_f);
        end


    end

    function out_N3 = fast_transp_mul_fft_y(z)
        z_up = up_samp(z);
        fz = fft2(z_up);
        xx_temp = ifft2(bsxfun(@times, conj(FXXy), fz));
        yy_temp = ifft2(bsxfun(@times, conj(FYYy), fz));
        zz_temp = ifft2(bsxfun(@times, conj(FZZy), fz));
        xy_temp = ifft2(bsxfun(@times, conj(FXYy), fz));
        xz_temp = ifft2(bsxfun(@times, conj(FXZy), fz));
        yz_temp = ifft2(bsxfun(@times, conj(FYZy), fz));

        indx = [1:n_boundry_p, size(xx_temp, 1) - n_boundry_p + 1:size(xx_temp, 1)];
        xx_temp(indx, :) = 0;
        xx_temp(:, indx) = 0;
        yy_temp(indx, :) = 0;
        yy_temp(:, indx) = 0;
        zz_temp(indx, :) = 0;
        zz_temp(:, indx) = 0;
        xy_temp(indx, :) = 0;
        xy_temp(:, indx) = 0;
        xz_temp(indx, :) = 0;
        xz_temp(:, indx) = 0;
        yz_temp(indx, :) = 0;
        yz_temp(:, indx) = 0;

        out_N3 = real([reshape(xx_temp, N, n_f); ...
            reshape(yy_temp, N, n_f); ...
            reshape(zz_temp, N, n_f); ...
            reshape(xy_temp, N, n_f); ...
            reshape(xz_temp, N, n_f); ...
            reshape(yz_temp, N, n_f)]);


        %         out_N3= real([reshape(padarray(xx_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(yy_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(zz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(xy_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(xz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f);...
        %             reshape(padarray(yz_temp(n_boundry_p+1:end-n_boundry_p,n_boundry_p+1:end-n_boundry_p,:),...
        %             [n_boundry_p,n_boundry_p]),N,n_f)]);
        %

        function out_N1_N3 = up_samp(x)
            out_N1_N3 = zeros(sqrt(N), sqrt(N), n_f);
            out_N1_N3(1:upsamplef:end, 1:upsamplef:end, :) = reshape(x, M, M, n_f);
        end


    end
end