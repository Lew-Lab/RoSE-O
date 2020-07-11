function l_max = upper_bound_Lipschitz_cnt(y, b, FPSFx, FPSFy, recovStruct)
%upper_bound_Lipschitz_cnt computes an upper bound on the Lipschitz
%constant of the Poisson negative log-likelihood  TODO: cite the relevent paper
%->---
%input
%->---
%y
%N
%n_f
%b
%M
%up_sample
%---->-
%output
%---->-
%l_max

%x_channel
FXXx = FPSFx.FXXx;
FYYx = FPSFx.FYYx;
FZZx = FPSFx.FZZx;
FXYx = FPSFx.FXYx;
FXZx = FPSFx.FXZx;
FYZx = FPSFx.FYZx;
%y_channel
FXXy = FPSFy.FXXy;
FYYy = FPSFy.FYYy;
FZZy = FPSFy.FZZy;
FXYy = FPSFy.FXYy;
FXZy = FPSFy.FXZy;
FYZy = FPSFy.FYZy;


N = recovStruct.n_grid_p;
n_f = recovStruct.num_frames;
M = recovStruct.img_size;
upsample = recovStruct.upsample_factor;

% b_t=min(b);
l_max = (max(bsxfun(@times, y, 1 ./ b)).^2) .* max(reshape((fast_mul_fft(ones(6 * N, n_f))), 2 * M^2, n_f));


    function out_N1 = down_samp(x)
        out_N1 = x(1:upsample:end, 1:upsample:end, :);
    end
    function out_N2 = fast_mul_fft(x)

        out_N2 = [down_samp(real(ifft2(bsxfun(@times, FXXx, fft2(xxgrid(x))) + bsxfun(@times, FYYx, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZx, fft2(zzgrid(x))) + bsxfun(@times, FXYx, fft2(xygrid(x))) + bsxfun(@times, FXZx, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZx, fft2(yzgrid(x)))))), ...
            down_samp(real(ifft2(bsxfun(@times, FXXy, fft2(xxgrid(x))) + bsxfun(@times, FYYy, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZy, fft2(zzgrid(x))) + bsxfun(@times, FXYy, fft2(xygrid(x))) + bsxfun(@times, FXZy, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZy, fft2(yzgrid(x))))))];
        %
        %    out_N2=down_samp(real(ifft2(bsxfun(@times,FXXx,fft2(xxgrid(x)))+bsxfun(@times,FYYx,fft2(yygrid(x)))+...
        %         bsxfun(@times,FZZx,fft2(zzgrid(x)))+bsxfun(@times,FXYx,fft2(xygrid(x)))+bsxfun(@times,FXZx,fft2(xzgrid(x)))+...
        %         bsxfun(@times,FYZx,fft2(yzgrid(x))))));
        %
        function out_N1_inN2 = xxgrid(x)
            out_N1_inN2_t = (reshape(x(1:N, :), sqrt(N), sqrt(N), n_f));
            out_N1_inN2 = (padarray(out_N1_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));

        end
        function out_N2_inN2 = yygrid(x)
            out_N2_inN2_t = (reshape(x(N + 1:2 * N, :), sqrt(N), sqrt(N), n_f));
            out_N2_inN2 = (padarray(out_N2_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = zzgrid(x)
            out_N3_inN2_t = (reshape(x(2 * N + 1:3 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = xygrid(x)
            out_N3_inN2_t = (reshape(x(3 * N + 1:4 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = xzgrid(x)
            out_N3_inN2_t = (reshape(x(4 * N + 1:5 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end

        function out_N3_inN2 = yzgrid(x)
            out_N3_inN2_t = (reshape(x(5 * N + 1:end, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end

    end
end