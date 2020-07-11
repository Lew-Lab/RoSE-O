%% background estimation using iterative wavelet transform

function backg = Wavelet_backg_est(SMLM_img, thresh, wavelet_level, wavename, iter)

%allocate space
backg = zeros(size(SMLM_img), 'single');
%get the number of frames
num_frame = size(SMLM_img, 3);

% loop over frames
for i = 1:num_frame

    %temorary variables
    img_t = SMLM_img(:, :, i);
    img_filt = img_t;
    ct = 0;
    %itterate
    while (ct < iter)
        % discrete wavelet transform
        [c, s] = wavedec2(img_filt, wavelet_level, wavename);
        c_t = zeros(size(c));
        c_t(1:s(1)*s(1)*1) = c(1:s(1)*s(1)*1);
        % inverse wavelet transform by only keeping low-freq components
        img_new = waverec2(c_t, s, wavename);

        if thresh > 0
            % filter values above current estimated background level.
            eps = sqrt(abs(img_filt)) / 2;
            ind = img_t > (img_new + eps);
            img_filt(ind) = img_new(ind) + eps(ind);

            % re-fine background
            [c, s] = wavedec2(img_filt, wavelet_level, wavename);
            c_t = zeros(size(c));
            c_t(1:s(1)*s(1)*1) = c(1:s(1)*s(1)*1);
            img_new = waverec2(c_t, s, wavename);
        end
        ct = ct + 1;
    end
    %end of iterations
    backg(:, :, i) = img_new;
end
end