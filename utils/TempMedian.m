function median_v = TempMedian(imgStack, w_size)
%TempMedian computes the median of intensity fluctations in imgStack per pixel
%over a certain window size (w_size)
%INPUT:
%imgStack:              array(m,m,n_f)- input image stack
%w_size:                scalar- number of frames for computation temporal
%median
%OUTPUT:
%median_val:            array(m,m,n_f)- estimated median values for all
%frames

% validate input
img_size = size(imgStack, 1);
num_frames = size(imgStack, 3);
median_v = zeros(img_size*size(imgStack, 2), num_frames);
if w_size > num_frames
    error('window size exceeds the number of frame in the input stack!')
end

num_w = floor(num_frames/w_size);

% rearrange
imgStack_t = reshape(imgStack, img_size*size(imgStack, 2), num_frames);

mean_t = mean(imgStack_t, 1);


for i = 1:num_w
    median_v_t = median(bsxfun(@times, imgStack_t(:, (i - 1) * w_size + 1:i * w_size), ...
        1 ./ mean_t((i-1) * w_size + 1:i * w_size)), 2);
    % de-normalize
    median_v(:, (i - 1)*w_size+1:i*w_size) = bsxfun(@times, median_v_t, mean_t((i-1) * w_size + 1:i * w_size));
end

if (num_frames / w_size) > num_w

    median_v_t = median(bsxfun(@times, imgStack_t(:, (num_w - 1) * w_size + w_size / 2:end), ...
        1 ./ mean_t((num_w-1) * w_size + w_size / 2:end)), 2);
    numLastSubFrames = rem(num_frames, w_size);
    median_v(:, num_w*w_size+1:end) = bsxfun(@times, median_v_t, ...
        mean_t(end -numLastSubFrames + 1:end));
end

median_v = reshape(median_v, size(imgStack));