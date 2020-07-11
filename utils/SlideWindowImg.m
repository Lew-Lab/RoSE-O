function img_croped = SlideWindowImg(img, initialPixel, finalPixel, varargin)
%SLIDEWINDOWIMG Summary of this function goes here
%   Detailed explanation goes here


s = opt2struct(varargin);
%image size
size_t = size(img);
imgSize = size_t(1);
num_frame = size_t(3:end);

%ROI size

HieghtROI = finalPixel(2) - initialPixel(2) + 1;
WidthROI = finalPixel(1) - initialPixel(1) + 1;

%size of the sliding window
windowSize = 121; %pixels
if imgSize < windowSize

    error(['input image is too small. Consider size of greater that ', num2str(windowSize), 'pixels'])
end

if isfield(s, 'slidevalue')
    slideing_value = s.slidevalue;
else
    slideing_value = 50; %pixels
end

if WidthROI < windowSize || HieghtROI < windowSize

    error('size of the ROI exceeds sliding window size')

end


%number of croped regions
num_croped_img_h = (floor((WidthROI-windowSize) / slideing_value) + 1);

num_croped_img_v = (floor((HieghtROI-windowSize) / slideing_value) + 1);


img_croped = zeros([windowSize, 2 * windowSize, num_croped_img_h * num_croped_img_v, num_frame(1, :)]);

%crop


for i = 1:num_croped_img_v

    indx_y = initialPixel(2) + (i - 1) * slideing_value: ...
        initialPixel(2) + (i - 1) * slideing_value + windowSize - 1;

    for j = 1:num_croped_img_h

        indx_x_L = initialPixel(1) + (j - 1) * slideing_value: ...
            initialPixel(1) + (j - 1) * slideing_value + windowSize - 1;

        indx_x_R = imgSize + initialPixel(1) + (j - 1) * slideing_value: ...
            imgSize + initialPixel(1) + (j - 1) * slideing_value + windowSize - 1;
        img_L = img(indx_y, indx_x_L, :);
        img_R = img(indx_y, indx_x_R, :);

        img_croped(:, :, (i - 1)*num_croped_img_h+j, :) = [img_L, img_R];

    end
end
end
