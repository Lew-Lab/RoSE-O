function rescale_phase_mask = config_scale_mask(phase_mask_para, imgPara)
% Generates a function, RescaleM to rescale phase mask to match
% experimental setup
% authors: Tianben Ding and Hesam Maizidi

configed_already = 0;
cashed_phase_mask = [];
cashed_mask_struc = [];
    function [pmask, MaskStruct] = RescaleM()

        if (~configed_already || isempty(cashed_mask_struc))


            maskName = cell2mat(phase_mask_para(1));
            %MaskSize = str2double(cell2mat(PhaseMaskPara(2)));
            rho_max = str2double(cell2mat(phase_mask_para(2)));
            xShift = str2double(cell2mat(phase_mask_para(3)));
            yShift = str2double(cell2mat(phase_mask_para(4)));
            maskRot = str2double(cell2mat(phase_mask_para(5)));
            newMaskName = strcat(maskName, '_biasAdjusted_stan_rescaled-resized_PupilRadius_', ...
                num2str(rho_max), '_Rot', num2str(maskRot), '_xyShift_', num2str(xShift), num2str(yShift));
            FileFormat = cell2mat(phase_mask_para(6));
            IsolPSFDist = str2double(cell2mat(phase_mask_para(7)));

            MaskStruct.maskName = maskName;
            MaskStruct.rho_max = rho_max;
            MaskStruct.xShift = xShift;
            MaskStruct.yShift = yShift;
            MaskStruct.maskRot = maskRot;
            MaskStruct.newMaskName = newMaskName;
            MaskStruct.FileFormat = FileFormat;
            MaskStruct.IsolPSFDist = IsolPSFDist;

            MaskSize = round((imgPara.pixel_size)^-1*imgPara.lambda*imgPara.Mag*rho_max/imgPara.NA);
            if mod(MaskSize, 2) ~= 0
                MaskSize = MaskSize + 1;
            end

            display('rescaling the phase mask to match the experimental parameters ...')


            folder_name = 'phasemask';
            addpath(folder_name);
            switch FileFormat
                case '.bmp'
                    mask = imread([maskName, '.bmp']);
                    % imagesc(mask);axis image
                    % title('Original Phase Mask','FontSize',16)
                    % set(gca,'FontSize',14)

                    bias = double(mask(1, 1));
                    mask = double(mask) - bias;
                    mask_size = size(mask);

                    mask = mask / bias * pi;

                case '.mat'

                    mask = angle(importdata([maskName, '.mat']));
                    % imagesc(mask);axis image
                    % title('Original Phase Mask','FontSize',16)
                    % set(gca,'FontSize',14)

                    bias = double(mask(1, 1));
                    mask = double(mask) - bias;
                    mask_size = size(mask);

            end


            %Pick up the aperture size of the mask
            mask_nonZeroStartV = find(mask(:, mask_size(2) / 2) ~= 0, 1, 'first');
            mask_nonZeroEndV = find(flipud(mask(:, mask_size(2) / 2)) ~= 0, 1, 'first');
            mask_nonZeroStartH = find(mask(mask_size(1) / 2, :) ~= 0, 1, 'first');
            mask_nonZeroEndH = find(flipud(mask(mask_size(2) / 2, :)) ~= 0, 1, 'first');
            mask_aper = mask_size(1) - min((mask_nonZeroStartV-1)+(mask_nonZeroEndV - 1), ...
                (mask_nonZeroStartH - 1)+(mask_nonZeroEndH - 1));

            scaleFactor = rho_max * 2 / mask_aper;
            newMaskSize = ceil(mask_size(1)*scaleFactor);

            if mod(newMaskSize, 2) ~= 0
                newMaskSize = newMaskSize + 1;
                % set the new mask size as even number
            end

            %adjust the phase mask size using nearest neighbor interpolation
            MaskResized = imresize(mask, [newMaskSize, NaN], 'nearest');

            if maskRot ~= 0
                MaskResized = rot90(MaskResized, maskRot);
            end

            if xShift > 0
                MaskResized = [zeros(size(MaskResized, 1), 2 * xShift), MaskResized];
            elseif xShift < 0
                MaskResized = [MaskResized, zeros(size(MaskResized, 1), 2 * -xShift)];
            end
            if yShift > 0
                MaskResized = [zeros(2 * yShift, size(MaskResized, 2)); MaskResized];
            elseif yShift < 0
                MaskResized = [MaskResized; zeros(2 * -yShift, size(MaskResized, 2))];
            end
            if size(MaskResized, 1) < MaskSize
                MaskResized = [zeros(floor((MaskSize-size(MaskResized, 1)) / 2), size(MaskResized, 2)); ...
                    MaskResized; zeros(floor((MaskSize-size(MaskResized, 1)) / 2), size(MaskResized, 2))];
            end
            if size(MaskResized, 2) < MaskSize
                MaskResized = [zeros(size(MaskResized, 1), floor((MaskSize-size(MaskResized, 2)) / 2)), ...
                    MaskResized, zeros(size(MaskResized, 1), floor((MaskSize-size(MaskResized, 2)) / 2))];
            end

            maskNew = MaskResized(size(MaskResized, 1)/2-floor(MaskSize / 2)+1:size(MaskResized, 1)/2+floor(MaskSize / 2), ...
                size(MaskResized, 2)/2-floor(MaskSize / 2)+1:size(MaskResized, 2)/2+floor(MaskSize / 2));

            configed_already = 1;
            cashed_phase_mask = exp(1i*maskNew);
            cashed_mask_struc = MaskStruct;

        else
            display('already scaled ...')
        end
        pmask = cashed_phase_mask;
        MaskStruct = cashed_mask_struc;
    end

rescale_phase_mask = @RescaleM;

end
