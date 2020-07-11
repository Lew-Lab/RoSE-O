function [grad_map_avg, gammaf] = grad_map2(gammahat, recovStruct)
%grad_map computes the gradient map from molecular parameter estimates
%(gammahat)
%->---
%input
%->---
%gammahat:     array(12*N,n_f) -molecular parameter estimates
%grid_p:       array(sqrt(N),1) -lateral grid points along x axis
%---->-
%output
%---->-
%gammaf:      array(12N,n_f) -initial estimate of molecular parameters
%grad_map_avg:   array(sqrt(N),sqrt(N),n_f) -gradient maps for of all frames

%------------------------------------------------------------

%% Summary

%------------------------------------------------------------
% In the first step,  gradmaps for each XX, YY and ZZ bases are constructed.
%In the second step, a joint gradmap for identifying molecules is obtained
%by averaging across X, YY and ZZ.
%In the third step, local maxima are extracted from the averaged gradmap
%and parameters of each basis plus the positions of  molecules are computed.
%Lastly, these parameters are mapped to joint parameter space represented
%by gammaf.

%%

%grid points passing through origin along x
grid_p = recovStruct.lateral_grid_p;

%number of grid points along x axis
grid_length = length(grid_p);

%grid points on a x-y plane corresponding to the input image
[gridx, gridy] = meshgrid(grid_p);

%number of grid points
n_grid_p = recovStruct.n_grid_p;

%number of frames
n_frames = recovStruct.num_frames;

%allocate space for output
gammaf = zeros(size(gammahat));
grad_map_avg = zeros(sqrt(n_grid_p), sqrt(n_grid_p), n_frames);

%% looping over the frames for computing gradmap

%==============================
for i = 1:n_frames

    %current molecular estimates
    z = gammahat(:, i);

    %molecular parameters related to XX basis
    z_xx = [z(1:n_grid_p); z(6 * n_grid_p + 1:7 * n_grid_p); z(7 * n_grid_p + 1:8 * n_grid_p)];

    %molecular parameters related to YY basis
    z_yy = [z(n_grid_p + 1:2 * n_grid_p); z(8 * n_grid_p + 1:9 * n_grid_p); z(9 * n_grid_p + 1:10 * n_grid_p)];

    %molecular parameters related to ZZ basis
    z_zz = [z(2 * n_grid_p + 1:3 * n_grid_p); z(10 * n_grid_p + 1:11 * n_grid_p); z(11 * n_grid_p + 1:12 * n_grid_p)];

    %molecular parameters related toXY basis
    z_xy = [z(3 * n_grid_p + 1:4 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %molecular parameters related to XZ basis
    z_xz = [z(4 * n_grid_p + 1:5 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %molecular parameters related to YZ basis
    z_yz = [z(5 * n_grid_p + 1:6 * n_grid_p); zeros(n_grid_p, 1); zeros(n_grid_p, 1)];

    %extract brightness-scaled second moments estimates (XX, YY and ZZ)
    z_secondMoments = reshape(z(1:3 * n_grid_p), sqrt(n_grid_p), sqrt(n_grid_p), 3);

    %determine brightness thresholds in each basis
    br_th = max(max(z_secondMoments)) * .1;

    %get the indices of pixels with low brightness
    indxx_low = (z_secondMoments < 1) + (z_secondMoments < repmat(br_th, sqrt(n_grid_p), sqrt(n_grid_p), 1));

    %set  indices with low brightness to zero
    z_secondMoments(indxx_low > 0) = 0;

    x_localized_xx_t = [];
    y_localized_xx_t = [];
    x_localized_yy_t = [];
    y_localized_yy_t = [];
    x_localized_zz_t = [];
    y_localized_zz_t = [];
    x_localized_xy_t = [];
    y_localized_xy_t = [];
    x_localized_xz_t = [];
    y_localized_xz_t = [];
    x_localized_yz_t = [];
    y_localized_yz_t = [];
    br_xx_t = [];
    br_yy_t = [];
    br_zz_t = [];
    br_xy_t = [];
    br_xz_t = [];
    br_yz_t = [];

    %checking if current frame contains molecular parameters
    if any(z_secondMoments(:) > 0)

        %find potential grid points for each basis
        %-----------------------------------------------------
        %XX basis
        [I, J] = find(z_secondMoments(:, :, 1) > 0);
        loc_xx = [J, I];
        %YY basis
        [I, J] = find(z_secondMoments(:, :, 2) > 0);
        loc_yy = [J, I];
        %ZZ basis
        [I, J] = find(z_secondMoments(:, :, 3) > 0);
        loc_zz = [J, I];

        %         [I,J]=find(z_2D>0);
        %         loc=[J,I];

        %compute gradmap for each basis XX, YY and ZZ
        %===========================
        grad_maps_xx = ComputeGradMap(loc_xx, z_xx);
        grad_maps_yy = ComputeGradMap(loc_yy, z_yy);
        grad_maps_zz = ComputeGradMap(loc_zz, z_zz);

        grad_maps_xx(grad_maps_xx < 0) = 0;
        grad_maps_yy(grad_maps_yy < 0) = 0;
        grad_maps_zz(grad_maps_zz < 0) = 0;

        %average
        grad_map_avg_t = (grad_maps_xx + grad_maps_yy + grad_maps_zz) / 3;
        grad_map_avg(:, :, i) = grad_map_avg_t;
        grad_map_avg_t = grad_map_avg_t / max(grad_map_avg_t(:)); % so that the most likely grid point has value of 1

        % get rid of edge pixels
        indx_edge = [1:3, size(grad_map_avg_t, 1) - 3:size(grad_map_avg_t, 1)];
        grad_map_avg_t(indx_edge, :, :) = 0;
        grad_map_avg_t(:, indx_edge, :) = 0;

        %allocate an array for position of local maxima

        %find local maxima
        %===========================
        if any(grad_map_avg_t(:) > 0.1) % filtering pixels with low likelihood of blinking

            [I, J] = find(grad_map_avg_t > 0.1);


            %  loop over indices in the current gradmap
            for ii = 1:length(I)

                % compute local gradmap with 2 neighbour pixels
                grad_map_avg_t_local = grad_map_avg_t(I(ii)-3:I(ii)+3, J(ii)-3:J(ii)+3);

                if (grad_map_avg_t(I(ii), J(ii)) == max(grad_map_avg_t_local(:))) % check if current index is a local max


                    loc_t = [J(ii), I(ii)];

                    %compute position and brightness in each basis
                    [~, x_localized_xx, y_localized_xx, br_localized_xx] = ComputeGradMap(loc_t, z_xx);
                    [~, x_localized_yy, y_localized_yy, br_localized_yy] = ComputeGradMap(loc_t, z_yy);
                    [~, x_localized_zz, y_localized_zz, br_localized_zz] = ComputeGradMap(loc_t, z_zz);


                    [~, x_localized_xy, y_localized_xy, br_localized_xy] = ComputeGradMap(loc_t, z_xy, 'crossTerms', true);
                    [~, x_localized_xz, y_localized_xz, br_localized_xz] = ComputeGradMap(loc_t, z_xz, 'crossTerms', true);
                    [~, x_localized_yz, y_localized_yz, br_localized_yz] = ComputeGradMap(loc_t, z_yz, 'crossTerms', true);


                    x_localized_xx_t = [x_localized_xx_t, x_localized_xx];
                    y_localized_xx_t = [y_localized_xx_t, y_localized_xx];
                    x_localized_yy_t = [x_localized_yy_t, x_localized_yy];
                    y_localized_yy_t = [y_localized_yy_t, y_localized_yy];
                    x_localized_zz_t = [x_localized_zz_t, x_localized_zz];
                    y_localized_zz_t = [y_localized_zz_t, y_localized_zz];
                    x_localized_xy_t = [x_localized_xy_t, x_localized_xy];
                    y_localized_xy_t = [y_localized_xy_t, y_localized_xy];
                    x_localized_xz_t = [x_localized_xz_t, x_localized_xz];
                    y_localized_xz_t = [y_localized_xz_t, y_localized_xz];
                    x_localized_yz_t = [x_localized_yz_t, x_localized_yz];
                    y_localized_yz_t = [y_localized_yz_t, y_localized_yz];


                    %get the brightness-scaled parameters in each channel
                    br_xx_t = [br_xx_t, br_localized_xx];
                    br_yy_t = [br_yy_t, br_localized_yy];
                    br_zz_t = [br_zz_t, br_localized_zz];
                    br_xy_t = [br_xy_t, br_localized_xy];
                    br_xz_t = [br_xz_t, br_localized_xz];
                    br_yz_t = [br_yz_t, br_localized_yz];

                end
            end
        end
    end

    %mapping estimated molecular parameters into gamma
    %=================================================
    if ~isempty(x_localized_xx_t) % if found a local max (source pixel) with high  likelihood

        for ll = 1:numel(x_localized_xx_t)

            br_t = br_xx_t(ll) + br_yy_t(ll) + br_zz_t(ll);

            x_localized = (br_xx_t(ll) * x_localized_xx_t(ll) + ...
                br_yy_t(ll) * x_localized_yy_t(ll) + br_zz_t(ll) * x_localized_zz_t(ll)) / br_t;

            y_localized = (br_xx_t(ll) * y_localized_xx_t(ll) + ...
                br_yy_t(ll) * y_localized_yy_t(ll) + br_zz_t(ll) * y_localized_zz_t(ll)) / br_t;

            distance = (gridx - x_localized).^2 + (gridy - y_localized).^2;
            %         distance_yy=(gridx-x_localized_yy_t(ll)).^2+(gridy-y_localized_yy_t(ll)).^2;
            %         distance_zz=(gridx-x_localized_zz_t(ll)).^2+(gridy-y_localized_zz_t(ll)).^2;
            %         distance_xy=(gridx-x_localized_xy_t(ll)).^2+(gridy-y_localized_xy_t(ll)).^2;
            %         distance_xz=(gridx-x_localized_xz_t(ll)).^2+(gridy-y_localized_xz_t(ll)).^2;
            %         distance_yz=(gridx-x_localized_yz_t(ll)).^2+(gridy-y_localized_yz_t(ll)).^2;
            %
            [gy, gx] = find(distance == min(distance(:)), 1);
            %         [gy_yy,gx_yy]=find(distance_yy==min(distance_yy(:)),1);
            %         [gy_zz,gx_zz]=find(distance_zz==min(distance_zz(:)),1);
            %         [gy_xy,gx_xy]=find(distance_xy==min(distance_xy(:)),1);
            %         [gy_xz,gx_xz]=find(distance_xz==min(distance_xz(:)),1);
            %         [gy_yz,gx_yz]=find(distance_yz==min(distance_yz(:)),1);
            %
            indx_lin = sub2ind(size(distance), gy, gx);
            %         indx_lin_yy=sub2ind(size(distance_yy),gy_yy,gx_yy);
            %         indx_lin_zz=sub2ind(size(distance_zz),gy_zz,gx_zz);
            %         indx_lin_xy=sub2ind(size(distance_xy),gy_xy,gx_xy);
            %         indx_lin_xz=sub2ind(size(distance_xz),gy_xz,gx_xz);
            %         indx_lin_yz=sub2ind(size(distance_xz),gy_yz,gx_yz);

            %indx_lin=find(distance==min(distance(:)),1);

            deltax_xx = x_localized_xx_t(ll) - gridx(1, gx);
            deltay_xx = y_localized_xx_t(ll) - gridy(gy, 1);

            gammaf(indx_lin, i) = br_xx_t(ll);
            gammaf(indx_lin+6*n_grid_p, i) = br_xx_t(ll) * deltax_xx * 10^-2;
            gammaf(indx_lin+7*n_grid_p, i) = br_xx_t(ll) * deltay_xx * 10^-2;

            deltax_yy = x_localized_yy_t(ll) - gridx(1, gx);
            deltay_yy = y_localized_yy_t(ll) - gridy(gy, 1);

            gammaf(indx_lin+n_grid_p, i) = br_yy_t(ll);
            gammaf(indx_lin+8*n_grid_p, i) = br_yy_t(ll) * deltax_yy * 10^-2;
            gammaf(indx_lin+9*n_grid_p, i) = br_yy_t(ll) * deltay_yy * 10^-2;


            deltax_zz = x_localized_zz_t(ll) - gridx(1, gx);
            deltay_zz = y_localized_zz_t(ll) - gridy(gy, 1);

            gammaf(indx_lin+2*n_grid_p, i) = br_zz_t(ll);
            gammaf(indx_lin+10*n_grid_p, i) = br_zz_t(ll) * deltax_zz * 10^-2;
            gammaf(indx_lin+11*n_grid_p, i) = br_zz_t(ll) * deltay_zz * 10^-2;

            %             deltax_xy=x_localized_xy(ll)-gridx(1,gx_xy);
            %             deltay_xy=y_localized_xy(ll)-gridy(gy_xy,1);
            gammaf(indx_lin+3*n_grid_p, i) = br_xy_t(ll);

            %             deltax_xz=x_localized_xz(ll)-gridx(1,gx_xz);
            %             deltay_xz=y_localized_xz(ll)-gridy(gy_xz,1);

            gammaf(indx_lin+4*n_grid_p, i) = br_xz_t(ll);

            %             deltax_yz=x_localized_yz(ll)-gridx(1,gx_yz);
            %             deltay_yz=y_localized_yz(ll)-gridy(gy_yz,1);

            gammaf(indx_lin+5*n_grid_p, i) = br_yz_t(ll);


            %         end
            %     end
            % end
        end
    end
end

%% local functions

%============================

%computing gradmap
%-----------------------------------------------------

    function [grad_map_t, x_localized, y_localized, br_localized] = ComputeGradMap(loc, z, varargin)

        %allocate space for GradMap
        grad_map_t = zeros(sqrt(n_grid_p), sqrt(n_grid_p)); % current gradient map

        %create localization parameters
        x_localized = [];
        y_localized = [];
        br_localized = [];
        %loop over indices
        indx_t = (loc(:, 1) < grid_length) .* (loc(:, 1) > 1) .* ...
            (loc(:, 2) < grid_length) .* (loc(:, 2) > 1);
        loc = loc(indx_t > 0, :);

        for k = 1:numel(loc(:, 1))


            %get the neighbour pixels
            %--------------------------------------------------
            % center point
            pm = [grid_p(loc(k, 1)), grid_p(loc(k, 2))]; % continous location
            im = (loc(k, 1) - 1) * grid_length + loc(k, 2); % index
            % left point
            pl = [grid_p(loc(k, 1) - 1), grid_p(loc(k, 2))];
            il = (loc(k, 1) - 2) * grid_length + loc(k, 2);
            % upper left point
            p_lup = [grid_p(loc(k, 1) - 1), grid_p(loc(k, 2) - 1)];
            i_lup = (loc(k, 1) - 2) * grid_length + loc(k, 2) - 1;
            % lower down point
            p_ld = [grid_p(loc(k, 1) - 1), grid_p(loc(k, 2) + 1)];
            i_ld = (loc(k, 1) - 2) * grid_length + loc(k, 2) + 1;
            % right point
            pr = [grid_p(loc(k, 1) + 1), grid_p(loc(k, 2))];
            ir = (loc(k, 1)) * grid_length + loc(k, 2);
            %upper right point
            p_rup = [grid_p(loc(k, 1) + 1), grid_p(loc(k, 2) - 1)];
            i_rup = (loc(k, 1)) * grid_length + loc(k, 2) - 1;
            % lower right point
            p_rd = [grid_p(loc(k, 1) + 1), grid_p(loc(k, 2) + 1)];
            i_rd = (loc(k, 1)) * grid_length + loc(k, 2) + 1;
            % lower  point
            pd = [grid_p(loc(k, 1)), grid_p(loc(k, 2) + 1)];
            id = (loc(k, 1) - 1) * grid_length + loc(k, 2) + 1;
            % upper point
            pu = [grid_p(loc(k, 1)), grid_p(loc(k, 2) - 1)];
            iu = (loc(k, 1) - 1) * grid_length + loc(k, 2) - 1;

            %compute convergence of neighbouring pixels to the center pixel
            %--------------------------------------------------
            %get the brightness and position of center pixel
            s = opt2struct(varargin);

            if isfield(s, 'crossterms') && s.crossterms
                positivityCond = false;
            else
                positivityCond = true;
            end

            if positivityCond

                if z(im) > 0

                    br_m = z(im);
                    x_c = pm(1, 1);
                    y_c = pm(1, 2);
                    imx = pm(1, 1) + (z(im + n_grid_p) / z(im)) * 10^2;
                    imy = pm(1, 2) + (z(im + 2 * n_grid_p) / z(im)) * 10^2;

                else
                    z(im) = 0;
                    x_c = pm(1, 1);
                    y_c = pm(1, 2);
                    br_m = 0;
                    imx = 0;
                    imy = 0;
                end

            else
                br_m = z(im);
                x_c = pm(1, 1);
                y_c = pm(1, 2);
                imx = 0;
                imy = 0;
            end


            % left point
            if positivityCond

                if (z(il) > 0)
                    dx = 10^2 * z(il+n_grid_p) / z(il); % in nm
                    dy = 10^2 * z(il+2*n_grid_p) / z(il); % in nm
                    ilx = pl(1, 1) + dx;
                    ily = pl(1, 2) + dy;
                    br_l = z(il);
                    x_i = pl(1);
                    y_i = pl(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];

                    % cosine of angle between the center grid point and
                    % the position vector of molecule at current grid point
                    theta_i = (G * r) / (norm(G) * norm(r));
                    % weight according to the brightness of molecule at current grid point
                    R_l = theta_i * br_l;
                else
                    br_l = 0;
                    R_l = 0;
                    ilx = 0;
                    ily = 0;
                end

            else
                br_l = z(il);
                R_l = 0;
                ilx = 0;
                ily = 0;
            end


            % upper left point

            if positivityCond

                if (z(i_lup) > 0)
                    dx = 10^2 * z(i_lup+n_grid_p) / z(i_lup);
                    dy = 10^2 * z(i_lup+2*n_grid_p) / z(i_lup);
                    i_lup_x = p_lup(1, 1) + dx;
                    i_lup_y = p_lup(1, 2) + dy;
                    br_lup = z(i_lup);
                    x_i = p_lup(1);
                    y_i = p_lup(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_lup = theta_i * br_lup;
                else
                    br_lup = 0;
                    R_lup = 0;
                    i_lup_x = 0;
                    i_lup_y = 0;
                end

            else
                br_lup = z(i_lup);
                R_lup = 0;
                i_lup_x = 0;
                i_lup_y = 0;
            end

            % lower  point

            if positivityCond

                if (z(i_ld) > 0)
                    dx = 10^2 * z(i_ld+n_grid_p) / z(i_ld);
                    dy = 10^2 * z(i_ld+2*n_grid_p) / z(i_ld);
                    i_ld_x = p_ld(1, 1) + dx;
                    i_ld_y = p_ld(1, 2) + dy;
                    br_ld = z(i_ld);
                    x_i = p_ld(1);
                    y_i = p_ld(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_ld = theta_i * br_ld;
                else

                    br_ld = 0;
                    R_ld = 0;
                    i_ld_x = 0;
                    i_ld_y = 0;
                end

            else
                br_ld = z(i_ld);
                R_ld = 0;
                i_ld_x = 0;
                i_ld_y = 0;
            end

            %upper right point

            if positivityCond

                if (z(i_rup) > 0)
                    dx = 10^2 * z(i_rup+n_grid_p) / z(i_rup);
                    dy = 10^2 * z(i_rup+2*n_grid_p) / z(i_rup);
                    i_rup_x = p_rup(1, 1) + dx;
                    i_rup_y = p_rup(1, 2) + dy;
                    br_rup = z(i_rup);
                    x_i = p_rup(1);
                    y_i = p_rup(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_rup = theta_i * br_rup;
                else
                    i_rup_x = 0;
                    i_rup_y = 0;
                    br_rup = 0;
                    R_rup = 0;
                end

            else
                i_rup_x = 0;
                i_rup_y = 0;
                br_rup = z(i_rup);
                R_rup = 0;
            end

            % lower right point
            if positivityCond

                if (z(i_rd) > 0)
                    dx = 10^2 * z(i_rd+n_grid_p) / z(i_rd);
                    dy = 10^2 * z(i_rd+2*n_grid_p) / z(i_rd);
                    i_rd_x = p_rd(1, 1) + dx;
                    i_rd_y = p_rd(1, 2) + dy;

                    br_rd = z(i_rd);
                    x_i = p_rd(1);
                    y_i = p_rd(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_rd = theta_i * br_rd;
                else
                    i_rd_x = 0;
                    i_rd_y = 0;
                    br_rd = 0;
                    R_rd = 0;
                end
            else
                i_rd_x = 0;
                i_rd_y = 0;
                br_rd = z(i_rd);
                R_rd = 0;
            end

            % right point
            if positivityCond

                if (z(ir) > 0)
                    dx = 10^2 * z(ir+n_grid_p) / z(ir);
                    dy = 10^2 * z(ir+2*n_grid_p) / z(ir);
                    irx = pr(1, 1) + dx;
                    iry = pr(1, 2) + dy;
                    br_r = z(ir);
                    x_i = pr(1);
                    y_i = pr(2);

                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_r = theta_i * br_r;
                else
                    irx = 0;
                    iry = 0;
                    br_r = 0;
                    R_r = 0;
                end
            else

                irx = 0;
                iry = 0;
                br_r = z(ir);
                R_r = 0;
            end

            % lower  point
            if positivityCond

                if (z(id) > 0)
                    dx = 10^2 * z(id+n_grid_p) / z(id);
                    dy = 10^2 * z(id+2*n_grid_p) / z(id);
                    idx = pd(1, 1) + dx;
                    idy = pd(1, 2) + dy;
                    br_d = z(id);
                    x_i = pd(1);
                    y_i = pd(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_d = theta_i * br_d;

                else
                    idx = 0;
                    idy = 0;
                    br_d = 0;
                    R_d = 0;
                end
            else
                idx = 0;
                idy = 0;
                br_d = z(id);
                R_d = 0;
            end

            % upper point
            if positivityCond

                if (z(iu) > 0)
                    dx = 10^2 * z(iu+n_grid_p) / z(iu);
                    dy = 10^2 * z(iu+2*n_grid_p) / z(iu);
                    iux = pu(1, 1) + dx;
                    iuy = pu(1, 2) + dy;
                    br_u = z(iu);
                    x_i = pu(1);
                    y_i = pu(2);
                    r = [x_c - x_i, y_c - y_i]';
                    G = [dx, dy];
                    theta_i = (G * r) / (norm(G) * norm(r));
                    R_u = theta_i * br_u;
                else
                    iux = 0;
                    iuy = 0;
                    br_u = 0;
                    R_u = 0;
                end

            else


                iux = 0;
                iuy = 0;
                br_u = z(iu);
                R_u = 0;
            end

            %compute scource coefficient of  center grid point (or pixel)
            %--------------------------------------------------
            R_v = [R_r, R_l, R_d, R_u, R_rup, R_rd, R_ld, R_lup];

            R = sum(R_v); %/(sum(br_v)+eps);% source coefficiet is the brightness-scaled average of

            %convergence values
            grad_map_t(loc(k, 2), loc(k, 1)) = R;

            cont_locx = [irx, ilx, idx, iux, i_rup_x, i_rd_x, i_ld_x, i_lup_x];
            cont_locy = [iry, ily, idy, iuy, i_rup_y, i_rd_y, i_ld_y, i_lup_y];
            br_v = [br_r, br_l, br_d, br_u, br_rup, br_rd, br_ld, br_lup];

            % compute  positin via weighted averaging
            %--------------------------------------------------
            cx = ((imx) * br_m + sum(1 .* br_v .* cont_locx)) / (eps + br_m + sum(1 .* br_v));
            cy = ((imy) * br_m + sum(1 .* br_v .* cont_locy)) / (eps + br_m + sum(1 .* br_v));

            %contious positions
            x_localized = [x_localized, cx];
            y_localized = [y_localized, cy];

            % molecule brightness estimate by aggregating neighbouring points
            br_localized = [br_localized, br_m + sum(1 .* br_v)];

        end

    end

end
