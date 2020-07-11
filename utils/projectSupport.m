function out = projectSupport(z, recovStruct, supportBases)

n_grid_p = recovStruct.n_grid_p;
subframe_l = recovStruct.subframe_l;

basisGrid = @(x, basisIndx)((reshape(x((basisIndx-1) * n_grid_p + 1:(basisIndx) * n_grid_p, :), ...
    sqrt(n_grid_p), sqrt(n_grid_p))));

for i = 1:subframe_l


    for ii = 1:6

        z_t = basisGrid(z(:, i), ii);

        z_t((supportBases{i, ii} == 0) > 0) = 0;

        if ii < 4
            indx_t = z_t <= 0;

            z_t(indx_t > 0) = 0;
        end
        out_t(:, :, ii) = z_t;

    end

    out(:, i) = reshape(out_t, n_grid_p*6, 1);
end
