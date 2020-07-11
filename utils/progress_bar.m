function num_char = progress_bar(p, num_char, prect_pos)
try

    cond = p > 1 || p < 0;
catch ME

    error('expected a scalar for precentage')
end


if cond

    error('precentage must lie between 0 and 1!')
end

if num_char == 0


    num_char = fprintf('%c', [repmat(' ', 1, prect_pos - 2), num2str(p * 100), ...
        '%[', repmat('=', 1, ceil(20 * p)), repmat('--', 1, floor(20 - 20 * p)), ']']);


else
    fprintf(repmat('\b', 1, num_char));

    if p == 1


        num_char = fprintf('%c', [repmat(' ', 1, prect_pos - 2), num2str(p * 100), ...
            '%[', repmat('=', 1, ceil(20 * p)), repmat('--', 1, floor(20 - 20 * p)), ']']);
        fprintf('\n')


    else


        num_char = fprintf('%c', [repmat(' ', 1, prect_pos - 2), num2str(p * 100), ...
            '%[', repmat('=', 1, ceil(20 * p)), repmat('--', 1, floor(20 - 20 * p)), ']']);


    end
end
