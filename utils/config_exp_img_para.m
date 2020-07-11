function get_exp_img_para = config_exp_img_para()
% Generates a function, Getpara to request for imaging system parameters
% and cashes the parameters
% author: Hesam Mazidi

ConfigedAlready = 0;
cashed_get_exp_img_para = [];

    function ImgPara = GetPara()
        if (~ConfigedAlready || isempty(cashed_get_exp_img_para))
            display('requesting experimental imaging system parameters ...')

            prompt1 = {; ...
                'maginification', ... % magnification  of the imaging system Mag=ftubelens/fobj/n;
                'numerical apertue', ... % NA of the objective lens
                'camera pixel size', ... % in nm
                'refractive index of the medium object is embedded', ... % e.g., oil
                'central emission wavelength', ...
                'camera pixel unit to photon', ...
                'upsampling ratio of object space to image space'};


            dlgTitle = 'Experimental parameters of the imaging system';
            num_lines = repmat([1, 80], size(prompt1, 2), 1);
            defaultans = {'111.1', ...
                '1.4', ...
                '6.5e-6', ...
                '1.518', ...
                '637e-9', ...
                '.49', ...
                '1'; ...
                };


            cashed_get_exp_img_para = inputdlg(prompt1, dlgTitle, num_lines, defaultans);
            ConfigedAlready = 1;
        else
            display('already have imaging system parameters ...')
        end

        exp_img_para = cashed_get_exp_img_para;


        ImgPara.Mag = str2double((exp_img_para(1)));
        ImgPara.NA = str2double((exp_img_para(2)));
        ImgPara.pixel_size = str2double((exp_img_para(3)));
        ImgPara.n = str2double(cell2mat(exp_img_para(4)));
        ImgPara.lambda = str2double(cell2mat(exp_img_para(5)));
        ImgPara.tophoton = str2double(cell2mat(exp_img_para(6)));
        ImgPara.up_sample = str2double(cell2mat(exp_img_para(7)));


    end

get_exp_img_para = @GetPara;
end
