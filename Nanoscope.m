classdef Nanoscope
    
    
    properties(Constant)
        pixelUpsample=1;% object-space pixel size upsampling factor
    end
    
    properties
        
        % imaging parameters
        %--------------------------------------------------------
        pixelSize=58.5;% object pixel size (nm)
        numApt=1.4;% numerical aperture
        emissWavelength=524;% central emission wavelength (nm)
        ADcount=.49;% photoelectron conversion
        offset=0;% baseline of camera (in digital units)
        refractiveIndx=1.51;% refractive index of the medium
        imageSize=151; % side length of the region of interest
        
        %===========================
        % phase mask parameters
        %--------------------------------------------------------
        
        
        % phase mask parameters
        %--------------------------------------------------------
        phaseMaskPara=struct('maskName','tri-spot',...%parameters of the phase mask mounted on SLM
            'pupilRadius',40,...
            'x_shift_center_position',0,...
            'y_shift_center_position',0,...
            'maskRotation',0);
        
        %===========================
        % phase mask mounted on SLM
        %---------------------------------------------------------
        
        phaseMask;%phase mask mounted on SLM
        
        %===========================
        %bases images of a dipole
        %---------------------------------------------------------
        
        XXxBasis%XX basis image at x_channel
        XXyBasis%XX basis image at y_channel
        YYxBasis%YY basis image at x_channel
        YYyBasis%YY basis image at y_channel
        ZZxBasis%ZZ basis image at x_channel
        ZZyBasis%ZZ basis image at y_channel
        XYxBasis%XY basis image at x_channel
        XYyBasis%XY basis image at y_channel
        XZxBasis%XZ basis image at x_channel
        XZyBasis%XZ basis image at y_channel
        YZxBasis%YZ basis image at x_channel
        YZyBasis%YZ basis image at y_channel
        brightnessScaling%Brightness scaling factor of a dipole
    end
    
    methods
        
        % constructor
        %--------------------------------------------------------
        function obj=Nanoscope(varargin)
            %Nanoscope initializes an instance of Nanoscope object
            %option-value pairs:
            %                                   pixelSize: scalar- object
            %                                   pixel size (nm)
            %                                   (default=58.5)
            %............................................................
            %                                    numApt:  scalar-  numerical
            %                                    aperture of the microscope
            %                                    (default=1.4)
            %............................................................
            %                                    emissWavelength: scalar-
            %                                    emission wavelength of
            %                                    light radiated from
            %                                    emitters (in nm) (default=637)
            %............................................................
            %                                    offset:     scalar or array
            %                                    (m*m) with m the image size
            %                                    (values of the offset are in the raw
            %                                    camera units) (default=0)
            %............................................................
            %                                    ADcount:  scalar-
            %                                    photoelectron conversion
            %                                    (default=.49)
            %............................................................
            %                                   refractiveIndx: scalar- the
            %                                   refractive index of the
            %                                   imaging medium
            %                                   (default=1.51)
            %............................................................
            %                                   imageSize: scalar- size of
            %                                   the sqaure region of interest to
            %                                   be analyzed
            %............................................................
            %                                   phaseMaskPara: structure
            %                                   with fields:
            %                                           maskname: string- name of
            %                                           the mask
            %............................................................
            %                                           pupilRadius:
            %                                           integer- radius of
            %                                           the pupil in units
            %                                           of SLM pixels
            %                                           (default=40)
            %............................................................
            %                                           x_shift_center_position:
            %                                           scalar- shift of
            %                                           the center pixel of
            %                                           the SLM w.r.t. the
            %                                           pupil in x
            %                                           direction (in units
            %                                           of SLM pixels)
            %                                           (default=0)
            %............................................................
            %                                           y_shift_center_position:
            %                                           scalar- shift of
            %                                           the center pixel of
            %                                           the SLM w.r.t. the
            %                                           pupil in y
            %                                           direction (in units
            %                                           of SLM pixels)
            %                                           (default=0)
            %............................................................
            %                                           maskRotation:
            %                                           integer- rotation of
            %                                           phase mask
            %                                           clockwise (in 90
            %                                           degree units)
            %............................................................
            %                                           zeroorder:
            %                                           array(1,2) in [0,1]
            %                                           models the pupil
            %                                           (for left and
            %                                           right)
            %                                           like
            %                                           (1-zeroorder)*pmask+zeroorder*abs(pmask)
            %                                           that is a weighted
            %                                           sum of pupil and clear apertur
            
            s=opt2struct(varargin);
            
            %=========================
            %intialize imaging parameters based on iput options
            %=========================
            if isfield(s,'pixelsize')
                if (isnumeric(s.pixelsize)&&s.pixelsize>0)
                    obj.pixelSize=s.pixelsize;
                else
                    msg=sprintf(['Expecting a positive numeric type for pixelSize.\n',...
                        'Setting pixelSize property to default %.2f'],obj.pixelSize);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'numapt')
                if (isnumeric(s.numapt)&&s.numapt>0)
                    obj.numApt=s.numapt;
                else
                    msg=sprintf(['Expecting a positive numeric type for numApt.\n',...
                        'Setting numApt property to default %.2f'],obj.numApt);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'emisswavelength')
                if (isnumeric(s.emisswavelength)&&s.emisswavelength>0)
                    obj.emissWavelength=s.emisswavelength;
                else
                    msg=sprintf(['Expecting a positive numeric type for emissWavelength.\n',...
                        'Setting emissWavelength property to default %.2f'],obj.emissWavelength);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'offset')
                if (isnumeric(s.offset)&&all(s.offset(:))>0)
                    obj.offset=s.offset;
                else
                    msg=sprintf(['Expecting a positive numeric type for offset.\n',...
                        'Setting offset property to default %.2f'],obj.offset);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'adcount')
                if (isnumeric(s.adcount)&&s.adcount>0)
                    obj.ADcount=s.adcount;
                else
                    msg=sprintf(['Expecting a positive numeric type for ADcount.\n',...
                        'Setting ADcount property to default %.2f'],obj.ADcount);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'refractiveindx')
                if (isnumeric(s.adcount)&&s.adcount>0)
                    obj.refractiveIndx=s.refractiveindx;
                else
                    msg=sprintf(['Expecting a positive numeric type for refractiveIndx.\n',...
                        'Setting refractiveIndx property to default %.2f'],obj.refractiveIndx);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            if isfield(s,'imagesize')
                if (isnumeric(s.imagesize)&&s.imagesize>0)
                    obj.imageSize=s.imagesize;
                else
                    msg=sprintf(['Expecting a positive numeric type for imageSize.\n',...
                        'Setting imageSize property to default %.2f'],obj.imageSize);
                    warning('Nanoscope:InconsistentInputType',msg)
                end
            end
            
            %display imaging parameters
            %--------------------------------------------------
            ImagingInfo=sprintf(['pixel size: %.2f\n','numerical aperture: %.2f\n',...
                'emission wavelength: %.2f\n','offset: %.2f\n',...
                'AD count: %.2f\n','medium refractive index:  %.2f\n',...
                'image size: %.2f\n'],  obj.pixelSize,  obj.numApt,obj.emissWavelength,...
                obj.offset, obj.ADcount, obj.refractiveIndx,obj.imageSize);
            
            fprintf('%s\n',repmat('=',1,20));
            fprintf('Imaging parameters\n');
            fprintf('%s\n',repmat('=',1,20));
            disp(ImagingInfo)
            
            %=========================
            %initialize phase mask parameters based on input options
            %=========================
            if isfield(s,'phasemaskpara')
                %                     if isfield(s.phasemaskpara,'maskname')
                %                         if (ischar(s.phasemaskpara.maskname))
                %                             obj.phaseMaskPara.maskName=s.phasemaskpara.maskname;
                %                         else
                %                             msg=sprintf(['Expecting a character array type for maskName.\n',...
                %                                 'Setting pupilRadius property to default %s'],obj.phaseMaskPara.maskName);
                %                             warning('Nanoscope:phaseMaskPara:InconsistentInputType',msg)
                %                         end
                %                     end
                %                     if isfield(s.phasemaskpara,'pupilradius')
                %                         if (isnumeric(s.phasemaskpara.pupilradius)&&s.phasemaskpara.pupilradius>0)
                %                             obj.phaseMaskPara.pupilRadius=s.phasemaskpara.pupilradius;
                %                         else
                %                             msg=sprintf(['Expecting a positive numeric type for pupilRadius.\n',...
                %                                 'Setting pupilRadius property to default %.2f', obj.phaseMaskPara.pupilRadius]);
                %                             warning('Nanoscope:phaseMaskPara:InconsistentInputType',msg)
                %                         end
                %                     end
                %                     if isfield(s.phasemaskpara,'x_shift_center_position')
                %                         if (isnumeric(s.phasemaskpara.x_shift_center_position))
                %                             obj.phaseMaskPara.x_shift_center_position=s.phasemaskpara.x_shift_center_position;
                %                         else
                %                             msg=sprintf(['Expecting a  numeric type for x_shift_center_position.\n',...
                %                                 'Setting x_shift_center_position property to default %.2f', obj.phaseMaskPara.x_shift_center_position]);
                %                             warning('Nanoscope:phaseMaskPara:InconsistentInputType',msg)
                %                         end
                %                     end
                %                     if isfield(s.phasemaskpara,'y_shift_center_position')
                %                         if (isnumeric(s.phasemaskpara.y_shift_center_position))
                %                             obj.phaseMaskPara.y_shift_center_position=s.phasemaskpara.y_shift_center_position;
                %                         else
                %                             msg=sprintf(['Expecting a  numeric type for y_shift_center_position.\n',...
                %                                 'Setting y_shift_center_position property to default %.2f', obj.phaseMaskPara.y_shift_center_position]);
                %                             warning('Nanoscope:phaseMaskPara:InconsistentInputType',msg)
                %                         end
                %                     end
                %
                %                      if isfield(s.phasemaskpara,'maskrotation')
                %                         if (isnumeric(s.phasemaskpara.maskrotation))
                %                             obj.phaseMaskPara.maskRotation=s.phasemaskpara.maskrotation;
                %                         else
                %                             msg=sprintf(['Expecting a  numeric type for maskRotation.\n',...
                %                                 'Setting maskRotation property to default %.2f', obj.phaseMaskPara.maskRotation]);
                %                             warning('Nanoscope:phaseMaskPara:InconsistentInputType',msg)
                %                         end
                %                      end
                %             end
                obj.phaseMaskPara=s.phasemaskpara;
            end
            %display phase mask parameters
            %--------------------------------------------------
            phaseMaskInfo=sprintf(['mask name: %s\n','pupil radius: %.2f\n'],...
                obj.phaseMaskPara.maskName,  obj.phaseMaskPara.pupilRadius);
            
            fprintf('%s\n',repmat('=',1,20));
            fprintf('Phase mask parameters\n');
            fprintf('%s\n',repmat('=',1,20));
            disp(phaseMaskInfo)
            %=========================
            
            %phase mask initialization
            
            % apply zero order term
            if isfield(s,'phasemaskpara') && isfield(s.phasemaskpara,'zeroorder')
                
                obj.phaseMask=Nanoscope.mountPhaseMask(obj,'zeroorder',...
                    s.phasemaskpara.zeroorder);
            elseif isfield(s,'phasemaskpara') && isfield(s.phasemaskpara,'zeroorder')&& isfield(s.phasemaskpara,'amplitudemask')
                
                obj.phaseMask=Nanoscope.mountPhaseMask(obj,'zeroorder',...
                    s.phasemaskpara.zeroorder,'amplitudemask',s.phasemaskpara.amplitudemask);
            elseif isfield(s,'phasemaskpara') && isfield(s.phasemaskpara,'amplitudemask')
                
                obj.phaseMask=Nanoscope.mountPhaseMask(obj,'amplitudemask',...
                    s.phasemaskpara.amplitudemask);
            else
                obj.phaseMask=Nanoscope.mountPhaseMask(obj);
            end
            %brightness scaling
            obj.brightnessScaling=obj.brightnessScalingCompute();
            
            %x_channel initialization
            obj.XXxBasis=obj.computeBasis(obj,'XX'...
                ,true,'x_channel',true,'crop',true);
            obj.YYxBasis=obj.computeBasis(obj,'YY'...
                ,true,'x_channel',true,'crop',true);
            obj.ZZxBasis=obj.computeBasis(obj,'ZZ'...
                ,true,'x_channel',true,'crop',true);
            obj.XYxBasis=obj.computeBasis(obj,'XY'...
                ,true,'x_channel',true,'crop',true);
            obj.XZxBasis=obj.computeBasis(obj,'XZ'...
                ,true,'x_channel',true,'crop',true);
            obj.YZxBasis=obj.computeBasis(obj,'YZ'...
                ,true,'x_channel',true,'crop',true);
            
            %y_channel initialization
            obj.XXyBasis=obj.computeBasis(obj,'XX'...
                ,true,'y_channel',true,'crop',true);
            obj.YYyBasis=obj.computeBasis(obj,'YY'...
                ,true,'y_channel',true,'crop',true);
            obj.ZZyBasis=obj.computeBasis(obj,'ZZ'...
                ,true,'y_channel',true,'crop',true);
            obj.XYyBasis=obj.computeBasis(obj,'XY'...
                ,true,'y_channel',true,'crop',true);
            obj.XZyBasis=obj.computeBasis(obj,'XZ'...
                ,true,'y_channel',true,'crop',true);
            obj.YZyBasis=obj.computeBasis(obj,'YZ'...
                ,true,'y_channel',true,'crop',true);
        end
        
        % properties set methods
        %---------------------------------------------------------
        function obj=set.pixelSize(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for pixel size')
                end
                if ~ (val>0)
                    error('Nanoscope:InconsistentInputValue','Expecting a positive value for pixel size')
                    
                end
                obj.pixelSize=val;
            end
            
        end
        
        function obj=set.numApt(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for numerical aperture')
                end
                if ~ (val>0)
                    error('Nanoscope:InconsistentInputValue','Expecting a positive value for numerical aperture')
                    
                end
                obj.numApt=val;
            end
            
            
        end
        
        
        
        function obj=set.emissWavelength(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for emission wavelength')
                end
                if ~ (val>0)
                    error('Nanoscope:InconsistentInputValue','Expecting a positive value for emission wavelength')
                    
                end
                obj.emissWavelength=val;
                
            end
        end
        
        
        function obj=set.refractiveIndx(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for refractive index')
                end
                if ~ (val>0)
                    error('Nanoscope:InconsistentInputValue','Expecting a positive value for refractive index')
                    
                end
                obj.refractiveIndx=val;
            end
            
            
        end
        
        
        function obj=set.ADcount(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for A/D count')
                end
                if ~ (val>0)
                    error('Nanoscope:InconsistentInputValue','Expecting a positive value for A/D count')
                    
                end
                obj.ADcount=val;
            end
            
        end
        
        function obj=set.offset(obj,val)
            if nargin>0
                if ~isnumeric(val)
                    error('Nanoscope:InconsistentInputType','Expecting a numeric type for offset')
                end
                if  any(val<0)
                    error('Nanoscope:InconsistentInputValue','Expecting  positive values for offset')
                    
                end
                obj.offset=val;
            end
            
        end
        
        
        %---------------------------------------------------------
        function obj=set.phaseMaskPara(obj,s)
            
            if nargin>0
                s_t=struct();
                
                if ~isstruct(s)
                    
                    error('Nanoscope:InconsistentInputType','Expecting a structure type as input')
                    
                end
                
                if isfield(s,'maskname')
                    
                    if ~(ischar(s.maskname))
                        
                        error('Nanoscope:InconsistentInputType','Expecting a character array for maskName')
                        
                    end
                    s_t.maskName=s.maskname;
                else
                    
                    s_t.maskName=obj.phaseMaskPara.maskName;
                end
                
                if isfield(s,'pupilradius')
                    
                    if ~((isnumeric(s.pupilradius) && (s.pupilradius>0)))
                        
                        error('Nanoscope:InconsistentInputType','Expecting a positive number for pupilRadius')
                        
                    end
                    s_t.pupilRadius=s.pupilradius;
                    
                else
                    s_t.pupilRadius=obj.phaseMaskPara.pupilRadius;
                end
                
                
                if isfield(s,'x_shift_center_position')
                    
                    if ~(isnumeric(s.x_shift_center_position))
                        
                        error('Nanoscope:InconsistentInputType','Expecting a numeric type for x_shift_center_position')
                    end
                    s_t.x_shift_center_position=s.x_shift_center_position;
                else
                    
                    s_t.x_shift_center_position=obj.phaseMaskPara.x_shift_center_position;
                end
                
                if isfield(s,'y_shift_center_position')
                    
                    if ~(isnumeric(s.y_shift_center_position))
                        
                        error('Nanoscope:InconsistentInputType','Expecting a  numeric type for y_shift_center_position')
                    end
                    s_t.y_shift_center_position=s.y_shift_center_position;
                    
                else
                    
                    s_t.y_shift_center_position=obj.phaseMaskPara.y_shift_center_position;
                    
                end
                
                if isfield(s,'maskrotation')
                    
                    if ~( floor(s.maskrotation)==s.maskrotation)
                        
                        error('Nanoscope:InconsistentInputType','Expecting a positive integer for maskRotation')
                        
                    end
                    s_t.maskRotation=s.maskrotation;
                    
                else
                    
                    s_t.maskRotation=obj.phaseMaskPara.maskRotation;
                    
                end
                
            end
            
            obj.phaseMaskPara=struct(s_t);
        end
        
    end
    
    %% methods for properties initialization
    %--------------------------------------------------------
    
    %%
    methods(Access = protected)
        
        function brightnessScaling=brightnessScalingCompute(obj)
            %brightnessScalingCompute computes the brightness scaling for
            %normalizing each basis (i.e., XX, YY , etc).
            %brightnessScaling corresponds to the YYy  basis image formed on the camera for a
            %dipole with \theta=pi/2 and \phi=pi/2
            
            %set Emitter properties to match YYy channel
            Emitter.polar_para.phiD=pi/2; Emitter.polar_para.thetaD=pi/2;
            Emitter.position_para.x=0;
            Emitter.position_para.y=0;
            Emitter.position_para.z=0;
            
            [~,brightnessScaling]=obj.simDipole_novotny(obj,Emitter);
        end
        
        %%
        function B=computeBasis(obj,varargin)
            %computeBasis computes the bases images (i.e., XX, YY, etc) corresponding to a
            %dipole.
            
            %get options
            opt=varargin(2:end);
            s=opt2struct(opt);
            
            Emitter.position_para.x=0;
            Emitter.position_para.y=0;
            Emitter.position_para.z=0;
            
            simDipole_novotny_h=@(Emitter)obj.simDipole_novotny(obj,Emitter);
            
            if isfield(s,'xx') ||isfield(s,'xy')||isfield(s,'xz')
                % XX
                Emitter.polar_para.phiD=0;
                Emitter.polar_para.thetaD=pi/2;
                
                
                [BXXx,BXXy]=simDipole_novotny_h(Emitter);
                
                
                if isfield(s,'x_channel')&& isfield(s,'xx')
                    
                    B=BXXx;
                    
                else
                    B=BXXy;
                    
                end
                
            end
            
            
            % YY
            Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.thetaD=pi/2;
            
            if isfield(s,'yy') ||isfield(s,'xy')||isfield(s,'yz')
                [BYYx,BYYy]=simDipole_novotny_h(Emitter);
                if isfield(s,'x_channel')&& isfield(s,'yy')
                    
                    B=BYYx;
                    
                else
                    B=BYYy;
                    
                end
                
            end
            
            
            % ZZ
            %             Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.phiD=0;
            Emitter.polar_para.thetaD=0;
            
            if isfield(s,'zz') ||isfield(s,'xz')||isfield(s,'yz')
                [BZZx,BZZy]=simDipole_novotny_h(Emitter);
                if isfield(s,'x_channel')&& isfield(s,'zz')
                    B=BZZx;
                    
                else
                    B=BZZy;
                    
                end
                
            end
            
            
            % XY
            Emitter.polar_para.phiD=pi/4;
            Emitter.polar_para.thetaD=pi/2;
            
            if isfield(s,'xy')
                [BXYxt,BXYyt]=simDipole_novotny_h(Emitter);
                if isfield(s,'x_channel')
                    B = 2*BXYxt - BXXx - BYYx;
                    
                else
                    B = 2*BXYyt - BXXy - BYYy;
                    
                end
                
            end
            
            % XZ
            Emitter.polar_para.phiD=0;
            Emitter.polar_para.thetaD=pi/4;
            
            if isfield(s,'xz')
                [BXZxt,BXZyt]=simDipole_novotny_h(Emitter);
                if isfield(s,'x_channel')
                    B = 2*BXZxt - BXXx - BZZx;
                else
                    B=2*BXZyt - BXXy - BZZy;
                end
                
            end
            % YZ
            Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.thetaD=pi/4;
            
            if isfield(s,'yz')
                [BYZxt,BYZyt]=simDipole_novotny_h(Emitter);
                if isfield(s,'x_channel')
                    
                    B=2*BYZxt - BYYx- BZZx;
                else
                    B = 2*BYZyt - BYYy - BZZy;
                end
            end
            
            % crop the basis images to match the desired image size
            %--------------------------------------------------------
            
            %accounting for photon loss
            brightness_scaling=obj.brightnessScaling;
            
            N_pupil=size(B,1);
            up_sample=obj.pixelUpsample;
            img_size=obj.imageSize;
            
            %handle for corping region of interest
            roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,....
                -up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,:);
            
            if isfield(s,'crop') && s.crop
                sumnorm=sum(sum(roi(brightness_scaling)));
                
                B=roi(B)/sumnorm;
                
            else
                sumnorm=sum(sum(brightness_scaling));
                B=B/sumnorm;
            end
            
        end
    end
    
    %%
    
    methods(Static)
        %---------------------------------------------------------
        function phaseMask_out=mountPhaseMask(obj,varargin)
            s=opt2struct(varargin);
            rho_max=obj.phaseMaskPara.pupilRadius;
            xShift=obj.phaseMaskPara.x_shift_center_position;
            yShift=obj.phaseMaskPara.y_shift_center_position;
            maskRot=obj.phaseMaskPara.maskRotation;
            maskName=obj.phaseMaskPara.maskName;
            
            maskSize=round((obj.pixelSize)^-1*(obj.emissWavelength*rho_max)/obj.numApt);
            if mod(maskSize,2)~=0
                maskSize=maskSize+1;
            end
            
            pmask = imread([maskName '.bmp']);
            bias = double(pmask(1,1));
            pmask = double(pmask) - bias;
            pmaskSize = size(pmask);
            
            pmask = pmask/bias*pi;
            
            
            %Pick up the aperture size of the mask
            mask_nonZeroStartV = find(pmask(:,pmaskSize(2)/2)~=0,1,'first');
            mask_nonZeroEndV = find(flipud(pmask(:,pmaskSize(2)/2))~=0,1,'first');
            mask_nonZeroStartH = find(pmask(pmaskSize(1)/2,:)~=0,1,'first');
            mask_nonZeroEndH = find(flipud(pmask(pmaskSize(2)/2,:))~=0,1,'first');
            mask_aper = pmaskSize(1)-min((mask_nonZeroStartV-1)+(mask_nonZeroEndV-1),...
                (mask_nonZeroStartH-1)+(mask_nonZeroEndH-1));
            
            scaleFactor = rho_max*2 / mask_aper;
            newMaskSize = ceil(pmaskSize(1)*scaleFactor);
            
            if mod(newMaskSize,2) ~= 0
                newMaskSize = newMaskSize + 1;
                % set the new mask size as even number
            end
            
            %adjust the phase mask size using nearest neighbor interpolation
            MaskResized = imresize(pmask,[newMaskSize NaN],'nearest');
            
            %amplitude modulation
            if isfield(s,'amplitudemask') 
                nameAmpMask=s.amplitudemask;
                amplitude_mask_struct=load(fullfile('phasemask',nameAmpMask));
                amplitude_mask = imresize(amplitude_mask_struct.xAmpMask_padded,[newMaskSize NaN],'nearest');
                
            else
                amplitude_mask=ones(size(MaskResized));
            end
            
            if maskRot ~= 0
                MaskResized = rot90(MaskResized,maskRot);
                amplitude_mask = rot90(amplitude_mask,maskRot);
            end
            
            if xShift > 0
                MaskResized = [zeros(size(MaskResized,1),2*xShift) MaskResized];
                amplitude_mask = [zeros(size(amplitude_mask,1),2*xShift) amplitude_mask];

            elseif xShift < 0
                MaskResized = [MaskResized zeros(size(MaskResized,1),2*-xShift)];
                amplitude_mask = [amplitude_mask zeros(size(amplitude_mask,1),2*-xShift)];
            end
            if yShift > 0
                MaskResized = [zeros(2*yShift,size(MaskResized,2)); MaskResized];
                amplitude_mask = [zeros(2*yShift,size(amplitude_mask,2)); amplitude_mask];
            elseif yShift < 0
                MaskResized = [MaskResized; zeros(2*-yShift,size(MaskResized,2))];
                amplitude_mask = [amplitude_mask; zeros(2*-yShift,size(amplitude_mask,2))];
            end
            if size(MaskResized,1) < maskSize
                MaskResized = [zeros(floor((maskSize-size(MaskResized,1))/2),size(MaskResized,2));...
                MaskResized;zeros(floor((maskSize-size(MaskResized,1))/2),size(MaskResized,2))];
            
                amplitude_mask = [zeros(floor((maskSize-size(amplitude_mask,1))/2),size(amplitude_mask,2));...
                amplitude_mask;zeros(floor((maskSize-size(amplitude_mask,1))/2),size(amplitude_mask,2))];            
            end
            if size(MaskResized,2) < maskSize
                MaskResized = [zeros(size(MaskResized,1),floor((maskSize-size(MaskResized,2))/2))...
                    MaskResized zeros(size(MaskResized,1),floor((maskSize-size(MaskResized,2))/2))];
                
                amplitude_mask = [zeros(size(amplitude_mask,1),floor((maskSize-size(amplitude_mask,2))/2))...
                    amplitude_mask zeros(size(amplitude_mask,1),floor((maskSize-size(amplitude_mask,2))/2))];                
            end
            
            maskNew = MaskResized(size(MaskResized,1)/2-floor(maskSize/2)+1:size(MaskResized,1)/2+floor(maskSize/2),...
                size(MaskResized,2)/2-floor(maskSize/2)+1:size(MaskResized,2)/2+floor(maskSize/2));
            
            amplitude_maskNew = amplitude_mask(size(amplitude_mask,1)/2-floor(maskSize/2)+1:size(amplitude_mask,1)/2+floor(maskSize/2),...
                size(amplitude_mask,2)/2-floor(maskSize/2)+1:size(amplitude_mask,2)/2+floor(maskSize/2));
            
            %the initil amplitude was inverse of the true amplitude!
            phaseMask=(amplitude_maskNew).*exp(1i*maskNew);
            
            if isfield(s,'zeroorder')
                phaseMask_x=(1-s.zeroorder(1))*phaseMask+s.zeroorder(1)*abs(phaseMask);
                phaseMask_y=(1-s.zeroorder(2))*phaseMask+s.zeroorder(2)*abs(phaseMask);
            else
                phaseMask_x=phaseMask;
                phaseMask_y=phaseMask;
            end
            
            
                
            phaseMask_out(:,:,1)=phaseMask_x;
            phaseMask_out(:,:,2)=phaseMask_y;
            
        end
    end
    
    
    %% methods for image formation
    %----------------------------------------------------
    methods(Static)
        
        %%
        function [imgx,imgy,Ex,Ey] = simDipole_novotny(Nanoscope,Emitter,varargin)
            %simDipole_novotny computes
            
            s=opt2struct(varargin);
            
            % get position parameters
            %--------------------------------------------------------
            fileds={'z','x','y'};
            position_para=Emitter.position_para; % in (nm)
            checkFileds(position_para,fileds); %validate fields
            z=position_para.z;
            deltax=position_para.x;
            deltay=position_para.y;
            
            % get molecular orientation parameters
            %--------------------------------------------------------
            fileds={'phiD','thetaD'};
            polar_para=Emitter.polar_para;
            checkFileds(polar_para,fileds); %validate fields
            phiD=polar_para.phiD;
            thetaD=polar_para.thetaD;
            
            % get emitter and imaging parameters
            %--------------------------------------------------------
            
            n1=Nanoscope.refractiveIndx;
            zh=0;% thickness of film
            z2=0;% distance from the emitter to the interface
            n2=1.33;% sample refractive index
            nh=n1;% thin film refractive index
            lambda = Nanoscope.emissWavelength; %wavelength (nm)
            NA =Nanoscope.numApt; %numerical aperture
            pmask=Nanoscope.phaseMask;
            N=size(pmask,1);
            pixel_size=Nanoscope.pixelSize;%object pixel size (nm)
            
            if isfield(s,'upsampling')
                upsampling=s.upsampling;
            else
                upsampling=1;%upsampling factor of image space
            end
            molecule_num=numel(deltax);
            %calculate both pupil and image plane sampling,
            %one will affect the other, so make sure not to introduce aliasing
            
            dx = n1*(pixel_size/(upsampling));%  in (nm) due to Abbe sine...
            % condition, scale by imaging medium r.i.
            dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
            
            % define pupil coordinates
            %--------------------------------------------------------
            [eta,xi] = meshgrid( ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(2*N*dx))...
                + (1/(2*dx))), ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(N*2*dx)) + (1/(2*dx))));
            x = lambda*(eta);
            y = lambda*(xi);
            [phi,rho] = cart2pol(x,y);
            rho_max = NA/n1;%pupil region of support determined by NA and imaging medium r.i.
            k1 = n1*(2*pi/lambda);
            kh = nh*(2*pi/lambda);
            k2 = n2*(2*pi/lambda);
            theta1 = asin(rho);%theta in matched medium
            thetah = asin((n1/nh)*sin(theta1));%theta in thin film
            theta2 = asin((n1/n2)*sin(theta1));%theta in mismatched medium
            
            % cache fixed terms
            costheta2=cos(theta2);
            costhetah=cos(thetah);
            costheta1=cos(theta1);
            sintheta1=sin(theta1);
            cosphi=cos(phi);
            sinphi=sin(phi);
            % Fresnel coefficients
            %--------------------------------------------------------
            tp_2h = 2*n2*costheta2./(n2*costhetah + nh*costheta2);
            ts_2h = 2*n2*costheta2./(nh*costhetah + n2*costheta2);
            tp_h1 = 2*nh*costhetah./(nh*costheta1 + n1*costhetah);
            ts_h1 = 2*nh*costhetah./(n1*costheta1 + nh*costhetah);
            
            rp_2h = (n2*costheta2 - nh*costhetah)./(n2*costheta2+ nh*costhetah);
            rs_2h = (nh*costheta2 - n2*costhetah)./(nh*costheta2+ n2*costhetah);
            rp_h1 = (nh*costhetah - n1*costheta1)./(nh*costhetah+ n1*costheta1);
            rs_h1 = (n1*costhetah - nh*costheta1)./(n1*costhetah+ nh*costheta1);
            
            % Axelrod's equations for E-fields at pupil plane
            %--------------------------------------------------------
            mux = reshape(sin(thetaD).*cos(phiD),1,1,molecule_num);
            muy = reshape(sin(thetaD).*sin(phiD),1,1,molecule_num);
            muz = reshape(cos(thetaD),1,1,molecule_num);
            tp = tp_2h.*tp_h1.*exp(1i*kh*costhetah*zh)./(1 + rp_2h.*rp_h1.*exp(2i*kh*zh*costhetah));
            ts = ts_2h.*ts_h1.*exp(1i*kh*costhetah*zh)./(1 + rs_2h.*rs_h1.*exp(2i*kh*zh*costhetah));
            
            Es = bsxfun(@times,ts.*(costheta1./costheta2).*(n1/n2),...
                (bsxfun(@times,muy,cosphi) - bsxfun(@times,mux,sinphi)));
            
            Ep = bsxfun(@times,tp,bsxfun(@times,(n1/n2).*costheta1,(bsxfun(@times,mux,cosphi) + ...
                bsxfun(@times,muy,sinphi))) - ...
                bsxfun(@times,bsxfun(@times,muz,sintheta1),(n1/n2)^2.*(costheta1./costheta2)));
            
            PupilFilt=(rho < rho_max).*1;%casting to numeric
            
            
            %computing E-fields
            %--------------------------------------------------------
            
            E_common=PupilFilt.*(1./sqrt(costheta1)).*exp(1i*kh*zh*costhetah).*...
                exp(1i*k2*z2*costheta2);
            Ey_common=bsxfun(@times,E_common,...
                (bsxfun(@times,cosphi,Es) + bsxfun(@times, sinphi,Ep)));
            Ex_common=bsxfun(@times,E_common,...
                (bsxfun(@times,cosphi,Ep) + bsxfun(@times, -sinphi,Es)));
            
            
            Ey_1=exp(1i*k1*(bsxfun(@times,reshape(z,1,1,molecule_num),costheta1))); %defocus term at pupil plane
            Ex_1=exp(1i*k1*(bsxfun(@times,reshape(z,1,1,molecule_num),costheta1))); %defocus term at pupil plane
            
            % aplly channel mismatch
            
            
            Ex_2=exp(1i*k1*(bsxfun(@times,reshape(-deltax,1,1,molecule_num),sintheta1.*cosphi))); %phase shift
            Ex_3=exp(1i*k1*(bsxfun(@times,reshape(deltay,1,1,molecule_num),sintheta1.*sinphi))); %phase shift
            
            if isfield(s,'channel_mismatch')
                deltax=deltax+s.channel_mismatch(1);
                deltay=deltay+s.channel_mismatch(2);
                
            end
            
            Ey_2=exp(1i*k1*(bsxfun(@times,reshape(deltax,1,1,molecule_num),sintheta1.*cosphi))); %phase shift
            Ey_3=exp(1i*k1*(bsxfun(@times,reshape(deltay,1,1,molecule_num),sintheta1.*sinphi))); %phase shift
            
            
            
            Ey_t=Ey_1.*Ey_2.*Ey_3;
            Ex_t=Ex_1.*Ex_2.*Ex_3;
            Ex=bsxfun(@times,Ex_t,Ex_common);
            Ey=bsxfun(@times,Ey_t,Ey_common);
            
            % for propagation from pupil plane E-field to image plane via tube-lens, paraxial
            % approximation is in force.
            %--------------------------------------------------------
            
            pmask_x=pmask(:,:,1);
            pmask_y=pmask(:,:,2);
            pmaskRot_y=rot90(pmask_y,2);
            pmaskRot_x=rot90(pmask_x,1);
            imgy = (fftshift(fft2(ifftshift(Ey.*repmat((pmaskRot_y),1,1,molecule_num)))));
            imgx= fliplr((fftshift(fft2(ifftshift(Ex.*repmat((pmaskRot_x),1,1,molecule_num))))));
            % image on the camera is the amplitude squared of the electric field
            %--------------------------------------------------------
            imgy = abs(imgy).^2;
            imgx = abs(imgx).^2;
        end
        %%
        function [imgx,imgy,Ex,Ey] = simDipole_novotny_costum(Nanoscope,Emitter,pmask,varargin)
            %simDipole_novotny computes
            
            s=opt2struct(varargin);
            
            % get position parameters
            %--------------------------------------------------------
            fileds={'z','x','y'};
            position_para=Emitter.position_para; % in (nm)
            checkFileds(position_para,fileds); %validate fields
            z=position_para.z;
            deltax=position_para.x;
            deltay=position_para.y;
            
            % get molecular orientation parameters
            %--------------------------------------------------------
            fileds={'phiD','thetaD'};
            polar_para=Emitter.polar_para;
            checkFileds(polar_para,fileds); %validate fields
            phiD=polar_para.phiD;
            thetaD=polar_para.thetaD;
            
            % get emitter and imaging parameters
            %--------------------------------------------------------
            
            n1=Nanoscope.refractiveIndx;
            zh=0;% thickness of film
            z2=0;% distance from the emitter to the interface
            n2=1.33;% sample refractive index
            nh=n1;% thin film refractive index
            lambda = Nanoscope.emissWavelength; %wavelength (nm)
            NA =Nanoscope.numApt; %numerical aperture
            N=size(pmask,1);
            pixel_size=Nanoscope.pixelSize;%object pixel size (nm)
            
            if isfield(s,'upsampling')
                upsampling=s.upsampling;
            else
                upsampling=1;%upsampling factor of image space
            end
            molecule_num=numel(deltax);
            %calculate both pupil and image plane sampling,
            %one will affect the other, so make sure not to introduce aliasing
            
            dx = n1*(pixel_size/(upsampling));%  in (nm) due to Abbe sine...
            % condition, scale by imaging medium r.i.
            dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
            
            % define pupil coordinates
            %--------------------------------------------------------
            [eta,xi] = meshgrid( ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(2*N*dx))...
                + (1/(2*dx))), ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(N*2*dx)) + (1/(2*dx))));
            x = lambda*(eta);
            y = lambda*(xi);
            [phi,rho] = cart2pol(x,y);
            rho_max = NA/n1;%pupil region of support determined by NA and imaging medium r.i.
            k1 = n1*(2*pi/lambda);
            kh = nh*(2*pi/lambda);
            k2 = n2*(2*pi/lambda);
            theta1 = asin(rho);%theta in matched medium
            thetah = asin((n1/nh)*sin(theta1));%theta in thin film
            theta2 = asin((n1/n2)*sin(theta1));%theta in mismatched medium
            
            % cache fixed terms
            costheta2=cos(theta2);
            costhetah=cos(thetah);
            costheta1=cos(theta1);
            sintheta1=sin(theta1);
            cosphi=cos(phi);
            sinphi=sin(phi);
            % Fresnel coefficients
            %--------------------------------------------------------
            tp_2h = 2*n2*costheta2./(n2*costhetah + nh*costheta2);
            ts_2h = 2*n2*costheta2./(nh*costhetah + n2*costheta2);
            tp_h1 = 2*nh*costhetah./(nh*costheta1 + n1*costhetah);
            ts_h1 = 2*nh*costhetah./(n1*costheta1 + nh*costhetah);
            
            rp_2h = (n2*costheta2 - nh*costhetah)./(n2*costheta2+ nh*costhetah);
            rs_2h = (nh*costheta2 - n2*costhetah)./(nh*costheta2+ n2*costhetah);
            rp_h1 = (nh*costhetah - n1*costheta1)./(nh*costhetah+ n1*costheta1);
            rs_h1 = (n1*costhetah - nh*costheta1)./(n1*costhetah+ nh*costheta1);
            
            % Axelrod's equations for E-fields at pupil plane
            %--------------------------------------------------------
            mux = reshape(sin(thetaD).*cos(phiD),1,1,molecule_num);
            muy = reshape(sin(thetaD).*sin(phiD),1,1,molecule_num);
            muz = reshape(cos(thetaD),1,1,molecule_num);
            tp = tp_2h.*tp_h1.*exp(1i*kh*costhetah*zh)./(1 + rp_2h.*rp_h1.*exp(2i*kh*zh*costhetah));
            ts = ts_2h.*ts_h1.*exp(1i*kh*costhetah*zh)./(1 + rs_2h.*rs_h1.*exp(2i*kh*zh*costhetah));
            
            Es = bsxfun(@times,ts.*(costheta1./costheta2).*(n1/n2),...
                (bsxfun(@times,muy,cosphi) - bsxfun(@times,mux,sinphi)));
            
            Ep = bsxfun(@times,tp,bsxfun(@times,(n1/n2).*costheta1,(bsxfun(@times,mux,cosphi) + ...
                bsxfun(@times,muy,sinphi))) - ...
                bsxfun(@times,bsxfun(@times,muz,sintheta1),(n1/n2)^2.*(costheta1./costheta2)));
            
            PupilFilt=(rho < rho_max).*1;%casting to numeric
            
            
            %computing E-fields
            %--------------------------------------------------------
            
            E_common=PupilFilt.*(1./sqrt(costheta1)).*exp(1i*kh*zh*costhetah).*...
                exp(1i*k2*z2*costheta2);
            Ey_common=bsxfun(@times,E_common,...
                (bsxfun(@times,cosphi,Es) + bsxfun(@times, sinphi,Ep)));
            Ex_common=bsxfun(@times,E_common,...
                (bsxfun(@times,cosphi,Ep) + bsxfun(@times, -sinphi,Es)));
            
            
            Ey_1=exp(1i*k1*(bsxfun(@times,reshape(z,1,1,molecule_num),costheta1))); %defocus term at pupil plane
            Ex_1=exp(1i*k1*(bsxfun(@times,reshape(z,1,1,molecule_num),costheta1))); %defocus term at pupil plane
            
            % aplly channel mismatch
            
            
            Ex_2=exp(1i*k1*(bsxfun(@times,reshape(-deltax,1,1,molecule_num),sintheta1.*cosphi))); %phase shift
            Ex_3=exp(1i*k1*(bsxfun(@times,reshape(deltay,1,1,molecule_num),sintheta1.*sinphi))); %phase shift
            
            if isfield(s,'channel_mismatch')
                deltax=deltax+s.channel_mismatch(1);
                deltay=deltay+s.channel_mismatch(2);
                
            end
            
            Ey_2=exp(1i*k1*(bsxfun(@times,reshape(deltax,1,1,molecule_num),sintheta1.*cosphi))); %phase shift
            Ey_3=exp(1i*k1*(bsxfun(@times,reshape(deltay,1,1,molecule_num),sintheta1.*sinphi))); %phase shift
            
            
            
            Ey_t=Ey_1.*Ey_2.*Ey_3;
            Ex_t=Ex_1.*Ex_2.*Ex_3;
            Ex=bsxfun(@times,Ex_t,Ex_common);
            Ey=bsxfun(@times,Ey_t,Ey_common);
            
            % for propagation from pupil plane E-field to image plane via tube-lens, paraxial
            % approximation is in force.
            %--------------------------------------------------------
            
            %apply a zero order effect
            %             amp=abs(pmask);
            %             phase_t=angle(pmask);
            %             pmask=.15*amp+.85*amp.*exp(1i*phase_t);
            pmaskRot=rot90(pmask,-1);
            imgy = (fftshift(fft2(ifftshift(Ey.*repmat((pmask),1,1,molecule_num)))));
            imgx= (fftshift(fft2(ifftshift(Ex.*repmat((pmaskRot),1,1,molecule_num)))));
            % image on the camera is the amplitude squared of the electric field
            %--------------------------------------------------------
            imgy = abs(imgy).^2;
            imgx = abs(imgx).^2;
        end
        
        
        %%
        
        function [BXXx,BYYx,BZZx,BXYx,BXZx,BYZx,...
                BXXy,BYYy,BZZy,BXYy,BXZy,BYZy]=computeBases(Nanoscope,Emitter,varargin)
            %computeBases computes
            
            num_molecules=numel(Emitter.position_para.x);
            %number of molecules
            
            simDipole_novotny_h=@(Emitter)Nanoscope.simDipole_novotny(Nanoscope,Emitter,varargin{:});
            % XX
            
            Emitter.polar_para.phiD=zeros(num_molecules,1);
            Emitter.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BXXx,BXXy]=simDipole_novotny_h(Emitter);
            
            % YY
            Emitter.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BYYx,BYYy]=simDipole_novotny_h(Emitter);
            
            % ZZ
            % Emitter.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter.polar_para.phiD=ones(num_molecules,1)*0;
            Emitter.polar_para.thetaD=ones(num_molecules,1)*0;
            
            [BZZx,BZZy]=simDipole_novotny_h(Emitter);
            
            % XY
            Emitter.polar_para.phiD=ones(num_molecules,1)*pi/4;
            Emitter.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BXYxt,BXYyt]=simDipole_novotny_h(Emitter);
            
            % XZ
            Emitter.polar_para.phiD=ones(num_molecules,1)*0;
            Emitter.polar_para.thetaD=ones(num_molecules,1)*pi/4;
            
            [BXZxt,BXZyt]=simDipole_novotny_h(Emitter);
            
            % YZ
            Emitter.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter.polar_para.thetaD=ones(num_molecules,1)*pi/4;
            
            
            [BYZxt,BYZyt]=simDipole_novotny_h(Emitter);
            
            BXYx = 2*BXYxt - BXXx - BYYx;
            BXZx = 2*BXZxt - BXXx - BZZx;
            BYZx = 2*BYZxt - BYYx- BZZx;
            
            BXYy = 2*BXYyt - BXXy - BYYy;
            BXZy = 2*BXZyt - BXXy - BZZy;
            BYZy = 2*BYZyt - BYYy - BZZy;
            
        end
        
        %%
        
        function [FPSFx,FPSFy,lateral_grid_p]=PSF_Fourier_tf(obj)
            
            % compute Fourier transforms
            % ---------------------------------------------------
            up_sample=obj.pixelUpsample;
            FPSFy=struct();
            FPSFx=struct();
            
            %y_channel
            FPSFy.FXXy=single(fft2(ifftshift(up_sample^2*obj.XXyBasis)));
            FPSFy.FYYy=single(fft2(ifftshift(up_sample^2*obj.YYyBasis)));
            FPSFy.FZZy=single(fft2(ifftshift(up_sample^2*obj.ZZyBasis)));
            FPSFy.FXYy=single(fft2(ifftshift(up_sample^2*obj.XYyBasis)));
            FPSFy.FXZy=single(fft2(ifftshift(up_sample^2*obj.XZyBasis)));
            FPSFy.FYZy=single(fft2(ifftshift(up_sample^2*obj.YZyBasis)));
            
            %x_channel
            FPSFx.FXXx=single(fft2(ifftshift(up_sample^2*(obj.XXxBasis))));
            FPSFx.FYYx=single(fft2(ifftshift(up_sample^2*(obj.YYxBasis))));
            FPSFx.FZZx=single(fft2(ifftshift(up_sample^2*(obj.ZZxBasis))));
            FPSFx.FXYx=single(fft2(ifftshift(up_sample^2*(obj.XYxBasis))));
            FPSFx.FXZx=single(fft2(ifftshift(up_sample^2*(obj.XZxBasis))));
            FPSFx.FYZx=single(fft2(ifftshift(up_sample^2*(obj.YZxBasis))));
            
            % compute  lateral grid points
            %----------------------------------------------------
            pixel_size=obj.pixelSize;
            img_size=obj.imageSize;
            lateral_grid_p=(-(img_size-1)/2:1/up_sample:(img_size-1)/2)*pixel_size;
            
        end
        
        
        function [FPSFx,FPSFy]=createPSFstruct(obj,varargin)
            
            %extract input options and set parameters
            %---------------------------------------------------------
            s=opt2struct(varargin);
            
            %determine upsampling
            if isfield(s,'upsampling')
                up_sample=s.upsampling;
            else
                up_sample=obj.pixelUpsample;
            end
            
            %determine image size
            if isfield(s,'imagesize')
                
                imageSize=s.imagesize;
            else
                imageSize=obj.imageSize;
                
            end
            
            %determine channel transmission ratio
            
            if isfield(s,'ytoxchanneltransratio')
                
                yToxChanTransRatio=s.ytoxchanneltransratio;
            else
                yToxChanTransRatio=1;
                
            end
            
            %deltax for computing gradients of PSFs
            deltax=10^-2; %in nm
            
            % output structures for PSFs
            FPSFy=struct();
            FPSFx=struct();
            
            %define a handle
            simDipole_novotny_h=@(Emitter) obj.simDipole_novotny(obj,...
                Emitter,'upsampling',up_sample);
            
            %XX
            %---------------------------------------------------------
            %set molecular orientation
            Emitter.polar_para.phiD=0;
            Emitter.polar_para.thetaD=pi/2;
            
            
            %set position parameters
            Emitter.position_para.x=0;
            Emitter.position_para.y=0;
            Emitter.position_para.z=0;
            
            
            [imgx,imgy] = simDipole_novotny_h(Emitter);
            
            %gradient along x
            Emitter.position_para.x=-0/2;
            [imgxdx1,imgydx1] =  simDipole_novotny_h(Emitter);
            
            Emitter.position_para.x=deltax;
            [imgxdx2,imgydx2] = simDipole_novotny_h(Emitter);
            Emitter.position_para.x=0;
            
            %gradient along y
            %                 Emitter.position_para.y=-deltax/2;
            [imgxdy1,imgydy1] =  simDipole_novotny_h(Emitter);
            
            Emitter.position_para.y=deltax;
            [imgxdy2,imgydy2] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.y=0;
            
            
            
            imgXXx=imgx;imgXXy=imgy;
            imgXXxdx=(imgxdx2-imgxdx1)/deltax;imgXXydx=(imgydx2-imgydx1)/deltax;
            imgXXxdy=(imgxdy2-imgxdy1)/deltax;imgXXydy=(imgydy2-imgydy1)/deltax;
            %YY
            %---------------------------------------------------------
            %set molecular orientation
            Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.thetaD=pi/2;
            
            %set position parameters
            Emitter.position_para.x=0;
            Emitter.position_para.y=0;
            
            [imgx,imgy] =  simDipole_novotny_h(Emitter);
            %                 Emitter.position_para.x=-deltax/2;
            [imgxdx1,imgydx1] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.x=+deltax;
            [imgxdx2,imgydx2] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.x=0;
            
            %                 Emitter.position_para.y=-deltax/2;
            
            [imgxdy1,imgydy1] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.y=deltax;
            [imgxdy2,imgydy2] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.y=0;
            
            imgYYx=imgx;imgYYy=imgy;
            imgYYxdx=(imgxdx2-imgxdx1)/deltax;imgYYydx=(imgydx2-imgydx1)/deltax;
            imgYYxdy=(imgxdy2-imgxdy1)/deltax;imgYYydy=(imgydy2-imgydy1)/deltax;
            
            % ZZ
            %---------------------------------------------------------
            %set molecular orientation
            %             Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.phiD=0;
            Emitter.polar_para.thetaD=0;
            
            
            [imgx,imgy] =  simDipole_novotny_h(Emitter);
            %                 Emitter.position_para.x=-deltax/2;
            [imgxdx1,imgydx1] = simDipole_novotny_h(Emitter);
            Emitter.position_para.x=deltax;
            [imgxdx2,imgydx2] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.x=0;
            
            %                 Emitter.position_para.y=-deltax/2;
            
            [imgxdy1,imgydy1] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.y=deltax;
            [imgxdy2,imgydy2] =  simDipole_novotny_h(Emitter);
            Emitter.position_para.y=0;
            
            imgZZx=imgx;imgZZy=imgy;
            imgZZxdx=(imgxdx2-imgxdx1)/deltax;imgZZydx=(imgydx2-imgydx1)/deltax;
            imgZZxdy=(imgxdy2-imgxdy1)/deltax;imgZZydy=(imgydy2-imgydy1)/deltax;
            
            % XY
            Emitter.polar_para.phiD=pi/4;
            Emitter.polar_para.thetaD=pi/2;
            
            [imgXYxt,imgXYyt]=simDipole_novotny_h(Emitter);
            
            imgXYx = 2*imgXYxt - imgXXx - imgYYx;
            imgXYy = 2*imgXYyt - imgXXy - imgYYy;
            
            % XZ
            Emitter.polar_para.phiD=0;
            Emitter.polar_para.thetaD=pi/4;
            
            
            [imgXZxt,imgXZyt]=simDipole_novotny_h(Emitter);
            
            imgXZx = 2*imgXZxt - imgXXx - imgZZx;
            imgXZy = 2*imgXZyt - imgXXy - imgZZy;
            
            
            % YZ
            Emitter.polar_para.phiD=pi/2;
            Emitter.polar_para.thetaD=pi/4;
            
            [imgYZxt,imgYZyt]=simDipole_novotny_h(Emitter);
            
            imgYZx=2*imgYZxt - imgYYx- imgZZx;
            imgYZy=2*imgYZyt - imgYYy- imgZZy;
            
            %crop and normalize
            %---------------------------------------------------------
            
            %handle for croping region of interest
            
            N_pupil=size(obj.phaseMask,1);
            
            roi=@(img)img(-up_sample*(imageSize-1)/2+N_pupil/2+1:1:up_sample*(imageSize-1)/2+N_pupil/2+1,....
                -up_sample*(imageSize-1)/2+N_pupil/2+1:1:up_sample*(imageSize-1)/2+N_pupil/2+1,:);
            
            sumNormalize=sum(sum(roi(obj.brightnessScaling)));
            %x_channel
            %---------------------------------------------------------
            imgXXx=bsxfun(@times,roi(imgXXx),1./sumNormalize);
            %gradient along x
            imgXXxdx=bsxfun(@times,roi(imgXXxdx),1./sumNormalize);
            %gradient along y
            imgXXxdy=bsxfun(@times,roi(imgXXxdy),1./sumNormalize);
            
            imgYYx=bsxfun(@times,roi(imgYYx),1./sumNormalize);
            %gradient along x
            imgYYxdx=bsxfun(@times,roi(imgYYxdx),1./sumNormalize);
            %gradient along y
            imgYYxdy=bsxfun(@times,roi(imgYYxdy),1./sumNormalize);
            
            imgZZx=bsxfun(@times,roi(imgZZx),1./sumNormalize);
            %gradient along x
            imgZZxdx=bsxfun(@times,roi(imgZZxdx),1./sumNormalize);
            %gradient along y
            imgZZxdy=bsxfun(@times,roi(imgZZxdy),1./sumNormalize);
            
            imgXYx = bsxfun(@times,roi(imgXYx),1./sumNormalize);
            imgXZx = bsxfun(@times,roi(imgXZx),1./sumNormalize);
            imgYZx = bsxfun(@times,roi(imgYZx),1./sumNormalize);
            
            %y_channel
            
            sumNormalize=(1/yToxChanTransRatio)*sumNormalize;
            %---------------------------------------------------------
            imgXXy=bsxfun(@times,roi(imgXXy),1./sumNormalize);
            %gradient along x
            imgXXydx=bsxfun(@times,roi(imgXXydx),1./sumNormalize);
            %gradient along y
            imgXXydy=bsxfun(@times,roi(imgXXydy),1./sumNormalize);
            
            imgYYy=bsxfun(@times,roi(imgYYy),1./sumNormalize);
            %gradient along x
            imgYYydx=bsxfun(@times,roi(imgYYydx),1./sumNormalize);
            %gradient along y
            imgYYydy=bsxfun(@times,roi(imgYYydy),1./sumNormalize);
            
            imgZZy=bsxfun(@times,roi(imgZZy),1./sumNormalize);
            %gradient along x
            imgZZydx=bsxfun(@times,roi(imgZZydx),1./sumNormalize);
            %gradient along y
            imgZZydy=bsxfun(@times,roi(imgZZydy),1./sumNormalize);
            
            imgXYy = bsxfun(@times,roi(imgXYy),1./sumNormalize);
            imgXZy = bsxfun(@times,roi(imgXZy),1./sumNormalize);
            imgYZy = bsxfun(@times,roi(imgYZy),1./sumNormalize);
            
            %compute FFTs (single precision)
            %---------------------------------------------------------
            
            %x_channel
            FPSFx.FXXx=single(fft2((fftshift(up_sample^2*imgXXx))));
            FPSFx.FYYx=single(fft2((fftshift(up_sample^2*imgYYx))));
            FPSFx.FZZx=single(fft2((fftshift(up_sample^2*imgZZx))));
            
            %gradients
            FPSFx.FXXxdx=single(fft2((fftshift(10^2*up_sample^2*imgXXxdx))));
            FPSFx.FXXxdy=single(fft2((fftshift(10^2*up_sample^2*imgXXxdy))));
            FPSFx.FYYxdx=single(fft2((fftshift(10^2*up_sample^2*imgYYxdx))));
            FPSFx.FYYxdy=single(fft2((fftshift(10^2*up_sample^2*imgYYxdy))));
            FPSFx.FZZxdx=single(fft2((fftshift(10^2*up_sample^2*imgZZxdx))));
            FPSFx.FZZxdy=single(fft2((fftshift(10^2*up_sample^2*imgZZxdy))));
            
            FPSFx.FXYx=single(fft2((fftshift(up_sample^2*imgXYx))));
            FPSFx.FXZx=single(fft2((fftshift(up_sample^2*imgXZx))));
            FPSFx.FYZx=single(fft2((fftshift(up_sample^2*imgYZx))));
            
            %y_channel
            FPSFy.FXXy=single(fft2((fftshift(up_sample^2*imgXXy))));
            FPSFy.FYYy=single(fft2((fftshift(up_sample^2*imgYYy))));
            FPSFy.FZZy=single(fft2((fftshift(up_sample^2*imgZZy))));
            
            %gradients
            FPSFy.FXXydx=single(fft2((fftshift(10^2*up_sample^2*imgXXydx))));
            FPSFy.FXXydy=single(fft2((fftshift(10^2*up_sample^2*imgXXydy))));
            FPSFy.FYYydx=single(fft2((fftshift(10^2*up_sample^2*imgYYydx))));
            FPSFy.FYYydy=single(fft2((fftshift(10^2*up_sample^2*imgYYydy))));
            FPSFy.FZZydx=single(fft2((fftshift(10^2*up_sample^2*imgZZydx))));
            FPSFy.FZZydy=single(fft2((fftshift(10^2*up_sample^2*imgZZydy))));
            
            FPSFy.FXYy=single(fft2((fftshift(up_sample^2*imgXYy))));
            FPSFy.FXZy=single(fft2((fftshift(up_sample^2*imgXZy))));
            FPSFy.FYZy=single(fft2((fftshift(up_sample^2*imgYZy))));
            
        end
        
    end
    
    methods
        
        %%
        function [imgx,imgy]=formImage(obj,Emitter,varargin)
            
            
            try
                molecule_num=numel(Emitter.position_para.x);
            catch
                error('Emitter does not have position_para.x field.')
            end
            
            if molecule_num>20
                error('number of molecules exceedes maximum allowable number (20)')
            end
            
            s=opt2struct(varargin);
            %object parameters
            pixel_size=obj.pixelSize;
            up_sample=obj.pixelUpsample;
            N_pupil=size(obj.phaseMask,1);
            
            %map orientation and rotational mobility to second moments
            
            
            rotMobil=Emitter.rotMobility;
            mux=sin(Emitter.theta).*cos(Emitter.phi);
            muy=sin(Emitter.theta).*sin(Emitter.phi);
            muz=cos(Emitter.theta);
            
            %consider only [-pi/2 pi/2] for mux muy
            
            Emitter.secondMoments.muxx=...
                rotMobil.*mux.^2+(1-rotMobil)/3;
            Emitter.secondMoments.muyy=...
                rotMobil.*muy.^2+(1-rotMobil)/3;
            Emitter.secondMoments.muzz=...
                rotMobil.*muz.^2+(1-rotMobil)/3;
            
            Emitter.secondMoments.muxy=...
                rotMobil.*mux.*muy;
            Emitter.secondMoments.muxz=...
                rotMobil.*mux.*muz;
            Emitter.secondMoments.muyz=...
                rotMobil.*muy.*muz;
            
            %molecules' positions
            xpos=Emitter.position_para.x; %in nm
            ypos=Emitter.position_para.y;% in nm
            
            %sub_pixel positions
            %             x_pixel = floor(xpos./(pixel_size));
            %             y_pixel = floor(ypos./(pixel_size));
            %             x_subpixel = xpos - x_pixel.*pixel_size;
            %             y_subpixel = ypos - y_pixel.*pixel_size;
            x_subpixel = xpos ;
            y_subpixel = ypos;
            if isfield(s,'displayinfo') &&s.displayinfo
                %display info
                %----------------------------------------------------
                display(strcat('muxx: ', num2str(Emitter.secondMoments.muxx)))
                display(strcat('muyy: ', num2str(Emitter.secondMoments.muyy)))
                display(strcat('muzz: ', num2str(Emitter.secondMoments.muzz)))
                %----------------------------------------------------
                %display info
                %----------------------------------------------------
                display(strcat('xpos: ', num2str(xpos)))
                display(strcat('ypos: ', num2str(ypos)))
                %----------------------------------------------------
                %display info
                %----------------------------------------------------
                display(strcat('x_pixel: ', num2str(x_pixel)))
                display(strcat('y_pixel: ', num2str(y_pixel)))
                display(strcat('x_subpixel: ', num2str(x_subpixel)))
                display(strcat('y_subpixel: ', num2str(y_subpixel)))
            end
            
            %define molecules with sub_pixel positions
            Emitter_t.position_para.x=x_subpixel;
            Emitter_t.position_para.y=y_subpixel;
            Emitter_t.position_para.z=0*y_subpixel;
            
            %allocate space
            
            img_size=obj.imageSize;
            imgy_shifted = zeros([[img_size,img_size]+N_pupil,molecule_num]);
            imgx_shifted = zeros([[img_size,img_size]+N_pupil,molecule_num]);
            imgy_corped = zeros([[img_size,img_size],molecule_num]);
            imgx_corped = zeros([[img_size,img_size],molecule_num]);
            
            [bx.XX,bx.YY,bx.ZZ, bx.XY,bx.XZ,bx.YZ,...
                by.XX,by.YY,by.ZZ,by.XY,by.XZ,by.YZ]=...
                Nanoscope.computeBases(obj,Emitter_t,varargin{:});
            
            %function handle for computing images formed on the camera
            img=@(bases,moments) bsxfun(@times,bases.XX,reshape(moments.muxx,1,1,molecule_num))...
                + bsxfun(@times,bases.YY,reshape(moments.muyy,1,1,molecule_num))+...
                + bsxfun(@times,bases.ZZ,reshape(moments.muzz,1,1,molecule_num)) +...
                bsxfun(@times,bases.XY, reshape(moments.muxy,1,1,molecule_num))+...
                bsxfun(@times,bases.XZ,reshape(moments.muxz,1,1,molecule_num))+...
                bsxfun(@times,bases.YZ,reshape(moments.muyz,1,1,molecule_num));
            
            imgx=img(bx,Emitter.secondMoments);
            imgy=img(by,Emitter.secondMoments);
            
            
            % shift images to the corresponding positions
            
            %             for nn=1:molecule_num
            %
            %                 start_indx_y=y_pixel(nn)+(img_size-1)/2;
            %                 end_indx_y=y_pixel(nn)+(img_size-1)/2+N_pupil-1;
            %                 start_indx_x=x_pixel(nn)+(img_size-1)/2;
            %                 end_indx_x=x_pixel(nn)+(img_size-1)/2+N_pupil-1;
            %
            %                 imgy_shifted(start_indx_y:end_indx_y,start_indx_x:end_indx_x,nn) = imgy(:,:,nn);
            %
            %                 imgx_shifted(start_indx_y:end_indx_y,start_indx_x:end_indx_x,nn) = imgx(:,:,nn);
            %             end
            
            
            
            % cropping the image to match the region of interest
            %             imgy_corped = imgy_shifted(N_pupil/2+1:img_size+...
            %                 N_pupil/2,N_pupil/2+1:img_size+N_pupil/2,:);
            %             imgx_corped = imgx_shifted(N_pupil/2+1:img_size+...
            %                 N_pupil/2,N_pupil/2+1:img_size+N_pupil/2,:);
            %accounting for photon loss and normalization
            
            Emitter_t.polar_para.phiD=pi/2;
            Emitter_t.polar_para.thetaD=pi/2;
            Emitter_t.position_para.x=0;
            Emitter_t.position_para.y=0;
            Emitter_t.position_para.z=0;
            
            [~,brightness_scaling]=obj.simDipole_novotny(obj,Emitter_t);
            
            %handle for croping region of interest
            roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,....
                -up_sample*(img_size-1)/2+N_pupil/2+2:1:up_sample*(img_size-1)/2+N_pupil/2+2,:);
            
            
            imgy_corped=roi(imgy);
            imgx_corped=roi(imgx);
            sumnorm=sum(sum(roi(brightness_scaling)));
            
            imgx=(imgx_corped)/sumnorm;
            imgy=(imgy_corped)/sumnorm;
            
        end
        
        function [CRLB_vector]=CRB_orinet(obj,Emitter,varargin)           
            
            
            s=opt2struct(varargin);
            
            if isfield(s,'brightness')
                br=s.brightness;
            else
                br=2000;
            end
            
            if isfield(s,'background')
                bg=s.background;
            else
                bg=5;
            end
            
            [imgx,imgy]=formImage(obj,Emitter);
            
            % dx
            Emitter_p=Emitter;
            Emitter_p.position_para.x= Emitter.position_para.x+.001;
            [imgxdx,imgydx]=formImage(obj,Emitter_p);
            
            gradx=br*([imgxdx,imgydx]-[imgx,imgy])/.001;
            % dy
            Emitter_p.position_para.x= Emitter.position_para.x;
            Emitter_p.position_para.y= Emitter.position_para.y+.001;
            [imgxdy,imgydy]=formImage(obj,Emitter_p);
            
            grady=br*([imgxdy,imgydy]-[imgx,imgy])/.001;
            % dtheta
            Emitter_p.position_para.y= Emitter.position_para.y;
            Emitter_p.theta= Emitter.theta+.001;
            [imgxdtheta,imgydtheta]=formImage(obj,Emitter_p);
            
            gradtheta=br*([imgxdtheta,imgydtheta]-[imgx,imgy])/.001;
            % dphi
             Emitter_p.theta= Emitter.theta;
             Emitter_p.phi= Emitter.phi+.001;
            [imgxdphi,imgydphi]=formImage(obj,Emitter_p);
            
             gradphi=br*([imgxdphi,imgydphi]-[imgx,imgy])/.001;

            % drotmob
             Emitter_p.phi= Emitter.phi;
             Emitter_p.rotMobility= Emitter.rotMobility+.001;
            [imgxdrotmob,imgydrotmob]=formImage(obj,Emitter_p);
            
            gradrotMob=br*([imgxdrotmob,imgydrotmob]-[imgx,imgy])/.001;
            
            FI=zeros(5,5);
            FI(1,2)=sum(sum((gradx.*grady)./([imgx,imgy]*br+bg)));
            FI(1,3)=sum(sum((gradx.*gradtheta)./([imgx,imgy]*br+bg)));
            FI(1,4)=sum(sum((gradx.*gradphi)./([imgx,imgy]*br+bg)));
            FI(1,5)=sum(sum((gradx.*gradrotMob)./([imgx,imgy]*br+bg)));
            
            FI(2,3)=sum(sum((grady.*gradtheta)./([imgx,imgy]*br+bg)));
            FI(2,4)=sum(sum((grady.*gradphi)./([imgx,imgy]*br+bg)));
            FI(2,5)=sum(sum((grady.*gradrotMob)./([imgx,imgy]*br+bg)));
            
            FI(3,4)=sum(sum((gradtheta.*gradphi)./([imgx,imgy]*br+bg)));
            FI(3,5)=sum(sum((gradtheta.*gradrotMob)./([imgx,imgy]*br+bg)));
            
            FI(4,5)=sum(sum((gradphi.*gradrotMob)./([imgx,imgy]*br+bg)));
            
            FI=FI'+FI;
            
            FI(1,1)=sum(sum(gradx.^2./([imgx,imgy]*br+bg)));
            FI(2,2)=sum(sum(grady.^2./([imgx,imgy]*br+bg)));
            FI(3,3)=sum(sum(gradtheta.^2./([imgx,imgy]*br+bg)));
            
            FI(4,4)=sum(sum(gradphi.^2./([imgx,imgy]*br+bg)));
            FI(5,5)=sum(sum(gradrotMob.^2./([imgx,imgy]*br+bg)));
            
            if Emitter.theta==0
                FI(4,:)=[];
                FI(:,4)=[];
                FI=FI+diag(ones(1,4)*eps);
                CRLB_vector=sqrt(diag(inv(FI)));
                CRLB_vector(5)=CRLB_vector(4);
                CRLB_vector(4)=intmax;
            else
                FI=FI+diag(ones(1,5)*eps);
                CRLB_vector=sqrt(diag(inv(FI)));
            end

            
        end
        
        %%
        function [imgx,imgy]=formImgIsotropic(obj,Emitter,varargin)
            
            s=opt2struct(varargin);
            
            %imaging parameters
            up_sample=obj.pixelUpsample;
            N_pupil=size(obj.phaseMask,1);
            img_size=obj.imageSize;
            num_molecules=numel(Emitter.position_para.x);
            
            
            Emitter_t.position_para.x=Emitter.position_para.x; %in nm
            Emitter_t.position_para.y=Emitter.position_para.y;% in nm
            Emitter_t.position_para.z=Emitter.position_para.z;% in nm;
            
            simDipole_novotny_h=@(Emitter)Nanoscope.simDipole_novotny(obj,Emitter,varargin{:});
            % XX
            
            Emitter_t.polar_para.phiD=zeros(num_molecules,1);
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BXXx,BXXy]=simDipole_novotny_h(Emitter_t);
            
            % YY
            Emitter_t.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BYYx,BYYy]=simDipole_novotny_h(Emitter_t);
            
            % ZZ
            Emitter_t.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*0;
            
            
            [BZZx,BZZy]=simDipole_novotny_h(Emitter_t);
            
            
            %sum images of XX ,YY, and  ZZ bases
            
            %handle for croping region of interest
            if isfield(s,'imgsize')
                img_size=s.imgsize;
            end
            roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,....
                -up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,:);
            imgx=BXXx+BYYx+BZZx;
            sumNormalize=sum(sum(roi(imgx)));
            imgx=bsxfun(@times,roi(imgx),1./sumNormalize);
            
            imgy=BXXy+BYYy+BZZy;
            sumNormalize=sum(sum(roi(imgy)));
            imgy=bsxfun(@times,roi(imgy),1./sumNormalize);
        end
        %%
        function [imgx,imgy]=formImgIsotropic_costum(obj,Emitter,pmask,varargin)
            
            s=opt2struct(varargin);
            
            %imaging parameters
            up_sample=obj.pixelUpsample;
            N_pupil=size(obj.phaseMask,1);
            img_size=obj.imageSize;
            num_molecules=numel(Emitter.position_para.x);
            
            
            Emitter_t.position_para.x=Emitter.position_para.x; %in nm
            Emitter_t.position_para.y=Emitter.position_para.y;% in nm
            Emitter_t.position_para.z=Emitter.position_para.z;% in nm;
            
            simDipole_novotny_h=@(Emitter)Nanoscope.simDipole_novotny_costum(obj,...
                Emitter,pmask,varargin{:});
            % XX
            
            Emitter_t.polar_para.phiD=zeros(num_molecules,1);
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BXXx,BXXy]=simDipole_novotny_h(Emitter_t);
            
            % YY
            Emitter_t.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*pi/2;
            
            
            [BYYx,BYYy]=simDipole_novotny_h(Emitter_t);
            
            % ZZ
            Emitter_t.polar_para.phiD=ones(num_molecules,1)*pi/2;
            Emitter_t.polar_para.thetaD=ones(num_molecules,1)*0;
            
            
            [BZZx,BZZy]=simDipole_novotny_h(Emitter_t);
            
            
            %sum images of XX ,YY, and  ZZ bases
            
            %handle for croping region of interest
            if isfield(s,'imgsize')
                img_size=s.imgsize;
            end
            roi=@(img)img(-up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,....
                -up_sample*(img_size-1)/2+N_pupil/2+1:1:up_sample*(img_size-1)/2+N_pupil/2+1,:);
            imgx=BXXx+BYYx+BZZx;
            sumNormalize=sum(sum(roi(imgx)));
            imgx=bsxfun(@times,roi(imgx),1./sumNormalize);
            
            imgy=BXXy+BYYy+BZZy;
            sumNormalize=sum(sum(roi(imgy)));
            imgy=bsxfun(@times,roi(imgy),1./sumNormalize);
        end
        %%
        
        function Emitter=SimLipidBilayerExp(obj,varargin)
            
            
            s=opt2struct(varargin);
            
            if isfield(s,'numberofmolecules')
                
                numMol=s.numberofmolecules;
                
            else
                numMol=1;
            end
            
            
            
            
            if isfield(s,'position')
                
                Emitter.position_para.x=s.position.x;
                Emitter.position_para.y=s.position.y;
                Emitter.position_para.z=s.position.z;
                
            else
                
                positionRangeVal=(obj.imageSize-1)/3*obj.pixelSize;
                Emitter.position_para.x=2*rand(1,numMol)*positionRangeVal-positionRangeVal ;
                Emitter.position_para.y=2*rand(1,numMol)*positionRangeVal-positionRangeVal ;
                Emitter.position_para.z=(2*rand(1,numMol)*positionRangeVal-positionRangeVal) *0;
                
            end
            
            if isfield(s,'theta')
                
                Emitter.theta=s.theta;
            else
                
                if isfield(s,'thetarange')
                    
                    try
                        rangeVal=(s.thetarange(2)-s.thetarange(1));
                        
                    catch ME
                        
                        error('expecting an array with 2 elements for range of theta')
                        
                    end
                    
                    if rangeVal<0
                        
                        error('expecting an array with 2 elements ([min_range,max_range]) for range of theta')
                    end
                else
                    
                    s.thetarange=[0 pi/2];
                end
                
                Emitter.theta=s.thetarange(1)+rand(1,numMol)*rangeVal;
            end
            
            if isfield(s,'phi')
                
                Emitter.theta=s.phi;
            else
                
                Emitter.phi=rand(1,numMol)*2*pi;
            end
            
            if isfield(s,'rotationalmobility')
                
                
                Emitter.rotMobility=s.rotationalmobility;
            else
                Emitter.rotMobility=rand(1,numMol);
            end






            
            
            
            %display emitter info
            
            if isfield(s,'displayemitterpara') && s.displayemitterpara
                %orientation
                %----------------------------------------------
                h=histogram2(Emitter.theta*(90/(pi/2)),Emitter.rotMobility,...
                    'DisplayStyle','tile','ShowEmptyBins','on', 'EdgeColor','none');
                h.XBinLimits=[0 90];
                h.YBinLimits=[0 1];
                h.BinWidth=[90*.05,.05];
                xlabel('\theta (degree)')
                ylabel('\gamma')
            end
            
        end
    end
    
    %% plotting methods
    %--------------------------------------------------------
    methods
        
        %%
        function visualizePhaseMask(obj,varargin)
            s_1=opt2struct(varargin);
            if (isfield(s_1,'reusefigure') && s_1.reusefigure)
                
                h=gcf;
            else
                h=figure;
            end
            
            [phaseMask_t,rho_max,sizePhaseMask]=circ(obj);
            
            
            ax1=subplot(1,2,1);
            imagesc(angle(phaseMask_t(sizePhaseMask/2-rho_max:sizePhaseMask/2+rho_max,...
                sizePhaseMask/2-rho_max:sizePhaseMask/2+rho_max)))
            axis image; axis off;
            title('angle ','FontSize',10)
            
            colorbar;
            ax2=subplot(1,2,2);
            imagesc(abs(phaseMask_t(sizePhaseMask/2-rho_max:sizePhaseMask/2+rho_max,...
                sizePhaseMask/2-rho_max:sizePhaseMask/2+rho_max)))
            axis image; axis off;
            title('amplitude','FontSize',10)
            colorbar;
            
            function [phaseMask_t,radius,sizePhaseMask]=circ(obj)
                
                
                phaseMask_t=obj.phaseMask(:,:,1);
                sizePhaseMask=size(phaseMask_t,1);
                radius=obj.phaseMaskPara.pupilRadius;
                
                eta=-sizePhaseMask/2:1:sizePhaseMask/2-1;
                [eta,zeta]=meshgrid(eta);
                
                indx=(eta.^2+zeta.^2)>radius^2;
                
                phaseMask_t(indx>0)=0;
                
            end
            
        end
        
        %%
        function visualizeBases(obj,varargin)
            
            s_1=opt2struct(varargin);
            if (isfield(s_1,'reusefigure') && s_1.reusefigure)
                h=gcf;
            else
                h=figure;
            end
            max_brightness=max(obj.YYyBasis(:));
            
            util_plot(h,obj.XXxBasis,obj.XXyBasis,[.01,.70],.25,.2,'basis','XX');
            util_plot(h,obj.YYxBasis,obj.YYyBasis,[.33,.70],.25,.2,'basis','YY');
            util_plot(h,obj.ZZxBasis,obj.ZZyBasis,[.67,.70],.25,.2,'basis','ZZ');
            util_plot(h,obj.XYxBasis,obj.XYyBasis,[.01,.25],.25,.2,'negativeVal',true,'basis','XY');
            util_plot(h,obj.XZxBasis,obj.XZyBasis,[.33,.25],.25,.2,'negativeVal',true,'basis','XZ');
            util_plot(h,obj.YZxBasis,obj.YZyBasis,[.67,.25],.25,.2,'negativeVal',true,'basis','YZ');
            
            function util_plot(handleIn,imgIn1,imgIn2,pos,w,h,varargin)
                s=opt2struct(varargin);
                ax1=axes(handleIn,'Position',[pos(1),pos(2),w,h]);
                imagesc(imgIn1/max_brightness);
                axis image;axis off
                
                %colorbar
                caxis manual
                cmax=max([max(imgIn1(:)),max(imgIn2(:))])/max_brightness;
                cmin=max([min(imgIn1(:)),min(imgIn2(:))])/max_brightness;
                colormap(ax1,parula);
                
                if (isfield(s,'negativeval') && s.negativeval)
                    caxis([-cmax,cmax])
                else
                    caxis([cmin,cmax])
                end
                c1=colorbar;
                
                set(c1.Label,'Rotation',90);
                c1.Position=[pos(1)+w+.005,pos(2)-h,.01,2*h];
                set(c1, 'YAxisLocation','left')
                
                %remove ticks since both pannels use the same range value
                c1.Ticks=[];
                
                % markup
                if (isfield(s,'basis') && (any(strcmp(s.basis,{'XX','XY'}))))
                    xLim=get(gca,'Xlim');
                    yLim=get(gca,'Ylim');
                    ht = text(0.9*xLim(1)-0.1*xLim(2),0.1*yLim(1)+0.9*yLim(2),...
                        'x-channel',...
                        'Color','k',...
                        'Rotation',90,...
                        'FontWeight','bold');
                end
                %title
                title(s.basis)
                ax2=axes(handleIn,'Position',[pos(1),pos(2)-h,w,h]);
                imagesc(imgIn2/max_brightness);
                axis image;axis off
                caxis manual
                if (isfield(s,'negativeval') && s.negativeval)
                    caxis([-cmax,cmax])
                else
                    caxis([cmin,cmax])
                end
                c2=colorbar;
                
                c2.Position=[pos(1)+w+.015,pos(2)-h,.01,2*h];
                
                % markup
                if (isfield(s,'basis') && (any(strcmp(s.basis,{'XX','XY'}))))
                    xLim=get(gca,'Xlim');
                    yLim=get(gca,'Ylim');
                    ht = text(0.9*xLim(1)-0.1*xLim(2),0.1*yLim(1)+0.9*yLim(2),...
                        'y-channel',...
                        'Color','k',...
                        'Rotation',90,...
                        'FontWeight','bold');
                end
            end
        end
        
        
        
        
        
        %% methods for loading TiFF data
        
        function [img,h,raw_img]=loadImg(obj,filename,varargin)
            
            
            if ~ ischar(filename)
                
                error('Nanoscope:loadImg:InconsistentInputType',...
                    'Expecting a character array for filename.')
            end
            
            s=opt2struct(varargin);
            
            
            if isfield(s,'fullpath') && s.fullpath
                
                FileTif=filename;
                
            else
                
                if 7==exist('dataset','dir')
                    
                    FileTif = fullfile('dataset',filename);
                else
                    error('Nanoscope:loadImg:DirNotFound',...
                        strcat('Expecting ',filename,' in a folder named as dataset'));
                end
                
            end
            
            
            
            % get image info
            
            info_image=imfinfo(FileTif);
            %remove unknown tags
            info_image(1).UnknownTags=[];
            number_images=length(info_image);
            img_sizex=info_image(1).Width;
            img_sizey=info_image(1).Height;
            raw_img=zeros(img_sizey,img_sizex,number_images);
            
            %check input and set appropriate parameter
            
            if isfield(s,'sizeroi') && isnumeric(s.sizeroi)
                
                %make sure image size is an odd number
                if mod(s.sizeroi,2)==0
                    sizeROI=s.sizeroi-1;
                else
                    sizeROI=s.sizeroi;
                end
            else
                sizeROI=obj.imageSize;
            end
            
            if isfield(s,'centerroi') && isnumeric(s.centerroi)
                x_c=s.centerroi(1);
                y_c=s.centerroi(2);
            else
                x_c=round(img_sizex/2);
                y_c=round(img_sizey/2);
            end
            
            %creat a Tiff object and read the data from TIFF file
            
            t=Tiff(FileTif,'r');
            
            for j=1:number_images
                setDirectory(t,j);
                raw_img(:,:,j)=t.read();
            end
            
            %close the Tiff object
            t.close();
            
            % crop images
            %---------------------------------------------------------
            if (isfield(s,'squareroi') && ~(s.squareroi))
                
                try
                    
                    if ((x_c+(sizeROI(1)-1)/2>img_sizex) || (x_c-(sizeROI(1)-1)/2)<=0 ...
                            || (y_c+(sizeROI(2)-1)/2>img_sizey)||(y_c-(sizeROI(2)-1)/2)<=0)
                        error('Nanoscope:loadImg:BadInput',...
                            'size of the ROI  exceeds the  size  of  input image.')
                    end
                catch ME
                    if(strcmp(ME.identifier,'MATLAB:badsubscript'))
                        error('Expecting sizeROI with two elements')
                    else
                        rethrow(ME)
                        
                    end
                end
                widht=sizeROI(1);
                height=sizeROI(2);
                % specify the width and height of the croped image
                %---------------------------------------------------------
                width_ROI=x_c-round((sizeROI(1)-1)/2):x_c+round((sizeROI(1)-1)/2);
                height_ROI=y_c-round((sizeROI(2)-1)/2):y_c+round((sizeROI(2)-1)/2);
                corped_img=raw_img(height_ROI,width_ROI,:);
            else
                
                %check the size of ROI
                if ((x_c+(sizeROI-1)/2>img_sizex) || (x_c-(sizeROI-1)/2)<=0 ...
                        || (y_c+(sizeROI-1)/2>img_sizey)||(y_c-(sizeROI-1)/2)<=0)
                    
                    error('Nanoscope:loadImg:BadInput',...
                        'size of the ROI  exceeds the  size  of  input image.')
                end
                widht=sizeROI(1);
                height=sizeROI(1);
                % specify the width and height of the croped image
                %---------------------------------------------------------
                width_ROI=x_c-(sizeROI-1)/2:x_c+(sizeROI-1)/2;
                height_ROI=y_c-(sizeROI-1)/2:y_c+(sizeROI-1)/2;
                corped_img=raw_img(height_ROI,width_ROI,:);
            end
            
            %subtract offset
            if isfield(s,'offset')
                offset_t=s.offset;
                dimOffset=size(s.offset);
                
                %check offset dimension to match croped image
                
                if ~ all(dimOffset==1) % if not a scalar
                    
                    if (dimOffset(1)~=widht || dimOffset(2)~=height)
                        
                        error('Nanoscope:loadImg:BadInput',...
                            'Expecting an offset image with the same size as the region of interest (sizeROI).')
                    end
                    
                    try
                        
                        % average offset over the stack of images
                        offset_t=sum(s.offset,3)/dimOffset(3);
                        
                    catch ME
                        if~(strcmp(ME.identifier,'MATLAB:badsubscript'))
                            rethrow(ME)
                        end
                    end
                    %
                    %                     if dimOffset(3)~=size(corped_img,3)
                    %
                    %                         % average offset over the stack of images
                    %                         offset_t=sum(s.offset,3)/dimOffset(3);
                    %                     end
                    
                end
                corped_img=bsxfun(@plus,corped_img,-offset_t);
            else
                corped_img=bsxfun(@plus,corped_img,-obj.offset);
            end
            
            %make sure all pixels take on positive values
            min_pixel_val=.001;
            non_pos_indx=corped_img<=0;
            corped_img(non_pos_indx)=min_pixel_val;
            
            %convert to photon counts
            if isfield(s,'adcount')
                pixel_val_to_photon=s.adcount;
            else
                pixel_val_to_photon=obj.ADcount;
            end
            
            img=single(corped_img.*pixel_val_to_photon);
            
            h=[];% empty handle
            if isfield(s,'visualize') &&s.visualize
                
                h=figure;
                imagesc(raw_img(:,:,ceil(number_images/2)))
                axis image
                hold on
                
                if isfield(s,'roi') &&s.roi
                    
                    % show ROI with a rectangle
                    rectangle('Position',[x_c-(widht-1)/2,y_c-(height-1)/2,widht,height],'EdgeColor','y')
                    
                    % display center of ROI with a marker
                    plot(x_c,y_c,'*y','MarkerSize',3)
                end
                
                % zoom in
                imagesc(img(:,:,ceil(number_images/2)))
                axis image
                drawnow
            end
        end
        
        %%
        function [SMLM_img,folder_path,posRect,centerROI]=loadImgInteractive(obj,varargin)
            %loadImgInteractive load raw images form a specificed path; a
            %ROI is selected by user;
            %it subtracts off-set and accounts for A/D count; it also applies a geometric
            %transform to find the corresponding ROI in the other channel
            %(optional)
            
            
            s=opt2struct(varargin);
            
            % specify the data folder
            if isfield(s,'datapath') && ischar(s.datapath)
                folder_path=s.datapath;
            else
                folder_path = uigetdir('','Folder of the desired set of camera images'); % specify the path of the desired set of camera images
            end
            if isfield(s,'offsetpath') && ischar(s.offsetpath)
                offset_folder_path=s.offsetpath;
            else
                offset_folder_path=uigetdir('','Folder of  the off-set image');% specify the path of the off-set images
            end
            
            
            addpath(folder_path);
            addpath(offset_folder_path);
            FilesInfo =dir(fullfile(folder_path,'*.tif'));% extracting  stack of images
            
            offset_files_info=dir(fullfile(offset_folder_path,'*.tif'));
            
            %catch possible errors
            if isempty(FilesInfo)
                error('Nanoscope:BadFileName',...
                    'The specified folder contains no file with .tif extension')
            end
            if isempty(offset_files_info)
                error('Nanoscope:BadFileName',...
                    'The specified folder for offset contains no file with .tif extension')
            end
            
            % load SMLM images
            %---------------------------------------------------------
            %get dimension info
            FileTif = fullfile(folder_path,FilesInfo(1).name);
            info_image=imfinfo(FileTif);
            number_axial_scans=length(FilesInfo);
            chip_img_width=info_image(1).Width; % width of the image captured on camera
            %launch an interactive ROI selection via a rectangle
            
            filename=fullfile(folder_path,FilesInfo(ceil(number_axial_scans/2)).name);
            
            % the image at focus is displayed
            [~,h_t,img_t]=obj.loadImg(filename,'fullpath',true,...
                'visualize',true);
            title('Select a region of interet','FontSize',11)
            
            if ~isfield(s,'nointeractive') || ~(s.nointeractive)
                
                %setup  an inteactive rectangle
                hRect=imrect(gca,[info_image(1).Width/3 ,info_image(1).Height/3,...
                    info_image(1).Height/3,info_image(1).Height/3]);
                
                %add an event listener to display chosen ROI at the top
                %left corner
                addID= addNewPositionCallback(hRect,...
                    @(positionRect) zoominfcn(positionRect,img_t));
                
                %get the position of the rectangle object
                posRect=round(wait(hRect));
                
                %get the rid of the listener
                removeNewPositionCallback(hRect,addID);
                
                %make sure the chosen ROI is valid
                if (posRect(3)<81 ||posRect(4)<81)
                    
                    error('Nanoscope:InconsistentInputValue',...
                        'Soryy! the ROI selected is small for 3D PSFs.')
                end
                
                
                %get ROI info
                if (isfield(s,'squareroi') && ~(s.squareroi))
                    %make sure ROI is large enough for sliding window algorithm
                    posRect(3:4)=posRect(3:4)+round(obj.imageSize/2);
                    sizeROI=[posRect(3),posRect(4)];
                    if mod(sizeROI(1),2)==0
                        sizeROI(1)=sizeROI(1)+1;
                    end
                    if mod(sizeROI(2),2)==0
                        sizeROI(2)=sizeROI(2)+1;
                    end
                    width=sizeROI(1);
                    height=sizeROI(2);
                else
                    s.squareroi=true;
                    sizeROI=min(posRect(3),posRect(4));
                    if mod(sizeROI,2)==0
                        sizeROI=sizeROI+1;
                    end
                    width=sizeROI(1);
                    height=sizeROI(1);
                end
                
                centerROI=[floor(posRect(1)+width/2),floor(posRect(2)+height/2)];
                
            else
                s.squareroi=true;
                posRect=[];
                try
                    sizeROI=s.sizeroi;
                catch ME
                    error('Expecting sizeROI option as input!')
                end
                
                if mod(sizeROI,2)==0
                    sizeROI=sizeROI+1;
                end
                width=sizeROI(1);
                height=sizeROI(1);
                
                try
                    centerROI=s.centerroi;
                catch ME
                    error('Expecting centerROI option as input!')
                end
            end
            
            
            
            close(h_t)
            
            number_frames=size(img_t,3);
            SMLM_raw_img=zeros(height,width,...
                number_frames,number_axial_scans);
            
            %set the flag to display image acquired at focus
            visulaizeFlag=zeros(1,number_axial_scans);
            visulaizeFlag(ceil((number_axial_scans+1)/2))=1;
            
            % read images
            for i=1:number_axial_scans
                filename = fullfile(folder_path,FilesInfo(i).name);
                
                SMLM_raw_img(:,:,:,i)=obj.loadImg(filename,'fullpath',true,...
                    'offset',0,...
                    'sizeROI',sizeROI,...
                    'centerROI',centerROI,...
                    'ADcount',1,...
                    'visualize',visulaizeFlag(i),...
                    'ROI',true,...
                    'squareroi',s.squareroi);
            end
            
            % load  offset image
            %---------------------------------------------------------
            offset_file_name = fullfile(offset_folder_path,offset_files_info(1).name);
            
            if  length(offset_files_info)>1
                
                error('Nanoscope:InconsistentInputValue',...
                    'offset (averaged) must contain a single stack of tif file.')
            end
            
            offfset_raw_img=obj.loadImg(offset_file_name,'fullpath',true,...
                'offset',0,...
                'sizeROI',sizeROI,...
                'centerROI',centerROI,...
                'ADcount',1,...
                'visualize',false,...
                'squareroi',s.squareroi);
            
            %average offset stack
            if isfield(s,'nooffset') && s.nooffset
                offfset_avg_img=0;
            else
                offfset_avg_img=sum(double(offfset_raw_img),3)/size(offfset_raw_img,3);
            end
            %subtract offset
            SMLM_img=bsxfun(@plus,double(SMLM_raw_img),-offfset_avg_img);
            
            %measurments shoud be positive
            SMLM_img((SMLM_img<=0)>0)=.001;
            
            %apply photoelecton conversion
            if ~isfield(s,'noadconversion') || ~s.noadconversion
                SMLM_img=SMLM_img.*obj.ADcount;
            end
            
            % select the appropriate region in other channel
            %--------------------------------------------------
            if isfield(s,'tform')
                if ~ isobject(s.tform)
                    error('Nanoscope:InconsistentInputType',...
                        'expecting an  images.geotrans.PolynomialTransformation2D object for input transform ');
                end
                %step 1- identify the current channel
                if centerROI(1)<chip_img_width/2
                    cur_channel='L';
                else
                    cur_channel='R';
                end
                
                %step 2- map the centerROI to an appropriate
                %coordinate
                ref_x_coordinate=chip_img_width/2; % the 0 coordinate for x-axis
                %used in obtaining the transformation
                if strcmp(cur_channel,'L')
                    
                    NewCenterROI(1)=ref_x_coordinate-centerROI(1);
                    NewCenterROI(2)=centerROI(2);
                else
                    NewCenterROI(1)=-ref_x_coordinate+centerROI(1);
                    NewCenterROI(2)=centerROI(2);
                    
                end
                
                %step 3- apply the transform on the center coordinate
                
                if isfield(s,'transformscale')
                    
                    scale= s.transformscale;
                    
                else
                    scale=1;
                end
                [transformed_centerROI(1),transformed_centerROI(2)]=...
                    transformPointsInverse(s.tform,NewCenterROI(1)*scale,NewCenterROI(2)*scale);
                
                %step 4- map the transformed coordinate to original
                %coordinate(pixel)
                transformed_centerROI=ceil(transformed_centerROI/scale);
                
                if isfield(s,'registershiftx')
                    deltax=s.registershiftx;
                else
                    deltax=0;
                    
                end
                
                if isfield(s,'registershifty')
                    
                    deltay=s.registershifty;
                else
                    deltay=0;
                end
                
                
                transformed_centerROI(1)=transformed_centerROI(1)-deltax;
                transformed_centerROI(2)=transformed_centerROI(2)-deltay;
                
                transformed_centerROI=ceil(transformed_centerROI);
                if strcmp(cur_channel,'L')
                    transformed_centerROI(1)=(transformed_centerROI(1))+ref_x_coordinate;
                else
                    transformed_centerROI(1)=-(transformed_centerROI(1))+ref_x_coordinate;
                end
                
                %step 5- extract the region on the other channel
                
                % read images
                for i=1:number_axial_scans
                    filename = fullfile(folder_path,FilesInfo(i).name);
                    
                    SMLM_raw_img_2(:,:,:,i)=obj.loadImg(filename,'fullpath',true,...
                        'offset',0,...
                        'sizeROI',sizeROI,...
                        'centerROI',transformed_centerROI,...
                        'ADcount',1,...
                        'visualize',visulaizeFlag(i),...
                        'ROI',true,...
                        'squareroi',s.squareroi);
                end
                
                offfset_raw_img_2=obj.loadImg(offset_file_name,'fullpath',true,...
                    'offset',0,...
                    'sizeROI',sizeROI,...
                    'centerROI',transformed_centerROI,...
                    'ADcount',1,...
                    'visualize',false,...
                    'squareroi',s.squareroi);
                %average offset stack
                if isfield(s,'nooffset') && s.nooffset
                    offfset_avg_img=0;
                else
                    offfset_avg_img=sum(double(offfset_raw_img_2),3)/size(offfset_raw_img_2,3);
                end
                %subtract offset
                SMLM_img_2=bsxfun(@plus,double(SMLM_raw_img_2),-offfset_avg_img);
                
                %measurments shoud be positive
                SMLM_img_2((SMLM_img_2<=0)>0)=.001;
                
                %apply photoelecton conversion
                if ~isfield(s,'noadconversion') || ~s.noadconversion
                    SMLM_img_2=SMLM_img_2.*obj.ADcount;
                end
                
                if strcmp(cur_channel,'L')
                    %              SMLM_img=flipud(SMLM_img);
                    %              SMLM_img_2=fliplr(SMLM_img_2);
                    %            SMLM_img_2=flipud(SMLM_img_2);
                    
                    SMLM_img=[SMLM_img_2,fliplr(SMLM_img)];
                else
                    SMLM_img_2=flipud(SMLM_img_2);
                    %              SMLM_img=fliplr(SMLM_img);
                    SMLM_img=flipud(SMLM_img);
                    
                    SMLM_img=[SMLM_img,SMLM_img_2];
                end
                
                
            end
            %local functions
            %--------------------------------------------------
            function zoominfcn(p,img)
                imagesc(img(round(p(2)):round(p(2))+round(p(4)),(round(p(1)):round(p(1))+round(p(3)))));
                axis image
            end
            
        end
            
    end
end
    