function projected_gamma = projection_cone2(gamma, r, recovStruct, varargin)
%projection_cone projects a vector onto a second-order cone
%constraint as well as a support constraint captured by the indx variable.
%->---
%input
%->---
%gamma:             array(3N,n_f) - molecular parameter estimates
%N:             scalar        - number of grid points
%r:             scalar        - constraint parameter associated with the secon-order cone
%n_f:           scalar        - number of frames
%optinal input parameters:
%indx:          logical array(3N,n_f) - indices of grid points with no
%molecules associated with them
%---->-
%output
%---->-
%projected_gamma:  array(3N,n_f) - projected molecular parameter estimates

numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:projectiong2:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

N = recovStruct.n_grid_p;

% set defaults for optional inputs
optargs = {[], [], []};
optargs(1:numvarargs) = varargin;

indx_xx = optargs{1};
indx_yy = optargs{2};
indx_zz = optargs{3};

%rearrange molecular parameters
%------------------------------
xx = gamma(1:N, :);
yy = gamma(N+1:2*N, :);
zz = gamma(2*N+1:3*N, :);
xxdx = gamma(6*N+1:7*N, :);
xxdy = gamma(7*N+1:8*N, :);
yydx = gamma(8*N+1:9*N, :);
yydy = gamma(9*N+1:10*N, :);
zzdx = gamma(10*N+1:11*N, :);
zzdy = gamma(11*N+1:12*N, :);
projected_gammaxy_t = gamma(3*N+1:4*N, :);
projected_gammaxz_t = gamma(4*N+1:5*N, :);
projected_gammayz_t = gamma(5*N+1:6*N, :);

% x1=gamma(1:N,:);
% x2=gamma(N+1:2*N,:);
% x3=gamma(2*N+1:3*N,:);
% x23=sqrt(x2.^2+x3.^2);

xxjoint = sqrt(xxdx.^2+xxdy.^2);
yyjoint = sqrt(yydx.^2+yydy.^2);
zzjoint = sqrt(zzdx.^2+zzdy.^2);

%compute the conditions of the projection operator
%-------------------------------------------------

%XX
cxx = (xx + r * xxjoint) / (1 + r^2);
cxx1 = xxjoint <= xx * r;
cxx2 = xxjoint <= -xx / r;
cxx3 = xxjoint > abs(xx*r);
cxx_t1 = (1 - cxx2) .* cxx1;
cxx_t2 = (1 - cxx2) .* cxx3 .* cxx;

projected_gammaxx_t = cxx_t1 .* xx + cxx_t2;
projected_gammaxx_tx = cxx_t1 .* xxdx + cxx_t2 .* r .* (xxdx) ./ (eps + xxjoint);
projected_gammaxx_ty = cxx_t1 .* xxdy + cxx_t2 .* r .* xxdy ./ (eps + xxjoint);
projected_gammaxx_t = projected_gammaxx_t .* (projected_gammaxx_t > 0);

%YY
cyy = (yy + r * yyjoint) / (1 + r^2);
cyy1 = yyjoint <= yy * r;
cyy2 = yyjoint <= -yy / r;
cyy3 = yyjoint > abs(yy*r);
cyy_t1 = (1 - cyy2) .* cyy1;
cyy_t2 = (1 - cyy2) .* cyy3 .* cyy;

projected_gammayy_t = cyy_t1 .* yy + cyy_t2;
projected_gammayy_tx = cyy_t1 .* yydx + cyy_t2 .* r .* (yydx) ./ (eps + yyjoint);
projected_gammayy_ty = cyy_t1 .* yydy + cyy_t2 .* r .* yydy ./ (eps + yyjoint);
projected_gammayy_t = projected_gammayy_t .* (projected_gammayy_t > 0);


%ZZ
czz = (zz + r * zzjoint) / (1 + r^2);
czz1 = zzjoint <= zz * r;
czz2 = zzjoint <= -zz / r;
czz3 = zzjoint > abs(zz*r);
czz_t1 = (1 - czz2) .* czz1;
czz_t2 = (1 - czz2) .* czz3 .* czz;

projected_gammazz_t = czz_t1 .* zz + czz_t2;
projected_gammazz_tx = czz_t1 .* zzdx + czz_t2 .* r .* (zzdx) ./ (eps + zzjoint);
projected_gammazz_ty = czz_t1 .* zzdy + czz_t2 .* r .* zzdy ./ (eps + zzjoint);
projected_gammazz_t = projected_gammazz_t .* (projected_gammazz_t > 0);


% c=(x1+r*x23)/(1+r^2);
% c1=x23<= x1*r;
% c2=x23<= -x1/r;
% c3=x23> abs(x1*r);
% c_t1=(1-c2).*c1;
% c_t2=(1-c2).*c3.*c;

% projected_gamma_t=c_t1.*x1+c_t2;
% projected_gamma_tx=c_t1.*x2+c_t2.*r.*(x2)./(eps+x23);
% projected_gamma_ty=c_t1.*x3+c_t2.*r.*x3./(eps+x23);
% projected_gamma_t=projected_gamma_t.*(projected_gamma_t>0);

%apply support constraints. Set indices in indx variable to zero
%---------------------------------------------------------------
% projected_gamma_t(indx)=0;
% projected_gamma_tx(indx)=0;
% projected_gamma_ty(indx)=0;

projected_gammaxx_t(indx_xx) = 0;
projected_gammaxx_tx(indx_xx) = 0;
projected_gammaxx_ty(indx_xx) = 0;

projected_gammayy_t(indx_yy) = 0;
projected_gammayy_tx(indx_yy) = 0;
projected_gammayy_ty(indx_yy) = 0;

projected_gammazz_t(indx_zz) = 0;
projected_gammazz_tx(indx_zz) = 0;
projected_gammazz_ty(indx_zz) = 0;


%set indices whos values are zero in both XX and YY bases to zero
indx_xy = (indx_xx .* indx_yy) > 0;

projected_gammaxy_t(indx_xy) = 0;

%set indices whos values are zero in both XX and ZZ bases to zero
indx_xz = (indx_xx .* indx_zz) > 0;
projected_gammaxz_t(indx_xz) = 0;

%set indices whos values are zero in both YY and ZZ bases to zero
indx_yz = (indx_yy .* indx_zz) > 0;
projected_gammayz_t(indx_yz) = 0;
%rearrange the molecular estimates
%---------------------------------
% projected_gamma=vertcat(projected_gamma_t,....
%     projected_gamma_tx,projected_gamma_ty);


projected_gamma(1:N, :) = projected_gammaxx_t;
projected_gamma(N+1:2*N, :) = projected_gammayy_t;
projected_gamma(2*N+1:3*N, :) = projected_gammazz_t;
projected_gamma(3*N+1:4*N, :) = projected_gammaxy_t;
projected_gamma(4*N+1:5*N, :) = projected_gammaxz_t;
projected_gamma(5*N+1:6*N, :) = projected_gammayz_t;
projected_gamma(6*N+1:7*N, :) = projected_gammaxx_tx;
projected_gamma(7*N+1:8*N, :) = projected_gammaxx_ty;
projected_gamma(8*N+1:9*N, :) = projected_gammayy_tx;
projected_gamma(9*N+1:10*N, :) = projected_gammayy_ty;
projected_gamma(10*N+1:11*N, :) = projected_gammazz_tx;
projected_gamma(11*N+1:12*N, :) = projected_gammazz_ty;
%projected_gamma(3*N+1:6*N,:)=gamma(3*N+1:6*N,:);
