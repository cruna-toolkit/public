%% Grid/Time-step estimator for 3D time-domain Euler simulations
% Estimates grid points and timestep from room size, c, CFL, PPW and f_max.
%
% Inputs:
%   L          : room dimensions [Lx Ly Lz] in meters
%   c          : speed of sound [m/s]
%   CFL        : CFL number (scheme dependent)
%   ppw        : points per wavelength
%   f_max      : maximum frequency to resolve [Hz]
%
% Optional:
%   u_max      : estimate of maximum flow speed magnitude [m/s] (default 0)
%   mult       : refinement factor (>=1), >1 -> finer grid (default 1)
%   dt_mode    : 'simple' or 'sum3d' timestep estimate (default 'sum3d')
%                simple: dt = CFL * min(dx,dy,dz) / (c + u_max)
%                sum3d : dt = CFL / ( (a/dx) + (a/dy) + (a/dz) ), a=c+u_max
%
% Outputs (printed):
%   Nx,Ny,Nz, total points
%   dx,dy,dz effective spacings
%   dt estimate, Nt for Tsim
%   Time-scale interpretation: crossings vs Lx and room diagonal

clear; clc;

%% ------------------- USER INPUTS -------------------
L     = [9 7 5];   % [m] room size: Lx Ly Lz
c     = 343;             % [m/s] speed of sound
CFL   = 0.75;            % [-] CFL number
ppw   = 12;              % [-] points per wavelength
f_max = 2000;            % [Hz] max frequency you want to resolve

u_max  = 0;              % [m/s] max flow speed estimate (0 for acoustics)
mult   = 1.0;            % refinement factor (>1 -> finer grid)
dt_mode = 'sum3d';       % 'simple' or 'sum3d'

Tsim  = 0.2;            % [s] simulation duration you want to cover
%% ---------------------------------------------------

%% Derived quantities (spatial)
lambda_min = c / f_max;                 % [m]
d_base     = lambda_min / ppw;          % [m]
d_target   = d_base / mult;             % [m]

% Use uniform spacing by default
dx = d_target; dy = d_target; dz = d_target;

% Number of grid points (include boundaries -> +1)
Nx = ceil(L(1)/dx) + 1;
Ny = ceil(L(2)/dy) + 1;
Nz = ceil(L(3)/dz) + 1;

% Effective spacings so domain fits exactly between points
dx_eff = L(1)/(Nx-1);
dy_eff = L(2)/(Ny-1);
dz_eff = L(3)/(Nz-1);

%% Timestep estimate
a_max = c + u_max; % max characteristic speed estimate

switch lower(dt_mode)
    case 'simple'
        dt = CFL * min([dx_eff, dy_eff, dz_eff]) / a_max;

    case 'sum3d'
        % More conservative / typical multi-D CFL estimate
        dt = CFL / ( (a_max/dx_eff) + (a_max/dy_eff) + (a_max/dz_eff) );

    otherwise
        error('dt_mode must be ''simple'' or ''sum3d''.');
end

% Time vector: ensure we cover at least Tsim
Nt = ceil(Tsim/dt) + 1;
t  = (0:Nt-1) * dt;

%% Time-scale interpretation (acoustic travel / "bulk convection" length)
Lx    = L(1);
Ldiag = norm(L); % sqrt(Lx^2 + Ly^2 + Lz^2)

% Acoustic travel times (speed = c)
t_ac_x    = Lx / c;
t_ac_diag = Ldiag / c;

% Crossings during Tsim
Ncross_ac_x    = Tsim / t_ac_x;       % = Tsim * c / Lx
Ncross_ac_diag = Tsim / t_ac_diag;    % = Tsim * c / Ldiag

% Using max signal speed a_max = c + u_max
t_a_x    = Lx / a_max;
t_a_diag = Ldiag / a_max;

Ncross_a_x    = Tsim / t_a_x;         % = Tsim * a_max / Lx
Ncross_a_diag = Tsim / t_a_diag;      % = Tsim * a_max / Ldiag

% "Bulk convection lengths" during Tsim
ell_ac = c * Tsim;        % [m] acoustic distance traveled in Tsim
ell_a  = a_max * Tsim;    % [m] max-signal distance traveled in Tsim

%% Print results
fprintf('--- Grid / time-step estimate ---\n');
fprintf('Room size [m]:      Lx=%.4g, Ly=%.4g, Lz=%.4g\n', L(1), L(2), L(3));
fprintf('c [m/s]=%.4g, f_max [Hz]=%.4g, PPW=%.4g -> lambda_min=%.4g m\n', c, f_max, ppw, lambda_min);
fprintf('Target spacing [m]: dx=dy=dz=%.6g (mult=%.3g)\n', d_target, mult);
fprintf('Grid points:        Nx=%d, Ny=%d, Nz=%d  (total=%.6g)\n', Nx, Ny, Nz, Nx*Ny*Nz);
fprintf('Effective spacing:  dx=%.6g, dy=%.6g, dz=%.6g\n', dx_eff, dy_eff, dz_eff);

fprintf('\n--- Timestep estimate ---\n');
fprintf('dt_mode:            %s\n', dt_mode);
fprintf('a_max=c+u_max [m/s]=%.6g (u_max=%.6g)\n', a_max, u_max);
fprintf('CFL:                %.6g\n', CFL);
fprintf('dt estimate [s]:    %.6g\n', dt);
fprintf('Tsim=%.6g s -> Nt=%d steps, t(end)=%.6g s\n', Tsim, Nt, t(end));

fprintf('\n--- Time-scale interpretation ---\n');
fprintf('Ldiag [m]:          %.6g\n', Ldiag);
fprintf('Acoustic travel time (x):        t_x    = %.6g s  (Lx/c)\n', t_ac_x);
fprintf('Acoustic travel time (diagonal): t_diag = %.6g s  (Ldiag/c)\n', t_ac_diag);
fprintf('Tsim crossings (acoustic) x:     Tsim/t_x    = %.6g\n', Ncross_ac_x);
fprintf('Tsim crossings (acoustic) diag:  Tsim/t_diag = %.6g\n', Ncross_ac_diag);

fprintf('\nUsing max signal speed a_max=c+u_max:\n');
fprintf('Travel time (x):        t_x(a_max)    = %.6g s\n', t_a_x);
fprintf('Travel time (diagonal): t_diag(a_max) = %.6g s\n', t_a_diag);
fprintf('Tsim crossings (a_max) x:        %.6g\n', Ncross_a_x);
fprintf('Tsim crossings (a_max) diag:     %.6g\n', Ncross_a_diag);

fprintf('\n--- Bulk distance traveled during Tsim ---\n');
fprintf('ell_ac = c*Tsim [m]:             %.6g\n', ell_ac);
fprintf('ell_a  = (c+u_max)*Tsim [m]:     %.6g\n', ell_a);
fprintf('ell_ac / Lx:                     %.6g\n', ell_ac / Lx);
fprintf('ell_ac / Ldiag:                  %.6g\n', ell_ac / Ldiag);
fprintf('ell_a  / Lx:                     %.6g\n', ell_a / Lx);
fprintf('ell_a  / Ldiag:                  %.6g\n', ell_a / Ldiag);