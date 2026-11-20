clear; clc;

% Critical response radius
R_um = 30;   % micrometers

% Minimum center-to-center spacing
d_um = 2 * R_um;


% Area occupied by one site in a hexagonal lattice
A_um2 = (sqrt(3) / 2) * d_um^2;

% Information density
density_bits_um2 = 1 / A_um2;

% Convert bits/um^2 to MB/cm^2
density_MB_cm2 = density_bits_um2 * 1e8 / (8e6);

fprintf('Radius = %.2f um\n', R_um);
fprintf('Spacing = %.2f um\n', d_um);
fprintf('Density = %.4e bits/um^2\n', density_bits_um2);
fprintf('Density = %.4e MB/cm^2\n', density_MB_cm2);