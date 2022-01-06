% Proracun motorskog nosaca MKE
% Autor: Milos D. Petrasinovic <mpetrasinovic@mas.bg.ac.rs>
% Proracun strukture letelica
% Masinski fakultet, Univerzitet u Beogradu
% Katedra za vazduhoplovstvo, Struktura letelica
% https://vazmfb.com
% Beograd, 2021
%
% ---------------
%
% Copyright (C) 2021 Milos Petrasinovic <info@vazmfb.com>
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%   
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%   
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ---------------
clear('all'), clc, close('all'), tic
disp([' --- ' mfilename ' --- ']);

% - Ulazni podaci

% Karaktetristike poprecnog preseka
Ds = 16; % [mm] spoljni precnik cevi
t = 1; % [mm] debljina zida cevi
% Karakteristike materijala
E = 210e9; % [Pa] Modul elasticnosti

% Koordinate cvorova
% nc(cvor, 1:3) = [x, y, z];
nc = [0, 0, 0; % cvor 1
      0, 700, 0; % cvor 2
      0, 0, 700; % cvor 3
      0, 700, 700; % cvor 4
      1200, 0, 350; % cvor 5
      1200, 700, 350; % cvor 6
      300, 37.5, 350; % cvor 7
      300, 662.5, 350; % cvor 8
      900, 350, 37.5; % cvor 9
      900, 350, 662.5;]; % cvor 10

% Cvorovi elemenata motorskog nosaca
% enn(element, :) = [cvor1, cvor2];
enn = [1, 7; % stap 1
      3, 7; % stap 2
      2, 8; % stap 3
      4, 8; % stap 4
      5, 10; % stap 5
      5, 9; % stap 6
      6, 10; % stap 7
      6, 9;];  % stap 8

% Cvorovi elemenata motora (beskonacno kruti elementi)
% enm(element, :) = [cvor1, cvor2];
enm = [7, 8; % stap motora 1
       8, 9; % stap motora 2
       9, 7; % stap motora 3
       7, 10; % stap motora 4
       10, 8; % stap motora 5
       10, 9;]; % stap motora 6

% Konturni uslovi
% bc(cvor, 1:3) = [Tx, Ty, Tz];
% 0 - slobodan
% 1 - ogranicen
bc(1, :) = [1, 1, 1]; % cvor 1
bc(2, :) = [1, 1, 1]; % cvor 2
bc(3, :) = [1, 1, 1]; % cvor 3
bc(4, :) = [1, 1, 1]; % cvor 4
bc(5, :) = [1, 1, 1]; % cvor 5
bc(6, :) = [1, 1, 1]; % cvor 6

% Spoljasnja opterecenja
% F(cvor, 1:3) = [Fx[N], Fy[N], Fz[N]];
F(7, :) = [3906.95, 0, -5036.15]; % cvor 7
F(8, :) = [3906.95, 0, -5036.15]; % cvor 8
F(9, :) = [3906.95, 0, -5036.15]; % cvor 9
F(10, :) = [3906.95, 0, -5036.15]; % cvor 10

% Dodatne promenljive
s = [1000, % [m] u [mm]
     200, % uvecanje prikaza pomeranja
     100, % duzina vektora za prikaz konturnih uslova
     1/50, % uvecanje vektora za prikaz spoljasnjih opterecenja
     100; % slobodan prostor oko modela
     0.2, % velicina strelice vektora
     -60]; % pomeranje oznaka cvorova

% - Definisanje matrice krutosti
d2r = 180/pi; % stepeni u radijane
en = [enn; enm]; % cvorovi svih elemenata
ne = size(en, 1); % broj elemenata
nen = size(enn, 1);
nem = size(enm, 1);
nn = size(nc, 1); % broj cvorova
dof = 3*nn; % broj stepeni slobode
D = zeros(dof, 1); % pomeranja
bc(end+1:nn, :) = zeros(length(size(bc, 1)+1:nn), 3); % konturni uslovi
F(end+1:nn, :) = zeros(length(size(F, 1)+1:nn), 3); % spoljasnja opterecenja
F = reshape(F.', [], 1); 
K = zeros(dof, dof); % matrica krutosti
sigma = zeros(nen, 1);

% Stepeni slobode cvorova
rdof = find(bc')'; % ograniceni stepeni slobode
fdof = find(~bc')'; % slobodni stepeni slobode

% Karakteristike poprecnog preseka (kruzna cev)
Ds = Ds/s(1); % [m]
t = t/s(1); % [m]
A = (Ds^2-(Ds-2*t)^2)*pi/4*ones(ne, 1); % [m^2] Povrsina poprecnog preseka
A(nen+1:end) = A(nen+1:end)*10^6; % povecavanje krutosti

% Odredjivanje matrice krutosti
nc = nc/s(1); % [m]
for i=1:ne 
  j = en(i, :);       
  edof = [3*j(1)-2, 3*j(1)-1, 3*j(1) ... % stepeni slobode elementa
         3*j(2)-2, 3*j(2)-1, 3*j(2)]; 
  L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % duzina elementa
      nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
      nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
  CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
  CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
  CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
  r = [CXx*CXx, CXx*CYx, CXx*CZx;
       CYx*CXx, CYx*CYx, CYx*CZx;
       CZx*CXx, CZx*CYx, CZx*CZx];
  
  % Matrica krutosti elementa u lokalnom ks elementa
  k_e = (E*A(i))/L_e;
  
  T = [r, -r; -r, r];  % Matrica transformacije    
  
  K(edof, edof) = K(edof, edof)+k_e*T; % Matrica krutosti elementa
end  

% - Resavanje MKE
% Resavanje redukovane jednacine konacnih elemenata
D1 = K(fdof, fdof)\F(fdof);
D(fdof) = D1;
F1 = K*D;
R = F1(rdof); % vektor reakcija veza

for i=1:nen
  j = en(i, :);       
  edof = [3*j(1)-2, 3*j(1)-1, 3*j(1) ... % stepeni slobode elementa
         3*j(2)-2, 3*j(2)-1, 3*j(2)]; 
  L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % duzina elementa
      nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
      nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
  CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
  CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
  CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
  sigma(i) = E/L_e*[-CXx, -CYx, -CZx, CXx, CYx, CZx]*D(edof); % [Pa] Normalni napon
end 

% Odredjivanje sila u stapovima
Q = sigma.*A(1:nen); % [N] Sila u stapu

% - Prikaz rezultata
% Prikaz vrednosti pomeranja i reakcija veza
disp(' -------------------- ');
disp(' Pomeranja')
i = 1:dof;
Dn = reshape(D, 3, []); % pomeranja cvorova
Dv = [reshape(repmat(1:nn, 3, 1), 1, []);...
  reshape(repmat(1:3, nn, 1).', [], 1).'; Dn(i)*s(1)];
disp(' Cvor | Komponenta | Pomeranje [mm]');
fprintf(' %3d | %3d | %14.10f\n', Dv);

% Prikaz vrednosti reakcija
disp(' -------------------- ');
disp(' Reakcije veza')
nrdof = mod(rdof, 3); % ograniceni stepen slobode u cvoru
nrdof(nrdof == 0) = ones(length(find(nrdof == 0)), 1)*3;
nr = (rdof-nrdof)/3+1; % cvor sa ogranicenjem
Fv = [nr.', nrdof.', R].';
disp(' Cvor | Komponenta | Reakcije [N]');
fprintf(' %3d | %3d | %14.10f\n', Fv);

% Prikaz sila u stapovima
disp(' -------------------- ');
disp(' Sile u stapovima')
Qv = [enn(:, 1), enn(:, 2), Q].';
disp(' Stap | Sila [N]');
fprintf(' %d-%d | %14.10f\n', Qv);

% - Prikaz proracunskog modela
% Prikaz polaznog modela sa opterecenjima
disp(' -------------------- ');
disp(' Prikaz proracunskog modela... ');
drawArrow = @(x, y, z, varargin) quiver3(x(1), y(1), z(1), ...
    x(2)+10^-5, y(2)+10^-5, z(2)+10^-5, 0, varargin{:});   
drawArrowMarker = @(x, y, z, m, varargin) [plot3([x(1); x(1)+x(2)], ...
  [y(1); y(1)+y(2)], [z(1); z(1)+z(2)], '-', varargin{:}), ...
  plot3(x(1)+x(2), y(1)+y(2), z(1)+z(2), m, varargin{:})];    

nc = nc*s(1);
rds = [[nc(en(:, 1), 1), nc(en(:, 2), 1)], ...
    [nc(en(:, 1), 2), nc(en(:, 2), 2)], ...
    [nc(en(:, 1), 3), nc(en(:, 2), 3)]]; % stapovi
  
figure(1);
box on, grid on, hold on
boje = get(gca, 'colororder');
plot3(rds(1:nen, 1:2).', rds(1:nen, 3:4).', rds(1:nen, 5:6).', ...
    'LineWidth', 2, 'Color', boje(2, :));
plot3(rds(nen+1:end, 1:2).', rds(nen+1:end, 3:4).', rds(nen+1:end, 5:6).', ...
    'LineWidth', 2, 'Color', boje(6, :));
plot3(nc(:, 1), nc(:, 2), nc(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', boje(2, :), 'Color', boje(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
Xlims = [min(nc(:, 1)), max(nc(:, 1))];
Ylims = [min(nc(:, 2)), max(nc(:, 2))];
Zlims = [min(nc(:, 3)), max(nc(:, 3))];
xlim([Xlims(1)-s(5), Xlims(2)+s(5)]);
ylim([Ylims(1)-s(5), Ylims(2)+s(5)]);
zlim([Zlims(1)-s(5), Zlims(2)+s(5)]);

% Prikaz konturnih uslova
for i=1:length(rdof)
  p = zeros(1, 3);
  p(nrdof(i)) = 1*s(3);
  drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
    [nc(nr(i), 3), p(3)], '+', 'Color', 'b', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b');
end

% Prikaz spoljasnjih opterecenja
for i=1:nn
  p = zeros(1, 3);
  if(norm(F((i-1)*3+1:(i-1)*3+3)) > 0)
    p = (F((i-1)*3+1:(i-1)*3+3))*s(4);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b', 'MaxHeadSize', s(6));
  end
end

% Oznake cvorova
text(nc(:, 1)+s(7), nc(:, 2)+s(7), nc(:, 3), ...
   num2str((1:nn).'), 'FontSize', 18, 'Color', boje(2, :));
print('fig-1.png', '-dpng', '-F:18');
print('fig-1.svg', '-dsvg', '-FCMU Serif:18');

% Prikaz deformisanog modela sa opterecenjima
disp(' Prikaz deformisanog modela... ');
ncd = nc+(Dn.'.*s(1)).*s(2);
rdsd = [[ncd(en(:, 1), 1), ncd(en(:, 2), 1)], ...
    [ncd(en(:, 1), 2), ncd(en(:, 2), 2)], ...
    [ncd(en(:, 1), 3), ncd(en(:, 2), 3)]]; % stapovi

figure(2);
box on, grid on, hold on
boje = get(gca, 'colororder');
plot3(rds(:, 1:2).', rds(:, 3:4).', rds(:, 5:6).', '--', ...
    'LineWidth', 1, 'Color', boje(1, :));
plot3(rdsd(1:nen, 1:2).', rdsd(1:nen, 3:4).', rdsd(1:nen, 5:6).', ...
    'LineWidth', 2, 'Color', boje(2, :));
plot3(rdsd(nen+1:end, 1:2).', rdsd(nen+1:end, 3:4).', rdsd(nen+1:end, 5:6).', ...
    'LineWidth', 2, 'Color', boje(6, :));
plot3(ncd(:, 1), ncd(:, 2), ncd(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', boje(2, :), 'Color', boje(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([Xlims(1)-s(5), Xlims(2)+s(5)]);
ylim([Ylims(1)-s(5), Ylims(2)+s(5)]);
zlim([Zlims(1)-s(5), Zlims(2)+s(5)]);
print('fig-2.png', '-dpng', '-F:18');
print('fig-2.svg', '-dsvg', '-FCMU Serif:18');

% - Kraj programa
disp(' -------------------- ');
disp(' Program je uspesno izvrsen... ');
disp([' Potrebno vreme: ' num2str(toc, '%.2f') ' sekundi']);
disp(' -------------------- ');