clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
% PERSAMAAN ADVEKSI-DIFUSI 2D
% METODE UPSTREAM
% SUMBER POLUTAN KONTINU
% GRAFIK TERHADAP WAKTU

identitas = 'Feby Syarief Al A`raaf - 26050124130087 - Oseanografi C';

% Parameter umum domain
Lx = 4000;
Ly = 4000;
dx = 100;
dy = 100;
T  = 3600;
dt = 10;

nx = Lx/dx + 1;
ny = Ly/dy + 1;
nt = T/dt;

x = 0:dx:Lx;
y = 0:dy:Ly;
waktu = (1:nt) * dt;

% Konsentrasi sumber kontinu
F = 100;

folder_output = 'output_Kontinu_Upstream_Adveksi_Difusi_2D_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

scenario(1).nama = 'Skenario 1';
scenario(1).Ad   = 2.5;
scenario(1).u    = 0.05;
scenario(1).v    = -0.5;
scenario(1).src  = [21, 21];

scenario(2).nama = 'Skenario 2';
scenario(2).Ad   = 25;
scenario(2).u    = 0.05;
scenario(2).v    = -0.5;
scenario(2).src  = [21, 21];

scenario(3).nama = 'Skenario 3';
scenario(3).Ad   = 2.5;
scenario(3).u    = -0.2;
scenario(3).v    = 0.01;
scenario(3).src  = [21, 21];

scenario(4).nama = 'Skenario 4';
scenario(4).Ad   = 10;
scenario(4).u    = 0.15;
scenario(4).v    = -0.10;
scenario(4).src  = [21, 21];

% 5 grid representatif
grid_pantau = [
    11 11
    16 16
    21 21
    26 26
    31 31
];
ngrid = size(grid_pantau, 1);

for s = 1:4
    
    Ad = scenario(s).Ad;
    u  = scenario(s).u;
    v  = scenario(s).v;
    
    isrc = scenario(s).src(1);
    jsrc = scenario(s).src(2);
    
    C = zeros(nx, ny);
    C(isrc, jsrc) = F;
    
    rx = Ad * dt / dx^2;
    ry = Ad * dt / dy^2;
    ax = u  * dt / dx;
    ay = v  * dt / dy;
    
    % Penyimpanan konsentrasi terhadap waktu
    Ct = zeros(ngrid, nt);
    
    for n = 1:nt
        Cn = C;
        
        for i = 2:nx-1
            for j = 2:ny-1
                
                if u >= 0
                    adv_x = ax * (Cn(i,j) - Cn(i-1,j));
                else
                    adv_x = ax * (Cn(i+1,j) - Cn(i,j));
                end
                
                if v >= 0
                    adv_y = ay * (Cn(i,j) - Cn(i,j-1));
                else
                    adv_y = ay * (Cn(i,j+1) - Cn(i,j));
                end
                
                C(i,j) = Cn(i,j) ...
                    - adv_x ...
                    - adv_y ...
                    + rx * (Cn(i+1,j) - 2*Cn(i,j) + Cn(i-1,j)) ...
                    + ry * (Cn(i,j+1) - 2*Cn(i,j) + Cn(i,j-1));
            end
        end
        
        C(1,:)   = C(2,:);
        C(end,:) = C(end-1,:);
        C(:,1)   = C(:,2);
        C(:,end) = C(:,end-1);
        
        C(isrc, jsrc) = F;
        
        % Simpan konsentrasi di 5 grid representatif
        for k = 1:ngrid
            i_plot = grid_pantau(k,1);
            j_plot = grid_pantau(k,2);
            Ct(k,n) = C(i_plot, j_plot);
        end
    end
    
    % Plot 5 subplot terhadap waktu
    fig = figure('Color','w','Position',[100 100 1400 800]);
    tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for k = 1:ngrid
        nexttile
        plot(waktu, Ct(k,:), 'LineWidth', 1.5)
        grid on
        xlabel('Waktu (detik)')
        ylabel('Konsentrasi (mg/L)')
        title(sprintf('%s - Grid (%d,%d)', ...
            scenario(s).nama, grid_pantau(k,1), grid_pantau(k,2)))
    end
    
    nexttile
    axis off
    text(0.1, 0.8, scenario(s).nama, 'FontSize', 14, 'FontWeight', 'bold')
    text(0.1, 0.6, sprintf('Ad = %.3f m^2/detik', scenario(s).Ad), 'FontSize', 12)
    text(0.1, 0.45, sprintf('u = %.3f m/detik', scenario(s).u), 'FontSize', 12)
    text(0.1, 0.30, sprintf('v = %.3f m/detik', scenario(s).v), 'FontSize', 12)
    text(0.1, 0.15, sprintf('Sumber = (%d, %d)', isrc, jsrc), 'FontSize', 12)
    
    sgtitle({
        sprintf('Persamaan Adveksi-Difusi 2D - Sumber Polutan Kontinu - Grafik terhadap Waktu - %s', scenario(s).nama), ...
        identitas
    })
    
    nama_file = sprintf('upstream_adveksi_difusi_2D_waktu_kontinu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);
end