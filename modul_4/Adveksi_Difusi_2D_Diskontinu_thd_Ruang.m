clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
% PERSAMAAN ADVEKSI-DIFUSI 2D
% SUMBER POLUTAN DISKONTINU
% GRAFIK TERHADAP RUANG

identitas = 'Feby Syarief Al A`raaf - 26050124130087 - Oseanografi C';

% Parameter umum domain
Lx = 4000;          % panjang perairan arah x (m)
Ly = 4000;          % lebar perairan arah y (m)
dx = 100;           % grid arah x (m)
dy = 100;           % grid arah y (m)
T  = 3600;          % lama simulasi (s)
dt = 10;            % langkah waktu (s)

nx = Lx/dx + 1;     % jumlah grid x
ny = Ly/dy + 1;     % jumlah grid y
nt = T/dt;          % jumlah langkah waktu

x = 0:dx:Lx;
y = 0:dy:Ly;

% Konsentrasi sumber diskontinu
nim = '26050124130087';
F = sum(nim - '0');     % jumlah digit NIM = 39 mg/L

folder_output = 'output_Diskontinu_Adveksi_Difusi_2D_Terhadap_Ruang';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

% 4 skenario
scenario(1).nama = 'Skenario 1';
scenario(1).Ad   = 2.5;
scenario(1).u    = 0.05;
scenario(1).v    = -0.5;

scenario(2).nama = 'Skenario 2';
scenario(2).Ad   = 25;
scenario(2).u    = 0.05;
scenario(2).v    = -0.5;

scenario(3).nama = 'Skenario 3';
scenario(3).Ad   = 2.5;
scenario(3).u    = -0.2;
scenario(3).v    = 0.01;

% Skenario 4 bebas
scenario(4).nama = 'Skenario 4';
scenario(4).Ad   = 10;
scenario(4).u    = 0.15;
scenario(4).v    = -0.10;

% Dua sumber diskontinu sesuai tabel
src1 = [29, 11];
src2 = [11, 29];

% Waktu representatif untuk 5 subplot
snapshot_steps = [1, round(nt/4), round(nt/2), round(3*nt/4), nt];
snapshot_steps = unique(snapshot_steps);
nsnap = length(snapshot_steps);

% Simulasi tiap skenario
for s = 1:4
    
    Ad = scenario(s).Ad;
    u  = scenario(s).u;
    v  = scenario(s).v;
    
    % Matriks konsentrasi awal
    C = zeros(nx, ny);
    
    % Sumber polutan diskontinu: hanya diberikan sekali di awal
    C(src1(1), src1(2)) = F;
    C(src2(1), src2(2)) = F;
    
    % Koefisien Upstream
    rx = Ad * dt / dx^2;
    ry = Ad * dt / dy^2;
    ax = u  * dt / dx;
    ay = v  * dt / dy;
    
    C_snap = cell(1, nsnap);
    idx_snap = 1;
    
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
        
        if idx_snap <= nsnap && n == snapshot_steps(idx_snap)
            C_snap{idx_snap} = C;
            idx_snap = idx_snap + 1;
        end
    end
    
    % Figure 1 skenario = 5 subplot
    fig = figure('Color','w','Position',[100 100 1400 800]);
    tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for k = 1:nsnap
        nexttile
        
        contourf(x/dx + 1, y/dy + 1, C_snap{k}', 20, 'LineStyle', 'none');
        hold on
        contour(x/dx + 1, y/dy + 1, C_snap{k}', 10, 'k');
        plot(src1(1), src1(2), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        plot(src2(1), src2(2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        hold off
        
        colorbar
        axis equal tight
        xlabel('grid x')
        ylabel('grid y')
        title(sprintf('%s - Waktu ke %d', scenario(s).nama, snapshot_steps(k)))
    end
    
    nexttile
    axis off
    text(0.1, 0.85, scenario(s).nama, 'FontSize', 14, 'FontWeight', 'bold')
    text(0.1, 0.68, sprintf('Ad = %.3f m^2/detik', scenario(s).Ad), 'FontSize', 12)
    text(0.1, 0.53, sprintf('u = %.3f m/detik', scenario(s).u), 'FontSize', 12)
    text(0.1, 0.38, sprintf('v = %.3f m/detik', scenario(s).v), 'FontSize', 12)
    text(0.1, 0.23, sprintf('F = %d mg/L', F), 'FontSize', 12)
    text(0.1, 0.08, sprintf('Sumber = (%d,%d) dan (%d,%d)', ...
        src1(1), src1(2), src2(1), src2(2)), 'FontSize', 12)
    
    sgtitle({
        sprintf('Persamaan Adveksi-Difusi 2D - Sumber Polutan Diskontinu - Grafik terhadap Ruang - %s', scenario(s).nama), ...
        identitas
    })
    
    nama_file = sprintf('adveksi_difusi_2D_ruang_diskontinu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);
end