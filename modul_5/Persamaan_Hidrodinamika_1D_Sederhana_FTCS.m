clear; clc; close all;

identitas = 'Feby Syarief Al A`raaf - 26050124130087 - Oseanografi C';

folder_u_waktu = '01_Grafik_Kecepatan_Terhadap_Waktu';
folder_u_ruang = '02_Grafik_Kecepatan_Terhadap_Ruang';
folder_z_waktu = '03_Grafik_Elevasi_Terhadap_Waktu';
folder_z_ruang = '04_Grafik_Elevasi_Terhadap_Ruang';

if ~exist(folder_u_waktu, 'dir'); mkdir(folder_u_waktu); end
if ~exist(folder_u_ruang, 'dir'); mkdir(folder_u_ruang); end
if ~exist(folder_z_waktu, 'dir'); mkdir(folder_z_waktu); end
if ~exist(folder_z_ruang, 'dir'); mkdir(folder_z_ruang); end

g = 9.81;

scenario(1).nama = 'Skenario 1';
scenario(1).L  = 6000;
scenario(1).dx = 30;
scenario(1).LS = 2250;
scenario(1).dt = 5;
scenario(1).A  = 1.780;
scenario(1).T  = 450;
scenario(1).H  = 3;

scenario(2).nama = 'Skenario 2';
scenario(2).L  = 6000;
scenario(2).dx = 60;
scenario(2).LS = 2250;
scenario(2).dt = 5;
scenario(2).A  = 1.780;
scenario(2).T  = 450;
scenario(2).H  = 3;

scenario(3).nama = 'Skenario 3';
scenario(3).L  = 6000;
scenario(3).dx = 15;
scenario(3).LS = 2250;
scenario(3).dt = 5;
scenario(3).A  = 1.780;
scenario(3).T  = 450;
scenario(3).H  = 3;

scenario(4).nama = 'Skenario 4';
scenario(4).L  = 6000;
scenario(4).dx = 30;
scenario(4).LS = 2250;
scenario(4).dt = 2;
scenario(4).A  = 1.780;
scenario(4).T  = 450;
scenario(4).H  = 3;

scenario(5).nama = 'Skenario 5';
scenario(5).L  = 6000;
scenario(5).dx = 30;
scenario(5).LS = 2250;
scenario(5).dt = 2;
scenario(5).A  = 1.080;
scenario(5).T  = 450;
scenario(5).H  = 3;

scenario(6).nama = 'Skenario 6';
scenario(6).L  = 6000;
scenario(6).dx = 30;
scenario(6).LS = 2250;
scenario(6).dt = 2;
scenario(6).A  = 2.780;
scenario(6).T  = 450;
scenario(6).H  = 3;

scenario(7).nama = 'Skenario 7';
scenario(7).L  = 6000;
scenario(7).dx = 30;
scenario(7).LS = 2250;
scenario(7).dt = 2;
scenario(7).A  = 1.780;
scenario(7).T  = 1350;
scenario(7).H  = 3;

scenario(8).nama = 'Skenario 8';
scenario(8).L  = 6000;
scenario(8).dx = 30;
scenario(8).LS = 2250;
scenario(8).dt = 2;
scenario(8).A  = 1.780;
scenario(8).T  = 450;
scenario(8).H  = 10;

for s = 1:length(scenario)

    L  = scenario(s).L;
    dx = scenario(s).dx;
    LS = scenario(s).LS;
    dt = scenario(s).dt;
    A  = scenario(s).A;
    T  = scenario(s).T;
    H  = scenario(s).H;

    x  = 0:dx:L;
    t  = 0:dt:LS;
    Nx = length(x);
    Nt = length(t);

    Z  = zeros(Nx, Nt);
    u  = zeros(Nx, Nt);

    Cs    = sqrt(g * H);
    lam   = Cs * T;
    k     = 2*pi / lam;
    sigma = 2*pi / T;
    CFL   = Cs * dt / dx;

    fprintf('%s | Cs=%.4f m/s | lambda=%.2f m | CFL=%.4f\n', ...
        scenario(s).nama, Cs, lam, CFL);

    for j = 1:Nx
        Z(j,1) = A * cos(k * x(j));
        u(j,1) = (A/H) * Cs * cos(k * (x(j) + 0.5*dx));
    end

    for n = 1:Nt-1

        tn1 = t(n+1);

        for j = 1:Nx-1
            u(j, n+1) = u(j, n) - g * dt/dx * (Z(j+1, n) - Z(j, n));
        end

        u(Nx, n+1) = (A/H) * Cs * cos(k * L - sigma * tn1);

        for j = 2:Nx
            Z(j, n+1) = Z(j, n) - H * dt/dx * (u(j, n+1) - u(j-1, n+1));
        end

        Z(1, n+1) = A * cos(k * x(1) - sigma * tn1);

    end

    idx_ruang = [1 21 41 61 81 101 121 141 161 181 201];
    idx_ruang = idx_ruang(idx_ruang <= Nx);
    if idx_ruang(end) ~= Nx
        idx_ruang = [idx_ruang, Nx];
    end

    idx_waktu = [1 11 21 31 41 51 101 201 301 401 451];
    idx_waktu = idx_waktu(idx_waktu <= Nt);
    if idx_waktu(end) ~= Nt
        idx_waktu = [idx_waktu, Nt];
    end

    fig = figure('Visible','off');
    hold on; grid on; box on;
    for kk = 1:length(idx_waktu)
        plot(x, Z(:, idx_waktu(kk)), 'LineWidth', 1.4);
    end
    xlabel('x (m)'); ylabel('Z (m)');
    title({sprintf('Elevasi vs Ruang - %s', scenario(s).nama), identitas});
    legenda = cell(1, length(idx_waktu));
    for kk = 1:length(idx_waktu)
        legenda{kk} = sprintf('t = %.0f s', t(idx_waktu(kk)));
    end
    legend(legenda, 'Location','eastoutside');
    saveas(fig, fullfile(folder_z_ruang, sprintf('elevasi_vs_ruang_skenario_%d.png',s)));
    close(fig);

    fig = figure('Visible','off');
    hold on; grid on; box on;
    for kk = 1:length(idx_ruang)
        plot(t, Z(idx_ruang(kk), :), 'LineWidth', 1.4);
    end
    xlabel('t (s)'); ylabel('Z (m)');
    title({sprintf('Elevasi vs Waktu - %s', scenario(s).nama), identitas});
    legenda = cell(1, length(idx_ruang));
    for kk = 1:length(idx_ruang)
        legenda{kk} = sprintf('j = %d', idx_ruang(kk));
    end
    legend(legenda, 'Location','eastoutside');
    saveas(fig, fullfile(folder_z_waktu, sprintf('elevasi_vs_waktu_skenario_%d.png',s)));
    close(fig);

    fig = figure('Visible','off');
    hold on; grid on; box on;
    for kk = 1:length(idx_waktu)
        plot(x, u(:, idx_waktu(kk)), 'LineWidth', 1.4);
    end
    xlabel('x (m)'); ylabel('u (m/s)');
    title({sprintf('Kecepatan vs Ruang - %s', scenario(s).nama), identitas});
    legenda = cell(1, length(idx_waktu));
    for kk = 1:length(idx_waktu)
        legenda{kk} = sprintf('t = %.0f s', t(idx_waktu(kk)));
    end
    legend(legenda, 'Location','eastoutside');
    saveas(fig, fullfile(folder_u_ruang, sprintf('kecepatan_vs_ruang_skenario_%d.png',s)));
    close(fig);

    fig = figure('Visible','off');
    hold on; grid on; box on;
    for kk = 1:length(idx_ruang)
        plot(t, u(idx_ruang(kk), :), 'LineWidth', 1.4);
    end
    xlabel('t (s)'); ylabel('u (m/s)');
    title({sprintf('Kecepatan vs Waktu - %s', scenario(s).nama), identitas});
    legenda = cell(1, length(idx_ruang));
    for kk = 1:length(idx_ruang)
        legenda{kk} = sprintf('j = %d', idx_ruang(kk));
    end
    legend(legenda, 'Location','eastoutside');
    saveas(fig, fullfile(folder_u_waktu, sprintf('kecepatan_vs_waktu_skenario_%d.png',s)));
    close(fig);

end

disp('Seluruh grafik selesai dibuat dan disimpan ke dalam 4 folder.');