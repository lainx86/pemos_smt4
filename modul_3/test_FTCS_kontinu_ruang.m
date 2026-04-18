clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi adveksi-difusi 1 dimensi dengan metode FTCS
% Sumber polutan kontinu
% Grafik: konsentrasi polutan terhadap ruang

x_digit = 0;
y_digit = 8;
z_digit = 7;

u_xyz    = 0.780;
C_source = 10 * z_digit;   % = 70 mg/L

% Waktu yang akan diplot (dalam indeks langkah waktu n)
n_plot = [0, ...
          10 + z_digit, 20 + z_digit, 30 + z_digit, 40 + z_digit, ...
          50 + z_digit, 100 + z_digit, 200 + z_digit, 300 + z_digit, ...
          500 + z_digit, 800 + z_digit, 1000 + z_digit];

% === PARAMETER DOMAIN (dari tabel) ===
L = 3000;   % panjang domain (m)
T = 7200;   % waktu simulasi (s)  <-- diubah dari 10800 ke 7200

% Grid sumber polutan kontinu = 15+z
m_kontinu = 15 + z_digit;   % = 22

% === SKENARIO (dari tabel) ===
scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 6;        % diubah dari 8 ke 6
scenario(1).u    = 0.5;      % diubah
scenario(1).D    = 4.0;      % diubah
scenario(1).m    = m_kontinu;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 30;
scenario(2).dt   = 6;
scenario(2).u    = 0.1;
scenario(2).D    = 2.0;
scenario(2).m    = m_kontinu;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 30;
scenario(3).dt   = 6;
scenario(3).u    = 0.5;
scenario(3).D    = 2.0;
scenario(3).m    = m_kontinu;

scenario(4).nama = 'Skenario 4';
scenario(4).dx   = 30;
scenario(4).dt   = 6;
scenario(4).u    = 0.1;
scenario(4).D    = 4.0;
scenario(4).m    = m_kontinu;

% === OUTPUT ===
folder_output = 'output_Kontinu_FTCS_Adveksi_Difusi_Terhadap_Ruang';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    u        = scenario(s).u;
    D        = scenario(s).D;
    m_source = scenario(s).m;

    x_space = 0:dx:L;
    t_grid  = 0:dt:T;

    Mmax = length(x_space);
    Nmax = length(t_grid) - 1;

    courant = u * dt / dx;
    alpha   = D * dt / dx^2;

    fprintf('dx = %.3f m, dt = %.3f s, u = %.3f m/s, D = %.3f m^2/s, m = %d\n', ...
            dx, dt, u, D, m_source);
    fprintf('Mmax = %d, Nmax = %d, Courant = %.4f, alpha = %.6f\n', ...
            Mmax, Nmax, courant, alpha);

    if m_source > Mmax
        fprintf('%s tidak dapat diproses karena grid sumber di luar domain.\n', scenario(s).nama);
        continue;
    end

    if any(n_plot > Nmax)
        n_tidak_valid = n_plot(n_plot > Nmax);
        fprintf('%s: n = %s melebihi Nmax = %d sehingga tidak digrafikkan.\n', ...
                scenario(s).nama, mat2str(n_tidak_valid), Nmax);
    end

    if alpha > 0.5
        fprintf('%s: peringatan, alpha > 0.5 sehingga skema FTCS berpotensi tidak stabil.\n', ...
                scenario(s).nama);
    end

    F0 = zeros(1, Mmax);
    F  = zeros(1, Mmax);

    snapshot = NaN(length(n_plot), Mmax);

    F0(:) = 0;
    F0(m_source) = C_source;

    idx0 = find(n_plot == 0, 1);
    if ~isempty(idx0)
        snapshot(idx0, :) = F0;
    end

    for J = 1:Nmax

        F0(1)    = F0(2);
        F0(Mmax) = F0(Mmax - 1);
        F0(m_source) = C_source;

        for i = 2:(Mmax - 1)
            if i == m_source
                F(i) = C_source;
            else
                F(i) = F0(i) ...
                     - (u * dt / (2 * dx)) * (F0(i + 1) - F0(i - 1)) ...
                     + (D * dt / dx^2) * (F0(i + 1) - 2 * F0(i) + F0(i - 1));
            end
        end

        F(1)    = F(2);
        F(Mmax) = F(Mmax - 1);
        F(m_source) = C_source;

        idx_snap = find(n_plot == J, 1);
        if ~isempty(idx_snap)
            snapshot(idx_snap, :) = F;
        end

        F0 = F;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(n_plot)
        if all(isnan(snapshot(k, :)))
            continue;
        end
        plot(x_space, snapshot(k, :), 'LineWidth', 1.2);
        legend_text{end + 1} = sprintf('n = %d', n_plot(k));
    end

    grid on;
    xlabel('Ruang x (m)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap ruang', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([x_space(1), x_space(end)]);

    y_max = max(snapshot(:), [], 'omitnan');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, max(C_source * 1.2, y_max * 1.1)]);
    end

    set(gca, 'FontSize', 11);
    legend(legend_text, 'Location', 'northwest', 'FontSize', 8, ...
           'NumColumns', 2, 'Box', 'on');
    hold off;

    nama_file = sprintf('ftcs_adveksi_difusi_ruang_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');