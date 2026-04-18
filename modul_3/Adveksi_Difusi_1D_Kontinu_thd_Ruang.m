clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi persamaan adveksi-difusi 1 dimensi dengan metode FTCS
% Konsentrasi polutan terhadap ruang pada beberapa waktu tertentu
%
% Catatan ilmiah:
% Pada script ini, "sumber polutan kontinu" dimodelkan sebagai
% titik grid internal dengan konsentrasi tetap (internal Dirichlet source),
% bukan sebagai source term eksplisit S(x,t).
%
% Persamaan yang disimulasikan:
% dC/dt + u dC/dx = D d2C/dx2
%
% Skema FTCS:
% C_i^(n+1) = C_i^n
%           - (u dt / (2 dx)) (C_(i+1)^n - C_(i-1)^n)
%           + (D dt / dx^2)   (C_(i+1)^n - 2C_i^n + C_(i-1)^n)

x_digit = 0;
y_digit = 8;
z_digit = 7;

C_source = 10 * z_digit;

n_plot = [0 300 600 900 1200];

L = 3000;
T = 7200;

m_kontinu = 15 + z_digit;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 6;
scenario(1).u    = 0.5;
scenario(1).D    = 4.0;
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

folder_output = 'output_Kontinu_FTCS_Adveksi_Difusi_Terhadap_Ruang';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\n==================================================\n');
    fprintf('Memproses %s\n', scenario(s).nama);

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
    fprintf('Mmax = %d, Nmax = %d\n', Mmax, Nmax);
    fprintf('Courant = %.6f, alpha = %.6f\n', courant, alpha);

    if m_source < 2 || m_source > Mmax - 1
        fprintf('%s tidak dapat diproses karena titik sumber terlalu dekat dengan batas/domain tidak valid.\n', ...
            scenario(s).nama);
        continue;
    end

    if any(n_plot > Nmax)
        n_tidak_valid = n_plot(n_plot > Nmax);
        fprintf('%s: n = %s melebihi Nmax = %d, sehingga tidak digrafikkan.\n', ...
            scenario(s).nama, mat2str(n_tidak_valid), Nmax);
    end

    stabil_diffusi = (alpha <= 0.5);
    stabil_gabungan = (courant^2 <= 2 * alpha);

    if ~stabil_diffusi
        fprintf(['PERINGATAN: alpha = %.6f > 0.5, sehingga skema eksplisit ', ...
                 'untuk komponen difusi berpotensi tidak stabil.\n'], alpha);
    end

    if ~stabil_gabungan
        fprintf(['PERINGATAN: Courant^2 = %.6f > 2*alpha = %.6f, ', ...
                 'sehingga kestabilan gabungan adveksi-difusi berpotensi tidak terpenuhi.\n'], ...
                 courant^2, 2*alpha);
    end

    if stabil_diffusi && stabil_gabungan
        fprintf('Syarat kestabilan numerik utama terpenuhi.\n');
    end

    C_old = zeros(1, Mmax);
    C_new = zeros(1, Mmax);

    snapshot = NaN(length(n_plot), Mmax);

    C_old(:) = 0;
    C_old(m_source) = C_source;

    idx0 = find(n_plot == 0, 1);
    if ~isempty(idx0)
        snapshot(idx0, :) = C_old;
    end

    for n = 1:Nmax

        C_old(1)    = C_old(2);
        C_old(Mmax) = C_old(Mmax - 1);

        for i = 2:(Mmax - 1)
            C_new(i) = C_old(i) ...
                     - (courant / 2) * (C_old(i + 1) - C_old(i - 1)) ...
                     + alpha * (C_old(i + 1) - 2 * C_old(i) + C_old(i - 1));
        end

        C_new(1)    = C_new(2);
        C_new(Mmax) = C_new(Mmax - 1);

        C_new(m_source) = C_source;

        C_new(C_new < 0) = 0;

        idx_snap = find(n_plot == n, 1);
        if ~isempty(idx_snap)
            snapshot(idx_snap, :) = C_new;
        end

        C_old = C_new;
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
        'FontSize', 14, 'FontWeight', 'bold');

    xlim([x_space(1), x_space(end)]);

    y_max = max(snapshot(:), [], 'omitnan');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, max(C_source * 1.2, y_max * 1.1)]);
    end

    set(gca, 'FontSize', 11);

    if ~isempty(legend_text)
        legend(legend_text, 'Location', 'northwest', ...
            'FontSize', 8, 'NumColumns', 2, 'Box', 'on');
    end

    hold off;

    nama_file = sprintf('ftcs_adveksi_difusi_ruang_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');
