clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C

% Simulasi adveksi 1 dimensi dengan metode Leapfrog (CTCS)
% Sumber polutan kontinu
% NIM ganjil, xyz = 087, sehingga z = 7 dan 0.zyx = 0.780

x_digit = 0;
y_digit = 8;
z_digit = 7;

u_xyz    = 0.780;
C_source = 10 * z_digit;

n_plot = [0, ...
          10 + z_digit, 20 + z_digit, 30 + z_digit, 40 + z_digit, ...
          50 + z_digit, 100 + z_digit, 200 + z_digit, 300 + z_digit, ...
          500 + z_digit, 800 + z_digit, 1000 + z_digit];

L = 3000;
T = 10800;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 30;
scenario(1).dt   = 8;
scenario(1).u    = u_xyz;
scenario(1).m    = 30 + z_digit;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 60;
scenario(2).dt   = 8;
scenario(2).u    = u_xyz;
scenario(2).m    = 60 + z_digit;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 15;
scenario(3).dt   = 8;
scenario(3).u    = u_xyz;
scenario(3).m    = 15 + z_digit;

scenario(4).nama = 'Skenario 4';
scenario(4).dx   = 30;
scenario(4).dt   = 8;
scenario(4).u    = 0.5 * u_xyz;
scenario(4).m    = 30 + z_digit;

scenario(5).nama = 'Skenario 5';
scenario(5).dx   = 30;
scenario(5).dt   = 8;
scenario(5).u    = 2 * u_xyz;
scenario(5).m    = 30 + z_digit;

scenario(6).nama = 'Skenario 6';
scenario(6).dx   = 30;
scenario(6).dt   = 8;
scenario(6).u    = -u_xyz;
scenario(6).m    = 30 + z_digit;

scenario(7).nama = 'Skenario 7';
scenario(7).dx   = 30;
scenario(7).dt   = 2;
scenario(7).u    = u_xyz;
scenario(7).m    = 30 + z_digit;

scenario(8).nama = 'Skenario 8';
scenario(8).dx   = 30;
scenario(8).dt   = 16;
scenario(8).u    = u_xyz;
scenario(8).m    = 30 + z_digit;

folder_output = 'output';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx       = scenario(s).dx;
    dt       = scenario(s).dt;
    u        = scenario(s).u;
    m_source = scenario(s).m;

    M_grid = round(L / dx);
    Mmax   = max(M_grid, m_source);
    Nmax   = round(T / dt);

    courant = u * dt / dx;

    fprintf('dx = %.3f m, dt = %.3f s, u = %.3f m/s, m = %d\n', dx, dt, u, m_source);
    fprintf('Mmax = %d, Nmax = %d, C = %.4f\n', Mmax, Nmax, courant);

    if any(n_plot > Nmax)
        n_tidak_valid = n_plot(n_plot > Nmax);
        fprintf('%s: n = %s melebihi Nmax = %d sehingga tidak digrafikkan.\n', ...
                scenario(s).nama, mat2str(n_tidak_valid), Nmax);
    end

    F00 = zeros(1, Mmax);
    F0  = zeros(1, Mmax);
    F   = zeros(1, Mmax);

    x_space = (1:Mmax) * dx;
    snapshot = NaN(length(n_plot), Mmax);

    F00(:) = 0;

    idx0 = find(n_plot == 0, 1);
    if ~isempty(idx0)
        snapshot(idx0, :) = F00;
    end

    for i = 2:(Mmax-1)
        F0(i) = F00(i) - (u * dt / (2 * dx)) * (F00(i+1) - F00(i-1));
    end

    F0(1)    = F0(2);
    F0(Mmax) = F0(Mmax-1);
    F0(m_source) = C_source;

    for J = 2:Nmax

        F00(1)    = F00(2);
        F00(Mmax) = F00(Mmax-1);

        F0(1)     = F0(2);
        F0(Mmax)  = F0(Mmax-1);

        F0(m_source) = C_source;

        for i = 2:(Mmax-1)
            F(i) = F00(i) - (u * dt / dx) * (F0(i+1) - F0(i-1));
        end

        F(1)    = F(2);
        F(Mmax) = F(Mmax-1);
        F(m_source) = C_source;

        idx_snap = find(n_plot == J, 1);
        if ~isempty(idx_snap)
            snapshot(idx_snap, :) = F;
        end

        F00 = F0;
        F0  = F;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(n_plot)

        if all(isnan(snapshot(k, :)))
            continue;
        end

        plot(x_space, snapshot(k, :), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('n = %d', n_plot(k));
    end

    grid on;
    xlabel('Ruang x (m)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap ruang', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C','FontSize', 14, 'FontWeight', 'bold');
    xlim([x_space(1), x_space(end)]);

    y_max = max(snapshot(:), [], 'omitnan');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, max(C_source * 1.2, y_max * 1.1)]);
    end

    set(gca, 'FontSize', 11);

    legend(legend_text, ...
           'Location', 'northwest', ...
           'FontSize', 8, ...
           'NumColumns', 2, ...
           'Box', 'on');

    hold off;

    nama_file = sprintf('ctcs_kontinu_ruang_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);
end

fprintf('\nSemua skenario telah selesai diproses.\n');