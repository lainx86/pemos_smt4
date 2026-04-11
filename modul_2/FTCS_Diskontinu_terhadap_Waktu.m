clear; clc; close all;

% Feby Syarief Al A`raaf
% 26050124130087
% Oseanografi C
%
% Simulasi difusi 1 dimensi metode FTCS
% Sumber polutan diskontinu
% Grafik: konsentrasi polutan terhadap waktu
% NIM = 087, sehingga:
% x = 0, y = 8, z = 7
% Konsentrasi sumber = 10*z = 70 mg/L

x = 0;
y = 8;
z = 7;

folder_output = 'output_Diskontinu_FTCS_Terhadap_Waktu';
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

Ad1 = 0;
L   = 3000;
T   = 7200;
C   = 10 * z;

scenario(1).nama = 'Skenario 1';
scenario(1).dx   = 15;
scenario(1).dt   = 6;
scenario(1).Ad2  = 0.5;
scenario(1).n    = 20 + z;
scenario(1).m    = 15 + z;

scenario(2).nama = 'Skenario 2';
scenario(2).dx   = 30;
scenario(2).dt   = 6;
scenario(2).Ad2  = 0.05;
scenario(2).n    = 20 + z;
scenario(2).m    = 30 + z;

scenario(3).nama = 'Skenario 3';
scenario(3).dx   = 15;
scenario(3).dt   = 6;
scenario(3).Ad2  = 1;
scenario(3).n    = 20 + z;
scenario(3).m    = 15 + z;

for s = 1:numel(scenario)

    fprintf('\nMemproses %s\n', scenario(s).nama);

    dx   = scenario(s).dx;
    dt   = scenario(s).dt;
    Ad2  = scenario(s).Ad2;
    nsrc = scenario(s).n;
    m    = scenario(s).m;

    n     = 0:dt:T;
    t     = 0:dx:L;
    nmax  = length(n);
    mmax  = length(t);

    alpha1 = (Ad1 * dt) / dx^2;
    alpha2 = (Ad2 * dt) / dx^2;

    fprintf('dx = %.3f m, dt = %.3f s, Ad2 = %.3f m^2/s, n = %d, m = %d\n', dx, dt, Ad2, nsrc, m);
    fprintf('nmax = %d, mmax = %d, alpha2 = %.6f\n', nmax, mmax, alpha2);

    if m > mmax
        fprintf('%s tidak dapat diproses karena grid sumber di luar domain.\n', scenario(s).nama);
        continue;
    end

    if nsrc > nmax
        fprintf('%s tidak dapat diproses karena waktu sumber di luar domain.\n', scenario(s).nama);
        continue;
    end

    if alpha2 > 0.5
        fprintf('%s: peringatan, alpha2 > 0.5 sehingga skema FTCS berpotensi tidak stabil.\n', scenario(s).nama);
    end

    F = zeros(nmax, mmax);

    rv = [10 12 14 16 19 21 23];
    rv_valid = rv(rv >= 1 & rv <= mmax);

    if length(rv_valid) < length(rv)
        fprintf('%s: beberapa grid rv berada di luar domain sehingga tidak digrafikkan.\n', scenario(s).nama);
    end

    history = NaN(nmax, length(rv_valid));

    F(nsrc, m) = C;

    for p = 1:length(rv_valid)
        history(1, p) = F(1, rv_valid(p));
    end

    for na = 1:(nmax - 1)

        if na == nsrc
            F(na, m) = C;
        end

        for i = 2:(mmax - 1)

            F(na + 1, i) = (1 - 2 * alpha2) * F(na, i) ...
                         + alpha2 * (F(na, i + 1) + F(na, i - 1));

        end

        F(na + 1, 1)    = F(na + 1, 2);
        F(na + 1, mmax) = F(na + 1, mmax - 1);

        for p = 1:length(rv_valid)
            history(na + 1, p) = F(na + 1, rv_valid(p));
        end

    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 600]);
    hold on;

    legend_text = {};

    for k = 1:length(rv_valid)
        plot(n, history(:, k), 'LineWidth', 1.2);
        legend_text{end+1} = sprintf('r = %d', rv_valid(k));
    end

    grid on;
    xlabel('Waktu t (detik)', 'FontSize', 12);
    ylabel('Konsentrasi polutan (mg/L)', 'FontSize', 12);
    title(sprintf('%s - Konsentrasi polutan terhadap waktu', scenario(s).nama), ...
          'Feby Syarief-0087-Ose C', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([n(1), n(end)]);

    y_max = max(history(:), [], 'omitnan');
    if isempty(y_max) || isnan(y_max) || y_max <= 0
        ylim([0 1]);
    else
        ylim([0, max(C * 1.2, y_max * 1.1)]);
    end

    set(gca, 'FontSize', 11);

    legend(legend_text, ...
           'Location', 'northwest', ...
           'FontSize', 8, ...
           'NumColumns', 2, ...
           'Box', 'on');

    hold off;

    nama_file = sprintf('ftcs_diskontinu_waktu_skenario_%d.png', s);
    saveas(fig, fullfile(folder_output, nama_file));
    close(fig);

    fprintf('%s selesai diproses.\n', scenario(s).nama);

end

fprintf('\nSemua skenario telah selesai diproses.\n');