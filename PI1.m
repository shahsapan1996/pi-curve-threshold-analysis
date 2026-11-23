%% PI_batch.m - threshold analysis for several PI files at once
clear; close all; clc

%% ---------- SETTINGS (edit these) ----------
colI      = 1;        % column with gain current
colP      = 2;        % column with power
splitCol  = 3;        % column that separates sweeps; [] = do not split
unitScale = 1000;     % A -> mA and W -> mW (use 1 if already in mA and mW)
nOffset   = 3;        % first N points used for dark-offset correction
lo = 70;  hi = 120;   % mA, linear-fit window (just above the knee)
lo_set = lo + [-5 0 5 10];      % window variations for the uncertainty
hi_set = hi + [-10 0 10];
kneeRange = [30 100]; % mA, region searched for the d2P/dI2 peak
sgWin     = 5;        % Savitzky-Golay smoothing window (odd number)
saveFigs  = true;     % save each figure as PNG next to the data
%% -------------------------------------------

[files, folder] = uigetfile({'*.csv;*.txt;*.xlsx;*.xls;*.dat','Data files'}, ...
    'Select PI data files', 'MultiSelect', 'on');
if isequal(files, 0), disp('No file selected'); return; end
if ischar(files), files = {files}; end

rows = {};

for n = 1:numel(files)
    fname = files{n};
    T = readtable(fullfile(folder, fname), 'VariableNamingRule', 'preserve');
    Iraw = T{:,colI};
    Praw = T{:,colP};

    if isempty(splitCol)
        grp = ones(height(T),1);
    else
        grp = round(T{:,splitCol}, 6);
        if numel(unique(grp)) > 5
            warning('%s: column %d has %d different values, not splitting', ...
                fname, splitCol, numel(unique(grp)));
            grp = ones(height(T),1);
        end
    end
    levels = unique(grp);
    fprintf('%s: %d rows, %d group(s)\n', fname, height(T), numel(levels));

    for k = 1:numel(levels)
        s = (grp == levels(k));
        I = Iraw(s)*unitScale;
        P = Praw(s)*unitScale;
        [I, iu] = unique(I);
        P = P(iu);
        if numel(I) < 8
            warning('%s (group %d): only %d points, skipped', fname, k, numel(I));
            continue
        end
        P = P - median(P(1:min(nOffset, numel(P))));

        % ---- Method B: second derivative ----
        Ps  = smoothdata(P, 'sgolay', sgWin, 'SamplePoints', I);
        d1  = gradient(Ps, I);
        d2  = gradient(d1, I);
        reg = I >= kneeRange(1) & I <= kneeRange(2);
        if nnz(reg) < 3
            warning('%s: too few points in knee range, skipped', fname);
            continue
        end
        Ireg = I(reg);  d2reg = d2(reg);
        [~, j] = max(d2reg);
        Ith_B = Ireg(j);

        % ---- Method A: linear fit above the knee ----
        m = I >= lo & I <= hi;
        if nnz(m) < 3
            warning('%s: too few points in fit window, skipped', fname);
            continue
        end
        pf = polyfit(I(m), P(m), 1);
        SE = pf(1);  b = pf(2);
        Ith_A = -b/SE;

        Ith_set = [];
        for a = lo_set
            for c = hi_set
                mm = I >= a & I <= c;
                if nnz(mm) >= 3
                    pp = polyfit(I(mm), P(mm), 1);
                    Ith_set(end+1) = -pp(2)/pp(1); %#ok<SAGROW>
                end
            end
        end
        Ith_A_err = std(Ith_set);

        % ---- Figure: PI curve (left) and d2P/dI2 (right) ----
        figure('Name', fname, 'Position', [100 100 1200 450]);
        tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
        title(tl, sprintf('%s   (group %d, value = %g)', fname, k, levels(k)), 'Interpreter', 'none');

        nexttile; hold on; grid on
        plot(I, P, 'o-');
        plot(I, SE*I + b, 'r-');
        xline(Ith_A, '--', sprintf('%.1f mA', Ith_A));
        xline(Ith_B, ':g', sprintf('%.1f mA', Ith_B));
        xlabel('Current (mA)'); ylabel('Power (mW)');
        ylim([min(P) max(P)*1.05]);
        legend('data', 'linear fit', 'Location', 'northwest');

        nexttile; hold on; grid on
        plot(I, d2, 'o-');
        xline(Ith_B, ':g');
        xlabel('Current (mA)'); ylabel('d^2P/dI^2 (mW/mA^2)');

        if saveFigs
            [~, base] = fileparts(fname);
            exportgraphics(gcf, fullfile(folder, sprintf('%s_group%d.png', base, k)), 'Resolution', 300);
        end

        fprintf('%s | group %d | I_th(lin) = %.1f +/- %.1f mA | SE = %.4f W/A | I_th(d2) = %.1f mA\n', ...
            fname, k, Ith_A, Ith_A_err, SE, Ith_B);
        rows(end+1,:) = {string(fname), k, Ith_A, Ith_A_err, SE, Ith_B}; %#ok<SAGROW>
    end
end

%% Summary table
if ~isempty(rows)
    results = cell2table(rows, 'VariableNames', ...
        {'File','Group','Ith_linear_mA','Ith_linear_err_mA','SE_W_per_A','Ith_d2_mA'});
    disp(results)
    writetable(results, fullfile(folder, 'PI_threshold_summary.csv'));
end