function line_ex(varargin)
%LINE_EX - Example of different RANSAC functions line estimation
%          Need NO arguments.

g_nSampLen = 2;
g_nIter = 50;
g_dThreshold = 2;
g_nTestLen = 10;

if nargin == 0

    ExPath = {fullfile(pwd,'SAC'); fullfile(pwd,'LINE')};

    % add all functions to path
    for i = 1:length(ExPath)
        addpath(ExPath{i});
    end

    f = figure('MenuBar', 'none', 'Name', 'RANSAC Example: LINE');
    vmenu = uicontextmenu;
    set(f, 'UIContextMenu', vmenu);

    axes;
    uimenu(vmenu, 'Label', 'View Legend', 'Checked', 'on', 'Callback', 'line_ex(''LegendTurn'',gcbo)');

    %clc;
    fprintf('///////////////////////////////////////////////////////////////////////////////\n');
    fprintf('                         SAC TOOL: LINE EXAMPLE                                \n');
    fprintf('///////////////////////////////////////////////////////////////////////////////\n');
    
    [Data, RealModel, Stat] = GenData;
    plot(Data(1, :), Data(2, :), '.', 'Tag', 'DATA');
    DrawLine(RealModel, Data, 'Tag', 'Real Model', 'Color', [0, 0, 0]);
    drawnow;
    hold on;
    
    fprintf('\nGENERATED MODEL:\n data_len = %d\n gen_inliers = %d\n', size(Data, 2), sum(Stat.gen_inliers));
    fprintf(' gen_mean = %f\n gen_abs_mean = %f\n', Stat.gen_mean, Stat.gen_abs_mean);
    fprintf(' res_inliers = %d\n', sum(Stat.res_inliers));
    fprintf(' res_mean = %f\n res_abs_mean = %f\n\n', Stat.res_mean, Stat.res_abs_mean);

    [RANSAC_mask, RANSAC_model] = RANSAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, g_dThreshold);
    DrawLine(RANSAC_model, Data, 'Tag', 'RANSAC', 'Color', [1, 0, 0]);
    drawnow;
    DIST = Dist(RANSAC_model, Data(:, find(RANSAC_mask)));
    fprintf('RANSAC MODEL:\n inliers = %d\n', sum(RANSAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = RANSAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = RANSAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));

    
    [MSAC_mask, MSAC_model] = MSAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, g_dThreshold);
    DrawLine(MSAC_model, Data, 'Tag', 'MSAC', 'Color', [0, 1, 0]);
    drawnow;
    DIST = Dist(MSAC_model, Data(:, find(MSAC_mask)));
    fprintf('MSAC MODEL:\n inliers = %d\n', sum(MSAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = MSAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = MSAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));

    
    [RRANSAC_mask, RRANSAC_model] = RRANSAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, g_dThreshold, g_nTestLen);
    DrawLine(RRANSAC_model, Data, 'Tag', 'RRANSAC', 'Color', [0, 1, 1]);
    drawnow;
    DIST = Dist(RRANSAC_model, Data(:, find(RRANSAC_mask)));
    fprintf('RRANSAC MODEL:\n inliers = %d\n', sum(RRANSAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = RRANSAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = RRANSAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));

    
    [NAPSAC_mask, NAPSAC_model] = NAPSAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, g_dThreshold, 2);
    DrawLine(NAPSAC_model, Data, 'Tag', 'NAPSAC', 'Color', [0.5, 0.5, 1]);
    drawnow;    
    DIST = Dist(NAPSAC_model, Data(:, find(NAPSAC_mask)));
    fprintf('NAPSAC MODEL:\n inliers = %d\n', sum(NAPSAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = NAPSAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = NAPSAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));

    [ZHANGSAC_mask, ZHANGSAC_model] = ZHANGSAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, g_dThreshold, 2);
    DrawLine(ZHANGSAC_model, Data, 'Tag', 'ZHANGSAC', 'Color', [1, 0.5, 0]);
    drawnow;    
    DIST = Dist(ZHANGSAC_model, Data(:, find(ZHANGSAC_mask)));
    fprintf('ZHANGSAC MODEL:\n inliers = %d\n', sum(ZHANGSAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = ZHANGSAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = ZHANGSAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    
    [MLESAC_mask, MLESAC_model] = MLESAC(Data, @TLS, g_nSampLen, @Dist, g_nIter, 1, 3);
    DrawLine(MLESAC_model, Data, 'Tag', 'MLESAC', 'Color', [1, 0.5, 0.5]);
    drawnow;    
    DIST = Dist(MLESAC_model, Data(:, find(MLESAC_mask)));
    fprintf('MLESAC MODEL:\n inliers = %d\n', sum(MLESAC_mask));
    fprintf(' mean = %f\n abs_mean = %f\n', mean(DIST), mean(abs(DIST)));
    tmp_mask = MLESAC_mask * 2 + Stat.gen_inliers;
    fprintf(' gen_in_miss = %d, gen_out_miss = %d\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));
    tmp_mask = MLESAC_mask * 2 + Stat.res_inliers;
    fprintf(' res_in_miss = %d, res_out_miss = %d\n\n', length(find(tmp_mask == 1)), length(find(tmp_mask == 2)));

    
    
    hold off;
    zoom on;

    ModelsName = get(sort(get(gca, 'Children')), 'Tag');
    legend(ModelsName, 0);

    for k = 1:length(ModelsName)
        uimenu(vmenu, 'Label', ModelsName{k}, 'Checked', 'on', 'Callback', 'line_ex(''Axes_Draw'',gcbo)');
    end

    % remove all function's paths
    for i = 1:length(ExPath)
        rmpath(ExPath{i});
    end
else
    feval(varargin{:});
end
return;

function LegendTurn(h)
    if strcmpi(get(h, 'Checked'), 'on') == 1
        set(h, 'Checked', 'off');
        legend('hide');
    else
        set(h, 'Checked', 'on');
        legend('show');
    end
    
return;

function Axes_Draw(h)
    
    Data = findobj('Parent', gca, 'Tag', get(h, 'Label'));
    
    if strcmpi(get(h, 'Checked'), 'on') == 1
        set(h, 'Checked', 'off');
        set(Data, 'Visible', 'off')
    else
        set(h, 'Checked', 'on');
        set(Data, 'Visible', 'on')
    end
    
return;