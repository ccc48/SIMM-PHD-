function [P_positive] = makeMatrixPositiveDefinite(P, epsilon)
% makeMatrixPositiveDefinite - 将非正定矩阵转换为正定矩阵
% 输入:
%   P - 待转换的非正定矩阵
%   epsilon - 替换负特征值的小正数，默认1e-6
% 输出:
%   P_positive - 转换后的正定矩阵

    % 设置默认的epsilon值
    if nargin < 2
        epsilon = 1e-6;
    end
    
    % 检查矩阵是否包含NaN或Inf
    if any(isnan(P(:))) || any(isinf(P(:)))
        % fprintf('警告: 输入矩阵包含NaN或Inf值，正在尝试修复...\n');
        
        % 替换NaN为0
        P(isnan(P)) = 0;
        
        % 替换Inf为极大值
        P(isinf(P)) = realmax;
        
        % 确保矩阵对称
        P = (P + P') / 2;
    end
    
    % 特征值分解
    try
        [V, D] = eig(P);
    catch ME
        % fprintf('错误: 特征值分解失败 (%s)，尝试正则化矩阵...\n', ME.message);
        
        % 添加小的正则化项
        reg_term = epsilon * eye(size(P));
        P = P + reg_term;
        
        % 重试特征值分解
        try
            [V, D] = eig(P);
        catch ME2
            error('无法对矩阵进行特征值分解: %s', ME2.message);
        end
    end
    
    % 提取特征值
    eigenvalues = diag(D);
    
    % 显示原始特征值
    % fprintf('原始矩阵的特征值:\n');
    for i = 1:length(eigenvalues)
        if isnan(eigenvalues(i)) || isinf(eigenvalues(i))
            % fprintf('λ%d = 无效值 (替换为epsilon)\n', i);
            eigenvalues(i) = epsilon;
        else
            % fprintf('λ%d = %.6e\n', i, eigenvalues(i));
        end
    end
    
    % 修正负特征值和无效值
    for i = 1:length(eigenvalues)
        if eigenvalues(i) <= 0 || isnan(eigenvalues(i)) || isinf(eigenvalues(i))
            eigenvalues(i) = epsilon;
        end
    end
    
    % 重构矩阵
    D_new = diag(eigenvalues);
    P_positive = V * D_new * V';
    
    % 确保结果矩阵对称
    P_positive = (P_positive + P_positive') / 2;
    
    % 验证特征值
    new_eigenvalues = eig(P_positive);
    
    % 显示修正后的特征值
    % fprintf('\n修正后矩阵的特征值:\n');
    % for i = 1:length(new_eigenvalues)
    %     fprintf('λ%d = %.6e\n', i, new_eigenvalues(i));
    % end
    
    % 检查正定性
    % if all(new_eigenvalues > 0)
    %     fprintf('\n转换后的矩阵是正定的!\n');
    % else
    %     fprintf('\n转换后的矩阵仍然不是正定的，请调整epsilon值!\n');
    % end
end    


%下面这个不能解决出现Nan特征值的情况
% function [P_positive] = makeMatrixPositiveDefinite(P, epsilon)
% % makeMatrixPositiveDefinite - 将非正定矩阵转换为正定矩阵
% % 输入:
% %   P - 待转换的非正定矩阵
% %   epsilon - 替换负特征值的小正数，默认1e-6
% % 输出:
% %   P_positive - 转换后的正定矩阵
% 
%     % 设置默认的epsilon值
%     if nargin < 2
%         epsilon = 1e-6;
%     end
% 
%     % 特征值分解
%     [V, D] = eig(P);
% 
%     % 提取特征值
%     eigenvalues = diag(D);
% 
%     % 显示原始特征值
%     % fprintf('原始矩阵的特征值:\n');
%     % for i = 1:length(eigenvalues)
%     %     fprintf('λ%d = %.6e\n', i, eigenvalues(i));
%     % end
% 
%     % 修正负特征值
%     for i = 1:length(eigenvalues)
%         if eigenvalues(i) <= 0
%             eigenvalues(i) = epsilon;
%         end
%     end
% 
%     % 重构矩阵
%     D_new = diag(eigenvalues);
%     P_positive = V * D_new * V';
% 
%     % 验证特征值
%     new_eigenvalues = eig(P_positive);
% 
%     % 显示修正后的特征值
%     % fprintf('\n修正后矩阵的特征值:\n');
%     % for i = 1:length(new_eigenvalues)
%     %     fprintf('λ%d = %.6e\n', i, new_eigenvalues(i));
%     % end
% 
%     % 检查正定性
%     % if all(new_eigenvalues > 0)
%     %     fprintf('\n转换后的矩阵是正定的!\n');
%     % else
%     %     fprintf('\n转换后的矩阵仍然不是正定的，请调整epsilon值!\n');
%     % end
% end    