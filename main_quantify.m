%% Main Virtual Clinical Trial Analysis and Metric Validation
% Script Purpose: 
% This script processes results from Virtual Clinical Trials (VCT) to evaluate 
% the predictive power of three immunity metrics: msImmunity (multiscale), 
% tsImmunity (tissue-scale), and molsImmunity (molecular-scale).
%
% Workflow:
% 1. Load therapy regimens and population-level simulation data for four cohorts.
% 2. Iteratively process 1,000 virtual patients (VPs) per cohort to calculate 
%    clinical response ($Percentage change$) and immunity scores.
% 3. Visualize the correlation between immunity metrics and tumor regression.
% 4. Perform ROC curve analysis to calculate AUC, quantifying the accuracy 
%    of metrics in predicting treatment responders.
%
% Data Dependencies: 
% Requires.mat files in 'VCT_Results' subfolders and therapy matrix files.

 clear
 Score ={};
        Therapy_martix = {};
  
        Percent_change = [];
        subfolder{4} = 'VP_Xelox_IR_antipd1\'; 
        Therapy_martix{4} = cell2mat(struct2cell(load('Therapy_xelox_antiPD1_IR.mat')));
        
        subfolder{3} = 'VP_Xelox_antipd1\';  
        Therapy_martix{3} = cell2mat(struct2cell(load('Therapy_xelox_antiPD1.mat')));
    
        subfolder{2} = 'VP_Xelox_IR\';
        Therapy_martix{2} = cell2mat(struct2cell(load('Therapy_xelox_IR.mat')));
    
        subfolder{1} = 'VP_Xelox\';
        Therapy_martix{1} = cell2mat(struct2cell(load('Therapy_xelox.mat')));

        print = 0;
  %% Data processing

 for score_index = 1:3
    for lll = 1:4
        kk =lll;
            f_path = fullfile('VCT_Results',subfolder{kk});
            S{kk} = cell2mat(struct2cell(load([f_path,'VP_LHS_Matrix.mat'])));
            Par_all = {};
            Score_matrix = [];
            Par_matrix = [];
            sample = 1000; 
            for i= 1:sample
                file=[f_path,'VP',num2str(i),'.mat'];
                A=cell2mat(struct2cell(load(file)));
                Par_all{i} = A.params;
                P_change(i) = (A.Vtumor_list(end)-A.Vtumor_list(1))/A.Vtumor_list(1);
                [Score_matrix(i,:),Par_matrix(i,:)] = Quantify(Par_all{i},Therapy_martix{kk});
            end
            [~, index] = sort(Score_matrix(:,score_index), 1, 'ascend');
            y_sort = [];
            score_sort = [];
            score_final = Score_matrix(:,score_index);
        switch score_index
            case 1
            msImmunity{lll} = score_final(index);
            P_final_ms{lll} = P_change(index);
            Index{lll} =index;
            case 2
            tsImmunity{lll} = score_final(index);
            P_final_ts{lll} = P_change(index);
            Index_tsm{lll} =index;
             case 3
            molsImmunity{lll} = score_final(index);
            P_final_mols{lll} = P_change(index);
            Index_mols{lll} =index;
        end  
        
         
    end

 end
 %%
% clear
colors_dot = lines(4)* 0.6 + 0.4;
colors_line = lines(4);
text1 = 'Xelox';
text2 = 'Xelox+RT';
text3 = 'Xelox+Anti-PD1';
text4 = 'Xelox+RT+Anti-PD1';
number = ["(A)","(B)","(C)","(D)"];

%% msImmunity vs Percentage change of tumor
Score_select = msImmunity;
P_select = P_final_ms;
for k =1:4
    figure(1)
    subplot(2,2,k)
    x = Score_select{k};
    y = P_select{k}';
        color_PD = [1.0, 0.96, 0.93];   
        color_CR = [0.94, 0.97, 1.0]; 
        alpha_val = 0.7; 
        x_lims = [0 ceil(max(x)*1.2)];
    patch([x_lims fliplr(x_lims)], [-30 -30 100 100], color_PD, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    patch([x_lims fliplr(x_lims)], [-100 -100 -30 -30], color_CR, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    yline(-70, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);
    
    yline(-30, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'Alpha', 0.8);
    yline(20, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);

    text(0.97, 0.8, 'PD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.45, 'SD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.25, 'PR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.08, 'CR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(1.01, 0.68, 'Non-Responders','Units', 'normalized','Rotation', 270,'Color', [0.6 0.3 0.3],'FontSize', 13, ...
        'FontWeight', 'bold','HorizontalAlignment', 'center','VerticalAlignment', 'bottom','Clipping', 'off');            
    text(1.01, 0.18, 'Responders', 'Units', 'normalized', 'Rotation', 270, 'Color', [0.2 0.4 0.7], 'FontSize', 13, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

    d(k)=scatter(x, y*100, 50, colors_dot(k,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none'); 
    hold on

    num_bins = 20;
    edges = quantile(x, linspace(0, 1, num_bins + 1));
    [~, bin_idx] = histc(x, edges);
    x_mean = accumarray(bin_idx(bin_idx > 0), x(bin_idx > 0), [], @mean);
    y_mean = accumarray(bin_idx(bin_idx > 0), y(bin_idx > 0), [], @mean);
    [x_mean, sort_idx] = sort(x_mean);
    y_mean = y_mean(sort_idx);
    l(k)= plot(x_mean,y_mean*100,LineWidth=3,Color=colors_line(k,:));
    hold on
    switch k
        case 1
            text_use = text1;
        case 2 
            text_use = text2;
        case 3
            text_use = text3;
        case 4
            text_use = text4;
    end
    legend([l(k),d(k)],{'Segmented average curve', 'Virtual patient data'},'Location','northeast', ...
           'Box','off','LineWidth',0.5,'FontSize',15,'Interpreter','latex'); 
    title(text_use,'FontSize',20, 'FontWeight', 'bold');

    ax = gca;
    ax.XAxis.FontSize = 16;
    ax.YAxis.FontSize = 16; 

    xl = xlabel("msImmunity",'FontSize',16, 'FontWeight', 'bold',Interpreter='tex');
    yl = ylabel('Percentage change of tumor (%)','FontSize',14, 'FontWeight', 'bold');
    set(gca,  'FontName', 'Arial', 'Box', 'on','LineWidth', 0.5);
    grid off; 
    
    xlim([0 ceil(min(10,max(x)))]);
end

figure(1)
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 36, 21]); 

%% tsImmunity vs Percentage change of tumor
Score_select = tsImmunity;
P_select = P_final_ts;
for k =1:4
    figure(2)
    subplot(2,2,k)
    x = Score_select{k};
    y = P_select{k}';
        color_PD = [1.0, 0.96, 0.93];   
        color_CR = [0.94, 0.97, 1.0]; 
        alpha_val = 0.7; 
        x_lims = [0 ceil(max(x)*1.2)];
    patch([x_lims fliplr(x_lims)], [-30 -30 100 100], color_PD, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    patch([x_lims fliplr(x_lims)], [-100 -100 -30 -30], color_CR, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    yline(-70, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);
    
    yline(-30, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'Alpha', 0.8);
    yline(20, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);

    text(0.97, 0.8, 'PD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.45, 'SD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.25, 'PR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.08, 'CR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(1.01, 0.68, 'Non-Responders','Units', 'normalized','Rotation', 270,'Color', [0.6 0.3 0.3],'FontSize', 13, ...
        'FontWeight', 'bold','HorizontalAlignment', 'center','VerticalAlignment', 'bottom','Clipping', 'off');            
    text(1.01, 0.18, 'Responders', 'Units', 'normalized', 'Rotation', 270, 'Color', [0.2 0.4 0.7], 'FontSize', 13, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

    d(k)=scatter(x, y*100, 50, colors_dot(k,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none'); 
    hold on

    num_bins = 20;
    edges = quantile(x, linspace(0, 1, num_bins + 1));
    [~, bin_idx] = histc(x, edges);
    x_mean = accumarray(bin_idx(bin_idx > 0), x(bin_idx > 0), [], @mean);
    y_mean = accumarray(bin_idx(bin_idx > 0), y(bin_idx > 0), [], @mean);
    [x_mean, sort_idx] = sort(x_mean);
    y_mean = y_mean(sort_idx);
    l(k)= plot(x_mean,y_mean*100,LineWidth=3,Color=colors_line(k,:));
    hold on
    switch k
        case 1
            text_use = text1;
        case 2 
            text_use = text2;
        case 3
            text_use = text3;
        case 4
            text_use = text4;
    end
    legend([l(k),d(k)],{'Segmented average curve', 'Virtual patient data'},'Location','northeast', ...
           'Box','off','LineWidth',0.5,'FontSize',15,'Interpreter','latex'); 
    title(text_use,'FontSize',20, 'FontWeight', 'bold');

    ax = gca;
    ax.XAxis.FontSize = 16;
    ax.YAxis.FontSize = 16; 

    xl = xlabel("tsImmunity",'FontSize',16, 'FontWeight', 'bold',Interpreter='tex');
    yl = ylabel('Percentage change of tumor (%)','FontSize',14, 'FontWeight', 'bold');
    set(gca,  'FontName', 'Arial', 'Box', 'on','LineWidth', 0.5);
    grid off; 
    xlim([0 ceil(max(x))*1.15]);
end

figure(2)
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 36, 21]); 


%% molsImmunity vs Percentage change of tumor
Score_select = molsImmunity;
P_select = P_final_mols;
for k =1:4
    figure(3)
    subplot(2,2,k)
    x = Score_select{k};
    y = P_select{k}';
        color_PD = [1.0, 0.96, 0.93];   
        color_CR = [0.94, 0.97, 1.0]; 
        alpha_val = 0.7; 
         x_lims = [0 (max(x)*1.2)];
    patch([x_lims fliplr(x_lims)], [-30 -30 100 100], color_PD, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    patch([x_lims fliplr(x_lims)], [-100 -100 -30 -30], color_CR, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    yline(-70, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);
    
    yline(-30, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'Alpha', 0.8);
    yline(20, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'Alpha', 0.6);

    text(0.97, 0.8, 'PD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.45, 'SD', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.25, 'PR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(0.97, 0.08, 'CR', 'Units', 'normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(1.01, 0.68, 'Non-Responders','Units', 'normalized','Rotation', 270,'Color', [0.6 0.3 0.3],'FontSize', 13, ...
        'FontWeight', 'bold','HorizontalAlignment', 'center','VerticalAlignment', 'bottom','Clipping', 'off');            
    text(1.01, 0.18, 'Responders', 'Units', 'normalized', 'Rotation', 270, 'Color', [0.2 0.4 0.7], 'FontSize', 13, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

    d(k)=scatter(x, y*100, 50, colors_dot(k,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none'); 
    hold on

    num_bins = 20;
    edges = quantile(x, linspace(0, 1, num_bins + 1));
    [~, bin_idx] = histc(x, edges);
    x_mean = accumarray(bin_idx(bin_idx > 0), x(bin_idx > 0), [], @mean);
    y_mean = accumarray(bin_idx(bin_idx > 0), y(bin_idx > 0), [], @mean);
    [x_mean, sort_idx] = sort(x_mean);
    y_mean = y_mean(sort_idx);
    l(k)= plot(x_mean,y_mean*100,LineWidth=3,Color=colors_line(k,:));
    hold on
    switch k
        case 1
            text_use = text1;
        case 2 
            text_use = text2;
        case 3
            text_use = text3;
        case 4
            text_use = text4;
    end
    legend([l(k),d(k)],{'Segmented average curve', 'Virtual patient data'},'Location','northeast', ...
           'Box','off','LineWidth',0.5,'FontSize',15,'Interpreter','latex'); 
    title(text_use,'FontSize',20, 'FontWeight', 'bold');

    ax = gca;
    ax.XAxis.FontSize = 16;
    ax.YAxis.FontSize = 16; 

    xl = xlabel("molsImmunity",'FontSize',16, 'FontWeight', 'bold',Interpreter='tex');
    yl = ylabel('Percentage change of tumor (%)','FontSize',14, 'FontWeight', 'bold');
    set(gca,  'FontName', 'Arial', 'Box', 'on','LineWidth', 0.5);
    grid off; 
    xlim([0 max(x)*1.15]);
end

figure(3)
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 36, 21]); 

%% AUC calculation ROC curve

response_thresh = -0.3; 
auc = zeros(3,4);  
figure(4);

Score_sets = {msImmunity, tsImmunity, molsImmunity};
P_sets     = {P_final_ms, P_final_ts, P_final_mols};
colors_set = lines(7);
lineStyles = {'-', '--', '-.'};   
Markers = {'^', 's', 'o'};
methodNames = {'msImmunity',  'tsImmunity', 'molsImmunity'};

for k = 1:4   
    subplot(2,2,k); hold on;
    
    for m = 1:3  
        Score_select = Score_sets{m};
        P_select    = P_sets{m};

        score_final = Score_select{k};
        P_change    = P_select{k};

        
        is_responder = P_change <= response_thresh;

     
        [X, Y, ~, AUC] = perfcurve(is_responder, score_final, true);
        auc(m,k) = AUC;

       
        h(m) = plot(X, Y, ...
             'LineWidth', 2.5, ...
             'Color', colors_set(k,:) , ...
             'LineStyle', lineStyles{m},...
             'Marker',Markers{m},...
             'MarkerSize', 9, ...
            'MarkerIndices', 1:110:length(X));
    end

    
    h_diag = plot([0 1], [0 1], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

 
    switch k
        case 1
            text_use = text1;
        case 2 
            text_use = text2;
        case 3
            text_use = text3;
        case 4
            text_use = text4;
    end
    title(text_use,'FontSize',20, 'FontWeight', 'bold');

    legend([h(1),h(2),h(3)], ...
      { sprintf('%s\nAUC=%.3f', methodNames{1}, auc(1,k)), ...
    sprintf('%s\nAUC=%.3f', methodNames{2}, auc(2,k)), ...
    sprintf('%s\nAUC=%.3f', methodNames{3}, auc(3,k))}, ...
       'Location','southeast', ...
       'Box','off','LineWidth',0.5,'FontSize',14,'Interpreter','tex');

    
    ax = gca; 
    ax.XAxis.FontSize = 16; % 
    ax.YAxis.FontSize = 16; % 
    xlabel('False Positive Rate','FontSize',16, 'FontWeight', 'bold');
    ylabel('True Positive Rate','FontSize',16, 'FontWeight', 'bold');
    set(gca,  'FontName', 'Arial', 'Box', 'on', ...
        'LineWidth', 0.5, 'TickDir','out');
    grid off;
    axis([-0.02 1.02 -0.02 1.02]);
grid off;
end

set(gcf, 'Units', 'centimeters', 'Position', [1, 1, 40, 26]);


