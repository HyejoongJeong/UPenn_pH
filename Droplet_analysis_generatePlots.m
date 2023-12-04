%% Generate figures and plots from stored analysis data -- Bleach correction for cy3 and cy5 channels built in here
% Emilia Leyes Porello

clear all
close all

% load data - user input required
import_folder = 'D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\Results\pH_2-NBDG_3hr'; % copy path of results folder you want to analyze!
load(strcat(import_folder,'\data_save.mat')) % load storage variables

for pos = 62%10:239 % USER INPUT - which positions to generate & save plots for
close all

pos_str = sprintf('%03d', pos); % string of "position" value
disp(strcat("Progress - position: ", pos_str))

% open droplet labeled image
% dropletIDfig = imread(strcat(import_folder,'\Droplet ID\', pos_str,'.png'));
% figure(1); imshow(dropletIDfig);

% referencing data cell corresponding to position, pos
idx = single_cell_droplets{pos}; % indices of droplets
idx2 = other_droplets{pos}; % indices of all other droplets
idx_empty = empty_droplets{pos}; % indices of zero cell droplets

% number of frames for movie "pos"
time = time_store{pos};

%% Manually excluding droplets
exclude = []; % input list of droplets to exclude - leave empty if none

if isempty(exclude) == 0 % delete droplets from list for plotting
    idx(ismember(idx, exclude)==1) = [];
    idx_empty(ismember(idx_empty, exclude) == 1) = [];
end

%% Plotting cell gfp (raw and normalized) and Cy5/cy3 ratio -- overlaid on same plot

% plot normalized cell gfp for single-cell droplets

if isempty(idx)==0
    figure(2);
    for i = idx
        j = find(idx==i);
        plot(smoothdata(gfp_cell_normalized_store{pos}(:,j),'omitnan'),'k-'); hold on; %smooth data, omitting NaN values
        xlim([1 time]);
%         ylim([1.15 1.5])
%         I = ~isnan(gfp_cell_normalized_store{pos}(:,j)); % identify frame index where gfp signal exists - for plot connecting
%         frac = sum(I)/time; % fraction of frames with GFP signal detected
%         if frac > 0.4 
%             fill = fillmissing(gfp_cell_normalized_store{pos}(:,j),'movmean',20); % fill missing values to enable proper smoothing 
% %             plot(smooth(find(I==1),gfp_cell_normalized_store{pos}(I,j)),'k-'); hold on;
%             plot(smooth(fill),'k-'); hold on; %smooth "filled" data

%         end
    end
    hold off; title(strcat('Single-cell GFP intensity (normalized) - position: ', pos_str));
    saveas(gcf,[import_folder,'\Plots\GFP\',pos_str,' GFP traces.png']);
end

% plot cy5/cy3 ratio for single-cell droplets
if isempty(idx)==0
    figure(3);
    for i = idx
        if i<= length(cy5_cy3_ratio_median_store{pos})
            j = find(idx==i);
            plot(smoothdata(cy5_cy3_ratio_median_store{pos}(:,i)),'g-'); plot(smooth(cy5_cy3_ratio_mean_store{pos}(:,i)),'b-'); hold on;
            xlim([1 time])
            ylim([0.75 1])
        end
    end
    legend('mean signals','median signals');
    hold off; title(strcat('Single-cell Cy5/Cy3 ratio - position: ', pos_str));
    saveas(gcf,[import_folder,'\Plots\cy5_cy3\',pos_str,' single cell - cy5_cy3 trace.png']);
end

% plot cy5/cy3 ratio for empty droplets
if isempty(idx_empty)==0
    figure(4); 
    for i = idx_empty
        if i<= length(cy5_cy3_ratio_median_store{pos})
            j = find(idx==i);
            plot(smoothdata(cy5_cy3_ratio_median_store{pos}(:,i)),'g-'); plot(smooth(cy5_cy3_ratio_mean_store{pos}(:,i)),'b-'); hold on;
            xlim([1 time])
            ylim([0.75 1])
        end
    end
    legend('mean signals','median signals');
    hold off; title(strcat('Empty cell Cy5/Cy3 ratio - position: ', pos_str));
    saveas(gcf,[import_folder,'\Plots\cy5_cy3\',pos_str,' empty - cy5_cy3 trace.png']);
end

%% Correcting for bleaching in cy3 and cy5 channels
empty_cy3_mean_list = cy3_mean_store{pos}(:,idx_empty);
mean_empty_cy3_mean_list = mean(empty_cy3_mean_list, 2);

empty_cy5_mean_list = cy5_mean_store{pos}(:,idx_empty);
mean_empty_cy5_mean_list = mean(empty_cy5_mean_list, 2);

% calculate percent drop in intensity of empty cell cy3 and cy5 signals
% over time (i.e., intensity of each time point compared to time 1)

bleaching_coeff_cy3_mean = ones(1,length(mean_empty_cy3_mean_list));
bleaching_coeff_cy5_mean = ones(1,length(mean_empty_cy3_mean_list));

for i = 2:length(mean_empty_cy3_mean_list)
    bleaching_coeff_cy3_mean(i) = ((mean_empty_cy3_mean_list(i)/mean_empty_cy3_mean_list(1)));
    bleaching_coeff_cy5_mean(i) = ((mean_empty_cy5_mean_list(i)/mean_empty_cy5_mean_list(1)));
end

figure(5); plot(bleaching_coeff_cy3_mean,'b','LineWidth',2), hold on, plot(bleaching_coeff_cy5_mean,'r','LineWidth',2), hold off;
legend('cy3', 'cy5', 'Location', 'northeast')
xlabel('time frames')
ylabel('fraction of bleaching over time')

% apply correction to all non-zero cell data
for j = 1:size(cy3_mean_store{pos},2)
    cy3_mean_store_corr{pos}(:,j) = cy3_mean_store{pos}(:,j).*(1+(1-bleaching_coeff_cy3_mean'));
    cy5_mean_store_corr{pos}(:,j) = cy5_mean_store{pos}(:,j).*(1+(1-bleaching_coeff_cy5_mean'));
end

cy5_cy3_ratio_mean_corrected{pos} = cy5_mean_store_corr{pos}./cy3_mean_store_corr{pos};

% plot CORRECTED cy5/cy3 ratio for single-cell droplets
if isempty(idx)==0
    figure(6);
    for i = idx
        if i<= length(cy5_cy3_ratio_mean_corrected{pos})
            j = find(idx==i);
            plot(smoothdata(cy5_cy3_ratio_mean_corrected{pos}(:,i)),'b-'); hold on;
            xlim([1 time])
            ylim([0.75 1])
        end
    end
    legend('mean signals');
    mean_single_cy5cy3_mean_bleachCorr = nanmean(cy5_cy3_ratio_mean_corrected{pos}(:,idx),2);
    plot(mean_single_cy5cy3_mean_bleachCorr,'-*r', 'LineWidth',2)
    hold off; title(strcat('Bleaching corrected - Single-cell Cy5/Cy3 ratio - position: ', pos_str));
    saveas(gcf,[import_folder,'\Plots\cy5_cy3_bleachCorrection\',pos_str,' single cell - cy5_cy3 trace.png']);
end

% plot CORRECTED cy5/cy3 ratio for empty droplets
if isempty(idx_empty)==0
    figure(7); 
    for i = idx_empty
        if i<= length(cy5_cy3_ratio_mean_corrected{pos})
            j = find(idx==i);
            plot(smoothdata(cy5_cy3_ratio_mean_corrected{pos}(:,i)),'b-'); hold on; 
            xlim([1 time])
            ylim([0.75 1])
        end
    end
    legend("mean signals")
    mean_empty_cy5cy3_mean_bleachCorr = nanmean(cy5_cy3_ratio_mean_corrected{pos}(:,idx_empty),2);
    plot(mean_empty_cy5cy3_mean_bleachCorr,'-*r', 'LineWidth',2)
    title(strcat('Bleaching corrected - Empty cell Cy5/Cy3 ratio - position: ', pos_str));
    saveas(gcf,[import_folder,'\Plots\cy5_cy3_bleachCorrection\',pos_str,' empty - cy5_cy3 trace.png']);
end
end

save(strcat(import_folder,'\data_save.mat'), "cy5_cy3_ratio_mean_corrected", "mean_empty_cy5cy3_mean_bleachCorr","mean_single_cy5cy3_mean_bleachCorr", "-append")