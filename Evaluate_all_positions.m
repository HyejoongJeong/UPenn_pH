% Combine data for all positions and plot together

% manually load data_save.mat file for the experiment you want to evaluate

first_pos = 1; % input index of first position with data
last_pos = 239; % input index of last position with data

% combining each cell entry (position) into single matrix of all droplets
% Note: we are using the mean cy5/cy3 ratio for the pH calculation but we
% can switch to using median by changing the variable below.
combined_cy5_cy3_ratio_mean_store = cat(2, cy5_cy3_ratio_mean_store{[first_pos:last_pos]});
combined_pH_calc = real(6.4-log(((1./combined_cy5_cy3_ratio_mean_store-0.15004)./(1.59979-1./combined_cy5_cy3_ratio_mean_store)).*6.45177));

% mapping single/empty cell indices to new combined matrix 
empty_binary = cell(1,last_pos);
for i = 1:length(empty_droplets)
    empty_binary{i} = zeros(1,size(cy5_cy3_ratio_mean_store{i}, 2));
    empty = empty_droplets{i}(empty_droplets{i} <= size(cy5_cy3_ratio_mean_store{i}, 2));
    empty_binary{i}(empty) = 1; % enter 1 at index of empty droplets (only for indices for which cy5/cy3 ratio is tracked)
end

single_binary = cell(1,last_pos);
for i = 1:length(single_cell_droplets)
    single_binary{i} = zeros(1,size(cy5_cy3_ratio_mean_store{i}, 2));
    single = single_cell_droplets{i}(single_cell_droplets{i} <= size(cy5_cy3_ratio_mean_store{i}, 2));
    single_binary{i}(single) = 1; % enter 1 at index of single droplets (only for indices for which cy5/cy3 ratio is tracked)
end

combined_empty_binary = cat(2, empty_binary{[first_pos:last_pos]});
combined_single_binary = cat(2, single_binary{[first_pos:last_pos]});

combined_empty_list = find(combined_empty_binary == 1);
combined_single_list = find(combined_single_binary==1);

% Generate plots of pH for single cell and empty droplets

% empty droplet pH plots
figure(1); plot(combined_pH_calc(:,combined_empty_list))
title("pH - Empty Droplets")
ylim([0 9])


% single cell droplet pH plots
figure(2); plot(combined_pH_calc(:,combined_single_list))
title("pH - Single Cell Droplets")
ylim([0 9])


% Calculate start-final pH

delta_pH_all = combined_pH_calc(1,:) - combined_pH_calc(end,:); 
delta_ph_empty = delta_pH_all(combined_empty_list);
delta_ph_single = delta_pH_all(combined_single_list);

% Box plot of delta_pH for empty and single cell droplets

pH_compare = nan(max([length(delta_ph_single) length(delta_ph_empty)]),2);
pH_compare(1:length(delta_ph_empty),1) = delta_ph_empty';
pH_compare(1:length(delta_ph_single),2) = delta_ph_single';

figure(3); 
boxplot(pH_compare,'Labels',{'empty','single cell'})
title('initial-final pH')
