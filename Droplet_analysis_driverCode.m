%% Driver code for droplet analysis
% Emilia Leyes Porello
%This code will call the function "fnx_droplet_analysis" using the image
%data loaded here. It can iterate over multiple videos/positions.

clear all

%% load initial variables

for vid = 62 % positions to iterate over
close all
 export = 'D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\Results\';
 results_save_name = 'pH_2-NBDG_3hr';  

%export = 'G:\Shared drives\Image analysis\Results\';
%results_save_name = '2023-02-20 2-NBDG';
    
    if vid < 100
        pos_3 = sprintf('%03d', vid);
        pos_4 = pos_3;
    else
        pos_3 = sprintf('%03d', vid);
        pos_4 = sprintf('%04d', vid);
    end

disp(['position: ', pos_3])

% load the particle analysis result (from FIJI)
 traj = readmatrix(['D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\',results_save_name,'\GFP_Trajectories\Position ', pos_3, '.csv']); % adapt name to call correct video trajectory

% load images 
total_pos = 240; % number of positions for this sample
droplet_diameter = 76; %87 % estimate diameter (pixels) of droplets conservatively -- can use imageJ measuring tool
path_Cy3 = ['D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\',results_save_name,'\SNARF-4F-Cy3\Position ', pos_3,'\SNARF-4F-Cy3_SNARF-4F-Cy3_10msec_image_'];
path_Cy5 = ['D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\',results_save_name,'\SNARF-4F-Cy5\Position ', pos_3,'\SNARF-4F-Cy5_SNARF-4F-Cy5_10msec_image_'];
path_GFP = ['D:\Upenn_since_122822\Microscope_data\Big data\2023-02-20_n_11_cell_damaged_2-NBDG\',results_save_name,'\GFP\Position ', pos_3,'\GFP_GFP_10msec_image_'];
%% call function
fxn_droplet_analysis_10_12_23_update(traj, path_Cy3, path_Cy5, path_GFP, export, pos_3, results_save_name, total_pos, droplet_diameter)

end
