function [] ...
    = fxn_droplet_analysis(traj, path_Cy3, path_Cy5, path_GFP, export, pos_3, results_save_name, total_pos, droplet_diameter)
% Process gfp, cy3, cy5 channels - 
% Output of this function is max, mean, and median droplet intensities
% (single-cell droplets only) as well as gfp intensity of single-cells in
% these same droplets
% Emilia Leyes Porello

%% initial variables

if exist(strcat(export,results_save_name,'\data_save','.mat')) == 2 % if results tracker already exists, load it
    load(strcat(export,results_save_name,'\data_save','.mat')) % load storage variables
end

time = max(traj(:,3))+1; % number of frames in the movie

x = zeros(time,traj(end,2)); % initialize the x,y coordinate (single cells)
y = zeros(time,traj(end,2));

pos = pos_3;
movie_idx=str2num(pos);

%% find droplets with one cell using particle analysis 

% organize the particle trajectories
for i=1:traj(end,2) % number of captured particles
    tmp = find(traj(:,2) == i); % for each particle index
    x(traj(tmp,3)+1,i) = traj(tmp,4); % x coordinates of each particle over time 
    y(traj(tmp,3)+1,i) = traj(tmp,5); % y 
    tmp2 = find(x(:,i) == 0);
    x(tmp2,i) = nan; % if the particle was not detected in a given frame, make "NaN"
    y(tmp2,i) = nan;
end
    
for i=1:time % time
     % load images
    I_gfp = imread(sprintf('%s%05d.tif',path_GFP,(movie_idx+1)+240*(i-1))); % load GFP images
    
    % segmentation of droplets - GFP
    mask = I_gfp; % use GFP images as a mask (can change to others)
%     mask = imadjust(mask,[0.001 0.003]); % ADDED THIS LINE FOR OLIGOMYCIN DATA - NOT NECESSARY FOR OTHERS 
    mask = imadjust(mask,[0.001 0.005]); % settings for 2023-02-20 2-NBDG data 
    thr = graythresh(mask); % graythresh
    if thr <=0
        bw = im2bw(mask,0.005); % if graythresh is too low, use manual value. Can change to other values
    else
        bw = im2bw(mask,thr); % use threshold to convert the image to a bw image
    end
    bw = imfill(bw,'holes'); 
    bw = bwareaopen(bw,600); % remove some small pixels 
    bw = imclearborder(bw,8); % remove any edge objects
    figure(10); imshow(bw),title(strcat('GFP - frame: ', int2str(i)));

    % region props - create a circle for each droplet, and resize it to 50%
    rpo = regionprops(bw,'Centroid','MajorAxisLength','MinorAxisLength');
    centers_temp = cat(1,rpo.Centroid);
    major_temp = cat(1,rpo.MajorAxisLength);
    minor_temp = cat(1,rpo.MinorAxisLength);
    
    % map objects identified by regionprops to correct droplet ID
    centers_idx = 1:length(rpo);
    if i == 1
       centroids_set_pos = centers_temp; % store position of ref wells (frame 1)
       centers=centers_temp;
       major = major_temp;
       minor = minor_temp;
    else
        for ref = 1:length(centroids_set_pos) % number of wells
            distance = sqrt((centers_temp(:,1)-centroids_set_pos(ref,1)).^2+(centers_temp(:,2)-centroids_set_pos(ref,2)).^2);
            [M,I]=min(distance); % value and index of min distance
            if M < 25
                centers(ref,:) = centers_temp(I,:); %re-assigning centroids indexing to match reference indices in frame 1
                major(ref,:) = major_temp(I); 
                minor(ref,:) = minor_temp(I);
                centers_idx(ref) = I;
            end
        end
    end

    droplet_x(i,:) = centers(:,1);
    droplet_y(i,:) = centers(:,2);
%     major = cat(1,rpo.MajorAxisLength);
%     minor = cat(1,rpo.MinorAxisLength);
    diameters = mean([major minor],2);
    radii = diameters/2;
    radii2 = radii*0.55; % can adjust to other values if needed - this is fraction of radius to consider for single-cell indentification
    
%     I3 = mat2gray(I_gfp);% for saving images
%     figure; imshow(I3); hold on; plot(x(i,:),y(i,:),'g.'); 
%     viscircles(centers,radii2); 

% particle's distance from the droplet's center
for j=1:length(centers) % droplet ID in order of objects in reference frame
    dist(j,:) = sqrt((x(i,:)-centroids_set_pos(j,1)).^2+(y(i,:)-centroids_set_pos(j,2)).^2);
    dum = find(dist(j,:) < radii2(j)); % find the particles within the circle (defined above)
%     plot(x(i,dum),y(i,dum),'bo'); % for saving images
    particle{i,j} = dum; % all particles within the x% inner radius
    num(i,j) = length(particle{i,j}); % number of particles in each droplet

end

% hold off;   

%     figure;
%     imshow(I_gfp*20);
%     title({'Droplet index - Frame: ',i})
%     for jjj = 1:length(centers)
%         x_plot = centers(jjj,1);
%         y_plot = centers(jjj,2);
%         text(x_plot-20, y_plot, int2str(jjj),'FontSize',10, 'Color' ,'w'); % annotate each obj with cell id
%         hold on;
%     end

end
    
for ii=1:length(centers)
    num2(ii) = max(num(:,ii)); % find the max number of particles 
end

idx = find(num2 == 1); % droplet index with one cell only

track_single = zeros(size(particle));
for ii = idx % convert "particle{}" data into matrix "track" for idx values
    for k = 1:time
        if isempty(particle{k,ii})==0
            track_single(k,ii) = particle{k,ii};
        end
    end
end

track_empty = zeros(size(particle));
for ii = 1:length(particle) % convert "particle{}" data into matrix "track" for idx values
    for k = 1:time
        if isempty(particle{k,ii})==1
            track_empty(k,ii) = 1; % enter 1 for droplets with zero cells in frame k
        end
    end
end

% eliminate droplets that have two different particles tracked as a single one
for j = 1:length(track_single)
    if mean(nonzeros(track_single(:,j))) ~= max(track_single(:,j))
        idx(idx==j) = []; % remove tracking of cells if two different cells were tracked as the same one in a single-cell droplet
    end
end

idx_empty = [];
% generate list of empty droplets
for j = 1:length(track_empty)
    if sum(track_empty(:,j)) == time
        idx_empty = [idx_empty j];
    end
end

member = ismember(1:length(centers),idx);
idx2 = find(member == 0); % index value of all droplets with zero or 2+ cells
idx_remove = []; % list to fill later for removal of "single cell" indices from idx

droplet_x(:,idx2) = [];
droplet_y(:,idx2) = [];

xf = zeros(time,length(idx)); % x,y coordinates of the droplets
yf = zeros(time,length(idx)); 
for i=1:length(idx) % number of droplets with one cell
    for j=1:time
        dum2 = particle{j,idx(i)};
        if isempty(dum2)==0
            xf(j,i) = x(j,dum2); % x,y coordinates of the cell
            yf(j,i) = y(j,dum2);
        else
            xf(j,i) = nan; % when the cell is not detected, NaN
            yf(j,i) = nan;
        end
    end
end

% % plot
% for i=1:time
%     I_gfp = imread(sprintf('%s%05d.tif',path_GFP,movie_idx+240*(i-1)));
%     I3 = mat2gray(I_gfp);
%     figure(3); imshow(I3); hold on; plot(xf(i,:),yf(i,:),'bo'); hold off;
%     saveas(gcf,sprintf('%s%s%03d.tif',export,'gfp_max/',i));
% end

%% measure the intensity! 
for i=1:time % time
    % load images
    I_cy3 = imread(sprintf('%s%05d.tif',path_Cy3,(movie_idx+1)+total_pos*(i-1)));
    I_cy5 = imread(sprintf('%s%05d.tif',path_Cy5,(movie_idx+1)+total_pos*(i-1)));
    I_gfp = imread(sprintf('%s%05d.tif',path_GFP,(movie_idx+1)+total_pos*(i-1)));
    
    % Segmentation and region props calcuated for each channel separately!   
    
    % CY5
    % segmentation of the droplets - Cy5 - only done in first frame
    if i == 1
        mask = I_cy5;
        mask = imadjust(mask,[0.0002 0.0019]); % ADDED THIS LINE FOR OLIGOMYCIN DATA - NOT NECESSARY FOR OTHERS 
        thr = graythresh(mask);
        if thr <=0
            bw = im2bw(mask,0.002);
        else
            bw = im2bw(mask,thr);
        end
        bw = imfill(bw,'holes');
        bw = bwareaopen(bw,20); % gfp_max
        bw = imclearborder(bw,8); % remove any edge objects
      figure(11); imshow(bw),title(strcat('cy5 - frame: ', int2str(i)));

        % region props - Cy5
        rpo_cy5 = regionprops(bw,I_cy5,'MaxIntensity','MeanIntensity','PixelList','Centroid','MajorAxisLength');

        % matching droplet IDs to frame 1 references - cy5
        centers_temp_cy5 = cat(1,rpo_cy5.Centroid);
        clear idx_cy5
        for ref = 1:length(centroids_set_pos) % number of ref wells
            distance = sqrt((centers_temp_cy5(:,1)-centroids_set_pos(ref,1)).^2+(centers_temp_cy5(:,2)-centroids_set_pos(ref,2)).^2);
            [M,I]=min(distance); % value and index of min distance
            if M < 25
                idx_cy5(ref) = I; % indices of frame i cells in same order as ref frame (1)
            else
                idx2 = [idx2, I];
                if ismember(I, idx)==1
                    idx_remove = [idx_remove find(idx==I)];
                    idx(idx==I)=[];
                end
            end
        end
    end
    
    disp(strcat("evaluating image:" , int2str((movie_idx+1)+total_pos*(i-1))));
        
    % extracting intensities - cy5
    for j=idx_cy5(idx_cy5~=0) % ID of droplets (in correct order)
        % Create a logical image of a circle with specified diameter, center, and image size.
        imageSizeX = size(I_cy5,2);
        imageSizeY = size(I_cy5,1);
        [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
        centerX = rpo_cy5(j).Centroid(1);
        centerY = rpo_cy5(j).Centroid(2);
        radius = droplet_diameter/2;
        circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
        shrink_pixList = regionprops(circlePixels,'PixelList');
        rpo_cy5(j).PixelList = shrink_pixList(1).PixelList;

        for k=1:length(rpo_cy5(j).PixelList)
            A(k) = I_cy5(rpo_cy5(j).PixelList(k,2),rpo_cy5(j).PixelList(k,1));
        end
        cy5_max_all(i,idx_cy5==j) = max(A);%rpo_cy5(j).MaxIntensity;
        cy5_median_all(i,idx_cy5==j) = median(A);
        cy5_mean_all(i,idx_cy5==j) = mean(A);%rpo_cy5(j).MeanIntensity;
        clear A
    end
    
    % CY3
    % segmentation of the droplets - Cy3 - only in frame 1
    if i == 1
        mask = I_cy3;
        mask = imadjust(mask,[0.0002 0.0019]); % ADDED THIS LINE FOR OLIGOMYCIN DATA - NOT NECESSARY FOR OTHERS 
        thr = graythresh(mask);
        if thr <=0
            bw = imbinarize(mask,0.002);
        else
            bw = imbinarize(mask,thr);
        end
        bw = imfill(bw,'holes');
        bw = bwareaopen(bw,20); 
        bw = imclearborder(bw,8); % remove any edge objects
        figure(13); imshow(bw),title(strcat('cy3 - frame: ', int2str(i)));
        
        % region props - Cy3
        rpo_cy3 = regionprops(bw,I_cy3,'MaxIntensity','MeanIntensity','PixelList','Centroid', 'MajorAxisLength');
    
        % matching droplet IDs to frame 1 references - cy3
        centers_temp_cy3 = cat(1,rpo_cy3.Centroid);
        clear idx_cy3
        for ref = 1:length(centroids_set_pos) % number of ref wells
            distance = sqrt((centers_temp_cy3(:,1)-centroids_set_pos(ref,1)).^2+(centers_temp_cy3(:,2)-centroids_set_pos(ref,2)).^2);
            [M,I]=min(distance); % value and index of min distance
            if M < 25
                idx_cy3(ref) = I; % indices of frame i cells in same order as ref frame (1)
            else
                idx2 = [idx2, I];
                if ismember(I, idx)==1
                    idx_remove = [idx_remove find(idx==I)];
                    idx(idx==I)=[];
                end
            end
        end
    end

    
    % extracting intensities - cy3
    for j=idx_cy3(idx_cy3~=0) % number of conserved droplets
        
        imageSizeX = size(I_cy3,2);
        imageSizeY = size(I_cy3,1);
        [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
        centerX = rpo_cy3(j).Centroid(1);
        centerY = rpo_cy3(j).Centroid(2);
        radius = droplet_diameter/2;
        circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
        shrink_pixList = regionprops(circlePixels,'PixelList');
        rpo_cy3(j).PixelList = shrink_pixList(1).PixelList;

        for k=1:length(rpo_cy3(j).PixelList)
            A(k) = I_cy3(rpo_cy3(j).PixelList(k,2),rpo_cy3(j).PixelList(k,1));
        end
        cy3_max_all(i,idx_cy3==j) = max(A); %rpo_cy3(j).MaxIntensity; % max intensity of each droplet
        cy3_median_all(i,idx_cy3==j) = median(A); % median 
        cy3_mean_all(i,idx_cy3==j) = mean(A); % rpo_cy3(j).MeanIntensity; % mean
        clear A
    end
    

    % segmentation of the droplets - GFP - done only in first frame
    if i == 1
        mask = I_gfp;
        mask = imadjust(mask,[0.001 0.005]); % ADDED THIS LINE FOR OLIGOMYCIN DATA - NOT NECESSARY FOR OTHERS 
        thr = graythresh(mask);
        if thr <=0
            bw = im2bw(mask,0.005);
        else
            bw = im2bw(mask,thr);
        end
        bw = imfill(bw,'holes');
        bw = bwareaopen(bw,20);
        bw = imclearborder(bw,8); % remove any edge objects
        figure(14); imshow(bw); title('GFP channel segmentation - first frame only');
        
        % region props - gfp
        rpo_gfp = regionprops(bw,I_gfp,'MaxIntensity','MeanIntensity','PixelList','Centroid');
        
        % matching droplet IDs to frame 1 references
        centers_temp_gfp = cat(1,rpo_gfp.Centroid);
    
        clear idx_gfp
        for ref = 1:length(centroids_set_pos) % number of ref wells
            distance = sqrt((centers_temp_gfp(:,1)-centroids_set_pos(ref,1)).^2+(centers_temp_gfp(:,2)-centroids_set_pos(ref,2)).^2);
            [M,I]=min(distance); % value and index of min distance
            if M < 25
                idx_gfp(ref) = I; % indices of frame i cells in same order as ref frame (1)
            else
                idx2 = [idx2, I];
                if ismember(I, idx)==1
                    idx_remove = [idx_remove find(idx==I)];
                    idx(idx==I)=[];
                end
            end
        end
    end

    % extracting intensities - gfp
    for j=idx_gfp(idx_gfp~=0) % number of droplets
        
        imageSizeX = size(I_gfp,2);
        imageSizeY = size(I_gfp,1);
        [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
        centerX = rpo_gfp(j).Centroid(1);
        centerY = rpo_gfp(j).Centroid(2);
        radius = droplet_diameter/2;
        circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
        shrink_pixList = regionprops(circlePixels,'PixelList');
        rpo_gfp(j).PixelList = shrink_pixList(1).PixelList;
        
        for k=1:length(rpo_gfp(j).PixelList)
            A(k) = I_gfp(rpo_gfp(j).PixelList(k,2),rpo_gfp(j).PixelList(k,1));
        end
        gfp_max_all(i,idx_gfp==j) = max(A); % rpo_gfp(j).MaxIntensity;
        gfp_median_all(i,idx_gfp==j) = median(A);
        gfp_mean_all(i,idx_gfp==j) = mean(A); %rpo_gfp(j).MeanIntensity;
        clear A
    end
end

% extract cell intensity from single-cell droplets
for i=1:time
    I_gfp = imread(sprintf('%s%05d.tif',path_GFP,(movie_idx+1)+total_pos*(i-1)));
    for j=1:length(idx)
        if isnan(xf(i,j)) == 0
%         intensity is calculated by average 9 pixels around the center of the cell
            gg = I_gfp(round(yf(i,j))-1:round(yf(i,j))+1, round(xf(i,j))-1:round(xf(i,j))+1);
            gfp_cell(i,j) = mean(mean(gg)); % cell intensity (in GFP channel)
%             gfp_cell(i,j) = max(gg, [], 'all');
        else
            gfp_cell(i,j) = nan; 
        end
    end
end

% for the droplets with a single cell only
% cy3_max = double(cy3_max_all(:,idx)); % max
% cy3_median = double(cy3_median_all(:,idx)); % median
% cy3_mean = double(cy3_mean_all(:,idx)); % mean
% cy5_max = double(cy5_max_all(:,idx));
% cy5_median = double(cy5_median_all(:,idx));
% cy5_mean = double(cy5_mean_all(:,idx));
% gfp_max = double(gfp_max_all(:,idx));
% gfp_median = double(gfp_median_all(:,idx));
gfp_mean = double(gfp_mean_all(:,idx));

% normalization for cell gfp (only if single-cell droplets exist)
if isempty(idx) == 0
    for i=1:length(idx)
        gfp_cell_normalized(:,i) = gfp_cell(:,i)./gfp_mean(:,i);
    end
end

%% ratio of cy5/cy3
dropCount = max(length(cy5_median_all), length(cy3_median_all));
check = 0;
if length(cy5_median_all) > length(cy3_median_all)
    diff = length(cy5_median_all) - length(cy3_median_all);
    cy5_median_all(:,dropCount-(diff-1):dropCount) = [];
    cy5_mean_all(:,dropCount-(diff-1):dropCount) = [];
    check = 1;
elseif length(cy3_median_all) > length(cy5_median_all)
    diff = length(cy3_median_all) - length(cy5_median_all);
    cy3_median_all(:,dropCount-(diff-1):dropCount) = [];
    cy3_mean_all(:,dropCount-(diff-1):dropCount) = [];
    check = 1;
end

% delete droplets from idx & idx_empty if either cy3 or cy5 signal doesn't exist for a given droplet
if check == 1 % only apply this if statement above applied
    I = ismember([dropCount-(diff-1):dropCount], idx);
    if I ~= 0
        idx(I) = [];
    end
    
    I = ismember([dropCount-(diff-1):dropCount], idx_empty);
    if I ~= 0
        idx_empty(I) = [];
    end
end

cy5_cy3_ratio_median = (double(cy5_median_all)./double(cy3_median_all));
cy5_cy3_ratio_mean = (double(cy5_mean_all)./double(cy3_mean_all));

%% Calculate initial and end values for single-cell GFP, and cy5/cy3 ratios
cy5_cy3_median_start_end(:,1) = cy5_cy3_ratio_median(1,:); % start median ratio
cy5_cy3_median_start_end(:,2) = cy5_cy3_ratio_median(end,:); % end mean ratio

cy5_cy3_mean_start_end(:,1) = cy5_cy3_ratio_mean(1,:); % start mean ratio
cy5_cy3_mean_start_end(:,2) = cy5_cy3_ratio_mean(end,:); % end mean ratio

if isempty(idx) == 0
    gfp_norm_start_end(:,1) = min(gfp_cell_normalized(1:4, :));
    gfp_norm_start_end(:,2) = max(gfp_cell_normalized(time-3:time,:));
end

%% Generate droplet labeling figure
xf(:,idx_remove) = []; % remove tracking of index idx_remove from single cell centroids
yf(:,idx_remove) = [];

figure; imshow(I_gfp*20);
title(['Droplet index labels - Position: ', pos])
for jjj = 1:length(centroids_set_pos)
    x_plot = centroids_set_pos(jjj,1);
    y_plot = centroids_set_pos(jjj,2);
    text(x_plot-30, y_plot, int2str(jjj),'FontSize',10, 'Color' ,'w'); % annotate each obj with cell id
    hold on;
end
for j = 1:length(xf(1,:))
    x_plot = nanmean(xf(:,j));
    y_plot = nanmean(yf(:,j));
    viscircles([x_plot, y_plot],60,'Color','r');% annotate single-cell droplet IDs
    hold on;
end
for kk = idx_empty
    xx = centroids_set_pos(kk,1);
    yy = centroids_set_pos(kk,2);
    viscircles([xx, yy],60, 'Color','b'); % annotate zero-cell droplet IDs
    hold on;
end

h = zeros(2, 1);
h(1) = plot(NaN,NaN,'or');
h(2) = plot(NaN,NaN,'ob');
legend(h, 'single-cell','zero-cell', 'Location', 'bestoutside');
exportgraphics(gcf,[export,results_save_name,'\Droplet ID\',pos,'.png'], 'ContentType','image');
% set(gcf, 'InvertHardcopy', 'off')
% saveas(gcf,strcat([export,results_save_name,'\',pos]),'png');

%% Saving data

if isempty(idx) == 1 % avoid crash if no single-cell droplets identified
    gfp_cell = [];
    gfp_cell_normalized = [];
    gfp_norm_start_end = [];
end 

cy3_max_store{movie_idx} = double(cy3_max_all); % max
cy3_median_store{movie_idx} = double(cy3_median_all); % median
cy3_mean_store{movie_idx} = double(cy3_mean_all); % mean
cy5_max_store{movie_idx} = double(cy5_max_all);
cy5_median_store{movie_idx} = double(cy5_median_all);
cy5_mean_store{movie_idx} = double(cy5_mean_all);
gfp_max_store{movie_idx} = double(gfp_max_all);
gfp_median_store{movie_idx} = double(gfp_median_all);
gfp_mean_store{movie_idx} = double(gfp_mean_all);
gfp_cell_store{movie_idx} = gfp_cell;
gfp_cell_normalized_store{movie_idx} = double(gfp_cell_normalized);
single_cell_droplets{movie_idx} = idx;
other_droplets{movie_idx} = idx2;
empty_droplets{movie_idx} = idx_empty;
time_store{movie_idx} = time;
cy5_cy3_ratio_mean_store{movie_idx} = cy5_cy3_ratio_mean;
cy5_cy3_ratio_median_store{movie_idx} = cy5_cy3_ratio_median;
cy5_cy3_mean_start_end_store{movie_idx} = cy5_cy3_mean_start_end;
cy5_cy3_median_start_end_store{movie_idx} = cy5_cy3_median_start_end;
gfp_norm_start_end_store{movie_idx} = gfp_norm_start_end;

save(strcat(export,results_save_name,'\data_save','.mat'),'cy3_max_store','cy3_median_store',...
    'cy3_mean_store','cy5_max_store','cy5_median_store','cy5_mean_store','gfp_max_store',...
    'gfp_median_store','gfp_mean_store','gfp_cell_store','gfp_cell_normalized_store',...
    'single_cell_droplets','other_droplets','empty_droplets','time_store',...
    'cy5_cy3_ratio_median_store', 'cy5_cy3_ratio_mean_store', 'cy5_cy3_mean_start_end_store',...
    'cy5_cy3_median_start_end_store', 'gfp_norm_start_end_store');
end

