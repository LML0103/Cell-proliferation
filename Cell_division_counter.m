close all
clear all
% set the folder of the .tif file here
input = '';
% set the filename of the .tif file inside the "input" folder here
filename='';   
     
% if desired set save_plots to true and specify the intended output folder
save_plots = false;
save_path = '';

% if desired set plot_ims to true (required to save the plots)
plot_ims = true;

cd(input)
info=imfinfo(filename);     
num_images=numel(info);  

% set the ground truth here
cell_divisions_gt=[0 4 7 2 3 6 2 7 4 1 2 4 19 16 15 8 14 8 18 7 11 1 3 13 14 11 7 11 10 9 3 10 10 3 4 15 8]'; %ground truth

sensitivity = 0.9465; % determind by hyperparameter tuning
threshold1_factor = 2.1375; % determind by hyperparameter tuning
threshold2_factor = 3.6; % determind by hyperparameter tuning
delete_circle_threshold_1 = 20; % determind by hyperparameter tuning
delete_circle_threshold_2 = 35; % determind by hyperparameter tuning
find_circle_threshold_low = 6; % determind by hyperparameter tuning
find_circle_threshold_high = 18; % determind by hyperparameter tuning
threshold_is_cell_division = 74; % determind by hyperparameter tuning
threshold_prev_cell_division = 15; % determind by hyperparameter tuning



cell_divisions=0; % initialization of the cell_division counter
cell_divisions_all = cell_divisions_gt;
prev_cell_divisions_centers = [];

for i=1:num_images
    im_BL=imread(filename,i,'Info',info);
    cell_divisions_image=0; % initialization of the cell division counter of the current image
    cell_divisions_centers=[];   
    if i~=1 % its not possible to detect cell divisions on the first image (unclear how long ago a cell division took place)
        diff = im_BL - im_BL_previous; 
        % the difference between the current and the previous image is
        % used for the division detection as thereby the movement of
        % cells gets visible (similar to the concept of an optical flow
        % algorithm). The main driving force for cell movement are cell
        % divisions and hence the usage of the difference of the
        % current and previous image leads to an increase in the
        % performance of the implemented algorithm

        threshold1=mean(im_BL,"all")*threshold1_factor; % later used to set background pixels to 0
        threshold2=mean(im_BL,"all")*threshold2_factor; % later used to set background pixels to 0

        % the following preprocessing has the main purpose of firstly
        % increasing the contrast between the pixel values of the cells
        % and the background and secondly making the cells more round
        % in order to increase the performance of the following
        % detection function
        diff(diff<=threshold1)=0; 
        
        diff = imgaussfilt(diff); % preprocessing: Gauss filter
        diff = imdilate(diff,offsetstrel('ball',3,3)); % preprocessing: dilation filter
        diff = imerode(diff,offsetstrel('ball',3,3)); % preprocessing: erosion filter
        diff = imgaussfilt(diff); % preprocessing: gauss filter
        diff = imdilate(diff,offsetstrel('ball',3,3)); % preprocessing: dilation filter
        diff = imerode(diff,offsetstrel('ball',3,3)); % preprocessing: erosion filter

        diff(diff<=threshold2)=0;
          
        % the folowing function finds circles in 2D data that are
        % brighter than the background
        [centers, radii]=imfindcircles(diff,[find_circle_threshold_low find_circle_threshold_high],'ObjectPolarity','bright','Sensitivity',sensitivity); %function that finds circels on the image
        
        % the following loop verifies that the above function did not 
        % put several circles around the same cell
        if length(centers) > 0
            for k=1:length(centers(:,1))
                diffs=centers - centers(k,:);
                diffs=abs(diffs(:,1))+abs(diffs(:,2));
                
                indi=find(diffs<=delete_circle_threshold_1);
                if length(indi)>1
                    centers(k,:)=[-100,-100];
                    
                end
            end
            % if several circles were found too close to each other,
            % all but one get deleted from the "cirlce" list
            indices = find(centers<0);
            centers(indices)=NaN;
        end

        % the following loop verifies that the above function did not 
        % put a big circle around a group (usually 2) cells together
        if length(centers) > 0
            for k=1:length(centers(:,1))
                diffs=centers - centers(k,:);
                diffs=abs(diffs(:,1))+abs(diffs(:,2));
                
                indi=find(diffs<=delete_circle_threshold_2);
                if length(indi)>2
                    centers(k,:)=[-100,-100];
                    
                end
            end
            % deleting the identified threshold transgressions
            indices = find(centers<0);
            centers(indices)=NaN;
        end
        centers_org = centers;
        
        % the following loop checks if there were cell divisions on the
        % same spots on the previous image and if so, deletes the circles
        % from the centers list.
        if length(centers) > 0
            if length(prev_cell_divisions_centers) > 0
                for k=1:length(centers(:,1))
                    if ~ isnan(centers(k,1))
                        distances_prev = prev_cell_divisions_centers - centers(k,:);
                        distances_sum_prev=abs(distances_prev(:,1))+abs(distances_prev(:,2));  
                        indices = find(distances_sum<threshold_prev_cell_division);
                        if min(distances_sum_prev) < threshold_prev_cell_division
                            centers(k,:) = [-100 -100];                                           
                        end
                    end
                end
                indices = find(centers<0);
                centers(indices)=NaN;
            end

            % the following loop finds the closest adjecing circles
            % (=cells) for each circle and if there are two in close
            % proximity, they get identified as a cell division.
            % As a reminder: the difference between to images is used,
            % so the "centers" mainly only consist of cells from cell
            % divisions.
            for k=1:length(centers(:,1))
                distances = centers - centers(k,:);
                distances_sum=abs(distances(:,1))+abs(distances(:,2));  
                indices = find(distances_sum<=threshold_is_cell_division);
                % if either 2 cells (one cell division) or 4 (two cell
                % divisions) are in close proximity goverend by the
                % optimized threshold_is_cell_division, it counts as
                % one/two cell divisions
                if length(indices) == 2 || length(indices) == 4
                    cell_divisions = cell_divisions + length(indices)/2;
                    for idx=1:length(indices) 
                        if length(cell_divisions_centers)==0
                            cell_divisions_centers=centers(indices(idx),:);
                        else
                            cell_divisions_centers(end+1,:)=centers(indices(idx),:);
                        end
                        % the detected cells are deleted from the list to
                        % prehibit the algorithm from counting a cell
                        % division multiple times
                        centers(indices(idx),:) = [NaN, NaN];
                    end
                end
            end

            for k=1:length(centers(:,1))
                distances = centers - centers(k,:);
                distances_sum=abs(distances(:,1))+abs(distances(:,2));  
                indices = find(distances_sum<=threshold_is_cell_division);
                % if in the remaining cells three are in close proximity it
                % is counted as one cell division
                if length(indices) == 3
                    cell_divisions = cell_divisions + 1;
                    for idx=1:length(indices) 
                        if length(cell_divisions_centers)==0
                            cell_divisions_centers=centers(indices(idx),:);
                        else
                            cell_divisions_centers(end+1,:)=centers(indices(idx),:);
                        end
                        centers(indices(idx),:) = [NaN, NaN];
                    end
                end
            end
        end

        % plotting of the results
        if plot_ims
            r=20*ones(length(cell_divisions_centers),1);
            figure_name=append('Image ',num2str(i));
            iptsetpref('ImshowBorder','tight'); %removal of figure border
            figure('Name',figure_name,'NumberTitle','off'); 
            imshow(im_BL,'InitialMagnification',100)
            h=viscircles(cell_divisions_centers,r);
            % saving of the plots
            if save_plots
                print(fullfile(save_path,[figure_name '.jpg']),"-dpng");
            end
        end
    end
    % the image data is saved to allow the division on the next iteration
    im_BL_previous = im_BL;
    prev_cell_divisions_centers = cell_divisions_centers;

    % saving the detected cell_divisions of the current image in a list
    if i == 1 
        cell_divisions_list=[cell_divisions]';
    else
        cell_divisions_list(end+1,:)=cell_divisions-sum(cell_divisions_list);

    end
    % adding the current cell divisions to an array of all cell divisions
    cell_divisions_all(i,2) = cell_divisions_list(i);
end

% suming up the total number of cell divisions
cell_divisions_all(end+1,1)=sum(cell_divisions_all(:,1));
cell_divisions_all(end,2) = sum(cell_divisions_all(:,2));

% error computation
error = (sum(abs(cell_divisions_all(1:end-1,2) - cell_divisions_all(1:end-1,1))));