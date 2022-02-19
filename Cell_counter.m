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

% initialization of the total number of cells
num_cells=0;

for i=1:num_images
    im_BL=imread(filename,i,'Info',info);
    % saving the original image for visualization purposes
    im_BL_org=im_BL;
    % applaying a Hough transform binarization filter
    im_BL = imbinarize(im_BL, graythresh(im_BL));
    % detecting the cells as circles
    [centers, radii]=imfindcircles(im_BL,[7 25],'ObjectPolarity','bright','Sensitivity',0.96);
    % deleting circles that are too close to each other
    for k=1:length(centers(:,1))
        diffs=centers - centers(k,:);
        diffs=abs(diffs(:,1))+abs(diffs(:,2));       
        indi=find(diffs<=30);
        if length(indi)>1
            centers(k,:)=[-100,-100];
        end
    end
    indices = find(centers<0);
    centers(indices)=NaN;

    % plotting the results
    if plot_ims
        figure_name=append('Image ', num2str(i));
        r=1*ones(length(centers),1);
        figure('Name',figure_name,'NumberTitle','off'); imshow(im_BL_org)
        h = viscircles(centers,r);  

        % saving the plot
        if save_plots
            print(fullfile(save_path,[figure_name '.jpg']),"-dpng");
        end
    end

    % saving the number of cells in the cell list and adding it to the
    % total number of cells
    num_cells_im=length(centers);
    num_cells=num_cells+num_cells_im;
    if i ==1
        num_cells_list=[num_cells_im];
    else
        num_cells_list(end+1)=num_cells_im;
      
    end
end

% printing the average number of detected cells
avg_num_cells=mean(num_cells_list)