input = 'C:\Users\marie\Desktop\CP2\Code\Existing Code\Input';
output = 'C:\Users\marie\Desktop\CP2\Code\Existing Code\Output\Tracking';
     
cd(input)

im=imread('as_p1.tif');  
filename='as_p1.tif';   
info=imfinfo(filename);     
num_images=numel(info);  
num_cells=0;


for i=1:num_images
    im_BL=imread(filename,i,'Info',info);
    im_BL_org=im_BL;
    im_BL = imbinarize(im_BL, graythresh(im_BL));
    [centers, radii]=imfindcircles(im_BL,[7 25],'ObjectPolarity','bright','Sensitivity',0.96);
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
    figure_name=append('Image ', num2str(i));
    r=1*ones(length(centers),1);
    figure('Name',figure_name,'NumberTitle','off'); imshow(im_BL_org)
    h = viscircles(centers,r);  
    num_cells_im=length(centers);
    num_cells=num_cells+num_cells_im;
    if i ==1
        num_cells_list=[num_cells_im];
    else
        num_cells_list(end+1)=num_cells_im;
      
    end
%     if i ==3
%         break
%     end
end

avg_num_cells=mean(num_cells_list)