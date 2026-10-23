%%load image
I=imread('two.jpg');
figure,imshow(I)
title('Original Image')
%%Binarize
Igray = rgb2gray(I);
level = 0.599;
binaryImage = im2bw(Igray,level);
figure,imshow(binaryImage)
title('Binarized image')
%%Sawp background and foreground pixels
[row col] = size(binaryImage);
binaryImage2=zeros(row,col);
  for i=1:row
      for j=1:col
          if binaryImage(i,j)==0
              binaryImage2(i,j)=1;
          else 
              binaryImage2(i,j)=0;
          end
    end
end
figure,imshow(binaryImage2)
title('Inverse binary Image')
%%Cleaning thining skel fill open close etc
binaryImage2 = bwmorph(binaryImage2,'clean',inf);
figure,imshow(binaryImage2)
title('Cleaned image')
binaryImage2 = imfill(binaryImage2, 'holes');
figure,imshow(binaryImage2)
title('Filled image')
% Get rid of anything less than 1% of the image
binaryImage2 = bwareaopen(binaryImage2, round(0.001*numel(binaryImage2)));
figure,imshow(binaryImage2)
%%Labeling image and regionprop for area
labeledImage = bwlabel(binaryImage2);
number_of_labeles = max(max(labeledImage))
for i= 1:number_of_labeles
      Im = (labeledImage==i)
    figure, imshow(Im)
end
measurements = regionprops(labeledImage, 'Area');
%%Boundaries
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);