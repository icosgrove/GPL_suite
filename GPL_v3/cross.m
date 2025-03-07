function cross = cross(index,n,radius,sz1)
% create a cross around a pixel n, with radius 'radius', and using values
% from the set of indices 'index'

top_right=[]; top_left=[]; bottom_right=[]; bottom_left=[];

middle = index(n);
top = (index(n)+1:index(n)+radius); %sp(top)=4; % top
bottom = (index(n)-radius:index(n)-1); %sp(bottom)=4; % bottom
right = (index(n)+sz1+2*radius:sz1+2*radius:index(n)+(sz1+2*radius)*radius); %sp(right)=4; % right
left = (index(n)-sz1-2*radius:-sz1-2*radius:index(n)-(sz1+2*radius)*radius); %sp(left)=4; % left
    
for m = 1:radius-1
    top_right = [top_right (top(radius-m)+sz1+2*radius:sz1+2*radius:top(radius-m)+(sz1+2*radius)*m)]; %sp(top_right)=3; % top right
    top_left = [top_left (top(radius-m)-sz1-2*radius:-sz1-2*radius:top(radius-m)-(sz1+2*radius)*m)]; %sp(top_left)=3; % top left
    bottom_right = [bottom_right (bottom(m+1)+sz1+2*radius:sz1+2*radius:bottom(m+1)+(sz1+2*radius)*m)]; %sp(bottom_right)=3; % bottom right
    bottom_left = [bottom_left (bottom(m+1)-sz1-2*radius:-sz1-2*radius:bottom(m+1)-(sz1+2*radius)*m)]; %sp(bottom_left)=3; % bottom left
end

cross = [middle top bottom right left top_right top_left bottom_right bottom_left];