function [chain] = radial_island(orig_index,sz1,sz2,parm)

%%%%%%%% To be re-done for optimization or omitted entirely for another
%%%%%%%% method



% load detections_GPL_v2_35_90_down90-30Hz_NUNAT_SB01_20210816_060302.x.x.mat
% orig_index = hyd.detection.calls(8).cm.index;
% sz1 = hyd.detection.calls(8).cm.size(1);
% sz2 = hyd.detection.calls(8).cm.size(2);
% parm = hyd.detection.parm;
sp = zeros(sz1,sz2);
sp(orig_index) = 1;

radius = parm.island_radius;
sp = [zeros(radius,sz2+2*radius); zeros(sz1,radius) sp zeros(sz1,radius); zeros(radius,sz2+2*radius)];
index = find(sp);


sp_copy = sp;
island = []; length_count = 0; c_length = []; island_count = 0; count=0;
s2 = size(sp); sz3 = s2(1); sz4 = s2(2); sz_max = sz3*sz4;
%% attempt 4



for w = 1:length(index) % loop over contour pixels
    
    if isfinite(sp(index(w))) == 1 % check if index is already icluded in an island
        
    next = index(w);
        
    while isempty(next) == 0 % loop until an island is created
        
        n = find(index==min(next));   
    
        
        cross1 = GPL_cross(index,n,radius,sz1);
%         sp(cross1) = 5; plot_sp(sp,'Cross',1);
        
        
        k = find(sp(cross1)~=0);
        sp(cross1(k)) = nan; % find contour pixels in radius around original pixel and get rid of them
    
    
    
        next = [];
    
        for m = 1:length(k) % loop over all NaN's
        
        
            cross2 = GPL_cross(cross1,m,radius,sz1);
            
            cross2(abs(cross2)~=cross2) = []; % kill negatives and indices higher than window max
            cross2(cross2 > sz_max) = [];
            
            
            % find 1'a around NaN's
            next = [next cross2(sp(cross2)==1)]; % save indices of the next set of 1's
            next = unique(next); % kill duplicates, most are overlaps
            
              
            
            if isempty(next) == 1 % double check in case of self islanding (search NaN's for nearby ones)
                border = find(isnan(sp));
                
                for b = 1:length(border)
                    check = [border(b)-(sz1+2*radius)-1:border(b)-(sz1+2*radius)+1 ...
                        border(b)-1:border(b)+1 border(b)+(sz1+2*radius)-1:border(b)+(sz1+2*radius)+1];
                    
                    if sum(~isnan(sp(check))) > 0 % skips nans with no possibility of 1 or 0 next to it
                        cross3 = GPL_cross(border,b,radius,sz1);
                    
                        passed = sp(cross3)==1;
                        next = [next cross3(passed)];
                        next = unique(next);
                        passed = [];
                        
                    end
                    check = [];
                end
                
            end
                      
        end % next set of ones is found
%         nans = find(isnan(sp));
%         sp_copy(nans)=3;
%         p = plot_sp(sp_copy,'Radius=1',0);
%         name = sprintf('fig%i.png',d);
%         saveas(p,name);
%         d = d+1;

        count = count+1;
    end
    
    pos_ind = find(isnan(sp))';
    c_length = [c_length length(pos_ind)-length_count*(sum(c_length(1:end)))];
    length_count = 1;
    
    isle = [c_length(end) pos_ind];
    island = [island isle];
    clear pos_ind
    island_count = island_count + 1;

    
    
    end
    
    island = unique(island,'stable');
end

for y = 1:island_count
    if island(1) >= length(island)
        island = [length(island) island];
    end
    if island(1) >= length(island)
        island(length(island):island(1)) = 1;
    end
    isl_ind{y} = island(2:island(1)+1);
    island(1:island(1)+1) = [];
    sp(isl_ind{y}) = y;
    if isempty(island) == 1
        break
    end
end

sp(:,1:radius) = []; 
sp(:,end-radius+1:end) = []; 
sp(1:radius,:) = []; 
sp(end-radius+1:end,:) = []; 

chain = [];
for f = 1:island_count
    ind = find(sp==f);
    i_ind = length(ind);
    chain = [chain i_ind ind'];
    clear ind l_ind
end
%plot_sp(sp,'Radial Test',1);


























