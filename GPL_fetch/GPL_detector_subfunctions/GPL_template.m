function [GPL_struct] = GPL_template(GPL_struct,sp,sp_whiten,start,finish,x,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_template will take each detection and determine where in time and
% frequency the energy spikes occured that originally flagged it as a call 
% during summation. This will idenitfy the contour of the call which is 
% used for later filtering. This function will also compress the contour 
% data and load it into GPL_struct. The original waveform of the call's
% single strongest contour can be computer here as well. 

% Written by Tyler Helble 
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numCalls,~] = size(start);

% Loop over each detection, extract contours and waveform.
for k = 1:numCalls 
    
    % Create whitened spectrogram with selected calls only
    sp_whiten_base = sp_whiten(:,start(k):finish(k)); 

    
    % Locate and create contours associated with the detection
    [cm,cm_max,cm_max2] = GPL_contour(sp_whiten_base,parm);

    
    % Load contour indices, values, and re-creation parameters into
    % GPL_struct for the given detection
    if(parm.cm_on == 1)
        GPL_struct = GPL_sparse(cm,'cm',k,GPL_struct);
    end
    
    if(parm.cm_max_on == 1)
        GPL_struct = GPL_sparse(cm_max,'cm_max',k,GPL_struct);
    end 
    
    if(parm.cm_max2_on == 1)
        GPL_struct = GPL_sparse(cm_max2,'cm_max2',k,GPL_struct);
    end

    %%% reconstruct with:
    %%%  cm = GPL_full('cm',k,GPL_struct);

    % Reconstruct the time-domain waveform of the strongest contour only, 
    % if requested.
    if parm.waveform_on == 1 
        
        % Locate the start/end of the detection in samples
        xv = x(1+(start(k)-1)*parm.fftOverlap:(finish(k)-1)*parm.fftOverlap);
        
        % Create time-domain waveform of cm_max and load into GPL_struct
        [waveform,~] = ww365(xv,cm_max,parm,sp(:,start(k):finish(k)));
        GPL_struct(k).cm_max_waveform = waveform;
    end

end % For: loop over all calls