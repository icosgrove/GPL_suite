function [GPL_struct] = GPL_measurements(GPL_struct,sp,quiet_fft,start,finish,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will perform simple measurements on a given detections
% contours. These measurements are used for binary filtration later.

% Originally written by Tyler Helble
% Tested, documented, cleaned, and heavily expanded by Ian, but not
% fundamentally modified.
% 04/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find a single RMS value f  or the sum over each column of singal-absent FFT
% slate.
quiet_fft = 2*quiet_fft.^2; 
quiet_fft = sum(quiet_fft,1); 
quiet_fft = sqrt(mean(quiet_fft));


% Loop over each call, reconstruct contours and take measurements.
for k = 1:length(start)
    
    %% Reconstruct Contours, if they exist. 
    % Use GPL-full to reconstruct contours from compressed state.
    if isfield(GPL_struct,'cm') 
        cm = GPL_full('cm',k,GPL_struct); 
    end
    if isfield(GPL_struct,'cm_max')
        cm_max =  GPL_full('cm_max',k,GPL_struct); 
    end
    if isfield(GPL_struct,'cm_max2')
        cm_max2 = GPL_full('cm_max2',k,GPL_struct); 
    end
      


    %% Duration Measurements
    % Measure the time duration of the contour in time bins and seconds.
    
    % Duration for cm: All islands
    if (parm.measure.cm_slope == 1) || (parm.measure.cm_duration == 1)
        if isfield(GPL_struct,'cm') 
        
            nonz = find(sum(cm)); % Locate nonzero columns
            num_Tbins = length(nonz); 
            GPL_struct(k).cm.duration_sec = num_Tbins * (parm.fftOverlap/parm.SampleFreq); % seconds
            GPL_struct(k).cm.duration_bin = num_Tbins; % bins
        
        end       
    end
    
    % Duration for cm_max: Single Strongest Contour
    if (parm.measure.cm_max_slope == 1) || (parm.measure.cm_max_duration == 1)
        if isfield(GPL_struct,'cm_max')
        
            nonz = find(sum(cm_max)); % Locate nonzero columns
            num_Tbins = length(nonz); 
            GPL_struct(k).cm_max.duration_sec = num_Tbins * (parm.fftOverlap/parm.SampleFreq); % seconds
            GPL_struct(k).cm_max.duration_bin = num_Tbins ; % bins
    
        end
    end

    % Duration for cm_max2: Single Second Strongest Contour
    if (parm.measure.cm_max2_slope == 1) || (parm.measure.cm_max2_duration == 1)
        if isfield(GPL_struct,'cm_max2')
        
            nonz = find(sum(cm_max2)); % Locate nonzero columns
            num_Tbins = length(nonz); 
            GPL_struct(k).cm_max2.duration_sec = num_Tbins * (parm.fftOverlap/parm.SampleFreq); % seconds
            GPL_struct(k).cm_max2.duration_bin = num_Tbins; % bins
        
        end
    end



    %% Slope Measurements

    % Slope for cm: all islands
    if parm.measure.cm_slope == 1 
        if isfield(GPL_struct,'cm')
            
            % Find linear fit params with GPL_slope
            [slope_ci, slope, slope_y_int] = GPL_slope(cm,parm);
            
            % Write outputs to GPL_struct - convert slope to Hz/s
            GPL_struct(k).cm.slope_ci = slope_ci;
            GPL_struct(k).cm.slope = slope*(parm.SampleFreq/parm.fftl)/(parm.fftOverlap/parm.SampleFreq);
            GPL_struct(k).cm.slope_y_int = slope_y_int;

        end
    end
    
    % Slope for cm_max: Single Strongest Contour
    if parm.measure.cm_max_slope == 1 
        if isfield(GPL_struct,'cm_max')
            
            % Find linear fit params with GPL_slope
            [slope_ci, slope, slope_y_int] = GPL_slope(cm_max,parm);
            
            % Write outputs to GPL_struct - convert slope to Hz/s
            GPL_struct(k).cm_max.slope_ci = slope_ci;
            GPL_struct(k).cm_max.slope = slope*(parm.SampleFreq/parm.fftl)/(parm.fftOverlap/parm.SampleFreq);
            GPL_struct(k).cm_max.slope_y_int = slope_y_int;

        end
    end
    
    % Slope for cm_max2: Single Second Strongest Contour
    if parm.measure.cm_max2_slope == 1 
        if isfield(GPL_struct,'cm_max2')
            
            % Find linear fit params with GPL_slope
            [slope_ci, slope, slope_y_int] = GPL_slope(cm_max2,parm);
            
            % Write outputs to GPL_struct - convert slope to Hz/s
            GPL_struct(k).cm_max2.slope_ci = slope_ci;
            GPL_struct(k).cm_max2.slope = slope*(parm.SampleFreq/parm.fftl)/(parm.fftOverlap/parm.SampleFreq);
            GPL_struct(k).cm_max2.slope_y_int = slope_y_int;

        end
    end



    
    %% Frequency Bandwidth Measurements 

    % Frequency Bandwidth for cm: All islands
    if parm.measure.cm_freq_bandwidth == 1
        if isfield(GPL_struct,'cm')
        
            num_Fbins = length(find(sum(cm,2))); % Find # of nonzero freq bins
            GPL_struct(k).cm.freq_bandwidth_hz = num_Fbins*(parm.SampleFreq/parm.fftl); % Hz
            GPL_struct(k).cm.freq_bandwidth_bin = num_Fbins; % bins
        
        end
    end
    
    % Frequency Bandwidth for cm_max: Single Strongest Contour
    if parm.measure.cm_max_freq_bandwidth == 1
        if isfield(GPL_struct,'cm_max')
   
            num_Fbins = length(find(sum(cm_max,2))); % Find # of nonzero freq bins
            GPL_struct(k).cm_max.freq_bandwidth_hz = num_Fbins*(parm.SampleFreq/parm.fftl); % Hz
            GPL_struct(k).cm_max.freq_bandwidth_bin = num_Fbins; % bins
        
        end
    end

    % Frequency Bandwidth for cm_max2: Single Second Strongest Contour
    if parm.measure.cm_max2_freq_bandwidth == 1
        if isfield(GPL_struct,'cm_max2')
        
            num_Fbins = length(find(sum(cm_max2,2))); % Find # of nonzero freq bins
            GPL_struct(k).cm_max2.freq_bandwidth_hz = num_Fbins*(parm.SampleFreq/parm.fftl); % Hz
            GPL_struct(k).cm_max2.freq_bandwidth_bin = num_Fbins; % bins
        
        end
    end




    %% Measure Recieved Level and Background Noise level

    % BG noise level same for each call, prev. calculated outisde for loop
    if parm.measure.spec_noise == 1
    
        GPL_struct(k).spec_noise = quiet_fft;

    end

    % Measure RL using all islands (cm):
    if parm.measure.spec_rl == 1
    
        if length(find(cm)) > parm.measure.RL_min_size % Contour minimum size
    
            % Extract spectrogram snippet of entire call window
            sp_subset = sp(:,start(k):finish(k)); 
            
            % Estimate RL of call vs. background snippet
            [spec_rl] = estimate_rl(cm,sp_subset);
        
        else                                
            spec_rl = nan;    
        end

        % Write output to GPL_struct
        GPL_struct(k).spec_rl = spec_rl; 
        
    end

    
    

    %% Measure Start/End & Min/Max Frequencies of the Contours
    
    % Frequency Limits for cm: All islands
    if parm.measure.cm_freq_limits == 1
        if isfield(GPL_struct,'cm')
            
            % Run Frequency measurements function
            [start_freq, end_freq, min_freq, max_freq, abs_bandwidth] = GPL_freq_measurements(cm,parm);

            % Load outputs into GPL_struct
            GPL_struct(k).cm.start_freq_hz = start_freq;
            GPL_struct(k).cm.end_freq_hz = end_freq;
            GPL_struct(k).cm.min_freq_hz = min_freq;
            GPL_struct(k).cm.max_freq_hz = max_freq;
            GPL_struct(k).cm.abs_bandwidth_hz = abs_bandwidth;
            
        end
    end
    
    % Frequency Limits for cm_max: Single strongest contour
    if parm.measure.cm_max_freq_limits == 1
        if isfield(GPL_struct,'cm_max')
            
            % Run Frequency measurements function
            [start_freq, end_freq, min_freq, max_freq, abs_bandwidth] = GPL_freq_measurements(cm_max,parm);

            % Load outputs into GPL_struct
            GPL_struct(k).cm_max.start_freq_hz = start_freq;
            GPL_struct(k).cm_max.end_freq_hz = end_freq;
            GPL_struct(k).cm_max.min_freq_hz = min_freq;
            GPL_struct(k).cm_max.max_freq_hz = max_freq;
            GPL_struct(k).cm_max.abs_bandwidth_hz = abs_bandwidth;
            
        end
    end
    
    % Frequency Limits for cm_max2: Single second strongest contour
    if parm.measure.cm_max2_freq_limits == 1
        if isfield(GPL_struct,'cm_max2')
            
            % Run Frequency measurements function
            [start_freq, end_freq, min_freq, max_freq, abs_bandwidth] = GPL_freq_measurements(cm_max2,parm);

            % Load outputs into GPL_struct
            GPL_struct(k).cm_max2.start_freq_hz = start_freq;
            GPL_struct(k).cm_max2.end_freq_hz = end_freq;
            GPL_struct(k).cm_max2.min_freq_hz = min_freq;
            GPL_struct(k).cm_max2.max_freq_hz = max_freq;
            GPL_struct(k).cm_max2.abs_bandwidth_hz = abs_bandwidth;
            
        end
    end


end % For: loop over all calls
     
end % Function: GPL_measurements





