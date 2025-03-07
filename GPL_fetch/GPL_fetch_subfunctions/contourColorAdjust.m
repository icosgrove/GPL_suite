function [fetched_call] = contourColorAdjust(fetched_call,parm,fig,time_range,freq_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will adjust the 0 values of the GPL contour window to be
% the floor of the colorbar used in the top spectrogram window. This
% should help show off the colors of the contours better.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty pixels
disp('Entering GPL Contour Window Settings.')
zero_idx = find(fetched_call.gpl_contour_window == 0);

% Set brightness and contrast
b = 0; c = 100;
ctrl = false;
while ~ctrl
    setType = input('Enter "b" for brightness, "c" for contrast, "d" to continue: ','s');
    switch setType
        case 'b'
            bstr = sprintf('Enter brightness (currently %d): ',b);
            b_copy = b;
            b = str2double(input(bstr,'s'));
            if isnan(b)
                dispErrorFlag(26)
                b = b_copy;
            else
                chgContour(fetched_call,parm,b,c,zero_idx,fig,time_range,freq_range)
            end
        case 'c'
            cstr = sprintf('Enter a contrast (currently %i): ',c);
            c_copy = c;
            c = str2double(input(cstr,'s'));
            if isnan(c)
                dispErrorFlag(26)
                c = c_copy;
            else
                chgContour(fetched_call,parm,b,c,zero_idx,fig,time_range,freq_range)
            end
        case 'd'
            ctrl = true;
        otherwise
            dispErrorFlag(25)
    end
end
end % Function: contourColorAdjust

