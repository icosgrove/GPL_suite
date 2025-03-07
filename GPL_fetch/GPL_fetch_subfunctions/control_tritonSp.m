% control_tritonSp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will plot the spectrogram created by Triton in GPL fetch,
% and allow the user to adjust the brightness and contrast
% Code is sourced from plot_specgram.m and mkspecgram.m from Triton version
% 1.95.20230315. All credit goes to the Scripps Whale Acoustics Laboratory.
% I am repurposing the code with minor changes.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User input to configure B/C

disp('Entering Spectrogram Settings')
subplot(3,1,1); hold on
proceed = false;
while ~proceed
    ch = input('Enter "b" for brightness, "c" for contrast, "d" when done to continue: ','s');
    switch ch
        case 'b'
            b_string = sprintf('New Brightness (currently %i): ',parm.fp1.brightness);
            b = input(b_string,'s');
            b = str2double(b);
            if isnan(b)
                dispErrorFlag(26)
                continue
            else
                parm.fp1.brightness = b;
            end
            cla
            [freq_f,time,sp_context,parm] = makeTritonSp(sub_data,parm);
            fetched_call.spectrogram = sp_context;
            plot_tritonSp
        case 'c'
            c_string = sprintf('New Contrast (currently %i): ',parm.fp1.contrast);
            c = input(c_string,'s');
            c = str2double(c);
            if isnan(c)
                dispErrorFlag(26)
                continue
            else
                parm.fp1.contrast = c;
            end
            cla
            [freq_f,time,sp_context,parm] = makeTritonSp(sub_data,parm);
            fetched_call.spectrogram = sp_context;
            plot_tritonSp
        case 'd'
            disp('Exiting Settings.')
            proceed = true;
        otherwise
            dispErrorFlag(25)
    end
end
hold off


