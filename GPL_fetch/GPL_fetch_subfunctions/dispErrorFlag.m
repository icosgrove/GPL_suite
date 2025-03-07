function dispErrorFlag(error_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function displays the issue that arose in GPL_fetch requiring the
% use re-do a call. 
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch error_flag
    case 0
        return % No error
    case 1
        disp('Error: Selected call is greater than number of calls in this window')
    case 2
        disp('Error: Already at the beginning, cannot go back to re-do')
    case 3
        disp('Manual Re-do of current call')
    case 4
        disp('Detection crosses over xwav file and is skipped')
    case 5
        disp('Error: Input must be an integer')
    case 6
        disp('Error: Empty Input')
    case 7
        disp('Error: Enter an integer >= 0')
    case 8
        disp('Error: Please pair detection(s) before flagging for reprocess')
    case 9
        disp('Error: Invalid special case input')
    case 10
        disp('Error: Invalid Input.')
    case 11
        disp('Error: Enter a positive integer or 0 for special cases')
    case 12
        disp('Error: Selected call number is greater than number of calls in window.')
    case 13
        disp('Error: Zero not allowed for this input.')
    case 14
        disp('Error: Please enter "c" or "n".')
    case 15
        disp('Warning: Possible incorrect xwav file selected, manual detection is out of time range')
    case 16
        disp('Error: Please enter "c" or "b".')
    case 17
        disp('User cancelled file selection.')
    case 18
        disp('Error: Please pair the primary call. (Adhocs & Window flags are saved)')
    case 19
        disp('Error: Pleae enter "c" or "new".')
    case 20
        disp('Error: Manual log has already been fully paired (last pick is final element.')
    case 21
        disp('Error: Please pair a case to the green manual detection before proceeding.')
    case 22
        disp('Error: Please enter "y" or "n".')
    case 23
        disp('Error: Please enter "c" or "a"')
    case 24
        disp('Error: No TF File entered in GPL_fetch.m')
    case 25
        disp('Error: Please enter "b": brightness, "c": contrast, "d" to continue')
    case 26
        disp('Error: Please enter a number.')
    case 27
        disp('Error: Please enter "sp" or "gpl".')
    case 28
        disp('Error: End & Enable Triton spectrogram to adjust Brightness/Contrast')
    case 29
        disp('Error: Please enter "1" or "2".')
end