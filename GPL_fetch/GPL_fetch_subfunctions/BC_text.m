function BC_text(parm,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays text for the brightness and contrast
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mode == 1 % Subplot 1
    b = parm.fp1.brightness;
    c = parm.fp1.contrast;
elseif mode == 2
    b = parm.fp1.gpl_brightness;
    c = parm.fp1.gpl_contrast;
end

% Check for decimals
if mod(b, 1) == 0
    b_str = sprintf('%d', b);
else
    b_str = sprintf('%.2f', b);
end

if mod(c, 1) == 0
    c_str = sprintf('%d', c);
else
    c_str = sprintf('%.2f', c);
end

% Add text
text('Units', 'normalized', 'Position', [1 -0.18 0], ...
    'String', sprintf('B: %s, C: %s', b_str, c_str), ...
    'HorizontalAlignment', 'right');