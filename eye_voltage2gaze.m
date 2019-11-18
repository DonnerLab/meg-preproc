function [x, y, p] = eye_voltage2gaze(raw, ranges, screen_x, screen_y, ch_mapping)  
    %Converts analog output of EyeLink 1000+ to gaze coordinates
    %Based in Niklas python script
    minvoltage = ranges(1);
    maxvoltage = ranges(2);
    minrange=0;
    maxrange = 1;
    screenright = screen_x(1);
    screenleft = screen_x(2);
    screenbottom = screen_y(1);
    screentop = screen_y(2);
    
%     %obtain the idx of the channels of interest
%     idx = ch_mapping(1);
%     idy = ch_mapping(2);
%     idp = ch_mapping(3);
        
    R = (raw(1,:)-minvoltage)/(maxvoltage-minvoltage);
    S = R*(maxrange-minrange)+minrange;
    x = S*(screenright-screenleft+1)+screenleft;
    
    R = (raw(2,:)-minvoltage)/(maxvoltage-minvoltage);
    S = R*(maxrange-minrange)+minrange;
    y = S*(screenbottom-screentop+1)+screentop;

    p = raw(3,:);  
end