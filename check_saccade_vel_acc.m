function [sacc_times] = check_saccade_vel_acc(xgaze, ygaze, Hz, threshold, acc_thresh, amp_thresh, ppd)
    
    %sacc = 0;  % initializing binary output flag that indicates whether this trial contains a saccade
    sacc_times = [];
    
    % get x and  y in degrees
    x = xgaze/ppd;
    y = ygaze/ppd;
    
    [velocity, acceleration] = get_velocity (x, y, Hz);
    saccades = double(velocity > threshold);  % getting vector of samples that violate velocity threshold
    
    borders = diff(saccades);   % start points of candidate saccades will be 1, end points will be -1, all others 0
    if saccades(1)>threshold, borders = [1 borders]; else borders = [0 borders]; end  % in case first sample violates threshold
    if saccades(end)>threshold, borders(end+1) = -1; end  % in case last sample violates threshold
    
    window_size = 3;
    win = ones(1,window_size)/double(window_size);
    x = conv(x, win, 'same'); y = conv(y, win, 'same');   % lightly smoothing gaze time series before applying amplitude threshold
    
    starts = find(borders==1); ends = find(borders==-1)-1;  % getting all start and end points of candidate saccades
    if length(starts)>length(ends), starts(end) = []; end   % if last saccade occurs right at trial end (and so doesn't have corresponding ends), discard
    for i = 1:length(starts)  % looping through all candidate saccades and only accepting if they also violate acceleration/amplitude thresholds
        if ~isempty(find(acceleration(starts(i):ends(i))>acc_thresh))  % applying acceleration threshold
            % applying amplitude threshold based on average gaze position over 3 samples before vs after saccade
            if i>1 && starts(i)-3>0  % a few exceptions to try to make sure preceding 3 samples are useable
                t1 = max([starts(i)-3,ends(i-1)]):starts(i)-1;
            elseif starts(i)<=2
                t1 = 1;
            else t1 = starts(i)-1;
            end
            if i<length(starts) && ends(i)+3<=length(x)  % a few exceptions to try to make sure following 3 samples are useable
                t2 = ends(i)+1:min([ends(i)+3,starts(i+1)]);
            elseif ends(i)>=length(x)-1
                t2 = length(x);
            else t2 = ends(i)+1;
            end
            p1 = [mean(x(t1)) mean(y(t1))];  % x,y coords of pre-saccade gaze position (in d.v.a., always positive)
            p2 = [mean(x(t2)) mean(y(t2))];  % x,y coords of post-saccade gaze position (in d.v.a., always positive)
            if ((p1(1)-p2(1)).^2 + (p1(2)-p2(2)).^2).^.5 > amp_thresh;  % applying amplitude threshold
                %sacc = 1;
                sacc_times(end+1,1:2) = [starts(i) ends(i)];
            end
        end
    end
    
%     figure(8)
%     plot(gaze_seg_test(1,:),'color',[0.5 0.5 0.5]), hold on
%     plot(gaze_seg_test(2,:),'color',[0.3 0.3 0.3]), hold on
%     scatter(idx_vel,-90*ones(length(idx_vel),1),'r'), hold on
%     scatter(idx_acc,-130*ones(length(idx_acc),1),'b'), hold on
%     scatter(idx_amp,-170*ones(length(idx_amp),1),'m')
%     legend('x','y','velocity exceeded','acceleration exceeded','amplitude exceeded','Location','bestoutside')
        
end

