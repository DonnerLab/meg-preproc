
function [velocity,acceleration] = get_velocity (x, y, Hz)

    % Compute velocity and acceleration of eye movements
    % Based on Niklas py script "The function asumes that the values in x,y are
    % sampled continuously at a rate specified by Hz"
    velocity_window_size = 3;
    Hz = double(Hz);
    distance = (diff(x).^2 + diff(y).^2).^.5;
    distance = [distance(1) distance];
    win = ones(1,velocity_window_size)/double(velocity_window_size);
    velocity = conv(distance, win, 'same');
    velocity = velocity / (velocity_window_size/Hz);
    acceleration = diff(velocity) / (1/Hz);
    acceleration = abs([acceleration(1) acceleration]);

end