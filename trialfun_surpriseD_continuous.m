function [trl, event] = trialfun_surpriseD_continuous(cfgin)
%TRIALFUN_SURPRISED_CONNTINUOUS trial definition function
%   for continuous data
%   cfgin.trialfun = 'trialfun_surpriseD'  

ID = cfgin.ID;

% TrgVal = 1 - Start of block
% TrgVal = 2 - End of block

% read header and event information
hdr = ft_read_header(cfgin.dataset);
event = ft_read_event(cfgin.dataset);

% Search for trigger events
trgVal = [event(find(strcmp('UPPT001',{event.type}))).value]';
trgSample = [event(find(strcmp('UPPT001',{event.type}))).sample]';

% Find specific triggers/samples for start and end of block
startSample = trgSample(find(trgVal==1));
endSample = trgSample(find(trgVal==2));


% Check whether start triggers are missing; add first start sample
if length(startSample) < length(endSample)
    startSample = [1+round(cfgin.trialdef.prestim*hdr.Fs); startSample];
end

% Check whether stop triggers are missing (interrupted block); add end
% sample
if length(startSample) > length(endSample)
    endSample = [endSample; trgSample(end)-round(cfgin.trialdef.poststim * hdr.Fs)];
end

% Check duration of block and discard if > 1 minutes
% In some cases EJG-1_01 there are 6 blocks identified because the task has started
% without eyelink, actually only 4 of them are usable
for i = 1:length(startSample)
    if (endSample(i)-startSample(i))/hdr.Fs < 60 % 60 seconds 
        startSample(i) = 0; % Set the start/end trigger samples that are not useable to zero to remove them afterwards
        endSample(i) = 0;
    end
end

startSample(startSample == 0) = [];
endSample(endSample == 0) = [];

if isfield(cfgin.trialdef, 'prestim')  
    trlOff = round(-cfgin.trialdef.prestim*hdr.Fs);
    startSample = max([ones(length(startSample),1),startSample + trlOff],[],2);
end


if isfield(cfgin.trialdef, 'poststim')
    endSample = min([event(end).sample*ones(length(endSample),1), endSample + round(cfgin.trialdef.poststim * hdr.Fs)],[],2);
end



% define trial matrix   trial x M (M = start sample; endSample; offset)
trl = [startSample, endSample, trlOff*ones(length(startSample),1)];
if strcmp('UDK-1_03',ID) && cfgin.block == 4 % Do not use 4th block, it was not started on purpose
   trl(4,:) = [];
end
trl = trl(cfgin.block,:);

end



