function [data,indNanTrl] = delete_artifact_Timescale( artifact_MJ,data)
% Takes artifacts which should be removed from data.trial and data.time
% and removes the samples.

% Insert NaNs for each artfifact
for iart = 1:size(artifact_MJ,1)
        
   % Insert NaNs for all the artifacts
   data.time{1}(artifact_MJ(iart,1):artifact_MJ(iart,2))   = NaN;
   
   data.trial{1}(:,artifact_MJ(iart,1):artifact_MJ(iart,2))   = NaN;
    
end

% Get NaN indeces
indNanTrl = isnan(data.trial{1}(1,:));

% Remove all the NaN values
data.time{1}(indNanTrl)=[];

data.trial{1}(:,indNanTrl)=[];

end


