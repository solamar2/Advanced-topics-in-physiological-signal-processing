function output_dataset = PSWE_v2(input_dataset)
%% Header
% PSWE_v2 - The PSWE algorithm. Consists of three modules:
%       1) PSD calculation
%       2) PSWE analysis for individual channels
%       3) Global PSWE analysis - looks for overlapping PSWE
%
% Inputs: input_dataset containing preprocessed data from import_edf
%
% Outputs:
%
%       PSD
%           PSD: Continous PSD
%           PSD by bandwidth - delta (1-3), theta (3-8), alpha (8-12)
%           beta (12-20) and gamma (20-50)
%
%       PSWE analysis for individual channels
%           eventData_allChans: PSWE data
%           num_of_events: Number of events
%           time_in_events: Percentage of time spent in events
%           meanFreq: Mean frequency of events
%           meanDur: Mean duration of events
%           medianDur: Median duration of events
%           medianFreq: Median frequency of events
%           medianDur: Median duration of events
%
%       Global PSWE
%           total_PSWE: total number of PSWEs detected
%           persistent_slowing: check variable that indicates if 99% of the
%               recording was detected as a single large event
%           PSWE_per_min: PSWE per minute
%           time_in_events: total time spent in PSWE as a percentage
%           meanDur: average duration of all PSWE
%           medianDur: median duration of all PSWE
%           meanMPF: average MPF of all PSWE
%           medianMPF: median MPF of all PSWE
%
%       Optional outputs:
%           eventPower: Contains raw PSWE segments from all channels

% Sol:
%     start_ind_global_pswe - start index 
%     end_ind_global_pswe- end index
%     num_of_chan -Number of channels involved in the event segment (in hamza code it's precent)
%     chansINevents - the chanels that detected the global PSWE

%
% Author: Hamza Imtiaz
% Date: 2023-06-26

%% Body
% PSWE parameters
freq_threshold = 6; % in hertz (Hz)
dur_threshold = 5; % in seconds (sec)
chan_threshold = 2; % minimum number of channels required to detect PSWE

% Initialize variables
output_dataset = input_dataset;  % Output same as input
fData = input_dataset.data;
srate = input_dataset.hdr.srate;
eventData_allChans = [];

% Check if data is present in input_dataset
if isempty(input_dataset.data)
    return;  % Exit the function
end

%% 1) PSD calculation
% Returns values in uV^2/Hz

% Initialize variables
% PSD = zeros(size(fData, 1), length(srate));
clear PSD
delta = zeros(size(fData, 1), 1);
theta = zeros(size(fData, 1), 1);
alpha = zeros(size(fData, 1), 1);
beta = zeros(size(fData, 1), 1);
gamma = zeros(size(fData, 1), 1);

for channelIdx = 1:size(fData, 1)
    % Continuous PSD
    [temp_PSD, f] = spectopo(fData(channelIdx, :), 0, srate, 'plot', 'off', 'freqfac', 10, 'verbose', 'off');
    % temp_PSD is in dB

    % BW PSD
    deltaIdx = f > 1 & f <= 3;
    thetaIdx = f > 3 & f <= 8;
    alphaIdx = f > 8 & f <= 12;
    betaIdx  = f > 12 & f <= 20;
    gammaIdx = f > 20 & f <= 50;

    % Convert dB to uV^2/Hz
    PSD(channelIdx, :) = 10.^(temp_PSD/10);
    delta(channelIdx) = mean(10.^((temp_PSD(deltaIdx))/10));
    theta(channelIdx) = mean(10.^((temp_PSD(thetaIdx))/10));
    alpha(channelIdx) = mean(10.^((temp_PSD(alphaIdx))/10));
    beta(channelIdx) = mean(10.^((temp_PSD(betaIdx))/10));
    gamma(channelIdx) = mean(10.^((temp_PSD(gammaIdx))/10));
end

clear channelIdx
%% The PSWE algorithm
% Calculate PSWE in individual channels
eventData_allChans = [];
medFreq_allChans = [];

for channelIdx = 1:size(fData, 1) % For all channels

    % Calculate median power frequency (MPF)
    [medFreq]=MPF(fData(channelIdx,:), srate); % Median power frequency

    % The patented PSWE Algorithm
    % Empty vectors
    eventData = zeros(1, length(medFreq));
    eventFreq = [];
    eventLen = [];

    % Find segments where medFreq is lower than threshold frequency for PSWE
    LessThanThresh = find(medFreq < freq_threshold & medFreq >= 1);
    LessThanThresh_diff = diff(LessThanThresh);

    % Initialize variables for tracking slow wave segments
    slow_wave_seg = [];   % Indices of slow wave segments
    k = 1;                % Index variable for current slow wave segment

    % Find the index of the first occurrence of '1' in LessThanThresh_diff
    i = find(LessThanThresh_diff == 1, 1);

    % Store the initial index value for reference
    prev_i = i;

    while ~isempty(LessThanThresh_diff) && ~isempty(i)
        % Start time of event
        ind_of_events(k,1,channelIdx) = LessThanThresh(prev_i) - 1; % Subtracting 1 to align with zero-based indexing

        % Update LessThanThresh_diff, the array containing PSWE segments to be evaluated
        LessThanThresh_diff = LessThanThresh_diff(i:end);
        slow_wave_seg = [slow_wave_seg find(LessThanThresh_diff~=1, 1)]; % Find number of slow wave segments

        % If slow wave segment is found, check if duration of segment is greater than threshold length
        if ~isempty(find(LessThanThresh_diff~=1, 1))
            if slow_wave_seg(end) >= dur_threshold
                % End time of event
                ind_of_events(k,2,channelIdx) = ind_of_events(k,1,channelIdx) + slow_wave_seg(end);

                % Calculate mean frequency and store in eventFreq
                eventFreq = [eventFreq mean(medFreq(ind_of_events(k,1,channelIdx)+1:ind_of_events(k,1,channelIdx)+slow_wave_seg(end)))];

                % Store duration of the segment in eventLen
                eventLen = [eventLen slow_wave_seg(end)];

                % Update eventData to mark the corresponding indices as 1
                eventData(ind_of_events(k,1,channelIdx)+1:ind_of_events(k,2,channelIdx)) = 1;

                % Increment index variable k
                k = k + 1;
            else
                ind_of_events(k,1,channelIdx) = 0; % Set start time of event to 0 if duration is less than threshold
            end
        else
            ind_of_events(k,1,channelIdx) = 0; % Set start time of event to 0 if no slow wave segment is found
        end

        % Update LessThanThresh_diff and i to omit the evaluated segment from the loop
        LessThanThresh_diff = LessThanThresh_diff(find(LessThanThresh_diff~=1, 1):end);
        i = find(LessThanThresh_diff == 1, 1);

        % Update prev_i to start at the end of the detected slow wave segment
        if ~isempty(slow_wave_seg) && ~isempty(i)
            prev_i = prev_i + slow_wave_seg(end) - 1 + i - 1;
        end
    end

    % Update LessThanThresh_diff to omit the evaluated segment from the loop
    LessThanThresh_diff = LessThanThresh_diff(find(LessThanThresh_diff ~= 1, 1):end);
    i = find(LessThanThresh_diff == 1, 1);

    if ~isempty(slow_wave_seg) && ~isempty(i)
        % Update prev_i to start at the end of the detected slow wave segment
        prev_i = prev_i + slow_wave_seg(end) - 1 + i - 1;
    end

    % Update output vectors with results after all segments have been analyzed
    num_of_events_per_chan(channelIdx) = length(find(slow_wave_seg >= dur_threshold)) / ((length(medFreq) - 1) / 60);
    time_in_events_per_chan(channelIdx) = sum(slow_wave_seg(slow_wave_seg >= dur_threshold)) / length(medFreq);
    meanFreq_per_chan(channelIdx) = mean(eventFreq, 'omitnan'); % will return NaN when there is no event in a channel
    meanDur_per_chan(channelIdx) = mean(eventLen, 'omitnan');

    medianFreq_per_chan(channelIdx) = median(eventFreq); % Median MPF
    medianDur_per_chan(channelIdx) = median(eventLen); % Median duration

    eventData_allChans = [eventData_allChans; eventData];
    medFreq_allChans = [medFreq_allChans; medFreq];

end

%% Global PSWE
% Calculate global PSWE; detected in multiple channels simultaneously

% Initialize variables
chansINevents = [];
event_seg = [];
k = 1;
ind_of_events = [];
maxPowerInd = [];
single_event = 0;

% Calculate the number of active channels at each time point
activeChans = sum(eventData_allChans);

% Find the indices where the number of active channels is greater than or equal to the threshold
eventIdx = find(activeChans >= chan_threshold);

% Calculate the duration of events that had more than chansN active channels
eventIdx_diff = diff(eventIdx);

% Check for persistent slowing or if no events were detected in any channel
if length(eventIdx) >= 0.99 * length(activeChans) % Persistent slowing
    flag_PSWE = 1;
    warning('Persistently slow recording detected');
elseif all(eventData_allChans(:) == 0) % No PSWEs
    flag_PSWE = 2;
else
    % PSWEs detected in at least 1 channel
    flag_PSWE = 0;
end

i = find(eventIdx_diff == 1, 1);  % Find the index of the first occurrence of 1 in eventIdx_diff
prev_i = i;  % Store the initial index

while ~isempty(eventIdx_diff) && ~isempty(i)  % Loop through all event segments

    ind_of_events(k, 1) = eventIdx(prev_i);  % Store the index of the event start
    eventIdx_diff = eventIdx_diff(i:end);  % Update eventIdx_diff by removing the processed segment

    event_seg = [event_seg find(eventIdx_diff ~= 1, 1)];  % Store the duration of the event segment

    % If only 1 PSWE is detected
    if isempty(event_seg) && ~isempty(eventIdx_diff)
        eventIdx_diff(end + 1) = 0;  % Append 0 to eventIdx_diff to handle single event case
        event_seg = [event_seg find(eventIdx_diff ~= 1, 1)];  % Update event_seg with the duration of the single event
        single_event = 1;  % Mark recording as single event
    end

    if ~isempty(find(eventIdx_diff ~= 1, 1))  % Check if there are more events in the segment
        if event_seg(end) >= dur_threshold  % Check if the event duration is longer than the length threshold
            ind_of_events(k, 2) = ind_of_events(k, 1) + event_seg(end) -1;  % Store the index of the event end

            % Segments of active channels during the PSWE event
            data = eventData_allChans(:, ind_of_events(k, 1):ind_of_events(k, 2));

            % Find the channel with the highest activity during the event
            [~, maxPowerInd] = max(sum(data, 2));

            % Determine which channels are involved
            [R, ~] = find(data > 0);

            startIdx = max(1, ind_of_events(k, 1));  % Start index for extracting raw data

            endIdx = min(length(fData), ind_of_events(k, 2));  % End index for extracting raw data

            % Store event-related data
            eventPower{k}.raw = fData(:, startIdx * srate:endIdx * srate);
            eventPower{k}.samplepoints = [startIdx * srate; endIdx * srate];
            eventPower{k}.duration = round(endIdx - startIdx);
            eventPower{k}.MaxPowerChan = maxPowerInd;

            chansINevents{k} = length(unique(R));  % Number of channels involved in the event segment
            involved_chan{k}=unique(R); % the channels numbers

            k = k + 1;
        else
            ind_of_events(k, 1) = 0;  % Set the event start index to 0 if the event duration is less than the threshold
        end
    else
        ind_of_events(k, 1) = 0;  % End of the segment
    end

    % Move to the next segment
    eventIdx_diff = eventIdx_diff(find(eventIdx_diff ~= 1, 1):end);  % Remove the processed segment from eventIdx_diff
    i = find(eventIdx_diff == 1, 1);  % Find the index of the next occurrence of 1 in eventIdx_diff

    if ~isempty(event_seg) && ~isempty(i)  % Update the previous index if event_seg and i are not empty
        prev_i = prev_i + event_seg(end) - 1 + i - 1;
    end
end

if isempty(chansINevents) % If no events were detected after checking for overlap
    flag_PSWE = 2;
end

if flag_PSWE == 0 % if signal has PSWE and is not persistently slow
    % Update fields
    num_of_events=length(chansINevents)/(length(eventData)/60);  % events per minute
    num_of_participating_chans=100*cellfun('length', chansINevents)/size(eventData_allChans,1); % percentage of channels in events
    time_in_events=sum(ind_of_events(:,2)-ind_of_events(:,1))/length(eventData);
    meanDur=mean(ind_of_events(:,2)-ind_of_events(:,1));
    medianDur=median(ind_of_events(:,2)-ind_of_events(:,1));
    persistent_slowing = 0

    %% Sol:
    start_ind_global_pswe=ind_of_events(:,1);
    end_ind_global_pswe=ind_of_events(:,2);

elseif flag_PSWE == 1 % if signal is persistently slow
    persistent_slowing = 1;
    % NaN values
    num_of_events = NaN; %Entire recording is PSWE, so events/min = 60
    num_of_participating_chans=NaN;
    time_in_events=0.99;
    meanDur=NaN;
    medianDur = NaN;
    chansINevents=NaN;
    mean_activeChans=NaN;
    max_activeChans=NaN;

    ind_of_events = NaN;
    eventPower = NaN;
    maxPowerInd = randi(size(fData,1));

    %% Sol:
    start_ind_global_pswe=NaN;
    end_ind_global_pswe=NaN;
    chansINevents=NaN;
    involved_chan=NaN;

elseif flag_PSWE == 2 % if signal has no PSWE and is not persistently slow
    persistent_slowing = 0;
    % Zero values
    num_of_events = 0;
    num_of_participating_chans=0;
    time_in_events=0;
    meanDur=0;
    medianDur = 0;
    chansINevents=0;
    mean_activeChans=0;
    max_activeChans=0;

    ind_of_events = 0;
    eventPower = 0;
    maxPowerInd = randi(size(fData,1));

        %% Sol:
    start_ind_global_pswe=0;
    end_ind_global_pswe=0;
    chansINevents=0;
    involved_chan=0;
end

%% Save results
% Continuous PSD
output_dataset.PSD.pxx=PSD;
output_dataset.PSD.f=f;

% PSD by bandwidth
output_dataset.PSD.BW.delta=delta;
output_dataset.PSD.BW.theta=theta;
output_dataset.PSD.BW.alpha=alpha;
output_dataset.PSD.BW.beta=beta;
output_dataset.PSD.BW.gamma=gamma;
output_dataset.PSD.BW.total=sum(PSD,2);

% PSWE per channel
output_dataset.PSWE.eventData_allChans_per_chan = eventData_allChans;
output_dataset.PSWE.num_of_events_per_chan = num_of_events_per_chan;
output_dataset.PSWE.time_in_events_per_chan = time_in_events_per_chan;
output_dataset.PSWE.meanFreq_per_chan = meanFreq_per_chan;
output_dataset.PSWE.meanDur_per_chan = meanDur_per_chan;
output_dataset.PSWE.medianFreq_per_chan = medianFreq_per_chan;
output_dataset.PSWE.medianDur_per_chan = medianDur_per_chan;
output_dataset.PSWE.medFreq_allChans = medFreq_allChans;

% Global PSWE
if num_of_events ~= 0
    output_dataset.PSWE.total_PSWE = size(ind_of_events, 1);
else
    output_dataset.PSWE.total_PSWE = 0;
end
output_dataset.PSWE.persistent_slowing = persistent_slowing;
output_dataset.PSWE.PSWE_per_min=num_of_events;
output_dataset.PSWE.time_in_events = time_in_events;

%
output_dataset.PSWE.activeChansN=num_of_participating_chans; 

output_dataset.PSWE.meanDur = meanDur;
output_dataset.PSWE.medianDur = medianDur;
output_dataset.PSWE.meanMPF = mean(output_dataset.PSWE.medianFreq_per_chan, 'omitnan');
output_dataset.PSWE.medianMPF = median(output_dataset.PSWE.medianFreq_per_chan, 'omitnan');
output_dataset.PSWE.rawData = eventPower;
%% Sol:
  output_dataset.PSWE.start_ind_global_pswe=start_ind_global_pswe;
  output_dataset.PSWE.end_ind_global_pswe=end_ind_global_pswe;
  output_dataset.PSWE.num_of_chan=chansINevents;
  output_dataset.PSWE.involved_chan=involved_chan;
  output_dataset.PSWE.flag_PSWE=flag_PSWE;

  for p = 1:length(involved_chan)
    % Get the vector from the current cell
    currentVector = involved_chan{p};
    % Calculate the length of the vector and add it to the sum
    totalLength(p) = numel(currentVector);
  end
  output_dataset.PSWE.mean_num_involved_chan=mean(totalLength);
  output_dataset.PSWE.max_num_involved_chan=max(totalLength);

% Raw PSWE segments for future studies
output_dataset.eventPower = eventPower;

end

%% Subfunctions

function medFreq = MPF(data, srate)
bufferLength = srate * 2;    % Buffer length in samples (2 seconds)
bufferOverlap = srate * 1;   % Overlap length in samples (1 second)

% Buffer the data
bdata = buffer(data, bufferLength, bufferOverlap);

% Initialize medFreq array
numSegments = size(bdata, 2);
medFreq = zeros(1, numSegments);

% Compute median power frequency for each buffer segment
for i = 1:numSegments
    % Load buffer segment
    y = bdata(:, i);

    % Compute median power frequency using built-in function
    [medFreq(i), ~] = medfreq(y, srate);
end
end