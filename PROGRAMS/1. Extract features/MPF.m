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