%% SMLM single emitter intensity analysis
%
% Link emitters when: d^2 = sum_axis (position_i-position_j)^2 /
% (sigma_i^2+sigma_j^2) <= chi^2.


close all;
clc;

inputCsv = "Locs_FOV8_345mbar_50ms.csv";

chi = 3.368;
intensityMultipleTolerance = 0.10;  % Relative tolerance around 2N, 3N, or 4N

tracksPerFigure = 20;  

%% 
inputCsv = string(inputCsv);
opts = detectImportOptions(inputCsv, 'VariableNamingRule', 'preserve');
linkedLocalizations = readtable(inputCsv, opts);

required = ["frame", "x[nm]", "y[nm]", "z[nm]", ...
    "UncertaintyX", "UncertaintyY", "UncertaintyZ", ...
    "intensity[photon]"];
missing = setdiff(required, string(linkedLocalizations.Properties.VariableNames));

if ~isempty(missing)
    error('Missing required CSV column(s): %s', strjoin(missing, ', '));
end

frames = linkedLocalizations.("frame");
position = [linkedLocalizations.("x[nm]"), ...
            linkedLocalizations.("y[nm]"), ...
            linkedLocalizations.("z[nm]")];
uncertainty = [linkedLocalizations.("UncertaintyX"), ...
               linkedLocalizations.("UncertaintyY"), ...
               linkedLocalizations.("UncertaintyZ")];
photons = linkedLocalizations.("intensity[photon]");

numberOfLocalizations = height(linkedLocalizations);
neighbors = cell(numberOfLocalizations, 1);
uniqueFrames = unique(frames).';

% Link emitters in consecutive frames.
for frame = uniqueFrames
    current = find(frames == frame);
    following = find(frames == frame + 1);
    if isempty(following)
        continue
    end
    for a = 1:numel(current)
        i = current(a);
        delta = position(following, :) - position(i, :);
        combinedVariance = uncertainty(following, :).^2 + uncertainty(i, :).^2;
        normalizedDistance = sqrt(sum(delta.^2 ./ combinedVariance, 2));
        accepted = following(normalizedDistance <= chi);
        for b = 1:numel(accepted)
            j = accepted(b);
            neighbors{i}(end + 1) = j; 
            neighbors{j}(end + 1) = i; 
        end
    end
end

% DBSCAN with min_samples=1: label connected epsilon neighborhoods.
trackId = zeros(numberOfLocalizations, 1);
numberOfTracks = 0;
for startIndex = 1:numberOfLocalizations
    if trackId(startIndex) ~= 0
        continue
    end
    numberOfTracks = numberOfTracks + 1;
    trackId(startIndex) = numberOfTracks;
    queue = startIndex;
    queuePosition = 1;
    while queuePosition <= numel(queue)
        currentIndex = queue(queuePosition);
        queuePosition = queuePosition + 1;
        for neighbor = neighbors{currentIndex}
            if trackId(neighbor) == 0
                trackId(neighbor) = numberOfTracks;
                queue(end + 1) = neighbor;
            end
        end
    end
end
linkedLocalizations.track_id = trackId;

firstFrame = zeros(numberOfTracks, 1);
lastFrame = zeros(numberOfTracks, 1);
onTimeFrames = zeros(numberOfTracks, 1);
numberOfDetections = zeros(numberOfTracks, 1);
numberOfUniqueFrames = zeros(numberOfTracks, 1);
duplicateFrameFlag = false(numberOfTracks, 1);
meanX_nm = zeros(numberOfTracks, 1);
meanY_nm = zeros(numberOfTracks, 1);
meanZ_nm = zeros(numberOfTracks, 1);

for id = 1:numberOfTracks
    members = find(trackId == id);
    memberFrames = frames(members);
    firstFrame(id) = min(memberFrames);
    lastFrame(id) = max(memberFrames);
    onTimeFrames(id) = lastFrame(id) - firstFrame(id) + 1;
    numberOfDetections(id) = numel(members);
    numberOfUniqueFrames(id) = numel(unique(memberFrames));
    duplicateFrameFlag(id) = numberOfDetections(id) ~= numberOfUniqueFrames(id);

    weights = 1 ./ sum(uncertainty(members, :).^2, 2);
    weights = weights ./ sum(weights);
    weightedPosition = sum(position(members, :) .* weights, 1);
    meanX_nm(id) = weightedPosition(1);
    meanY_nm(id) = weightedPosition(2);
    meanZ_nm(id) = weightedPosition(3);
end

tracks = table((1:numberOfTracks).', firstFrame, lastFrame, onTimeFrames, ...
    numberOfDetections, numberOfUniqueFrames, duplicateFrameFlag, ...
    meanX_nm, meanY_nm, meanZ_nm, ...
    'VariableNames', {'track_id', 'first_frame', 'last_frame', ...
    'on_time_frames', 'n_localizations', 'n_unique_frames', ...
    'duplicate_frame_flag', 'mean_x_nm', 'mean_y_nm', 'mean_z_nm'});

% Select the track-duration range, currently (5,10]
longTracks = tracks(tracks.on_time_frames > 5 & ...
    tracks.on_time_frames <= 10, :);
longTrackIds = longTracks.track_id;
longEmitterLocalizations = linkedLocalizations( ...
    ismember(linkedLocalizations.track_id, longTrackIds), :);

%%
% Find candidate 2N, 3N, and 4N intensity events within each selected
% track by calculating deltaF/F. First and last frames won't be N. 
validateattributes(intensityMultipleTolerance, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative', '<', 1}, ...
    mfilename, 'intensityMultipleTolerance');
eventTrackId = zeros(0, 1);
referenceFrame = zeros(0, 1);
referenceIntensityPhotons = zeros(0, 1);
eventFrame = zeros(0, 1);
eventIntensityPhotons = zeros(0, 1);
assignedMultiple = zeros(0, 1);
observedRatio = zeros(0, 1);
relativeDeviation = zeros(0, 1);

for longIndex = 1:height(longTracks)
    id = longTracks.track_id(longIndex);
    members = find(linkedLocalizations.track_id == id);
    [trackFrames, order] = sort(frames(members));
    trackPhotons = photons(members);
    trackPhotons = trackPhotons(order);

    interiorReferenceMask = trackFrames > min(trackFrames) & ...
        trackFrames < max(trackFrames);
    referenceIndices = find(interiorReferenceMask);

    for targetIndex = 1:numel(trackFrames)
        bestDeviation = inf;
        bestReferenceIndex = NaN;
        bestMultiple = NaN;
        bestRatio = NaN;

        for candidateReferenceIndex = referenceIndices.'
            referenceValue = trackPhotons(candidateReferenceIndex);
            targetValue = trackPhotons(targetIndex);
            if candidateReferenceIndex == targetIndex || ...
                    referenceValue <= 0 || targetValue <= referenceValue
                continue
            end

            ratio = targetValue / referenceValue;
            for multiple = 2:4
                deviation = abs(ratio - multiple) / multiple;
                if deviation <= intensityMultipleTolerance && ...
                        deviation < bestDeviation
                    bestDeviation = deviation;
                    bestReferenceIndex = candidateReferenceIndex;
                    bestMultiple = multiple;
                    bestRatio = ratio;
                end
            end
        end

        if isfinite(bestDeviation)
            eventTrackId(end + 1, 1) = id; 
            referenceFrame(end + 1, 1) = ...
                trackFrames(bestReferenceIndex); 
            referenceIntensityPhotons(end + 1, 1) = ...
                trackPhotons(bestReferenceIndex); 
            eventFrame(end + 1, 1) = trackFrames(targetIndex); 
            eventIntensityPhotons(end + 1, 1) = ...
                trackPhotons(targetIndex); 
            assignedMultiple(end + 1, 1) = bestMultiple; 
            observedRatio(end + 1, 1) = bestRatio; 
            relativeDeviation(end + 1, 1) = bestDeviation; 
        end
    end
end

intensityEvents = table(eventTrackId, referenceFrame, ...
    referenceIntensityPhotons, eventFrame, eventIntensityPhotons, ...
    assignedMultiple, observedRatio, relativeDeviation, ...
    'VariableNames', {'track_id', 'reference_frame', ...
    'reference_intensity_photons', 'event_frame', ...
    'event_intensity_photons', 'assigned_multiple', 'observed_ratio', ...
    'relative_deviation'});

maximumOnTime = max(onTimeFrames);
frameBins = (1:maximumOnTime).';
counts = accumarray(onTimeFrames, 1, [maximumOnTime, 1]);
onTimeCounts = table(frameBins, counts, ...
    'VariableNames', {'on_time_frames', 'number_of_linked_emitters'});

[inputFolder, inputStem] = fileparts(inputCsv);
if strlength(inputFolder) == 0
    inputFolder = ".";
end
outputStem = fullfile(inputFolder, inputStem + "_matlab");

% Plot exactly one emitter per subplot. Figures are paginated so that
% a large number of long tracks does not produce unreadably small axes.
validateattributes(tracksPerFigure, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'tracksPerFigure');
numberOfLongTracks = height(longTracks);
numberOfPages = ceil(numberOfLongTracks / tracksPerFigure);

for page = 1:numberOfPages
    firstTrackIndex = (page - 1) * tracksPerFigure + 1;
    lastTrackIndex = min(page * tracksPerFigure, numberOfLongTracks);
    tracksOnPage = lastTrackIndex - firstTrackIndex + 1;
    numberOfColumns = min(4, tracksOnPage);
    numberOfRows = ceil(tracksOnPage / numberOfColumns);

    intensityFigure = figure('Color', 'white', ...
        'Name', sprintf('Intensity for selected tracks, page %d', page), ...
        'Position', [80, 60, 1400, 850]);
    layout = tiledlayout(intensityFigure, numberOfRows, numberOfColumns, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for longIndex = firstTrackIndex:lastTrackIndex
        id = longTracks.track_id(longIndex);
        members = find(linkedLocalizations.track_id == id);
        [trackFrames, order] = sort(frames(members));
        trackPhotons = photons(members);
        trackPhotons = trackPhotons(order);

        trackAxes = nexttile(layout);
        hold(trackAxes, 'on');
        plot(trackAxes, trackFrames, trackPhotons, '-o', ...
            'Color', [0.16, 0.47, 0.71], 'LineWidth', 1.2, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.16, 0.47, 0.71]);
        ylim([100 7000]);

        trackEvents = intensityEvents(intensityEvents.track_id == id, :);

        % Label every distinct interior reference N once, even if it
        % supports more than one 2N/3N/4N assignment.
        if ~isempty(trackEvents)
            referencePairs = [trackEvents.reference_frame, ...
                trackEvents.reference_intensity_photons];
            [~, uniqueReferenceRows] = unique(referencePairs, 'rows', 'stable');
            for referenceIndex = uniqueReferenceRows.'
                referenceX = trackEvents.reference_frame(referenceIndex);
                referenceY = ...
                    trackEvents.reference_intensity_photons(referenceIndex);
                plot(trackAxes, referenceX, referenceY, 's', ...
                    'Color', [0.05, 0.55, 0.25], 'MarkerSize', 7, ...
                    'LineWidth', 1.3);
                text(trackAxes, referenceX, referenceY, ...
                    sprintf('N: %.0f  ', referenceY), ...
                    'Color', [0.03, 0.42, 0.18], 'FontSize', 8, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'Interpreter', 'none');
            end
        end

        for eventIndex = 1:height(trackEvents)
            eventX = trackEvents.event_frame(eventIndex);
            eventY = trackEvents.event_intensity_photons(eventIndex);
            multiple = trackEvents.assigned_multiple(eventIndex);
            plot(trackAxes, eventX, eventY, 'o', ...
                'Color', [0.80, 0.10, 0.10], 'MarkerSize', 7, ...
                'LineWidth', 1.3);
            text(trackAxes, eventX, eventY, ...
                sprintf('  %dN: %.0f', multiple, eventY), ...
                'Color', [0.70, 0.05, 0.05], 'FontSize', 8, ...
                'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
                'Interpreter', 'none');
        end
        hold(trackAxes, 'off');

        xlabel(trackAxes, 'Frame');
        ylabel(trackAxes, 'Photons / localization');
        title(trackAxes, sprintf('Emitter %d (%d frames)', ...
            id, longTracks.on_time_frames(longIndex)));
%             grid(trackAxes, 'on');
        box(trackAxes, 'off');
        % if numel(trackFrames) <= 20
        %     xticks(trackAxes, unique(trackFrames));
        % end
    end

    title(layout, sprintf(['Intensity traces for selected tracks: page %d of %d ', ...
        '(one emitter per subplot)'], page, numberOfPages));
    pageStem = outputStem + "_emitters_4_to_10_frames_intensity_page_" + ...
        compose("%02d", page);

end

fprintf('Input localizations: %d\n', numberOfLocalizations);
fprintf('Linked tracks: %d\n', numberOfTracks);
fprintf('Tracks lasting more than one frame: %d\n', sum(onTimeFrames > 1));
fprintf('Tracks selected for intensity analysis: %d\n', height(longTracks));
fprintf('Candidate 2N/3N/4N intensity points: %d\n', height(intensityEvents));
fprintf('Selected tracks containing candidate intensity events: %d\n', ...
    numel(unique(intensityEvents.track_id)));
fprintf('Tracks flagged with duplicate detections in one frame: %d\n', ...
    sum(duplicateFrameFlag));
fprintf('Maximum on-time: %d frames\n', maximumOnTime);


% %% Mean intensity of every selected long emitter
% 
% % Calculate one arithmetic mean intensity per selected emitter using all
% % localization intensities in that emitter's track.
% meanIntensityPhotons = zeros(height(longTracks), 1);
% 
% for longIndex = 1:height(longTracks)
%     id = longTracks.track_id(longIndex);
% 
%     meanIntensityPhotons(longIndex) = mean( ...
%         photons(linkedLocalizations.track_id == id));
% end
% 
% % Add the mean intensity to the existing longTracks table.
% longTracks.mean_intensity_photons = meanIntensityPhotons;
% 
% % Create a compact output table containing one row per selected emitter.
% longEmitterMeanIntensities = table( ...
%     longTracks.track_id, ...
%     meanIntensityPhotons, ...
%     'VariableNames', { ...
%         'track_id', ...
%         'mean_intensity_photons'});
% 
% % %% Mean intensity of emitters exhibiting intensity stepping
% % 
% % % A stepping emitter has at least one detected 2N, 3N, or 4N event.
% % steppingTrackIds = unique(intensityEvents.track_id, 'stable');
% % 
% % steppingMeanIntensityPhotons = zeros(numel(steppingTrackIds), 1);
% % numberOfSteppingEvents = zeros(numel(steppingTrackIds), 1);
% % 
% % for steppingIndex = 1:numel(steppingTrackIds)
% %     id = steppingTrackIds(steppingIndex);
% % 
% %     % Mean intensity calculated from every localization in this track.
% %     steppingMeanIntensityPhotons(steppingIndex) = mean( ...
% %         photons(linkedLocalizations.track_id == id));
% % 
% %     % Number of detected 2N/3N/4N points in this track.
% %     numberOfSteppingEvents(steppingIndex) = sum( ...
% %         intensityEvents.track_id == id);
% % end
% % 
% % % Create one output row per stepping emitter.
% % steppingEmitterMeanIntensities = table( ...
% %     steppingTrackIds, ...
% %     steppingMeanIntensityPhotons, ...
% %     numberOfSteppingEvents, ...
% %     'VariableNames', { ...
% %         'track_id', ...
% %         'mean_intensity_photons', ...
% %         'number_of_detected_stepping_events'});
% % 
% % writetable(longEmitterMeanIntensities, ...
% %     outputStem + "_selected_emitter_mean_intensities.csv");
% % 
% % writetable(steppingEmitterMeanIntensities, ...
% %     outputStem + "_stepping_emitter_mean_intensities.csv");
