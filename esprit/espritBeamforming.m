function espritBeamforming()

    % ------------------------------------------------
    %             ESPRIT BEAMFORMING DEMO
    % Simulation of several sources around the array
    %
    %          Jose Ignacio Dominguez Simon
    %
    %            Array Signal Processing
    %           Aalborg University - 2015
    % ------------------------------------------------
    %
    % Edit the setup parameters in this file before
    % running it.
    %
    % Joe.

    % ---------------------- SIMULATION PARAMETERS -----------------------

    AMOUNT_OF_SENSORS = 35;
    SENSORS_POSITION_ORIGIN = [0, 0];
    SENSORS_DISPLACEMENT_VECTOR = [0.05, 0]; % <-- if the separation between sensors is too big, you get spatial aliasing unless you reduce the frequency of the sources
    SENSOR_NOISE = -60;

    PROPAGATION_SPEED = 340;
    SAMPLE_RATE = 44100;

    AMOUNT_OF_SOURCES = 3;
    inputSignal_length = 20000;

    sourcesAngles_origin = [0 45 90];
    sourcesAngles_speed = [1 1 1];

    SHOW_SOURCE_SIGNALS = 0;

%     freq = (PROPAGATION_SPEED/(2 * norm(SENSORS_DISPLACEMENT_VECTOR))) /1.7;
    freq = 1133.3333;
    inputSignal_freqs = [freq, freq/1.47, freq/1.99];
    inputSignal_amplitudes = [1, 1, 1];


    % --------------------------- MAIN PROGRAM ---------------------------

    inputSignals = zeros(AMOUNT_OF_SOURCES, inputSignal_length);
    sensorCoordinates = zeros(AMOUNT_OF_SENSORS, 2);

    disp(sprintf(' ------------------------------------------------\n             ESPRIT BEAMFORMING DEMO\n\n          Jose Ignacio Dominguez Simon\n                Tobias Van Baarsel\n\n            Array Signal Processing\n           Aalborg University - 2015\n ------------------------------------------------\n'));

    figureHandler = figure('name', 'ESPRIT algorithm interactive demo - I. D. Simon, T. V. Baarsel.', 'NumberTitle','off', 'menubar', 'none');
    frameCounter = 1;

    updateSourceSignals();
    updateSensorsPlacement();
    currentAngles = sourcesAngles_origin(1 : AMOUNT_OF_SOURCES);

    % Find amount of needed iterations (based on the sources' initial positions)
    amountOfIterations = min((180 - currentAngles(1 : AMOUNT_OF_SOURCES)) ./ sourcesAngles_speed(1 : AMOUNT_OF_SOURCES));

    % Initialize all data variables
    estimatedAngles = zeros(AMOUNT_OF_SOURCES, amountOfIterations);
    estimationErrors = zeros(2, 1, AMOUNT_OF_SOURCES);
    eigenvectors = zeros(AMOUNT_OF_SENSORS, AMOUNT_OF_SOURCES, amountOfIterations);
    C0s = zeros(AMOUNT_OF_SOURCES, AMOUNT_OF_SOURCES, amountOfIterations);
    covariances = zeros(AMOUNT_OF_SENSORS, AMOUNT_OF_SENSORS, amountOfIterations);

    % Initialize structure for all origin data
    dataToProcess = struct('AMOUNT_OF_SENSORS', AMOUNT_OF_SENSORS, ...
                           'SENSORS_POSITION_ORIGIN', SENSORS_POSITION_ORIGIN, ...
                           'SENSORS_DISPLACEMENT_VECTOR', SENSORS_DISPLACEMENT_VECTOR, ...
                           'SENSOR_NOISE', SENSOR_NOISE, ...
                           'PROPAGATION_SPEED', PROPAGATION_SPEED, ...
                           'SAMPLE_RATE', SAMPLE_RATE, ...
                           'AMOUNT_OF_SOURCES', AMOUNT_OF_SOURCES, ...
                           'SOURCE_ANGLES', currentAngles(1 : AMOUNT_OF_SOURCES), ...
                           'inputSignals', inputSignals(1 : AMOUNT_OF_SOURCES, :), ...
                           'inputSignal_freqs', inputSignal_freqs(1 : AMOUNT_OF_SOURCES));

    disp(sprintf('Processing data from setup. Please, wait ...'));
    for currentIteration = 1 : amountOfIterations

        % Update new angle(s) to process
        currentAngles = currentAngles + sourcesAngles_speed(1 : AMOUNT_OF_SOURCES);
        dataToProcess.SOURCE_ANGLES = currentAngles(1 : AMOUNT_OF_SOURCES);

        % Process the current angle(s)
        processedData = processData(dataToProcess);

        % Accumulate results from processing current angle(s)
        estimatedAngles(:, currentIteration) = processedData.estimatedAngles;
        estimationErrors(1, currentIteration, :) = currentAngles;
        estimationErrors(2, currentIteration, :) = processedData.currentErrors;
        eigenvectors(:, :, currentIteration) = processedData.U;
        C0s(:, :, currentIteration) = processedData.C0;
        covariances(:, :, currentIteration) = processedData.covarianceMatrix;

        % Adds additional fields to the structure, needed to plot the data
        processedData.sensorCoordinates = sensorCoordinates;
        processedData.AMOUNT_OF_SOURCES = AMOUNT_OF_SOURCES;
        processedData.AMOUNT_OF_SENSORS = AMOUNT_OF_SENSORS;
        processedData.SOURCE_ANGLES = currentAngles;
        processedData.SENSORS_DISPLACEMENT_VECTOR = SENSORS_DISPLACEMENT_VECTOR;
        processedData.SENSOR_NOISE = SENSOR_NOISE;
        processedData.inputSignal_freqs = inputSignal_freqs;
        processedData.estimationError = estimationErrors;
        processedData.U1 = processedData.U(1 : size(processedData.U, 1) - 1, :);
        processedData.U2 = processedData.U(2 : size(processedData.U, 1), :);
        processedData.inputSignals = inputSignals;
        processedData.SAMPLE_RATE = SAMPLE_RATE;
        processedData.PROPAGATION_SPEED = PROPAGATION_SPEED;
        processedData.plotLimit = currentIteration;
        processedData.sourcesAngles_speed = sourcesAngles_speed;

        % Show current iteration plots
        if ~updatePlot(figureHandler, processedData)
            return
        end

        % Save current frame for the animation
        currentFrame = getframe(figureHandler);
        if frameCounter == 1 
            [imageToSave, colorsMap] = rgb2ind(currentFrame.cdata, 256, 'nodither');
        else
            [imageToSave(:, :, 1, frameCounter), colorsMap] = rgb2ind(currentFrame.cdata, colorsMap, 'nodither');
        end
        frameCounter = frameCounter + 1;
    end

    % Saves animation to file
    disp('Saving GIF ...')
    fileNumber = 1;
    filenameToSave = sprintf('simulation_animation_%d.gif', fileNumber);
    while exist(filenameToSave, 'file')
        fileNumber = fileNumber + 1;
        filenameToSave = sprintf('simulation_animation_%d.gif', fileNumber);
    end
    imwrite(imageToSave, colorsMap, filenameToSave, 'DelayTime', 0, 'LoopCount', 0);
    disp(sprintf('Saved to "%s"', filenameToSave));
    
    disp(sprintf('Initial processing done! Use left and right arrows to navigate ...\n'))

    % Wires the callback for the keyboard presses
    if ~isempty(figureHandler) && ~ishandle(figureHandler)
        disp('[main function] Closed window. Terminating.')
        return
    end
    set(figureHandler,'KeyPressFcn',{@lineCallback});
    function lineCallback(~, second)
        shouldUpdatePlot = 0;
        if strcmp(second.Key, 'leftarrow')
            if processedData.plotLimit > 1
                processedData.plotLimit = processedData.plotLimit -1;
                shouldUpdatePlot = 1;
            end
        end
        if strcmp(second.Key, 'rightarrow')
            if processedData.plotLimit < amountOfIterations
                processedData.plotLimit = processedData.plotLimit +1;
                shouldUpdatePlot = 1;
            end
        end

        if shouldUpdatePlot
            % Update data for the plot
            processedData.SOURCE_ANGLES = sourcesAngles_origin(1 : AMOUNT_OF_SOURCES) + processedData.plotLimit * sourcesAngles_speed(1 : AMOUNT_OF_SOURCES);
            processedData.estimatedAngles = estimatedAngles(:, processedData.plotLimit);
            processedData.covarianceMatrix = covariances(:, :, processedData.plotLimit);
            processedData.U = eigenvectors(:, :, processedData.plotLimit);
            processedData.U1 = processedData.U(1 : size(processedData.U, 1) - 1, :);
            processedData.U2 = processedData.U(2 : size(processedData.U, 1), :);
            processedData.C0 = C0s(:, :, processedData.plotLimit);

            disp(sprintf('Loading data from iteration %d of %d ...', processedData.plotLimit, amountOfIterations));

            % Refresh plot window
            updatePlot(figureHandler, processedData);
        end
    end

    % ----------------------- INTERNAL FUNCTIONS -----------------------

    function updateSourceSignals()

        inputSignals = zeros(AMOUNT_OF_SOURCES, inputSignal_length);
        % Generate a sine signal with a hamming window        
        rng(0);
        inputSignals(1,:) = sin(2 * pi * ...
            inputSignal_freqs(1) * [1 : inputSignal_length] / SAMPLE_RATE) ...
            .* (hamming(inputSignal_length)'.^2) ...
            * inputSignal_amplitudes(1);
        rng(0);
        inputSignals(2,:) = sin(2 * pi * ...
            inputSignal_freqs(2) * [1 : inputSignal_length] / SAMPLE_RATE) ...
            .* (hamming(inputSignal_length)'.^2) ...
            * inputSignal_amplitudes(2);
        rng(0);
        inputSignals(3,:) = sin(2 * pi * ...
            inputSignal_freqs(3) * [1 : inputSignal_length] / SAMPLE_RATE) ...
            .* (hamming(inputSignal_length)'.^2) ...
            * inputSignal_amplitudes(3);

        % ------- PLOT (AND SAVE IMAGE) OF SOURCE SIGNALS WITH NOISE -------

        if SHOW_SOURCE_SIGNALS
            SENSOR_NOISE = 3;

            figureHandler = figure();
            set(figureHandler, 'Position', [500 100 1000 800])

            timeAxis = [1 : length(inputSignals(1,:))] / SAMPLE_RATE;
            for i = 1 : 3
                subplot(3,2,(2*(i-1) + 1))
                plot(timeAxis, inputSignals(i,:) + randn(1, inputSignal_length) * 10^(SENSOR_NOISE/20));
                ylim([-3 3])
                xlim([0 max(timeAxis)])
                grid
                xlabel('Time [s]', 'FontSize', 14);
                ylabel('Amplitude [.]', 'FontSize', 14);
                set(gca, 'FontSize', 14);
                legend({sprintf('Signal %d: %.1f [Hz]', i, inputSignal_freqs(i))})
            end
            subplot(3,2,1)
            title(sprintf('Source signals.\nAdditive white noise at %.1f [dB]', SENSOR_NOISE), 'FontSize', 14)

            timeRange = [2500 : 4000];
            timeAxis = [1 : length(inputSignals(1,:))] / SAMPLE_RATE;
            for i = 1 : 3
                subplot(3, 2, (2 * (i)))
                plot(timeAxis(timeRange), inputSignals(i, timeRange) + randn(1, length(timeRange)) * 10^(SENSOR_NOISE/20));
                ylim([-3 3])
                xlim([min(timeRange) max(timeRange)]/SAMPLE_RATE)
                grid
                xlabel('Time [s]', 'FontSize', 14);
                ylabel('Amplitude [.]', 'FontSize', 14);
                set(gca, 'FontSize', 14);
                legend({sprintf('Signal %d: %.1f [Hz]', i, inputSignal_freqs(i))})
            end
            subplot(3,2,2)
            title(sprintf('Source signals (zoom).\nAdditive white noise at %.1f [dB]', SENSOR_NOISE), 'FontSize', 14)

            currentFrame = getframe(figureHandler);
            [imageToSave, colorsMap] = rgb2ind(currentFrame.cdata, 256, 'nodither');
            imwrite(imageToSave, colorsMap, sprintf('source_signals_noise_%d.gif', SENSOR_NOISE), 'DelayTime', 0, 'LoopCount', 0);
        end
    end

    % ------------------------ SENSOR PLACEMENT ------------------------

    function updateSensorsPlacement()

        % Places the sensors according the array geometry
        for index = 1 : AMOUNT_OF_SENSORS

            % Uniform linear array
            currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
                (index - 1) * SENSORS_DISPLACEMENT_VECTOR;

            sensorCoordinates(index, :) = currentSensorCoordinates;
        end
    end

end

function returnedResult = updatePlot(figureHandler, processedData)

    VERTICAL_ERROR_RANGE = [-100 100];

    signalColors = [0.5 0.1 0.2
                    0.1 0.2 0.5
                    0.2 0.5 0.1];

    % Checks if the plot window still exists
    if ~isempty(figureHandler) && ~ishandle(figureHandler)
        disp('[updatePlot() function] Closed window. Terminating.')
        returnedResult = 0;
        return
    end

    % Cleans window and sets size
    clf
    windowPosition = get(figureHandler, 'Position');
    set(figureHandler, 'Position', [windowPosition(1) windowPosition(2) 1400 600])

    % --------------- SOURCES POSITION AND ESTIMATIONS ---------------

    subplot(2,3,1)

    % Trick to display the legend as desired
    plot([10 11], [10 11], '--', 'LineWidth', 4, 'Color', signalColors(1, :))
    hold on
    plot([10 11], [10 11], 'LineWidth', 12, 'Color', signalColors(1, :) + 0.5)
    legend({'Estimated direction', 'Real direction'}, 'Location', 'southeast', 'FontSize', 12)

    % Plot the half circumference
    circleFunction = [cos([0 : 180] / 180 * pi);
                      sin([0 : 180] / 180 * pi)];
    plot(circleFunction(1,:), circleFunction(2,:), 'LineWidth', 2);

    xlim([-1.3 1.3])
    ylim([-0.4 1.2])
    grid

    % Plot the array sensors
    for i = 1 : processedData.AMOUNT_OF_SENSORS
        plot(processedData.sensorCoordinates(i,1) / 2, processedData.sensorCoordinates(i,2) / 2, 'o', 'Color', [0.2 0.5 0.8]);
    end

    % Plot the source(s) line(s)
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        x_actual = cos(processedData.SOURCE_ANGLES(i) / 180 * pi) * 0.95;
        y_actual = sin(processedData.SOURCE_ANGLES(i) / 180 * pi) * 0.95;
        plot([0 x_actual], [0 y_actual], 'LineWidth', 12, 'Color', signalColors(i, :) + 0.5)

        plot([x_actual / 0.95], [y_actual / 0.95], 'x', 'MarkerSize', 12, 'LineWidth', 8, 'Color', [0.2 .2 .8])
        plot([x_actual / 0.95], [y_actual / 0.95], 'x', 'MarkerSize', 8, 'LineWidth', 3, 'Color', [0.8 .2 .3])

        text(x_actual * 1.2, y_actual * 1.2, sprintf('%.0f', processedData.SOURCE_ANGLES(i)), 'FontSize', 14, 'HorizontalAlign', 'center');
    end

    % Plot the estimated line(s)
    for i = 1 : processedData.AMOUNT_OF_SOURCES

        if imag(processedData.estimatedAngles(i)) == 0
            x_estimated = cos(processedData.estimatedAngles(i)) * 0.95;
            y_estimated = sin(processedData.estimatedAngles(i)) * 0.95;
            plot([0 x_estimated], [0 y_estimated], '--', 'LineWidth', 4, 'Color', signalColors(i, :))
        else
            x_estimated = cos(real(processedData.estimatedAngles(i))) * 0.95;
            y_estimated = sin(real(processedData.estimatedAngles(i))) * 0.95;
            plot([0 x_estimated], [0 y_estimated], '--', 'LineWidth', 4, 'Color', signalColors(i, :))

            text(x_estimated * 0.5 + 0.05, y_estimated * 0.5 + 0.05, '  Wrong!', 'FontSize', 14, 'Color', [0.8 .2 .3])
%             pause(1)
        end
    end

%     text(1.04, 1.1, sprintf('Amount of sensors: %d', processedData.AMOUNT_OF_SENSORS), 'FontSize', 12)
%     text(1.04, 1, sprintf('Noise: %.1f [dB]', processedData.SENSOR_NOISE), 'FontSize', 12)
%     for i = 1 : processedData.AMOUNT_OF_SOURCES
%         text(1.04, 1-(i *0.1), sprintf('Freq. #%d: %.1f [Hz]', i, processedData.inputSignal_freqs(i)), 'FontSize', 12);
%     end

    set(gca, 'FontSize', 11)
    set(gca,'XTickLabel', {});
    set(gca,'YTickLabel', {});
    title(sprintf('Rotating source detected angle'), 'FontSize', 14)


    % -------------------- ANGLE ESTIMATION ERROR --------------------

    subplot(2,3,4)

    for i = 1 : processedData.AMOUNT_OF_SOURCES
        plot(processedData.estimationError(1, :, i), real(processedData.estimationError(2, :, i)), 'LineWidth', 2, 'Color', signalColors(i, :));
        hold on
    end
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        plot(processedData.estimationError(1, processedData.plotLimit, i), real(processedData.estimationError(2, processedData.plotLimit, i)), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', signalColors(i, :));
    end

    xlim([0 180])
    ylim(VERTICAL_ERROR_RANGE)
    set(gca, 'FontSize', 11)        
    xlabel('Actual source angle [degrees]', 'FontSize', 14)
    ylabel('Estimation error [degrees]', 'FontSize', 14)
    grid
    title('Angle estimation error', 'FontSize', 14)


    % ------------------- ARRAY COVARIANCE MATRIX --------------------

    subplot(2,3,2)

    realPart = real(processedData.covarianceMatrix);
    imagPart = imag(processedData.covarianceMatrix);
    anglePart = angle(processedData.covarianceMatrix);

    realPart = (realPart - min(realPart(:)));
    realPart = realPart / max(realPart(:));

    imagPart = (imagPart - min(imagPart(:)));
    imagPart = imagPart / max(imagPart(:));

    anglePart = (anglePart - min(anglePart(:)));
    anglePart = anglePart / max(anglePart(:));

    imagesc([realPart imagPart anglePart ])
    set(gca,'YDir','normal')

    axisValues = [[1 : processedData.AMOUNT_OF_SENSORS] [1 : processedData.AMOUNT_OF_SENSORS] [1 : processedData.AMOUNT_OF_SENSORS]];
    valuesXaxis = num2cell(axisValues);
    ax = gca;
    ax.XTickLabel = valuesXaxis(ax.XTick);
    hold on

    plot([0 processedData.AMOUNT_OF_SENSORS], [0 processedData.AMOUNT_OF_SENSORS], 'LineWidth', 1.5, 'Color', [1 0 0])
    plot([processedData.AMOUNT_OF_SENSORS processedData.AMOUNT_OF_SENSORS*2], [0 processedData.AMOUNT_OF_SENSORS], 'LineWidth', 1.5, 'Color', [1 0 0])
    plot([processedData.AMOUNT_OF_SENSORS*2 processedData.AMOUNT_OF_SENSORS*3], [0 processedData.AMOUNT_OF_SENSORS], 'LineWidth', 1.5, 'Color', [1 0 0])

    % Plot separation lines
    plot([processedData.AMOUNT_OF_SENSORS processedData.AMOUNT_OF_SENSORS], [0 processedData.AMOUNT_OF_SENSORS], 'LineWidth', 2, 'Color', [1 1 1]);
    plot([processedData.AMOUNT_OF_SENSORS processedData.AMOUNT_OF_SENSORS] * 2, [0 processedData.AMOUNT_OF_SENSORS], 'LineWidth', 2, 'Color', [1 1 1]);

    set(gca, 'FontSize', 11)
    xlabel('Sensor number', 'FontSize', 14)
    ylabel('Sensor number', 'FontSize', 14)
    title(sprintf('Array covariance matrix\nReal part                    Imag part                       Angle'), 'FontSize', 14)


    % -------------- EIGENVECTORS OF COVARIANCE MATRIX ---------------

    subplot(2,3,5)

    eigenvectorsMatrix = angle([processedData.U2, fliplr(processedData.U1 * processedData.C0)]);
    imagesc(eigenvectorsMatrix)
    set(gca,'YDir','normal')
    hold on

    cells = get(gca,'XTickLabel');
    axisNumbers = zeros(length(cells), 1);
    for i = 1 : length(cells)
        axisNumbers(i) = str2double(cell2mat(cells(i)));
    end
    center_x_pos = (max(axisNumbers) + min(axisNumbers))/2;
    plot([center_x_pos center_x_pos], [0 size(eigenvectorsMatrix, 1)], 'Color', [1 1 1]);

    set(gca, 'FontSize', 11)
    ylabel('Eigenvector length = Sources count', 'FontSize', 14)
    set(gca,'XTickLabel', {});
    title(sprintf('Computed eigenvectors checking\n   |<------------------ U -----------------> | <----------------- U'' ----------------->|  '), 'FontSize', 14)


    % ----------------------- SOURCE SIGNALS -------------------------
    subplot(2,3,6)

    timeAxis = [1 : length(processedData.inputSignals(1,:))] / processedData.SAMPLE_RATE;
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        rng(0);
        plot(timeAxis, processedData.inputSignals(i,:) + randn(1, length(processedData.inputSignals)) * 10^(processedData.SENSOR_NOISE/20) - 2*(i - 1));
        hold on
    end
    ylim([-1.5-2*(i-1) 1.5])
    xlim([0 max(timeAxis)])
    grid
    xlabel('Time [s]', 'FontSize', 14);
    ylabel('Amplitude [.]', 'FontSize', 14);
    set(gca, 'FontSize', 11);

    legendCells = cell(processedData.AMOUNT_OF_SOURCES, 1);
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        legendCells{i} = sprintf('Signal %d: %.1f [Hz]', i, processedData.inputSignal_freqs(i));
    end
    legend(legendCells, 'FontSize', 12)

    % Places plots in their desired position
    subplot(2,3,1)
    set(gca,'position', [0.03    0.55    0.3    0.4])
    subplot(2,3,2)
    set(gca,'position', [0.37    0.55    0.3    0.37])
    subplot(2,3,4)
    set(gca,'position', [0.04    0.06    0.29    0.42])
    subplot(2,3,5)
    set(gca,'position', [0.37    0.06    0.3    0.355])
    subplot(2,3,6)
    set(gca,'position', [0.71    0.06    0.28    0.355])

    % Displays text with details of setup
    textToShow = sprintf('                                     Current setup details:');
    textToShow = strcat(textToShow, sprintf('\n\n Amount of sources:     %d', processedData.AMOUNT_OF_SOURCES));
    textToShow = sprintf('%s\n Amount of sensors:     %d\n', textToShow, processedData.AMOUNT_OF_SENSORS);
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        textToShow = sprintf('%s\n Source #%d angle:         %d [deg]', textToShow, i, processedData.SOURCE_ANGLES(i));
        textToShow = sprintf('%s\n Estimated angle #%d:     %.0f [deg]', textToShow, i, processedData.estimatedAngles(i) * 180 / pi);
    end
    textToShow = sprintf('%s\n\n Sensors distance:        %.2f [m]', textToShow, norm(processedData.SENSORS_DISPLACEMENT_VECTOR));
    textToShow = strcat(textToShow, sprintf('\n Critical frequency:        %.1f [Hz]', processedData.PROPAGATION_SPEED / norm(processedData.SENSORS_DISPLACEMENT_VECTOR) /2));
    textToShow = strcat(textToShow, sprintf('\n Propagation speed:      %.f [m/s]', processedData.PROPAGATION_SPEED));
    textToShow = strcat(textToShow, sprintf('\n Sample rate:                %.f [S/s]', processedData.SAMPLE_RATE));
    textToShow = sprintf('%s\n Additive noise level:     %.0f [dB]\n', textToShow, processedData.SENSOR_NOISE);
%     textToShow = sprintf('%s\n Covariance matrix rank:                    %d', textToShow, rank(processedData.covarianceMatrix));
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        textToShow = sprintf('%s\n Source %d ang. speed:   %.0f [deg/iteration]', textToShow, i, processedData.sourcesAngles_speed(i));
    end
    textToShow = sprintf('%s\n', textToShow);
    for i = 1 : processedData.AMOUNT_OF_SOURCES
        textToShow = sprintf('%s\n Signal %d frequency:      %.1f [Hz]', textToShow, i, processedData.inputSignal_freqs(i));
    end
    uicontrol('Style','text', 'Position',[975 280 400 300], 'String', textToShow, 'FontSize', 13, 'HorizontalAlignment', 'left');

    % Forces drawing of window
    drawnow
    returnedResult = 1;

end

function returnedData = processData(dataToProcess)

    % Create sensors' "recordings"
    sensorSignals = zeros(dataToProcess.AMOUNT_OF_SENSORS, length(dataToProcess.inputSignals));

    % Recreate the Vandermonde vector with associated delays
    for currentSource = 1 : dataToProcess.AMOUNT_OF_SOURCES

        % Construct the new Vandermonde vector
        if dataToProcess.SOURCE_ANGLES(currentSource) == 90
            theta = 0;
        else
            theta = dataToProcess.SOURCE_ANGLES(currentSource) / 180 * pi;
        end
        lambda = dataToProcess.PROPAGATION_SPEED / dataToProcess.inputSignal_freqs(currentSource);
        xsi = norm(dataToProcess.SENSORS_DISPLACEMENT_VECTOR) / lambda * cos(theta);
        vandermonde = exp(1i * 2 * pi * xsi * [0 : dataToProcess.AMOUNT_OF_SENSORS-1]');

        for currentSensor = 1 : dataToProcess.AMOUNT_OF_SENSORS
            sensorSignals(currentSensor, :) = sensorSignals(currentSensor, :) + ...
                dataToProcess.inputSignals(currentSource, :) * vandermonde(currentSensor) + ...
                randn(1, length(dataToProcess.inputSignals)) * 10^(dataToProcess.SENSOR_NOISE/20);

        end
    end

    % Estimation of covariance matrix from sensed signals
    covarianceMatrix = zeros(dataToProcess.AMOUNT_OF_SENSORS);
    timeRange = [1, length(dataToProcess.inputSignals)];

    for t = timeRange(1) : timeRange(2)

        covarianceMatrix = covarianceMatrix + ...
            sensorSignals(:, t) * sensorSignals(:, t)';
    end
    covarianceMatrix = covarianceMatrix / (timeRange(2) - timeRange(1) + 1);

    % Obtain eigenvectors from the covariance matrix -> U
    [eigenvectors, ~, flags] = eigs(rot90(covarianceMatrix), dataToProcess.AMOUNT_OF_SOURCES);    
    if flags
        disp('[Warning] Not all eigenvectors of the covariance matrix converged!');
    end

    % Form U1 and U2 from U
    U = eigenvectors;
    U1 = U(1:size(U,1)-1, :);
    U2 = U(2:size(U,1), :);

    % Solve C for U2 = U1 * C
    C0 = inv(corr(U1)) * corr(U1, U2); % <-- this seems to be the right order for the: corr(U1, U2)

    if ~any(isnan(C0(:))) && ~any(isinf(C0(:)))
        % Compute eta, the eigenvalues of C (as many as sources)
        [~, eigenvalues] = eigs(C0);
        eta = diag(eigenvalues);

        % Find the estimated angles
        xsi = angle(eta) / (2 * pi);
        lambdas = dataToProcess.PROPAGATION_SPEED ./ dataToProcess.inputSignal_freqs(1 : dataToProcess.AMOUNT_OF_SOURCES)';
        estimatedAngles = sort(acos(lambdas ./ norm(dataToProcess.SENSORS_DISPLACEMENT_VECTOR) .* xsi));

        % Compute the error for the current source(s) angle(s)
        currentErrors = (estimatedAngles - (dataToProcess.SOURCE_ANGLES(1 : dataToProcess.AMOUNT_OF_SOURCES) / 180 * pi)') / pi * 180;

    else
        disp('TODO: finish this here ...');
    end

    returnedData = struct('estimatedAngles', estimatedAngles, ...
                          'currentErrors', currentErrors, ...
                          'U', U, ...
                          'C0', C0, ...
                          'covarianceMatrix', covarianceMatrix);
end
