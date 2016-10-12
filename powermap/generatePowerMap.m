function generatePowerMap()
    % Generates the amplitude / power map for a set of sources with given characteristics.
    % Saves results to "amplitudeMap.mat".
    %
    % Joe.

    clc
    disp(sprintf('\n |   Array Signal Processing - ESPRIT algorithm demo   |'));
    disp(sprintf(' |       Acoustics and Audio Technology - 2015         |'));
    disp(sprintf(' |                    (Power map)                      |'));
    disp(sprintf(' |                                                     |'));
    disp(sprintf(' |  Jose Ignacio Dominguez Simon - Tobias Van Baarsel  |\n\n'));

    % ----------- Obtain excitation signal(s) and source positions -----------

    SCALING_FACTOR = 1;
    SAMPLE_RATE = 48000 * SCALING_FACTOR;
    PROPAGATION_SPEED = 340 * SCALING_FACTOR;

    amountOfSignals = 3;
    inputSignal_length = 10000 * SCALING_FACTOR;
    inputSignal_freqs = [1000, 1000, 1000] * SCALING_FACTOR;
    inputSignal_coordinates = [2, 5
                               8, 5
                               2, 11] * SCALING_FACTOR;
    inputSignal_amplitudes = [1, 0.5, 0.5];

    boundaries = [ 5, 1, 5, 8;
                  -5, 8, 5, 8] * SCALING_FACTOR;

    x_limits = [-5, 10] * SCALING_FACTOR;
    y_limits = [-5, 13] * SCALING_FACTOR;
    GRID_RESOLUTION = 0.5 * SCALING_FACTOR;

    % Generate sinus with hamming window
    inputSignals = zeros(amountOfSignals, inputSignal_length);
    for i = 1 : amountOfSignals
        inputSignals(i, 1 : inputSignal_length) = sin(2 * pi * ...
            inputSignal_freqs(i) * [1 : inputSignal_length] / SAMPLE_RATE) ...
            .* (hamming(length(inputSignals))'.^2) ...
            * inputSignal_amplitudes(i);
    end

    % Map input signals with coordinates to the "sources-map"
    sourcesMap = zeros(amountOfSignals, inputSignal_length + 2);
    for i = 1 : amountOfSignals
        sourcesMap(i, 1) = inputSignal_coordinates(i, 1);
        sourcesMap(i, 2) = inputSignal_coordinates(i, 2);
        sourcesMap(i, 3:length(sourcesMap)) = inputSignals(i, :);
    end

    % -- Generated sources map table:
    % A row for each source, with its x,y coordinates and the signal on time
    %
    %   x   y   signal_transmitted
    %   _   _   __________________ ...
    %   _   _   __________________ ...
    %   _   _   __________________ ...


    % ------------- Test: generation of the amplitude map -------------

tic

    x_coordinates = [x_limits(1) : GRID_RESOLUTION : x_limits(2)];
    y_coordinates = [y_limits(1) : GRID_RESOLUTION : y_limits(2)];

    amplitudeMap = zeros(length(y_coordinates), length(x_coordinates));

    disp(sprintf(' > Creating map\n    Processing row: %3d of %3d.', 0, length(x_coordinates)));
    for x_index = 1 : length(x_coordinates)

        disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b%3d of %3d.', x_index, length(x_coordinates)));

        for y_index = 1 : length(y_coordinates)

            accumulatedSignal = zeros(1, inputSignal_length);

            for currentSource = 1 : amountOfSignals
                % Obtain distance from current sensor to each source
                currentSensorCoordinates = [x_coordinates(x_index), ...
                                            y_coordinates(y_index)];
                currentDistance = sqrt( ...
                    (currentSensorCoordinates(1) - sourcesMap(currentSource, 1))^2 + ...
                    (currentSensorCoordinates(2) - sourcesMap(currentSource, 2))^2 );

                % Delay signal from that transmitter according to that distance
                currentDelay = currentDistance / PROPAGATION_SPEED;
                fractionalDelay = generateFractionalDelayIR(currentDelay, SAMPLE_RATE);
                newReceivedSignal = conv(...
                    sourcesMap(currentSource, 3 : length(sourcesMap)), ...
                    fractionalDelay);
                newReceivedSignal = newReceivedSignal(1 : inputSignal_length);

                % Amplitude compensation
                currentAmplitude = 1 / (GRID_RESOLUTION + currentDistance);
                newReceivedSignal = newReceivedSignal * currentAmplitude;
                
                % Accumulate received signal from current transmitting source
                accumulatedSignal = accumulatedSignal + newReceivedSignal;
            end

            amplitudeMap(y_index, x_index) = 20*log10(rms(accumulatedSignal));
        end
    end

    elapsedTime = toc;
    disp(sprintf('\n    Took: %d:%.1f [mm:ss]', round(elapsedTime/60), elapsedTime));

    % ------------- Sensors positions and received signals -------------

    amountOfSensors = 3;
    firstSensorPosition = [0, 0];
    sensorDisplacementVector = [0.1, 0];

    sensorSignals = zeros(amountOfSensors, inputSignal_length);
    sensorCoordinates = zeros(amountOfSensors, 2);

    for currentSensor = 1 : amountOfSensors
        for currentSource = 1 : amountOfSignals
            % Obtain distance from current sensor to each source
            currentSensorCoordinates = firstSensorPosition + ...
                (currentSensor-1) * sensorDisplacementVector;
            currentDistance = sqrt( ...
                (currentSensorCoordinates(1) - sourcesMap(currentSource, 1))^2 + ...
                (currentSensorCoordinates(2) - sourcesMap(currentSource, 2))^2 );

            % Delay signal from that transmitter according to that distance
            currentDelay = currentDistance / PROPAGATION_SPEED;
            fractionalDelay = generateFractionalDelayIR(currentDelay, SAMPLE_RATE);
            newReceivedSignal = conv(...
                sourcesMap(currentSource, 3 : inputSignal_length), ...
                fractionalDelay);
            newReceivedSignal = newReceivedSignal(1 : inputSignal_length);
 
            % Amplitude compensation
            currentAmplitude = 1 / currentDistance;
            newReceivedSignal = newReceivedSignal * currentAmplitude;

            % Accumulate received signal from current transmitting source
            sensorSignals(currentSensor, :) = ...
                sensorSignals(currentSensor, :) + newReceivedSignal;
        end

        % Saves the position of each sensor
        sensorCoordinates(currentSensor, :) = [currentSensorCoordinates(1) currentSensorCoordinates(2)];

    end

    save('amplitudeMap.mat', 'amplitudeMap', 'x_coordinates', 'y_coordinates', ...
        'GRID_RESOLUTION', 'sourcesMap', 'amountOfSignals', 'boundaries', ...
        'inputSignal_freqs', 'PROPAGATION_SPEED', 'SAMPLE_RATE', ...
        'sensorCoordinates', 'sensorSignals');

    plotAmplitudeMap_image('amplitudeMap.mat');

end
