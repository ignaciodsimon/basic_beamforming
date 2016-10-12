function bartlettBeamforming()
    % ------------------------------------------------
    %      GENERALISED BARTLETT BEAMFORMING DEMO
    % Single source rotating around array simulation
    %
    %          Jose Ignacio Dominguez Simon
    %
    %            Array Signal Processing
    %           Aalborg University - 2015
    % ------------------------------------------------
    %
    % This program simulates an array under anechoic
    % conditions, where a source is playing describing
    % a circular movement around the array. The signals
    % are weighted using a modified Vandermonde vector
    % and the power for all points in the space (of a
    % predefined grid) is presented as an image.
    %
    % Joe.

    % -==- SIMULATION CONSTANTS

    AMOUNT_OF_SENSORS = 25;
    SENSORS_POSITION_ORIGIN = [0, 0];
    SENSORS_DISPLACEMENT_VECTOR = [0.10, 0];
    SENSORS_DISPLACEMENT_FACTOR = [0.01 0];      % <-- implement this to make the array non-uniform
    SENSOR_NOISE = -60;

    SOURCE_CIRCLE_RADIUS = 8;
    SOURCE_POSITION = [5, 7];
    PROPAGATION_SPEED = 340;
    SAMPLE_RATE = 44100;

    circularSourceAngles = [-180 : 1 : -170];
    x_range = [-10 10];
    y_range = [-10 10];
    grid_resolution = 0.5;


    % -==- OBTAIN SOURCE AUDIO

    inputSignal_length = 10000;

    loadedAudio = load('audio_fragments.mat');
    startPoint = 1;

%     inputSignal = conv(ones(4, 1)/ 4, loadedAudio.fragments(1,:));
%     inputSignal = loadedAudio.fragments(1, startPoint : startPoint + inputSignal_length -1);

%     inputSignal_freq = PROPAGATION_SPEED / norm(SENSORS_DISPLACEMENT_ORIGIN) / 2;
    inputSignal_freq = 100;
    inputSignal_amplitude = 1;

    % Generate a sine signal with a hamming window        
    inputSignal = sin(2 * pi * ...
        inputSignal_freq * [1 : inputSignal_length] / SAMPLE_RATE) ...
        .* (hamming(inputSignal_length)'.^2) ...
        * inputSignal_amplitude;

    

    % Obtain distance from current sensor to each source
    sensorCoordinates = zeros(AMOUNT_OF_SENSORS, 2);
    for index = 1 : AMOUNT_OF_SENSORS

        % Uniform linear array
        currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
            (index - 1) * SENSORS_DISPLACEMENT_VECTOR;

        % Spiral array
        currentAngle = 10/180*pi * index.^1.3;
        currentRadius = 0.25 + index/15;
        currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
            [cos(currentAngle) sin(currentAngle)] * currentRadius;

        % Non uniform linear array
        currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
            (index - 1) * SENSORS_DISPLACEMENT_VECTOR/1.3 + (index^2 * SENSORS_DISPLACEMENT_FACTOR);

        % Conical array
        if index < AMOUNT_OF_SENSORS/2
            currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
                (index - 1) * [0.1 0.1] .* [1 1+index^2/80];
        else
            currentSensorCoordinates = SENSORS_POSITION_ORIGIN + ...
                (index -AMOUNT_OF_SENSORS/2 - 1) * [0.1 -0.1] .* [1 (1+(index-AMOUNT_OF_SENSORS/2)^2/100)] + [0 -1];            
        end


        sensorCoordinates(index, :) = currentSensorCoordinates;
    end

%     figure
%     plot(sensorCoordinates(:, 1), sensorCoordinates(:, 2), 'o')
%     xlim([-2 5])
%     ylim([-5 5])
%     grid
%     return
%     




    x_axis = min(x_range) : grid_resolution : max(x_range);
    y_axis = min(y_range) : grid_resolution : max(y_range);

    outputPower = zeros(length(y_axis), length(x_axis));
    outputPowers = zeros(length(y_axis), length(x_axis), length(circularSourceAngles));

    disp(sprintf('\n------------------------------------------------\n     GENERALISED BARTLETT BEAMFORMING DEMO\n Single source rotating around array simulation\n\n         Jose Ignacio Dominguez Simon\n               Tobias Van Barsel\n\n           Array Signal Processing\n          Aalborg University - 2015\n------------------------------------------------\n'))

    disp(sprintf('Sensors count: %d', AMOUNT_OF_SENSORS));
    disp(sprintf('Sensors noise: %.1f [dB]', SENSOR_NOISE));
    disp(sprintf('Source circumference radius: %.1f [m]', SOURCE_CIRCLE_RADIUS));
    disp(sprintf('Propagation speed: %.1f [m/s]', PROPAGATION_SPEED));
    disp(sprintf('Sample rate: %.1f [Hz]', SAMPLE_RATE));
    disp(sprintf('Signal frequency: %.1f [Hz]', inputSignal_freq));
    disp(sprintf('Angles range: %.1f to %.1f [degrees]', min(circularSourceAngles), max(circularSourceAngles)));
    disp(sprintf('Grid resolution: %.3f [m]', grid_resolution));
    
    disp(sprintf('Expected aprox. output data size: %.2f [MB]\n', size(outputPowers,1) * size(outputPowers,2) * size(outputPowers,3) * 8 / 1000000));

    disp(sprintf('Scanning started ...\n'));
    for currentCircularSourceAngleIndex = 1 : length(circularSourceAngles)

        currentCircularSourceAngle = circularSourceAngles(currentCircularSourceAngleIndex);
        SOURCE_POSITION = [cos(currentCircularSourceAngle / 180 * pi) sin(currentCircularSourceAngle / 180 * pi)] * SOURCE_CIRCLE_RADIUS;
        disp(sprintf(' > Scanning angle %.2f [degrees] (%d of %d)', currentCircularSourceAngle, currentCircularSourceAngleIndex, length(circularSourceAngles)));

        % -==- CREATE SENSOR "RECORDINGS"

        lambda = inputSignal_freq / PROPAGATION_SPEED;
        sensorSignals = zeros(AMOUNT_OF_SENSORS, inputSignal_length);

        % Recreate the Vandermonde vector with associated delays
        for i = 1 : AMOUNT_OF_SENSORS

            % Obtain distance from current sensor to each source
            currentSensorCoordinates = sensorCoordinates(i, :);
            currentDistanceVector = SOURCE_POSITION - currentSensorCoordinates;

            % Update angle to current source
    %         theta = atan(currentDistanceVector(2) / currentDistanceVector(1));

            % Update the xsi value for each angle tested
    %         xsi = norm(SENSORS_DISPLACEMENT_ORIGIN) / lambda * cos(theta);
            xsi = norm(currentDistanceVector) / lambda;

            currentAmplitude = 1 / (1 + norm(currentDistanceVector));

            % "Delay" source signal as necessary
            vandermonde = currentAmplitude * exp(1i * 2 * pi * xsi);% * (i-1));
            sensorSignals(i, :) = inputSignal * vandermonde + ...
                randn(1, inputSignal_length) * 10^(SENSOR_NOISE/20);

            % Save the coordinates of the sensor
            sensorCoordinates(i, :) = currentSensorCoordinates;

        end


        % -==- ESTIMATION OF COVARIANCE MATRIX FROM SENSED SIGNALS

    %     covariances = zeros(AMOUNT_OF_SENSORS, AMOUNT_OF_SENSORS, length(sensorSignals));

    %     figureHandler = figure();
    %     set(figureHandler, 'Position', [500 100 400 600])    

        covarianceMatrix = zeros(AMOUNT_OF_SENSORS);
        timeRange = [1, inputSignal_length];
        counter = 1;
        for t = timeRange(1) : timeRange(2)
            covarianceMatrix = covarianceMatrix + ...
                sensorSignals(:, t) * sensorSignals(:, t)';        

    %         clf
    %         subplot(2,1,1)
    %         plot(sensorSignals(:, t), 'x-')
    %         xlabel('Input sensor number')
    %         ylabel('Sensor amplitude [.]')
    %         title('Captured signals')
    %         grid
    % %         ylim([-0.5 0.5])
    %         subplot(2,1,2)
    %         imagesc(covarianceMatrix)
    %         title(sprintf('Estimated array covariance matrix - Time averages: %d', counter));
    %         xlabel('Sensor number')
    %         ylabel('Sensor number')
    % 
    %         f = getframe(figureHandler);
    %         if counter == 1 
    %             [im, map] = rgb2ind(f.cdata, 65535, 'nodither');
    %         else
    %             [im(:,:,1,counter), map] = rgb2ind(f.cdata, map, 'nodither');
    %         end
    %         counter = counter + 1;
    % 
    %         disp(sprintf('frame: %d / %d', t, length(sensorSignals)))

        end

        covarianceMatrix = covarianceMatrix / (timeRange(2) - timeRange(1) + 1);


        % -==- BARTLETT OUTPUT POWER SCANNING

        % Not sure about this: what frequency to use?    <-- INVESTIGATE THIS

        for x_pos_index = 1 : length(x_axis)
            x_pos = x_axis(x_pos_index);
            for y_pos_index = 1 : length(y_axis)
                y_pos = y_axis(y_pos_index);

                % Recreate the Vandermonde vector with associated delays
                vandermonde = zeros(AMOUNT_OF_SENSORS, 1);
                for sensorIndex = 1 : AMOUNT_OF_SENSORS

                    % Obtain distance from current sensor to each source
                    currentSensorCoordinates = sensorCoordinates(sensorIndex, :);
                    currentDistanceVector = [x_pos y_pos] - currentSensorCoordinates;

                    % Update angle to current source
    %                 theta = atan(currentDistanceVector(2) / currentDistanceVector(1));

                    % Update the xsi value for each angle tested
                    xsi = norm(currentDistanceVector) / lambda;

                    % Use distance to current scan point to compensate amplitude
                    currentAmplitude = (1 + norm(currentDistanceVector));

                    % Obtain current entry for the Vandermonde vector
                    vandermonde(sensorIndex) = currentAmplitude * exp(1i * 2 * pi * xsi);% * (sensorIndex-1));
                end

                % Compute the output power for current tested angle
                outputPower(y_pos_index, x_pos_index) = (ctranspose(vandermonde) * covarianceMatrix * vandermonde) / (norm(vandermonde)^2);        
            end
        end
        outputPower = abs(outputPower);

        outputPowers(:, :, currentCircularSourceAngleIndex) = outputPower;
    end

    disp(' > Saving generated data ...')
    fileNumber = 1;
    filenameToSave = sprintf('bartlett_beamformed_data_%d.mat', fileNumber);
    while exist(filenameToSave, 'file')
        fileNumber = fileNumber + 1;
        filenameToSave = sprintf('bartlett_beamformed_data_%d.mat', fileNumber);
    end
    save(filenameToSave, 'outputPowers', 'x_axis', 'y_axis', 'SOURCE_POSITION', 'sensorCoordinates', 'AMOUNT_OF_SENSORS', 'circularSourceAngles');
    disp(sprintf(' > Saved to "%s"\n', filenameToSave));

    disp('All done. Bye!')

return
