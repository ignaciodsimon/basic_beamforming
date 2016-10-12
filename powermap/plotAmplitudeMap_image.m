function plotAmplitudeMap_image(filename)
    % Displays the amplitude / power map generated using the function "generatePowerMap.m"
    % It includes the boundaries (if any) and array elements in the plot.
    %
    % Joe.

    loadedData = load(filename);

    figureHandler = figure();

    % Power map
    imagesc(loadedData.x_coordinates, loadedData.y_coordinates, loadedData.amplitudeMap);
    set(gca,'YDir','normal')
    hold on

    % Crosses on top of sources
    for i = 1 : loadedData.amountOfSignals
        
        x_pos = loadedData.sourcesMap(i, 1);
        y_pos = loadedData.sourcesMap(i, 2);

        plot(x_pos, y_pos, 'x', 'MarkerSize', 20, 'LineWidth', 6, 'Color', [1 1 1]);
        plot(x_pos, y_pos, 'x', 'MarkerSize', 10, 'LineWidth', 4, 'Color', [1 .2 .3]);
    end

    % Boundaries 
    for i = 1 : size(loadedData.boundaries, 1)
        
        plot([loadedData.boundaries(i, 1) loadedData.boundaries(i, 3)], ...
              [loadedData.boundaries(i, 2) loadedData.boundaries(i, 4)], ...
              'LineWidth', 6, 'Color', [1 1 1]);

        plot([loadedData.boundaries(i, 1) loadedData.boundaries(i, 3)], ...
              [loadedData.boundaries(i, 2) loadedData.boundaries(i, 4)], ...
              'LineWidth', 4, 'Color', [1 .2 .3]);
    end

    % Sensors positions
    for i = 1 : size(loadedData.sensorCoordinates, 1);
        plot(loadedData.sensorCoordinates(i, 1), loadedData.sensorCoordinates(i, 2), ...
            'x', 'MarkerSize', 20, 'LineWidth', 6, 'Color', [1 1 1]);
        plot(loadedData.sensorCoordinates(i, 1), loadedData.sensorCoordinates(i, 2), ...
            'x', 'MarkerSize', 10, 'LineWidth', 4, 'Color', [.2 .3 1]);
    end

    % Rest of the furniture of the plot
    xlabel('X axis [m]', 'FontSize', 14);
    ylabel('Y axis [m]', 'FontSize', 14);

    set(figureHandler, 'Position', [500 100 1200 800])    

    colorBarHandler = colorbar('eastoutside', 'FontSize', 14);
	colorBarHandler.Label.String = 'Power [dB]';

    title(sprintf('Freq: %.1f [Hz] - Propagation speed: %.1f [m/s] - Sample rate: %.1f [Hz]', ...
        loadedData.inputSignal_freqs(1), ...
        loadedData.PROPAGATION_SPEED, ...
        loadedData.SAMPLE_RATE), ...
        'FontSize', 14)

end
