function generateGIFfromBeamformedData(filename)
    % Plots all frames included in the generated file from the beamforming simulation.
    % Also, generates a GIF file at the end.
    %
    % Joe.

    disp(sprintf('Loading data from file "%s" ...', filename));
    loadedBeamformedData = load(filename);
    figureHandler = figure();
    frameCounter = 1;

    disp(sprintf('Generating plots ...'));
    for currentCircularSourceAngleIndex = 1 : length(loadedBeamformedData.circularSourceAngles)

        currentCircularSourceAngle = loadedBeamformedData.circularSourceAngles(currentCircularSourceAngleIndex);      
        SOURCE_POSITION = [cos(currentCircularSourceAngle / 180 * pi) sin(currentCircularSourceAngle / 180 * pi)] * 8;

        clf
        set(figureHandler, 'Position', [500 100 1200 800])    
        imagesc(loadedBeamformedData.x_axis, loadedBeamformedData.y_axis, loadedBeamformedData.outputPowers(:, :, currentCircularSourceAngleIndex));%,'EdgeColor', 'none', 'LineStyle', 'none', 'FaceLighting', 'phong');
        set(gca,'YDir','normal')
        hold on
        title(sprintf('Source angle: %.1f   (%d/%d)', currentCircularSourceAngle, currentCircularSourceAngleIndex, length(loadedBeamformedData.circularSourceAngles)), 'FontSize', 14);

        % Plotted first to create a legend with only two entries
        plot(SOURCE_POSITION(1), SOURCE_POSITION(2), 'x', 'MarkerSize', 7, 'LineWidth', 5, 'Color', [ 1 .2 .3]);
        plot(loadedBeamformedData.sensorCoordinates(1, 1), loadedBeamformedData.sensorCoordinates(1, 2), 'o', 'MarkerSize', 4, 'Color', [.2 .3 1]);
        legend({'Source', 'Sensors'}, 'FontSize', 14);

        plot(SOURCE_POSITION(1), SOURCE_POSITION(2), 'x', 'MarkerSize', 15, 'LineWidth', 10, 'Color', [ 1 1 1]);
        plot(SOURCE_POSITION(1), SOURCE_POSITION(2), 'x', 'MarkerSize', 7, 'LineWidth', 5, 'Color', [ 1 .2 .3]);

        for i = 1 : loadedBeamformedData.AMOUNT_OF_SENSORS
            currentSensorCoordinates = loadedBeamformedData.sensorCoordinates(i, :);
            plot(currentSensorCoordinates(1), currentSensorCoordinates(2), 'o', 'MarkerSize', 7, 'Color', [1 1 1]);
            plot(currentSensorCoordinates(1), currentSensorCoordinates(2), 'o', 'MarkerSize', 4, 'Color', [.2 .3 1]);
        end

        xlabel('X [m]', 'FontSize', 14)
        ylabel('Y [m]', 'FontSize', 14)
        set(gca,'FontSize', 14)

    %     drawnow

        currentFrame = getframe(figureHandler);
        if frameCounter == 1 
            [imageToSave, colorsMap] = rgb2ind(currentFrame.cdata, 256, 'nodither');
        else
            [imageToSave(:, :, 1, frameCounter), colorsMap] = rgb2ind(currentFrame.cdata, colorsMap, 'nodither');
        end
        frameCounter = frameCounter + 1;

%        pause(0.05)

    end

    disp('Saving GIF ...')
    fileNumber = 1;
    filenameToSave = sprintf('%s_animation_%d.gif', filename, fileNumber);
    while exist(filenameToSave, 'file')
        fileNumber = fileNumber + 1;
        filenameToSave = sprintf('%s_animation_%d.gif', filename, fileNumber);
    end
%     filenameToSave = 'animation123.gif';
    imwrite(imageToSave, colorsMap, filenameToSave, 'DelayTime', 0, 'LoopCount', 0);
    disp(sprintf('Saved to "%s"', filenameToSave));

end
