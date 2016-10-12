function sincIR = generateFractionalDelayIR(fractionalDelay, sampleRate)
    % Generates the impulse response of a pure delay, allowing for
    % fractional values.
    %
    % Joe.

    if fractionalDelay * sampleRate > 100
        sincIR = zeros(round(fractionalDelay * 2 * sampleRate), 1);
    else
        sincIR = zeros(200, 1);
    end

    if fractionalDelay == round(fractionalDelay)
        sincIR(fractionalDelay + 1) = 1;
    else
        for i = 1 : length(sincIR)
            sincArgument = pi * (i - fractionalDelay * sampleRate);
            if sincArgument ~= 0
                sincIR(i) = sin(sincArgument) / sincArgument;
            else
                sincIR(i) = 1;
            end
        end
    end
end
