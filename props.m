classdef props
    properties
        spacing     % m
        sound_speed % m/s
        frequency   % Hz
        amplitude   % dB
        dt          % s
    end

    methods
        function res = amp(obj)
            pr = 20e-6;
            res = pr * 10^(obj.amplitude/20);
        end

        function res = wavelength(obj)
            res = obj.sound_speed / obj.frequency;
        end

        function res = wavenumber(obj)
            res = 2 * pi / obj.wavelength();
        end
    end
end
