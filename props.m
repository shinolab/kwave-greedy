% Copyright (c) 2022 Shun Suzuki. All rights reserved.
% File              : props.m
% License           : MIT
% Author            : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>
% Date              : 15/07/2022
% Last Modified Date: 15/07/2022
% Last Modified By  : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>

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
