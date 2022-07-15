% Copyright (c) 2022 Shun Suzuki. All rights reserved.
% File              : main.m
% License           : MIT
% Author            : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>
% Date              : 15/07/2022
% Last Modified Date: 15/07/2022
% Last Modified By  : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>

clearvars;

trans_spacing = 10e-3; % m
trans_x = 10;
trans_y = 10;

% center of transducers is [0 0 0]
trans_positions = gen_transducers(trans_x, trans_y, trans_spacing);

z = 100e-3;
foci = [
    0 30e-3 z;
    0 -30e-3 z;
    30e-3 0 z;
    -30e-3 0 z;
    ]; % m

n = size(trans_positions, 1);
m = size(foci, 1);

props = props();
props.spacing = 1e-3;    % m
props.sound_speed = 340; % m/s
props.frequency = 40e3;  % Hz
props.amplitude = 117;   % dB
props.dt = 25e-6 / 32;   % s (I'm not sure this is good)

% foci is required to determine how much range to calculate
c = calculator(trans_positions, foci, props);

phase_div = 16;
phase_candidates = zeros(phase_div, 1);
for i = 1:phase_div
    phase_candidates(i) = 2 * pi * i / phase_div;
end

trans_phases = zeros(1, n);
select = randperm(n);
cache = zeros(m, c.NtRms);
objective = @(x) sum(rms(x, 2));
for i = select
    max_phase = 0;
    max_p     = 0;
    max_field = zeros(m, c.NtRms);
    for j = 1:phase_div
        phase = phase_candidates(j);
        tmp = c.calc_field(i, foci, phase);
        p = objective(tmp + cache);
        if p > max_p
            max_p = p;
            max_phase = phase;
            max_field = tmp;
        end
    end
    cache = cache + max_field;
    trans_phases(i) = max_phase;
end

% plot
observe_x = [-50e-3 50e-3]; % m
observe_y = [-50e-3 50e-3]; % m
observe_z = [z z];  % m
p_rms = c.calc_rms(trans_phases, observe_x, observe_y, observe_z);

figure;
imagesc(p_rms(:, :));
xticks(1:10:size(p_rms, 2));
xlabels = 1:10:size(p_rms, 2);
xlabels = xlabels - ones(1, length(xlabels)) * (size(p_rms, 2)+1)/2;
xticklabels(string(xlabels));
xlabel('y-position [mm]');
yticks(1:10:size(p_rms, 1));
ylabels = 1:10:size(p_rms, 1);
ylabels = ylabels - ones(1, length(ylabels)) * (size(p_rms, 1)+1)/2;
yticklabels(string(ylabels));
ylabel('x-position [mm]');
colormap(jet(256));
ylabel(colorbar, 'Pressure [Pa]');
