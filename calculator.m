% Copyright (c) 2022 Shun Suzuki. All rights reserved.
% File              : calculator.m
% License           : MIT
% Author            : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>
% Date              : 15/07/2022
% Last Modified Date: 15/07/2022
% Last Modified By  : Shun Suzuki <suzuki@hapis.k.u-tokyo.ac.jp>

classdef calculator
    properties
        single_source_sensor_data
        trans_positions
        p
        rx
        ry
        rz
        rms_idx
        NtRms
    end

    methods(Static)
        function res = get_convex_hull(p)
            r_min = min(p);
            r_max = max(p);
            res = [r_min r_max];
        end

        function idx = posToIdx(pos, rx, ry, rz, spacing)
            ix = round((pos(1) - rx(1)) / spacing) + 1;
            iy = round((pos(2) - ry(1)) / spacing) + 1;
            iz = round((pos(3) - rz(1)) / spacing) + 1;
            idx = [ix iy iz];
        end

        function grid = gridNum(rx, ry, rz, p)
            grid_x = round((rx(2)- rx(1)) / p.spacing) + 1;
            grid_y = round((ry(2)- ry(1)) / p.spacing) + 1;
            grid_z = round((rz(2)- rz(1)) / p.spacing) + 1;
            grid = [grid_x grid_y grid_z];
        end

        function res = maxDist(rx, ry, rz)
            res = norm([(rx(2) - rx(1)) (ry(2) - ry(1)) (rz(2) - rz(1))]);
        end
    end

    methods
        function obj = calculator(trans_positions, foci, p)
            obj.trans_positions = trans_positions;
            obj.p = p;

            rx = calculator.get_convex_hull([trans_positions(:, 1)', foci(:, 1)']);
            ry = calculator.get_convex_hull([trans_positions(:, 2)', foci(:, 2)']);
            rz = calculator.get_convex_hull([trans_positions(:, 3)', foci(:, 3)']);
            rx = rx + [min(trans_positions(:, 1)) max(trans_positions(:, 1))];
            ry = ry + [min(trans_positions(:, 2)) max(trans_positions(:, 2))];
            obj.rx = rx;
            obj.ry = ry;
            obj.rz = rz;

            grid = calculator.gridNum(rx, ry, rz, p);
            kgrid = kWaveGrid(grid(1), p.spacing, grid(2), p.spacing, grid(3), p.spacing);
            
            rms_idx = obj.rmsIdx(rx, ry, rz, p.dt);
            kgrid.t_array = (0:rms_idx(2) + round(2 / (p.dt * p.frequency))) * p.dt;

            obj.rms_idx = rms_idx;
            obj.NtRms = rms_idx(2) - rms_idx(1) + 1;

            gridc = {grid(:)'};
            source.p_mask = zeros(gridc{:});
            source_idx = calculator.posToIdx([0 0 0], rx, ry, rz, p.spacing);
            source.p_mask(source_idx(1), source_idx(2), source_idx(3)) = 1;
            source.p = p.amp() * sin(2 * pi * p.frequency * kgrid.t_array);

            sensor.mask = zeros(gridc{:});
            sensor.mask(:, :, end) = 1;
            sensor.record = {'p'};

            medium.sound_speed = p.sound_speed;

            input_args = {'PMLInside', false, 'PlotPML', false, 'DataCast', 'single', 'CartInterp', 'nearest'};
            single_source_sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

            obj.single_source_sensor_data = reshape(single_source_sensor_data.p, [grid(1) grid(2) kgrid.Nt]);
        end

        function res = rmsIdx(obj, rx, ry, rz, dt)
            max_dist = calculator.maxDist(rx, ry, rz);
            reach_time = max_dist / obj.p.sound_speed;
            rms_start = round(reach_time / dt);
            rms_end  = rms_start + round(1 / (dt * obj.p.frequency));
            res = [rms_start rms_end];
        end

        function res = calc_field(obj, trans_idx, foci, phase)
            m = size(foci, 1);

            tmp = zeros(m, obj.NtRms);

            tr = obj.trans_positions(trans_idx, :);
            shift = round(phase / (2 * pi) / obj.p.frequency / obj.p.dt); 
            for i = 1:m
                focus = foci(i, :);
                pos = focus - tr;
                sensor_idx = calculator.posToIdx(pos, obj.rx, obj.ry, obj.rz, obj.p.spacing);
                tmp(i, :) = obj.single_source_sensor_data(sensor_idx(1), sensor_idx(2), (obj.rms_idx(1) + shift):(obj.rms_idx(2) + shift));
            end
            
            res = tmp;
        end

        function res = calc_rms(obj, trans_phases, observe_x, observe_y, observe_z)
            rx_ = calculator.get_convex_hull([obj.trans_positions(:, 1)', observe_x(1), observe_x(2)]);
            ry_ = calculator.get_convex_hull([obj.trans_positions(:, 2)', observe_y(1), observe_y(2)]);
            rz_ = calculator.get_convex_hull([obj.trans_positions(:, 3)', observe_z(1), observe_z(2)]);

            grid = calculator.gridNum(rx_, ry_, rz_, obj.p);

            kgrid = kWaveGrid(grid(1), obj.p.spacing, grid(2), obj.p.spacing, grid(3), obj.p.spacing);
            
            rms_idx_ = obj.rmsIdx(rx_, ry_, rz_, obj.p.dt);
            kgrid.t_array = (0:rms_idx_(2)) * obj.p.dt;

            gridc = {grid(:)'};
            source.p_mask = zeros(gridc{:});
            source.p = zeros(size(obj.trans_positions, 1), kgrid.Nt);
            for i = 1:size(obj.trans_positions, 1)
                tr = obj.trans_positions(i, :);
                phase = trans_phases(i);
                tr_idx = calculator.posToIdx(tr, rx_, ry_, rz_, obj.p.spacing);
                source.p_mask(tr_idx(1), tr_idx(2), tr_idx(3)) = 1;
                source.p(i, :) = obj.p.amp() * sin(2 * pi * obj.p.frequency * kgrid.t_array + phase * ones(1, kgrid.Nt));
            end

            sensor.mask = zeros(gridc{:});
            obs_st = calculator.posToIdx([observe_x(1) observe_y(1) observe_z(1)], rx_, ry_, rz_, obj.p.spacing);
            obs_ed = calculator.posToIdx([observe_x(2) observe_y(2) observe_z(2)], rx_, ry_, rz_, obj.p.spacing);
            sensor.mask(obs_st(1):obs_ed(1), obs_st(2):obs_ed(2), obs_st(3):obs_ed(3)) = 1;
            sensor.record = {'p'};

            medium.sound_speed = obj.p.sound_speed;

            input_args = {'PMLInside', false, 'PlotPML', false, 'DataCast', 'single', 'CartInterp', 'nearest'};
            sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

            p_rms = rms(sensor_data.p(:, rms_idx_(1):rms_idx_(2)), 2);
            res = reshape(p_rms, [obs_ed(1) - obs_st(1) + 1, obs_ed(2) - obs_st(2) + 1]);
        end
    end
end

