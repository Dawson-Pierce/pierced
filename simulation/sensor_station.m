classdef sensor_station < handle

    properties
        location
        view_vector
        max_angle
        grid_slope
        grid_base
        internal_tracker
    end

    methods
        function obj = sensor_station(varargin) 
            p = inputParser;  
            addParameter(p, 'location', [0 0 0]); 
            addParameter(p, 'view_vector', [1 0 0]); 
            addParameter(p, 'max_angle', pi/3);  
            addParameter(p, 'grid_slope', 0.05);  
            addParameter(p, 'grid_base', 0.01);  
            addParameter(p, 'internal_tracker', []);  
            
            % Parse known arguments
            parse(p, varargin{:});

            obj.location = p.Results.location;
            obj.view_vector = p.Results.view_vector / norm(p.Results.view_vector);
            obj.max_angle = p.Results.max_angle;
            obj.grid_slope = p.Results.grid_slope;
            obj.grid_base = p.Results.grid_base;
            obj.internal_tracker = p.Results.internal_tracker;
        end 

        function plot(obj)
        
            range = 10;
            half_ang = obj.max_angle;

            fwd = obj.view_vector(:) / norm(obj.view_vector);

            tmp = [0;0;1];
            if abs(dot(tmp, fwd)) > 0.9
                tmp = [0;1;0];
            end

            right = cross(fwd, tmp);   right = right / norm(right);
            up    = cross(right, fwd); up = up / norm(up);

            R_up_pos    = axang2rotm([right'  half_ang]);
            R_up_neg    = axang2rotm([right' -half_ang]);
            R_right_pos = axang2rotm([up'     half_ang]);
            R_right_neg = axang2rotm([up'    -half_ang]);

            v1_dir = R_up_pos   * (R_right_pos * fwd);
            v2_dir = R_up_pos   * (R_right_neg * fwd);
            v3_dir = R_up_neg   * (R_right_neg * fwd);
            v4_dir = R_up_neg   * (R_right_pos * fwd);

            O  = obj.location(:);
            V1 = O + range * v1_dir;
            V2 = O + range * v2_dir;
            V3 = O + range * v3_dir;
            V4 = O + range * v4_dir;

            verts = [O'; V1'; V2'; V3'; V4'];

            faces = [
                1 2 3; 
                1 3 4; 
                1 4 5; 
                1 5 2; 
            ];

            patch('Vertices', verts, 'Faces', faces, ...
                  'FaceColor', [0 0.6 1], 'FaceAlpha', 0.2, ...
                  'EdgeColor', 'b','EdgeAlpha', 0, 'LineWidth', 0.5, ...
                  'BackFaceLighting', 'reverselit');
        
        end


        function pnts = get_points(obj,TR,T)

            if isa(TR,'cell') & isa(T,'cell')
                for k = 1:length(TR)
                    old_points = TR{k}.Points';
        
                    new_points = (T{k} * [old_points; ones(1, size(old_points, 2))])';
                    
                    TR_transformed = triangulation(TR{k}.ConnectivityList,new_points(:,1:3));
                    
                    ptCloud = mesh2pc(TR_transformed,"SamplingMethod","PoissonDiskSampling");
                
                    loc = T{k}(1:3,4) - obj.location';

                    [ptCloud, ~] = removeHiddenPoints(ptCloud,loc');
                
                    gridStep = obj.grid_slope * norm(loc) + obj.grid_base;
                
                    ptCloudPreMerge = pcdownsample(ptCloud,"gridAverage",gridStep);

                    if k == 1
                        ptCloudOut = ptCloudPreMerge;
                    else
                        ptCloudOut = pcmerge(ptCloudOut, ptCloudPreMerge, 0.01);
                    end
                end

                [ptCloudOut, ~] = removeHiddenPoints(ptCloudOut,obj.view_vector);
            else
                old_points = TR.Points';
            
                new_points = (T * [old_points; ones(1, size(old_points, 2))])';
                
                TR_transformed = triangulation(TR.ConnectivityList,new_points(:,1:3));
                
                ptCloud = mesh2pc(TR_transformed,"SamplingMethod","PoissonDiskSampling");
            
                loc = T(1:3,4) - obj.location';

                [ptCloud, ~] = removeHiddenPoints(ptCloud,loc');
            
                gridStep = obj.grid_slope * norm(loc) + obj.grid_base;
            
                ptCloudOut = pcdownsample(ptCloud,"gridAverage",gridStep); 

            end
            
            pts = ptCloudOut.Location; 
            vecs = pts - obj.location;
            vecs_norm = vecs ./ vecnorm(vecs,2,2);
            cosTheta = vecs_norm * obj.view_vector(:);
            angles = acos(cosTheta); 
            inside_idx = (-obj.max_angle <= angles) & (angles <= obj.max_angle);
        
            pnts = ptCloudOut.Location(inside_idx,:);
        end
    end
end