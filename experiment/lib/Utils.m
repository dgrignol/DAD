classdef Utils

    properties(Constant)
        %stimuli type enum
        likelihood = 0;
        path_duration = 1;
        path_duration_norm = 2;
    end

    methods(Static)
        function inputs = GetUserDoubles(title, dims, prompts, defaults)
            
            notValid = true;
            while notValid

                inputs = [];
                userInput = inputdlg(prompts, title, dims, defaults);

                for index = 1:length(userInput)
                    input = str2double(userInput{index});
                    if isnan(input) || input <= 0
                        disp('Invalid input. Please enter positive numeric values.');
                        notValid = true;
                        break;
                    end

                    inputs = [inputs input];
                    notValid = false;
                end
            end

        end

        function angle = ComputeVisualAngleDegrees(viewDistance, monitorWidth)
            %formula of visual angle in rad
            visualAngle = 2*atan(monitorWidth/(2*viewDistance));
            %return in degree
            angle = rad2deg(visualAngle);
        end

        function sign = RandPosNeg()
            signs = [-1 1];
            
            sign = signs(randi([1 2], 1));
        end

        %generates a random angle in degrees
        function angle = RandAngleDegree()
            radAngle = rand(1, 1)*2*pi;

            angle = rad2deg(radAngle);
        end

        %computes curvyness factor with respec to the config parameters, it
        %can be of both random valence and magnitude, if randValence then
        %the factor sign is randomly assigned else it will always be
        %negative to swap the previous one; if randMagnitude the factor is
        %multiplied by a rand 0.xxx ratio else remains 1 to mantain
        %previous value
        function curvFactor = ComputeCurvyness(randValence, randMagnitude)
            
            if randValence
                valence = Utils.RandPosNeg();
            else %the valence alternates between +/- instead of random
                valence = -1; 
            end

            if randMagnitude
                magnitude = rand(1, 1);
            else %the factor remains the same
                magnitude = 1;
            end

            curvFactor = valence*magnitude;
        end

        %clean the comments and the last function after you finish, might
        %not need the check of boundaries, see the comment that i wrote in
        %marissas


        %computes the offset for the starting point coordinates that would
        %be added to the absolute center, and makes sure that its respects the
        %dotRectSize boundaries; it loops until randomly getting and accepted
        %starting point
        function coord = ComputeOffsetCoordinates(aroundCenterfactor, rectSize)            
            %get the rect mid point coord
            midPoint = rectSize/2;

            %get the offset, which should be +/- rands <1
            offset = rand(1, 2).*[Utils.RandPosNeg() Utils.RandPosNeg()]*aroundCenterfactor;
            coord = offset.*midPoint; %this would be a scaled wrt the midpoint offset

            %make sure you are still within the rectSize boundary --> so
            %this loop seems redundant since you already take half the
            %recSize and multiply by a num < 1
            
            % while (coord(1) < -rectSize(1)/2) || (coord(1) > rectSize(1)/2) ...
            %         || (coord(2) < -rectSize(2)/2) || (coord(2) > rectSize(2)/2)
            % 
            %     %get the offset, which should be +/- rands <1
            %     offset = rand(1, 2).*[Utils.RandPosNeg() Utils.RandPosNeg()]*aroundCenterfactor;
            %     coord = offset.*midPoint; %this would be a scaled wrt the midpoint offset
            % end
        end

        %takes a degree angle as input and return a direction unitary
        %vector by computing its sin and cos
        function vector = GetDirectionVector(angle, dim)
            radAngle = deg2rad(angle);
            if dim == 2
                vector = [cos(radAngle(1)), sin(radAngle(1)), cos(radAngle(2)), sin(radAngle(2))];
            else
                vector = [cos(radAngle), sin(radAngle)];  
            end
        end

        function vector = GetSumDirectionVector(initialAngle, increment, N)
            radInitialAngle = deg2rad(initialAngle);
            radIncrement = deg2rad(increment);

            factor = sin((N*radIncrement)/2);
            factor = factor / sin(radIncrement/2);

            componentAngle = ((N-1)*radIncrement)/2;
            componentAngle = componentAngle + radInitialAngle;

            vector = [cos(componentAngle), sin(componentAngle)];
            vector = factor*vector;
        end

        %return a logical grid in which the dot position is set in its
        %respective grid
        function grid = GetFrameInGrid(numXgrids, gridXsize, numYgrids, gridYsize, dotXY)
            grid = zeros(numXgrids, numYgrids);
            for indexXGrid  = 1:numXgrids
                if (gridXsize*(indexXGrid-1) < dotXY(1)) ...
                    && (dotXY(1) <= gridXsize*(indexXGrid))
                    for indexYGrid = 1:numYgrids
                        if (gridYsize*(indexYGrid-1) < dotXY(2)) ...
                                && (dotXY(2) <= gridYsize*(indexYGrid))
                            grid(indexYGrid, indexXGrid) = grid(indexYGrid, indexXGrid) + 1;
                            break;
                        end
                    end
                    break;
                end
            end
        end

    end
end