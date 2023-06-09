%% Power Prediction Algorithm
% The algorithm identifies the highest performing cells and saves the mean
% MPP. Depending on the type and length of the microcrack, power
% predictions are generated and saved into a double array.

tic
% number of solar cells
N=10;
% read measured MPP values for validation set
mpp_cell = xlsread("measured_mpp_val.xlsx");

% used to store the MPP of the highest performing cells
peak_perf = zeros(N,1);
for i = 1:N 
    % highest performing cells have no cracks detected and the percent of
    % black pixels is less than 5%
    if isempty(line_data(i).lengths) == 1 && percent_black(i) < 5
        peak_perf(i) = mpp_cell(i);
    end
end

% takes the average of highest performing solar cells
peak_performance = mean(nonzeros(peak_perf));

% Originally, the training set was used to find the highest performing
% cells using the following code:

% read training data
% mpp_cell = xlsread("training_mpp.xlsx");

% Solar cells 10 through 20 were determined to be the highest performing
% using the performance factor.
% peak_performance = mean(mpp_cell([10:20],1));

% preallocate memory for predicted power array
predicted_pw = zeros(N,1);

% iterate through each solar cell and make power predictions
for i=1:N
    % check for non defective cells
    if isempty(line_data(i).lengths) == 1 && percent_black(i) < 5
        predicted_pw(i) = peak_performance;
    else 
        % check for cell inactive areas with cracks
        if percent_black(i) > 5 && length(line_data(i).angles) >= 1
            % store all the angles for that solar cell
            angles = [line_data(i).angles];

            for k=1:length(line_data(i).angles)
                % check for vertical lines
                if line_data(i).angles(k) > -5 && line_data(i).angles(k) < 5
                    % find index of solar cells with vertical cracks
                    vert_index = find(angles > -5 & angles < 5);
                    % apply linear regression model for vertically cracked
                    % solar cells with inactive areas
                    predicted_pw(i) = -0.006 * percent_black(i) + -.00095 * line_data(i).lengths{vert_index} + peak_performance;
                    break
                % check for horizontal lines
                elseif abs(line_data(i).angles(k)) > 80 && abs(line_data(i).angles(k)) < 100
                    % find index of horizontal cracks
                    hor_index = find(abs(angles) > 80 & abs(angles) < 100);
                    % apply linear regression moel for horizontally cracked
                    % solar cells with inactive areas
                    predicted_pw(i) = -0.0055 * percent_black(i) + -.00028 * line_data(i).lengths{hor_index} + peak_performance;
                    break
                % check for diagonal cracks
                elseif abs(line_data(i).angles(k)) > 35 && abs(line_data(i).angles(k)) < 50
                    % find index of diagonal cracks
                    diag_index = find(abs(angles) > 35 & abs(angles) < 50);
                    % apply linear regression moel for diagonally cracked
                    % solar cells with inactive areas
                        predicted_pw(i) = -0.0031 * percent_black(i) + -.00025 * line_data(i).lengths{diag_index} + peak_performance;
                        break
                end
            end
        
        % no cracks just inactive areas
        elseif  percent_black(i) > 5
            % apply linear regression model
           predicted_pw(i) = -0.0031 * percent_black(i) + peak_performance;
       
        % identify different types of cracks  
        else
            for k=1:length(line_data(1,i).angles) 
                angles = [line_data(i).angles];

                % check for horizontal cracks
                if abs(line_data(i).angles(k)) > 80 && abs(line_data(i).angles(k)) < 100
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                            % check for horizontal cracks and find index 
                            hor_index = find(abs(angles) > 80 & abs(angles) < 100);
                            % take the mean of the lengths and model as one crack
                            agg_crack = sum([line_data(i).lengths{hor_index}]);
                            % apply linear regression model
                            predicted_pw(i) = -0.00028 * agg_crack + peak_performance ;
                            break
                    % only one horizontal crack was detected
                    else
                        % apply linear regression model
                        predicted_pw(i) = -0.00028 * line_data(i).lengths{1} + peak_performance ;
                    end

                % check for vertical lines
                elseif line_data(i).angles(k) > -5 && line_data(i).angles(k) < 5
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                        % find index
                        vert_index = find(angles > -5 & angles < 5);
                        agg_crack = sum([line_data(i).lengths{vert_index}]);
                        % apply linear regression model
                        predicted_pw(i) = -0.0001 * agg_crack + peak_performance ;
                        break
                    % only one vertical crack was detected
                    else
                        % apply linear regression model
                        predicted_pw(i) = -0.0001 * line_data(i).lengths{1} + peak_performance ;
                    end
                
                % check for diagonal cracks
                elseif abs(line_data(i).angles(k)) > 35 && abs(line_data(i).angles(k)) < 50
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                        % find index of diagonal cracks
                        diag_index = find(abs(angles) > 35 & abs(angles) < 50);
                        agg_crack = sum([line_data(i).lengths{diag_index}]);
                        % apply linear regression model
                        predicted_pw(i) = -0.00025 * agg_crack + peak_performance ;
                        break
                    % only one diagonal crack was detected
                    else
                        % apply linear regression model
                        predicted_pw(i) = -0.00025 * line_data(i).lengths{1} + peak_performance;
                    end
                else
                    % apply linear regession model for remaining cracks
                    predicted_pw(i) = -0.0002 * line_data(i).lengths{1} + peak_performance;
                end
            end
        end
    end
end
%% Plots

% approximation error
figure(2)
scatter((1:10)', abs(((mpp_cell - predicted_pw) ./ mpp_cell) .* 100), "filled", "MarkerFaceColor", [0, 0.2, .5])

grid on 
xlim([0 11])
xlabel("Control Number");
ylabel("Approximation Error (%)");
title("Approximation Error [Third Iteration]") 

% residual plot
figure(3)
hold on
scatter((1:10)', mpp_cell - predicted_pw, "filled", "MarkerFaceColor", [0, 0.2, .5])
yline(0, "--", "color",[.89, 0.11, .14], "LineWidth",2)

grid on 
xlim([0 11])
xlabel("Control Number");
ylabel("Residuals");
ylim([-.04 .04])
title("Residual Plot [Third Iteration]") 

