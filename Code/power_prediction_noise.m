%% Power Prediction Algorithm with Noise
% This is the first version of the predictive algorithm. It adds random
% noise to the linear regression model by using the standard deviation of
% the solar cells that are the highest performing.

tic
% Power Prediction
mpp_cell = xlsread("measured_mpp_val.xlsx");
N=10;
peak_perf = zeros(N,1);
for i = 1:N 
    if isempty(line_data(i).lengths) == 1 && percent_black(i) < 5
        peak_perf(i) = mpp_cell(i);
    end
end
    

% takes the average of highest performing solar cells
peak_performance = mean(nonzeros(peak_perf));
% finds standard deviation of highest performing cells
peak_std = std(nonzeros(peak_perf));

% preallocate memory for predicted power array
predicted_pw = zeros(N,1);

for i=1:N
    % check for non defective cells
    if isempty(line_data(i).lengths) == 1 && percent_black(i) < 5
        noise = rand*((peak_performance + peak_std) -  (peak_performance - peak_std));
        r = randi([0 1],1,1);
        if r == 0
            predicted_pw(i) = peak_performance - peak_std;
         % subtraction
        else
            predicted_pw(i) = peak_performance + peak_std;
        end
    % there are inactive areas or cracks
    else 
        noise = rand*((peak_performance + peak_std) -  (peak_performance - peak_std));
        r = randi([0 1],1,1);
        % check for cell inactive areas with cracks
        if percent_black(i) > 5 && length(line_data(i).angles) >= 1
            % store all the angles for that solar cell
            angles = [line_data(i).angles];
            for k=1:length(line_data(i).angles)
                % check for vertical lines
                if line_data(i).angles(k) > -5 && line_data(i).angles(k) < 5
                    % find index of solar cells with vertical cracks
                    vert_index = find(angles > -5 & angles < 5);
                    if  r == 0
                        % linear regression model for cracks and inactive area + noise
                        predicted_pw(i) = -0.002 * percent_black(i) + -.0001 * line_data(i).lengths{vert_index} + peak_performance + noise;
                        break
                    else
                        predicted_pw(i) = -0.002 * percent_black(i) + -.0001 * line_data(i).lengths{vert_index} + peak_performance - noise;
                        break
                    end
                % check for horizontal lines
                elseif abs(line_data(i).angles(k)) > 80 && abs(line_data(i).angles(k)) < 100
                    % find index of horizontal cracks
                    hor_index = find(abs(angles) > 80 & abs(angles) < 100);
                    if  r == 0
                        % linear regression model for cracks and inactive area + noise
                        predicted_pw(i) = -0.0031 * percent_black(i) + -.00015 * line_data(i).lengths{hor_index} + peak_performance + noise;
                        break
                    else
                        predicted_pw(i) = -0.0031 * percent_black(i) + -.00015 * line_data(i).lengths{hor_index} + peak_performance - noise;
                        break
                    end
                % check for diagonal cracks
                elseif abs(line_data(i).angles(k)) > 35 && abs(line_data(i).angles(k)) < 50
                    % find index of diagonal cracks
                    diag_index = find(abs(angles) > 35 & abs(angles) < 50);
                    if  r == 0
                        % linear regression model for cracks and inactive area + noise
                        predicted_pw(i) = -0.0031 * percent_black(i) + -.00025 * line_data(i).lengths{diag_index} + peak_performance + noise;
                        break
                    else
                        predicted_pw(i) = -0.0031 * percent_black(i) + -.00025 * line_data(i).lengths{diag_index} + peak_performance - noise;
                        break
                    end
                end
            end
        
            
        % no cracks just inactive areas
        elseif  percent_black(i) > 5
            if  r == 0
               % linear regression model + noise
               predicted_pw(i) = -0.0031 * percent_black(i) + peak_performance + noise;
            else
                predicted_pw(i) = -0.0031 * percent_black(i) + peak_performance - noise;
            end
       
        % identify different types of cracks  
        else
            % check for diagonal cracks (+-45 deg w/ 10 deg upper and lower
            % bound
            for k=1:length(line_data(1,i).angles) 
                angles = [line_data(i).angles];

                % check for horizontal cracks
                if abs(line_data(i).angles(k)) > 80 && abs(line_data(i).angles(k)) < 100
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                        % horizontal cracks
                        
                            % check for horizontal cracks and find index 
                            hor_index = find(abs(angles) > 80 & abs(angles) < 100);
                            % take the mean of the lengths and model as one crack
                            mean_crack = mean([line_data(i).lengths{hor_index}]);
                            if r == 0
                                predicted_pw(i) = -0.00015 * mean_crack + peak_performance + noise;
                                break
                            else
                                predicted_pw(i) = -0.00015 * mean_crack + peak_performance - noise;
                                break
                            end
                    else
                        if r == 0
                            predicted_pw(i) = -0.00015 * line_data(i).lengths{1} + peak_performance + noise;
                        else
                            predicted_pw(i) = -0.00015 * line_data(i).lengths{1} + peak_performance - noise;
                        end
                    end

                elseif line_data(i).angles(k) > -5 && line_data(i).angles(k) < 5
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                        % find index
                        vert_index = find(angles > -5 & angles < 5);
                        mean_crack = mean([line_data(i).lengths{vert_index}]);
                        if r == 0
                                predicted_pw(i) = -0.0001 * mean_crack + peak_performance + noise;
                                break
                            else
                                predicted_pw(i) = -0.0001 * mean_crack + peak_performance - noise;
                                break
                        end
                    else
                        if r == 0
                            predicted_pw(i) = -0.0001 * line_data(i).lengths{1} + peak_performance + noise;
                        else
                            predicted_pw(i) = -0.0001 * line_data(i).lengths{1} + peak_performance - noise;
                        end
                    end
                
                % check for diagonal cracks
                elseif abs(line_data(i).angles(k)) > 35 && abs(line_data(i).angles(k)) < 50
                    % check if multiple cracks are detected for one cell
                    if length(line_data(i).angles) > 1
                        % find index of diagonal cracks
                        diag_index = find(abs(angles) > 35 & abs(angles) < 50);
                        mean_crack = mean([line_data(i).lengths{diag_index}]);
                        if r == 0
                                predicted_pw(i) = -0.00025 * mean_crack + peak_performance + noise;
                                break
                            else
                                predicted_pw(i) = -0.00025 * mean_crack + peak_performance - noise;
                                break
                        end
                    else
                        if r == 0
                            predicted_pw(i) = -0.00025 * line_data(i).lengths{1} + peak_performance + noise;
                        else
                            predicted_pw(i) = -0.00025 * line_data(i).lengths{1} + peak_performance - noise;
                        end
                    end

                else
                    if r == 0
                        predicted_pw(i) = -0.0002 * line_data(i).lengths{1} + peak_performance + noise;
                    else
                        predicted_pw(i) = -0.0002 * line_data(i).lengths{1} + peak_performance - noise;
                    end
                end
            end
        end
    end
end
toc
%% Plots
% approximation error
hold off
figure(2)
hold on
scatter((1:10)', abs(((mpp_cell - predicted_pw) ./ mpp_cell) .* 100), "filled", "MarkerFaceColor", [0, 0.2, .5])
% % plot((1:35)',0,'--', "DisplayName",txt, "color", [.89, 0.11, .14]) 
% % scatter((1:35)', predicted_pw, "filled", "MarkerFaceColor", [.89, 0.11, .14], "DisplayName","Predicted MPP")
grid on 
xlim([0 11])

xlabel("Control Number");
ylabel("Approximation Error (%)");

title("Approximation Error") %% vertical crack detection

% residual plot
figure(3)
hold on
scatter((1:10)', mpp_cell - predicted_pw, "filled", "MarkerFaceColor", [0, 0.2, .5])
yline(0, "--", "color",[.89, 0.11, .14], "LineWidth",2)
% % plot((1:35)',0,'--', "DisplayName",txt, "color", [.89, 0.11, .14]) 
% % scatter((1:35)', predicted_pw, "filled", "MarkerFaceColor", [.89, 0.11, .14], "DisplayName","Predicted MPP")
grid on 
xlim([0 11])

xlabel("Control Number");
ylabel("Residuals");

title("Residual Plot")