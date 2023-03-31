
%% Plot IV Curves
figure(1);
hold on
    
count = 1;
for i = 6:25
    
    filename = strcat('01_10_ter_',num2str(i),'.xlsx');
    data = xlsread(filename);
    txt = ['Control #', num2str(i)];
    plot(data(:,1), data(:,2),'DisplayName', txt, "LineWidth",2);
    count= count +1;   
end
grid on;
xlabel("Voltage (V)");
ylabel("Current (A)");
legend(Location="northwest");
title("IV Curves ");

hold off
%% %% Extracting MPP from each curve
MPP = zeros(20,1);
ISC = zeros(20,1);
VOC = zeros(20,1);
mpp_idx = 1;
for i = 6:25
    filename = strcat('01_10_ter_',num2str(i),'.xlsx');
    data = xlsread(filename);
    [pw_max, idx] = max(data(:,2).*data(:,1));
    VOC(mpp_idx) = max(data(:,1));
    ISC(mpp_idx) = data(1,2);
    MPP(mpp_idx) = pw_max;
    mpp_idx = mpp_idx + 1;
    
end

%% Plot PV Curves w/ MPP
figure(2);
hold on
   
count = 1;
for i = 6:25
    
    filename = strcat('01_10_ter_',num2str(i),'.xlsx');
    data = xlsread(filename);
    txt = ['Control #', num2str(i), ', ', (num2str(MPP(i-5),3)), 'W'];
    plot(data(:,1), data(:,2).*data(:,1),'DisplayName', txt, "LineWidth",2);

    count= count + 1;   
end
grid on;
xlabel("Voltage (V)");
ylabel("Power (W)");
legend(Location="northwest");
title("Power Curves");

hold off

%% Scatter Plot
figure(3)
scatter((6:25)', MPP, "filled", "MarkerFaceColor", [0, 0.2, .5])
grid on 
xlabel("Control");
ylabel("MPP (W)");
title("MPP Values")


%% Residual Plot

file = xlsread("mpp_trials.xlsx");
figure(4)
scatter((6:25)', file(:,4)-file(:,3), "filled", "MarkerFaceColor", [0, 0.2, .5])
grid on 
ylim([-0 .4])
xlabel("Control");
ylabel("MPP Difference from Jan. 8 and Jan. 10 (W)");
title("Residual Plot of MPP Values")