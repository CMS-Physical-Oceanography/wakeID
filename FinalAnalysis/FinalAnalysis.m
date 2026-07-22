data = load("C:\Users\cls7049\OneDrive - UNC-Wilmington\WakeProjData\Summer2026\SpectrogramTrainingData\final_boat_parameters.mat");

%% Figure 1: Distance from camera Vs. Distance from pressure sensor:

% 1. Extract the nested values into standard numeric arrays using arrayfun
slopes = arrayfun(@(s) s.Pressure.WakeSlope_dt_df_s_per_Hz, data.finalBoatParameters);
cameraDistances = arrayfun(@(s) s.Camera.DistanceFromPressureFt, data.finalBoatParameters);

% 2. Calculate the distance measured by the pressure sensor
dbypressure = slopes .* 2.562;

% 3. Calculate standard deviation and mean for both values. Use these to
% compute max and min range
STDcutoff = 2;

Cdsig = std(cameraDistances);
Pdsig = std(dbypressure);
Cdmean = mean(cameraDistances);
Pdmean = mean(dbypressure);

Cdmax = Cdmean + STDcutoff*Cdsig;
Cdmin = Cdmean - STDcutoff*Cdsig;
Pdmax = Pdmean + STDcutoff*Pdsig;
Pdmin = Pdmean - STDcutoff*Pdsig;

% 4. Filter outliers by removing data that is outside of 2 standard
% deviations

dummy = ones(size(cameraDistances));
for i = 1:length(cameraDistances)
    if cameraDistances(i) <= Cdmin || cameraDistances(i) >= Cdmax ...
            || dbypressure(i) <= Pdmin || dbypressure(i) >= Pdmax 
        dummy(i) = 0;
    end
end 
cameraDistances = cameraDistances(dummy == 1);
dbypressure = dbypressure(dummy == 1);

% 5. Generate and the fit line
m = dbypressure\cameraDistances;

xfit = linspace(min(dbypressure), max(dbypressure), 100);
yfit = m * xfit; % Calculate y-values using only the slope (intercept is 0)

% 6. Plot the data
figure
scatter(dbypressure, cameraDistances, 8, 'filled', 'o')
hold on

plot(xfit, yfit, 'r-', 'LineWidth', 1)
xlim([0 max(dbypressure)])
ylim([0 max(cameraDistances)])
axis equal
grid on

xlabel("Distance by pressure sensor (ft)", "FontSize",15)
ylabel("Distance by camera (ft)","FontSize",15)
title("Camera Distance Vs. Pressure Sensor Distance","FontSize",20)
