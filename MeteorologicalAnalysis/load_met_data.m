% load WB pier met data
fin = '../data/WrightsvilleBeachPier_2022_0922_1011.csv';
% header:
% "Date","Time (GMT)","Wind Speed (kn)","Wind Dir (deg)","Wind Gust (kn)","Air Temp (°F)","Baro (mb)","Humidity (%)","Visibility (nm)"
%
% open data-file and create pointer on file system:
fid    = fopen(fin);
header = fgetl(fid);
% scan file for data:
data= textscan(fid,'%s %s %f %f %f %f %f %f %f','delimiter',{',','"'},'treatasempty','-','multipleDelimsAsOne',true);
%
% create time vector
WBn       = size(data{1},1);
WBtString = cat(2,char(data{1}),repmat(' ',WBn,1),char(data{2}));
WBt       = datenum(tString)-4/24;% GMT to EDT
%
% get wind speed/direction
WBwindSpeed = data{3};% knots:
WBwindDir   = data{4};
%
WBwindSpeedmph = 1.15078*windSpeed;% miles per hour
WBwindSpeedmps = 0.51444*windSpeed;% meters per second
%
%
% load MBNB-North met data
fin = '../data/MasonboroNorth_MetStation.csv';
% open data-file and create pointer on file system:
fid    = fopen(fin);
% scan file for data:
data= textscan(fid,'%s %f %s %f %s','headerlines',10,'delimiter',{','});
%
% create time vector
MSNBtString = char(data{1});
MSNBt       = datenum(tString);% EDT
%
% get wind speed/direction
MSNBwindSpeed = data{2};% knots:
MSNBwindDir   = data{4};
%
MSNBwindSpeedmph = 1.15078*windSpeed;% miles per hour
MSNBwindSpeedmps = 0.51444*windSpeed;% meters per second
%
%
fig1 = figure;
ax1 = subplot(2,1,1);
plot(WBt,WBwindSpeedmps,'-r',MSNBt,MSNBwindSpeedmps,'-b')
ylabel(ax1,'Speed [m/s]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','tickdir','out','fontsize',15)
datetick('x','mm/dd')
ax2 = subplot(2,1,2);
plot(WBt,WBwindDir,'-r',MSNBt,MSNBwindDir,'-b')
datetick('x','mm/dd')
ylabel(ax2,'Direction [$^\circ$]','interpreter','latex')
xlabel(ax2,'date [mm/dd]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out','fontsize',15)
%
fig2 = figure;
ax3 = subplot(2,1,1);
plot(WBt,WBwindSpeedmps.*sind(WBwindDir),'-r',MSNBt,MSNBwindSpeedmps.*sind(MSNBwindDir),'-b')
ylabel(ax3,'East [m/s]','interpreter','latex')
set(ax3,'ticklabelinterpreter','latex','tickdir','out','fontsize',15)
datetick('x','mm/dd')
ax4 = subplot(2,1,2);
plot(WBt,-WBwindSpeedmps.*cosd(WBwindDir),'-r',MSNBt,-MSNBwindSpeedmps.*cosd(MSNBwindDir),'-b')
datetick('x','mm/dd')
ylabel(ax4,'North [m/s]','interpreter','latex')
xlabel(ax4,'date [mm/dd]','interpreter','latex')
set(ax4,'ticklabelinterpreter','latex','tickdir','out','fontsize',15)