% This script was written by Qicheng Tang to process Shale Hills czo sap flux data
% Unfortunately, current open-source softwares, including
% Baseliner(Matlab), TREX (R)
% and Aquaflux (R), failed to calculate sap flux on a daily basis. Thus, we
% still need to calculate sap flux 'by hand'. 
% WARNING: THIS CODE CANNOT SOLVE DAILY DRIFT
% The code 1. read the excel data; 2. QAQC, plot the data; 3. Generate the
% daily total sap flux, maximum flux and export them to an excel sheet; 
% Set directory

cd("F:\Projects\Chapter 3\Data");

% Read excel spreadsheet
Shale_Hills_OTT = readtable('Master Data Sheet Chapter 3.xlsx','Sheet','RTH1_ten_min');


% Set time period that you are interested of

% SDateString = '01-June-2016';
% formatIn = 'dd-mmm-yyyy';
% SD = datenum(SDateString,formatIn);
% 
% EDateString = '01-January-2017';
% formatIn = 'dd-mmm-yyyy';
% ED = datenum(EDateString,formatIn);

%% Caculate vapor pressure deficitclc
a = Shale_Hills_OTT{:,8};
b = Shale_Hills_OTT{:,9};
e = Shale_Hills_OTT{:,10};
c = Shale_Hills_OTT{:,15};
c= str2double(c);
d = Shale_Hills_OTT{:,16};
d= str2double(d);
f = Shale_Hills_OTT{:,17};
f = str2double(f);

[m,n] = find(isnan(a));
[p q] = size(a);
Shift = m(:,1);
Row = p;
clear m n p q;
Air_Temp = NaN(Row,1);
Air_RH = NaN(Row,1);
Wind_Speed = NaN(Row,1);


for i = 1:1:Row

    if i<Shift
    Air_Temp(i,1) = a(i,1);
    Air_RH(i,1) = b(i,1);
    Wind_Speed(i,1) = e(i,1);
    else
    Air_Temp(i,1) = c(i,1);
    Air_RH(i,1) = d(i,1);
    Wind_Speed(i,1) = f(i,1);
    end
    
end


OTT_Date = Shale_Hills_OTT{:,1};
OTT_Date = OTT_Date - 4/24; % Daylight saving, 4 hours, from March to November
VPD = 0.61078.*exp(17.269*Air_Temp./(273.3+Air_Temp)).*(1-Air_RH/100);
plot(OTT_Date,VPD);
saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 1 VPD.png']);
clear a b c d i Shift;
%% Sap flux calculation

Shale_Hills_sap = readtable('Master Data Sheet Chapter 3.xlsx','Sheet','SPMS_Sap_Flow');
Garner_Run_sap = readtable('Master Data Sheet Chapter 3.xlsx','Sheet','LRMS_Alt_SapFlow');

a = Shale_Hills_sap{:,3};
b = Shale_Hills_sap{:,4};
c = Shale_Hills_sap{:,5};
d = Shale_Hills_sap{:,6};
e = Shale_Hills_sap{:,11};
e = str2double(e);
f = Shale_Hills_sap{:,12};
f = str2double(f);
g = Shale_Hills_sap{:,13};
g = str2double(g);
h = Shale_Hills_sap{:,14};
h = str2double(h);

[m,n] = find(isnan(a));
Shift = m(:,1);
[p q] = size(a);
Row = p;
clear m n p q;
SH_Tree_1 = NaN(Row,1);
SH_Tree_2 = NaN(Row,1);
SH_Tree_3 = NaN(Row,1);
SH_Tree_4 = NaN(Row,1);

% Import data
for i = 1:1:Row
    if i<Shift
    
        SH_Tree_1(i,1) = a(i,1);
        SH_Tree_2(i,1) = b(i,1);
        SH_Tree_3(i,1) = c(i,1);
        SH_Tree_4(i,1) = d(i,1);
        
        
    else
    
       SH_Tree_1(i,1) = e(i,1);
        SH_Tree_2(i,1) = f(i,1);
        SH_Tree_3(i,1) = g(i,1);
        SH_Tree_4(i,1) = h(i,1);
        
        
    end  
end

% QAQC
for i =1:1:Row
    
    if SH_Tree_1(i,1)<1 || SH_Tree_1(i,1)>30 
        SH_Tree_1(i,1) = NaN;
    end
    
    if SH_Tree_2(i,1)<1 || SH_Tree_2(i,1)>30 
        SH_Tree_2(i,1) = NaN;
    end
    
    if SH_Tree_3(i,1)<1 || SH_Tree_3(i,1)>30 
        SH_Tree_3(i,1) = NaN;
    end
    
    if SH_Tree_4(i,1)<1 || SH_Tree_4(i,1)>30 
        SH_Tree_4(i,1) = NaN;
    end
    
end

SH_Sap_Date = Shale_Hills_sap{:,1};
SH_Sap_Date = SH_Sap_Date - 4/24;
plot(SH_Sap_Date,SH_Tree_1,SH_Sap_Date,SH_Tree_2,SH_Sap_Date,SH_Tree_3,SH_Sap_Date,SH_Tree_4);
xlabel('Date');
ylabel('sap flux');
saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 2 Shale Hills sap flux all year plot.png']);
legend('Tree 1','Tree 2','Tree 3','Tree 4');
clear a b c d e f g h i Shift Row;


GR_Tree_1 = Garner_Run_sap{:,2};
GR_Tree_2 = Garner_Run_sap{:,3};
GR_Tree_3 = Garner_Run_sap{:,4};
GR_Tree_4 = Garner_Run_sap{:,5};

[p q] = size(GR_Tree_1);
Row = p;
clear p q;

for i =1:1:Row
    
    if GR_Tree_1(i,1)<1 || GR_Tree_1(i,1)>30 
        GR_Tree_1(i,1) = NaN;
    end
    
    if GR_Tree_2(i,1)<1 || GR_Tree_2(i,1)>30 
        GR_Tree_2(i,1) = NaN;
    end
    
    if GR_Tree_3(i,1)<1 || GR_Tree_3(i,1)>30 
        GR_Tree_3(i,1) = NaN;
    end
    
    if GR_Tree_4(i,1)<1 || GR_Tree_4(i,1)>30 
        GR_Tree_4(i,1) = NaN;
    end
    
end


GR_Sap_Date = Garner_Run_sap{:,1};
GR_Sap_Date = GR_Sap_Date - 4/24;

plot(GR_Sap_Date,GR_Tree_1,GR_Sap_Date,GR_Tree_2,GR_Sap_Date,GR_Tree_3,GR_Sap_Date,GR_Tree_4);
xlabel('Date');
ylabel('sap flux');
saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 3 Garner Run sap flux all year plot.png']);
legend('Tree 1','Tree 2','Tree 3','Tree 4');
clear Row i;
%% Calculate hourly VPD
SDateString = '01-June-2019';
formatIn = 'dd-mmm-yyyy';
SD = datenum(SDateString,formatIn);

EDateString = '01-January-2020';
formatIn = 'dd-mmm-yyyy';
ED = datenum(EDateString,formatIn);


OTT_Date = datenum(OTT_Date);
SSD = find(OTT_Date(:,1)==SD);
SED = find(OTT_Date(:,1)==ED);
leng = SED - SSD;
OTT_Hour_ID = NaN(leng,2);
Residual = NaN(leng,1);


% First sort out hourly data

for i = 1:1:leng

Residual(i,1) = OTT_Date(SSD+i-1,1)-floor(OTT_Date(SSD+i-1,1));
if mod(round(Residual(i,1)*24,4),1) == 0
    OTT_Hour_ID(i,1) = SSD+i-1;      
    if abs((OTT_Date(SSD+i-1,1)+5/144-OTT_Date(SSD+i+4,1))*1000)<1
        OTT_Hour_ID(i,2) = 1;
    end
end
end

OTT_Hour_ID = OTT_Hour_ID(~isnan(OTT_Hour_ID(:,2)));
m = size(OTT_Hour_ID);
m = m(1,1);
VPD_hourly = NaN(m,2);

for i = 1:1:m

VPD_hourly(i,1) = OTT_Date(OTT_Hour_ID(i,1),1);
VPD_hourly(i,2) = VPD(OTT_Hour_ID(i,1),1);

end
%% Calculate Days with all VPD data


[p q] = size(OTT_Hour_ID);
VPD_Daily_Date = NaN(p,1);
n = 0;

for i = 1:1:(p-23)
temp_1 = OTT_Hour_ID(i,1);
temp_2 = OTT_Hour_ID(i+23,1);
if mod(round(OTT_Date(temp_1,1),4),1) == 0
    n = n+1;
    if abs((OTT_Date(temp_1,1)+23/24-OTT_Date(temp_2,1))*1000) <1
       VPD_Daily_Date(i,1) = OTT_Date(temp_1,1);
    end
end

end

VPD_Daily_Date = VPD_Daily_Date(~isnan(VPD_Daily_Date));
%% Calculate Daily Shale Hills and Garner Run data

SH_Sap_Date = datenum(SH_Sap_Date);
GR_Sap_Date = datenum(GR_Sap_Date);
SH_SSD = find(SH_Sap_Date(:,1)==SD);
SH_SED = find(SH_Sap_Date(:,1)==ED);
GR_SSD = find(GR_Sap_Date(:,1)==SD);
GR_SED = find(GR_Sap_Date(:,1)==ED);

SH_length = SH_SED - SH_SSD;
GR_length = GR_SED - GR_SSD;
SH_Daily_Date = NaN(SH_length,1);
GR_Daily_Date = NaN(GR_length,1);

for i = 1:1:SH_length-23

temp_1 = SH_Sap_Date(SH_SSD+i-1,1);
temp_2 = SH_Sap_Date(SH_SSD+i+22,1);
    
if mod(round(temp_1,4),1) == 0    
    if abs((temp_1+23/24-temp_2)*1000) <1
       SH_Daily_Date(i,1) = temp_1;
    end
end

end


for i = 1:1:GR_length-23

temp_1 = GR_Sap_Date(GR_SSD+i-1,1);
temp_2 = GR_Sap_Date(GR_SSD+i+22,1);
    
if mod(round(temp_1,4),1) == 0    
    if abs((temp_1+23/24-temp_2)*1000) <1
       GR_Daily_Date(i,1) = temp_1;
    end
end

end

SH_Daily_Date = SH_Daily_Date(~isnan(SH_Daily_Date));
GR_Daily_Date = GR_Daily_Date(~isnan(GR_Daily_Date));
%% Find overlap dates

leng = size(SH_Daily_Date);
leng = leng(1,1);
Store = NaN(leng,3);

n = 0; % Number of overlapping dates

Common_dates = NaN(leng,1);

for i = 1:1:leng

a = SH_Daily_Date(i,1);

% a_1 = find(SH_Daily_Date == a);
% a_2 = find(GR_Daily_Date == a);
% a_3 = find(VPD_Daily_Date == a);
   
if find(SH_Daily_Date == a) & find(GR_Daily_Date == a) & find(VPD_Daily_Date == a)
    Common_dates(i,1) = SH_Daily_Date(i,1);
    n = n+1;
end

end


Common_dates = Common_dates(~isnan(Common_dates));

% Create a matrix to store all selected data, vpd, sap shale hills, sap
% garner run, in hourly time format
%% Calculate sap flux, modified from SSHCZOsap_2020.m

SH_T_max = NaN(ED-SD+1,4);
SH_Dsap = NaN(ED-SD+1,4);
SH_Hsap = NaN(SH_SED-SH_SSD+1,4);
SH_inter = NaN(ED-SD+1,4);
SH_Day = NaN(ED-SD+1,1);

GR_T_max = NaN(ED-SD+1,4);
GR_Dsap = NaN(ED-SD+1,4);
GR_Hsap = NaN(GR_SED-GR_SSD+1,4);
GR_inter = NaN(ED-SD+1,4);
GR_Day = NaN(ED-SD+1,1);

m = 1:1:8;

for i = 1:1:(ED-SD+1) % Number of days for Shale Hills
   
   
   % Store data inside cells
   
   temp_SS = find(SH_Sap_Date(:,1)==(SD+i-1));
   temp_SE = find(SH_Sap_Date(:,1)==(SD+i));
   temp_VS = find(VPD_hourly(:,1)==(SD+i-1));
   temp_VE = find(VPD_hourly(:,1)==(SD+i));
   
   
   if (temp_VE - temp_VS) == (temp_SE - temp_SS)
       temp = [SH_Tree_1(temp_SS:temp_SE,1) ...
           SH_Tree_2(temp_SS:temp_SE,1) ...
           SH_Tree_3(temp_SS:temp_SE,1) ...
           SH_Tree_4(temp_SS:temp_SE,1) ...
           VPD_hourly(temp_VS:temp_VE,2)]; 
       SH_Day(i,1) = SD+i-1;
   else
       continue;
   end
   
   
   for n = 1:1:4
       % Sensor order
           % From 0:00 am to 7:00 am 
           coefs = polyfit(temp(m,5),temp(m,n),1);
           SH_inter(i,n) = coefs(2);
           SH_T_max(i,n) = max(SH_inter(i,n),max(temp(:,n)));
           SH_Dsap(i,n) = sum(0.0119*(SH_T_max(i,n)./temp(:,n)-1).^1.23);
           for j = temp_SS:1:temp_SE
           SH_Hsap(j,n) = 0.0119*(SH_T_max(i,n)/temp(j-temp_SS+1,n)-1)^1.23;
           end
   end
    
end

plot(SH_Day,SH_Dsap);
datetick('x','mmm')
xlabel('Date');
ylabel('Sap flux');
legend('Tree 1','Tree 2','Tree 3','Tree 4');
saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 4 Shale Hills daily sap flux.png']);




for i = 1:1:(ED-SD+1) % Number of days for Shale Hills
   
   
   % Store data inside cells
   
   temp_SS = find(GR_Sap_Date(:,1)==(SD+i-1));
   temp_SE = find(GR_Sap_Date(:,1)==(SD+i));
   temp_VS = find(VPD_hourly(:,1)==(SD+i-1));
   temp_VE = find(VPD_hourly(:,1)==(SD+i));
   
   
   if (temp_VE - temp_VS) == (temp_SE - temp_SS)
       temp = [GR_Tree_1(temp_SS:temp_SE,1) ...
           GR_Tree_2(temp_SS:temp_SE,1) ...
           GR_Tree_3(temp_SS:temp_SE,1) ...
           GR_Tree_4(temp_SS:temp_SE,1) ...
           VPD_hourly(temp_VS:temp_VE,2)]; 
       GR_Day(i,1) = SD+i-1;
   else
       continue;
   end
   
   
   for n = 1:1:4
       % Sensor order
           % From 0:00 am to 7:00 am 
           coefs = polyfit(temp(m,5),temp(m,n),1);
           GR_inter(i,n) = coefs(2);
           GR_T_max(i,n) = max(GR_inter(i,n),max(temp(:,n)));
           GR_Dsap(i,n) = sum(0.0119*(GR_T_max(i,n)./temp(:,n)-1).^1.23);
           for j = temp_SS:1:temp_SE
           GR_Hsap(j,n) = 0.0119*(GR_T_max(i,n)/temp(j-temp_SS+1,n)-1)^1.23;
           end
   end
    
end

plot(GR_Day,GR_Dsap);
datetick('x','mmm')
xlabel('Date');
ylabel('Sap flux');
legend('Tree 1','Tree 2','Tree 3','Tree 4');
saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 5 Garner Run daily sap flux.png']);
%% Sort out RAD data that overlaps with the above results, store in Common_dates

Shale_Hills_RAD = readtable('Master Data Sheet Chapter 3.xlsx','Sheet','EC_Rad');
RAD_Date = Shale_Hills_RAD{:,1};
RAD = Shale_Hills_RAD{:,7};

RAD_Date = datenum(RAD_Date);
RSD = find(RAD_Date(:,1) == SD);
RED = find(RAD_Date(:,1) == ED);
leng = RED - RSD;
RAD_Hour_ID = NaN(leng,2);
Residual = NaN(leng,1);


for i = 1:1:leng
Residual(i,1) = RAD_Date(RSD+i-1,1)-floor(RAD_Date(RSD+i-1,1));
if mod(round(Residual(i,1)*24,4),1) == 0
    RAD_Hour_ID(i,1) = RSD+i-1;      
    if abs((RAD_Date(RSD+i-1,1)+5/144-RAD_Date(RSD+i+4,1))*1000)<1
        RAD_Hour_ID(i,2) = 1;
    end
end
end

RAD_Hour_ID = RAD_Hour_ID(~isnan(RAD_Hour_ID(:,2)));
m = size(RAD_Hour_ID);
m = m(1,1);
RAD_hourly = NaN(m,2);

for i = 1:1:m

RAD_hourly(i,1) = RAD_Date(RAD_Hour_ID(i,1),1);
RAD_hourly(i,2) = VPD(RAD_Hour_ID(i,1),1);

end

[p q] = size(RAD_Hour_ID);
RAD_Daily_Date = NaN(p,1);
n = 0;

for i = 1:1:(p-23)
temp_1 = RAD_Hour_ID(i,1);
temp_2 = RAD_Hour_ID(i+23,1);
if mod(round(RAD_Date(temp_1,1),4),1) == 0
    n = n+1;
    if abs((RAD_Date(temp_1,1)+23/24-RAD_Date(temp_2,1))*1000) <1
       RAD_Daily_Date(i,1) = RAD_Date(temp_1,1);
    end
end

end

RAD_Daily_Date = RAD_Daily_Date(~isnan(RAD_Daily_Date));

leng = size(Common_dates);
leng = leng(1,1);
% Store = NaN(leng,3);

for i = 1:1:leng

a = Common_dates(i,1);
   
if find(RAD_Daily_Date(:,1) == a)
     Common_dates(i,1) = a;
else
    Common_dates(i,1) = NaN;
end

end


Common_dates = Common_dates(~isnan(Common_dates));
%% Caculate PET based on overlapped dates
% T_stamp_S = All_S{:,1};
% T_stamp_S = datenum(T_stamp_S);

leng = size(Common_dates);
leng = leng(1,1);
All_S = NaN(leng*144,5);

for i = 1:1:leng
    b = Common_dates(i,1);
    for j = 1:1:144
    a = find(OTT_Date(:,1)== b);
    a_2 = find(RAD_Date(:,1)==b);
    All_S(144*(i-1)+j,1) = OTT_Date(a+j-1,1);
    All_S(144*(i-1)+j,2) = Wind_Speed(a+j-1,1);
    All_S(144*(i-1)+j,3) = Air_Temp(a+j-1,1);
    All_S(144*(i-1)+j,4) = Air_RH(a+j-1,1);
    All_S(144*(i-1)+j,5) = RAD(a_2+j-1,1);
    end
end


T_stamp_S = All_S(:,1);
u_2_S = All_S(:,2);
Air_T_S = All_S(:,3);
RH_Max_S = All_S(:,4);
% RH_Min = All{:,11};
S_R_S = All_S(:,5);
% RH = (RH_Max+RH_Min)/2;
%This has to be replaced by daily air temperature

A_P_S = 101.3*((293-0.0065*500)/293)^5.26;
Psy_C_S = 0.000665*A_P_S;

N_S = size(Air_T_S);
N_S = N_S(1,1);
N_S = N_S/144;


Max_T_S = NaN(N_S,1);
Min_T_S = NaN(N_S,1);
Ave_T_S = NaN(N_S,1);
M_SR_S = NaN(N_S,1);
SSVP_S = NaN(N_S,1); % Slope of saturation vapor pressure
M_WS_S = NaN(N_S,1);
e_s_S = NaN(N_S,1);
e_a_S = NaN(N_S,1);

ET_S = NaN(N_S,1);
Day_S = NaN(N_S,1);


for i = 1:1:N_S
   
       Day_S(i,1) = T_stamp_S(144*(i-1)+1,1);    
       
       M_SR_S(i,1) = mean(S_R_S(144*(i-1)+1:144*(i-1)+144,1));
    
       Max_T_S(i,1) = max(Air_T_S(144*(i-1)+1:144*(i-1)+144,1));
       Min_T_S(i,1) = min(Air_T_S(144*(i-1)+1:144*(i-1)+144,1));
       Ave_T_S(i,1) = (Max_T_S(i,1)+Min_T_S(i,1))/2;
       
       SSVP_S(i,1) = 4098*0.61078*exp(17.269*Ave_T_S(i,1)/(273.3+Ave_T_S(i,1)))/(Ave_T_S(i,1)+237.3)^2;
       
       M_WS_S(i,1) = mean(u_2_S(144*(i-1)+1:144*(i-1)+144,1));
       
       e_max = 0.6108*exp(17.27*(Max_T_S(i,1))/(Max_T_S(i,1)+237.3));
       e_min = 0.6108*exp(17.27*(Min_T_S(i,1))/(Min_T_S(i,1)+237.3));
       
       e_s_S(i,1) = (e_max+e_min)/2;
       e_a_S(i,1) = e_min*(max(RH_Max_S(144*(i-1)+1:144*(i-1)+144,1)/100));
       
       ET_S(i,1) = (0.408*SSVP_S(i,1)*M_SR_S(i,1)+Psy_C_S*(900/(Ave_T_S(i,1)+273))...
           *M_WS_S(i,1)*(e_s_S(i,1)-e_a_S(i,1)))/(SSVP_S(i,1)+Psy_C_S*(1+0.34*M_WS_S(i,1))); 
      
end

% plot(Day_S,ET_S);
% hold;
%% Plot PET vs. Sap flux

Shale_Hills_Sap_Mean = NaN(length(Common_dates(:,1)),1);
Garner_Run_Sap_Mean = NaN(length(Common_dates(:,1)),1);

%% Select trees that you'd like to plot

j = [1 2 4]; % 1 2 3 4 are tree numbers

for i = 1:1:length(Common_dates(:,1))
b = Common_dates(i,1);
a = find(SH_Day(:,1)==b);
if a
Shale_Hills_Sap_Mean(i,1) = mean(SH_Dsap(a,j));
Garner_Run_Sap_Mean(i,1) = mean(GR_Dsap(a,j));
end

end



figure()
yyaxis left

plot(Common_dates,Shale_Hills_Sap_Mean,'r-',Common_dates,Garner_Run_Sap_Mean,'b-');

set(gca,'YColor','k');

% ylim([0 0.2])

ylabel('Daily sap flux (cm/s)');

yyaxis right

plot(Common_dates,ET_S,'k');

set(gca,'YColor','k')

datetick('x','dd-mm','keepticks');

legend('Shale Hills sap','Garner Run sap','PET')

xlabel('Month');

ylabel('Potential Evapotranspiration (mm/day)')

set(gca,'FontSize',14);

saveas(gcf,['F:/Projects/Chapter 3/Outputs/Figure 6 ET vs. sap flux.png']);
% legend('Tree 1','Tree 2','Tree 3','Tree 4');