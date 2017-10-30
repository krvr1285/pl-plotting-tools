% Plot PL and UV-vis data
close all

cd '/Users/Kristina/Google Drive/Dukovic group/Data/PL/171026'

%% Importing all files
%Sample 1 Emission
    S11=dlmread('1_1.txt','',2,0);
    S12=dlmread('1_2.txt','',2,0);
    S13=dlmread('1_3.txt','',2,0);
    %Excitation
    S11ex=dlmread('1_1ex.txt','',2,0);
    S12ex=dlmread('1_2ex.txt','',2,0);
    S13ex=dlmread('1_3ex.txt','',2,0);
    
%Background scans
    S01=dlmread('0_1.txt','',2,0);
    S02=dlmread('0_2.txt','',2,0);
    S03=dlmread('0_3.txt','',2,0);
    %Excitation scan
%     S01ex=dlmread('0_1ex','',3,0);
%     S02ex=dlmread('0_2ex.txt','',3,0);
%     S03ex=dlmread('0_3ex.txt','',3,0);

    
%UV-vis data
cd '/Users/Kristina/Google Drive/Dukovic group/Data/UV-vis/171026'
S11uv= csvread('1_1.csv',3,0);
S12uv = csvread('1_2.csv',3,0);
S13uv= csvread('1_3.csv',3,0);

%S01uv=csvread('0_1.csv',3,0);
S02uv=csvread('0_2.csv',3,0);
S03uv=csvread('0_3.csv',3,0);
cd '/Users/Kristina/Google Drive/Dukovic group/Data/PL/171026'

%% Background correction
%Sample 1 Emission
    S11corr=S11(:,2)-S01(:,2);
    S12corr=S12(:,2)-S02(:,2);
    S13corr=S13(:,2)-S03(:,2);
    %Excitation
    S11excorr=S11ex(:,2); %-S01ex(:,2);
    S12excorr=S12ex(:,2); %- S02ex(:,2);
    S13excorr=S13ex(:,2);%-S03ex(:,2);

%% UV vis correction
%S1
S11uvcorr=S11uv(:,2);%-S01uv(:,2);
S12uvcorr=S12uv(:,2)-S02uv(:,2);
S13uvcorr=S13uv(:,2)-S03uv(:,2);


%% Fig 1 UV vis data
figure(1)
hold on

p1=plot(S11uv(:,1),S11uvcorr,S12uv(:,1),S12uvcorr,S13uv(:,1),S13uvcorr);

%Change all lines to be Line Width 2
k=1;
for k=1:length(p1)
    p1(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3')
xlabel('Wavelength (nm)')
ylabel('Absorbance')
axis([250 650 -0.02 0.6])
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig1.tiff')

hold off

%% Fig 2 Raw S1 data
figure(2)
hold on

p2=plot(S11(:,1),S11(:,2),S12(:,1),S12(:,2),S13(:,1),S13(:,2));
k=1;
for k=1:length(p2)
    p2(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig2.tiff')

hold off


%% Fig 3 Both sets of Data, bkg subtracted
figure(3)
hold on

p2=plot(S11(:,1),S11corr,S12(:,1),S12corr,S13(:,1),S13corr);


for k=1:length(p2)
    %p2(k).Marker='.';
    p2(k).LineWidth=2;
    
    k=k+1;
end

legend('S1_1','S1_2','S1_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig4.tiff')

hold off

%% Fig 4 Excitation Data
figure(4)
hold on

p2=plot(S11ex(:,1),S11excorr,S12ex(:,1),S12excorr,S13ex(:,1),S13excorr);
k=1;
for k=1:length(p2)
    p2(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig5.tiff')

hold off

%% Fig 5 Peak Fits
figure(5)
hold on

p2=plot(S11(:,1),S11corr,S12(:,1),S12corr,S13(:,1),S13corr);
%Find peaks, raw data
[~,locsS11]=findpeaks(S11corr,'MinPeakDistance',250);
[~,locsS12]=findpeaks(S12corr,'MinPeakDistance',250);
[~,locsS13]=findpeaks(S13corr,'MinPeakDistance',250);

p2b = plot(S11(locsS11,1),S11corr(locsS11),S12(locsS12,1),S12corr(locsS12),S13(locsS13,1),S13corr(locsS13));
%Change plotting parameters of plot2 and plot2b (figure data and peak
%finding)

for k=1:length(p2)
    p2(k).LineWidth=2;
    p2b(k).Marker='*';
    p2b(k).MarkerSize=3;
    k=k+1;
end
%Fit peaks to a gaussian (gauss2) and put fits into a cell
fits={fit(S11(27:127,1),S11corr(27:127),'gauss2') fit(S12(27:127,1),S12corr(27:127),'gauss2'),...
    fit(S13(27:127,1),S13corr(27:127),'gauss2')};
% plot fits

for k=1:length(fits)
    p3(k)=plot(fits{k});
   
    p3(k).LineWidth=0.5;
    p3(k).Color ='k';
    
end   

p1=legend('S1_1','S1_2','S1_3');
axis([500 675 0 14e6])
xlabel('Wavelength (nm)')
ylabel('CPS')

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig6.tiff')

hold off

%% Fig 6 Plot I0/I vs AA concentration
%Use results of "findpeaks" for I 
figure(6)
hold on

S1I0=S11corr(locsS11);

I0 = [S1I0/S11corr(locsS11) S1I0/S12corr(locsS12) S1I0/S13corr(locsS13)];
AAconc = [0 0.05 0.1];

for k=1:3
    p4(k)= plot(AAconc(k),I0(1,k),'Marker','x','MarkerSize',5);
    k=k+1;
end

%Fit a line to data sets
for k=1:1
    coeffs(k,:) = polyfit(AAconc, I0(k,:), 1);
    % Get fitted values
    fittedX = linspace(min(AAconc), max(AAconc), 200);
    fittedY = polyval(coeffs(k,:), fittedX);
    % Plot the fitted line
    plot(fittedX, fittedY, 'k', 'LineWidth', 0.5);
    
end
labels={'Sample 1' 'slope' num2str(round(coeffs(1,1),2)) 'yint' num2str(round(coeffs(1,2),2))};
labelS1=strjoin(labels(1,:));

legend(labelS1);
xlabel('[Ascorbate], M');
ylabel('I_0/I');

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig7.tiff')

hold off