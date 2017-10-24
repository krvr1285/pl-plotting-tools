% Plot PL and UV-vis data
close all

cd '/Users/Kristina/Google Drive/Dukovic group/Data/PL/171018'

%% Importing all files
%Sample 1 Emission
    S11=dlmread('1_1.txt','',3,0);
    S12=dlmread('1_2.txt','',3,0);
    S13=dlmread('1_3.txt','',3,0);
    %Excitation
    S11ex=dlmread('1_1ex.txt','',3,0);
    S12ex=dlmread('1_2ex.txt','',3,0);
    S13ex=dlmread('1_3ex.txt','',3,0);
%Sample 2 Emission
    S21=dlmread('2_1.txt','',3,0);
    S22=dlmread('2_2.txt','',3,0);
    S23=dlmread('2_3.txt','',3,0);
    %Excitation
    S21ex=dlmread('2_1ex.txt','',3,0);
    S22ex=dlmread('2_2ex.txt','',3,0);
    S23ex=dlmread('2_3ex.txt','',3,0);   
%Background scans
    S02=dlmread('0_2.txt','',3,0);
    S03=dlmread('0_3.txt','',3,0);
    %612 n m
    S02ex1=dlmread('0_2ex1.txt','',3,0);
    S03ex1=dlmread('0_3ex1.txt','',3,0);
    %609 nm
    S02ex2=dlmread('0_2ex2.txt','',3,0);
    S03ex2=dlmread('0_3ex2.txt','',3,0);
    %this background scan was copied from data from 171012
    S01 = dlmread('bufferem450.txt','',3,0);
    S01ex2 = dlmread('bufferex609.txt','',3,0);
    
%UV-vis data
cd '/Users/Kristina/Google Drive/Dukovic group/Data/UV-vis/171018'
S11uv= csvread('1_1.csv',3,0);
S12uv = csvread('1_2.csv',3,0);
S13uv= csvread('1_3.csv',3,0);
S21uv=csvread('2_1.csv',3,0);
S22uv=csvread('2_2.csv',3,0);
S23uv=csvread('2_3.csv',3,0);
S01uv=csvread('0_1.csv',3,0);
S02uv=csvread('0_2.csv',3,0);
S03uv=csvread('0_3.csv',3,0);
cd '/Users/Kristina/Google Drive/Dukovic group/Data/PL/171018'

%% Background correction
%Sample 1 Emission
    S11corr=S11(:,2)-S01(:,2);
    S12corr=S12(:,2)-S02(:,2);
    S13corr=S13(:,2)-S03(:,2);
    %Excitation
    S11excorr=S11ex(:,2) -S01ex2(:,2);
    S12excorr=S12ex(:,2) - S02ex2(:,2);
    S13excorr=S13ex(:,2)-S03ex2(:,2);
%Sample 2 Emission
    S21corr=S21(:,2)-S01(:,2);
    S22corr=S22(:,2)-S02(:,2);
    S23corr=S23(:,2)-S03(:,2);
    %Excitation
    S21excorr=S21ex(:,2)-S01ex2(:,2); %Don't have a 612 nm scan for the buffer ex
    S22excorr=S22ex(:,2)-S02ex1(:,2);
    S23excorr=S23ex(:,2)-S03ex1(:,2);

%% UV vis correction
%S1
S11uvcorr=S11uv(:,2)-S01uv(:,2);
S12uvcorr=S12uv(:,2)-S02uv(:,2);
S13uvcorr=S13uv(:,2)-S03uv(:,2);
%S2
S21uvcorr=S21uv(:,2)-S01uv(:,2);
S22uvcorr=S22uv(:,2)-S02uv(:,2);
S23uvcorr=S23uv(:,2)-S03uv(:,2);

%% Fig 1 UV vis data
figure(1)
hold on

p1=plot(S11uv(:,1),S11uvcorr,S12uv(:,1),S12uvcorr,S13uv(:,1),S13uvcorr,...
    S21uv(:,1),S21uvcorr,S22uv(:,1),S22uvcorr,S23uv(:,1),S23uvcorr);
%Change all lines to be Line Width 2
k=1;
for k=1:length(p1)
    p1(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3','S2_1','S2_2','S2_3')
xlabel('Wavelength (nm)')
ylabel('Absorbance')
axis([250 650 -0.02 0.6])
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig1.tiff')

hold off

%% Fig 2 Raw S1 and S2 data
figure(2)
hold on

p2=plot(S11(:,1),S11(:,2),S12(:,1),S12(:,2),S13(:,1),S13(:,2),S21(:,1),S21(:,2),...
    S22(:,1),S22(:,2),S23(:,1),S23(:,2));
k=1;
for k=1:length(p2)
    p2(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3','S2_1','S2_2','S2_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig2.tiff')

hold off

%% Fig 3 CdSe/CdS-MPA Data (raw), AA backgrounds
figure(3)
hold on

p3=plot(S21(:,1),S21(:,2),S22(:,1),S22(:,2),S23(:,1),S23(:,2),S01(:,1),S01(:,2),S02(:,1),S02(:,2),S03(:,1),S03(:,2));

k=1;
for k=1:length(p3)
    p3(k).LineWidth=2;
    k=k+1;
end

legend('S2_1','S2_2','S2_3','S0_1','S0_2','S0_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig3.tiff')

hold off

%% Fig 4 Both sets of Data, bkg subtracted
figure(4)
hold on

p2=plot(S11(:,1),S11corr,S12(:,1),S12corr,S13(:,1),S13corr,S21(:,1),S21corr,...
    S22(:,1),S22corr,S23(:,1),S23corr);


for k=1:length(p2)
    %p2(k).Marker='.';
    p2(k).LineWidth=2;
    
    k=k+1;
end

legend('S1_1','S1_2','S1_3','S2_1','S2_2','S2_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig4.tiff')

hold off

%% Fig 5 Excitation Data
figure(5)
hold on

p2=plot(S11ex(:,1),S11excorr,S12ex(:,1),S12excorr,S13ex(:,1),S13excorr,S21ex(:,1),S21excorr,S22ex(:,1),S22excorr,S23ex(:,1),S23excorr);

k=1;
for k=1:length(p2)
    p2(k).LineWidth=2;
    k=k+1;
end

legend('S1_1','S1_2','S1_3','S2_1','S2_2','S2_3')
xlabel('Wavelength (nm)')
ylabel('CPS')
%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig5.tiff')

hold off

%% Fig 6 Peak Fits
figure(6)
hold on

p2=plot(S11(:,1),S11corr,S12(:,1),S12corr,S13(:,1),S13corr,S21(:,1),S21corr,...
    S22(:,1),S22corr,S23(:,1),S23corr);
%Find peaks, raw data
[~,locsS11]=findpeaks(S11corr,'MinPeakDistance',300);
[~,locsS12]=findpeaks(S12corr,'MinPeakDistance',300);
[~,locsS13]=findpeaks(S13corr,'MinPeakDistance',300);
[~,locsS21]=findpeaks(S21corr,'MinPeakDistance',300);
[~,locsS22]=findpeaks(S22corr,'MinPeakDistance',300);
[~,locsS23]=findpeaks(S23corr,'MinPeakDistance',300);

p2b = plot(S11(locsS11,1),S11corr(locsS11),S12(locsS12,1),S12corr(locsS12),S13(locsS13,1),...
S13corr(locsS13),S12(locsS21,1),S21corr(locsS21),S21(locsS22,1),S22corr(locsS22),S23(locsS23,1),S23corr(locsS23));
%Change plotting parameters of plot2 and plot2b (figure data and peak
%finding)

for k=1:length(p2)
    p2(k).Marker='x';
    p2(k).MarkerSize=2;
    p2b(k).Marker='*';
    p2b(k).MarkerSize=3;
    k=k+1;
end
%Fit peaks to a gaussian (gauss2) and put fits into a cell
fits={fit(S11(27:127,1),S11corr(27:127),'gauss2') fit(S12(27:127,1),S12corr(27:127),'gauss2'),...
    fit(S13(27:127,1),S13corr(27:127),'gauss2'),fit(S21(27:127,1),S21corr(27:127),'gauss2')...
    fit(S22(27:127,1),S22corr(27:127),'gauss2'),fit(S23(27:127,1),S23corr(27:127),'gauss2')};
% plot fits

for k=1:length(fits)
    p3(k)=plot(fits{k});
   
    p3(k).LineWidth=0.5;
    p3(k).Color ='k';
    
end   
% for k=1:length(fits)
%        p4(k)=plot(((fits{k}.b1*fits{k}.c1+fits{k}.b2*fits{k}.c2)/(fits{k}.c1+fits{k}.c2)),((fits{k}.a1*fits{k}.c1+fits{k}.a2*fits{k}.c2)/(fits{k}.c1+fits{k}.c2)),'Marker','x','MarkerSize',10);
%     k=k+1;
% end

% S11fit=fit(S11(27:127,1),S11corr(27:127),'gauss2');
% S12fit=fit(S12(27:127,1),S12corr(27:127),'gauss2');
% S13fit=fit(S13(27:127,1),S13corr(27:127),'gauss2');
% S21fit=fit(S21(27:127,1),S21corr(27:127),'gauss2');
% S22fit=fit(S22(27:127,1),S22corr(27:127),'gauss2');
% S23fit=fit(S23(27:127,1),S23corr(27:127),'gauss2');

% m1=plot(S11fit);
% m2=plot(S12fit);
% m3=plot(S13fit);
% m4=plot(S21fit);
% m5=plot(S22fit);
% m6=plot(S23fit);
% 
% m1.LineWidth=1;
% m2.LineWidth=1;
% m3.LineWidth=1;
% m4.LineWidth=1;
% m5.LineWidth=1;
% m6.LineWidth=1;
% 
% m1.Color=[0 0 0];
% m2.Color=[0 0 0];
% m3.Color=[0 0 0];
% m4.Color=[0 0 0];
% m5.Color=[0 0 0];
% m6.Color=[0 0 0];



p1=legend('S1_1','S1_2','S1_3','S2_1','S2_2','S2_3');
axis([550 700 0 14e5])
xlabel('Wavelength (nm)')
ylabel('CPS')

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig6.tiff')

hold off

%% Fig 7 Plot I0/I vs AA concentration
%Use results of "findpeaks" for I 
figure(7)
hold on
locsS21=85;
S1I0=S11corr(locsS11);
S2I0=S21corr(locsS21);

I0 = [S1I0/S11corr(locsS11) S1I0/S12corr(locsS12) S1I0/S13corr(locsS13);...
    S2I0/S21corr(locsS21) S2I0/S22corr(locsS22) S2I0/S23corr(locsS23)];
AAconc = [0 0.05 0.1];

for k=1:3
    p4(k)= plot(AAconc(k),I0(1,k),'Marker','x','MarkerSize',5,'Color','r');
    p5(k)=plot(AAconc(k),I0(2,k),'Marker','*','MarkerSize',5,'Color','b');
    k=k+1;
end

%Fit a line to data sets
for k=1:2
    coeffs(k,:) = polyfit(AAconc, I0(k,:), 1);
    % Get fitted values
    fittedX = linspace(min(AAconc), max(AAconc), 200);
    fittedY = polyval(coeffs(k,:), fittedX);
    % Plot the fitted line
    plot(fittedX, fittedY, 'k', 'LineWidth', 0.5);
    
end
labels={'Sample 1' 'slope' num2str(round(coeffs(1,1),2)) 'yint' num2str(round(coeffs(1,2),2));...
    'Sample 2' 'slope' num2str(round(coeffs(2,1),2)) 'yint' num2str(round(coeffs(2,2),2))};
labelS1=strjoin(labels(1,:));
labelS2=strjoin(labels(2,:));
%legend('Sample 1','Sample 2');
legend(labelS1,labelS2);
xlabel('[Ascorbate], M');
ylabel('I_0/I');

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig7.tiff')

hold off