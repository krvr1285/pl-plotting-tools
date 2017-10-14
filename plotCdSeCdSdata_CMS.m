% Plot CdSe/CdS PL and UV-vis data
close all
%delete(findall(0,'Type','figure'))

%% Importing all files
%CdSe/CdS in toluene
    %Em 1 (450 nm)
    CdSetolem450=dlmread('CdSetolem450.txt','',3,0);
    tolem450=dlmread('tolueneem450.txt','',3,0);
    %Em 2 (450-2 nm)
    CdSetolem4502=dlmread('CdSetolem450_2.txt','',3,0);
    tolem4502=dlmread('tolueneem450_2.txt','',3,0);
    %Em 3 (350 nm)
    CdSetolem350=dlmread('CdSetolem350.txt','',3,0);
    tolem350=dlmread('tolueneem350.txt','',3,0);

    %Ex 1 (612 nm)
    CdSetolex612=dlmread('CdSetolex612.txt','',3,0);
    tolex612=dlmread('tolueneex612.txt','',3,0);

%CdSe/CdS in buffer
    CdSebuffem450=dlmread('CdSebuffem450.txt','',3,0);
    buffem450=dlmread('bufferem450.txt','',3,0);

    CdSebuffem4502=dlmread('CdSebuffem450_2.txt','',3,0);
    buffem4502=dlmread('bufferem450_2.txt','',3,0);

    CdSebuffem350=dlmread('CdSebuffem350.txt','',3,0);
    buffem350=dlmread('bufferem350.txt','',3,0);

    CdSebuffex609=dlmread('CdSebuffex609.txt','',3,0);
    buffex609=dlmread('bufferex609.txt','',3,0);

%CdSe/CdS in buffer, 100 mM AA
    CdSeAAem450=dlmread('CdSeAAem450.txt','',3,0);
    buffem450=dlmread('bufferem450.txt','',3,0);

    CdSeAAem4502=dlmread('CdSeAAem450_2.txt','',3,0);
    buffem4502=dlmread('bufferem450_2.txt','',3,0);

    CdSeAAem350=dlmread('CdSeAAem350.txt','',3,0);
    buffem350=dlmread('bufferem350.txt','',3,0);

    CdSeAAex609=dlmread('CdSeAAex609.txt','',3,0);
    buffex609=dlmread('bufferex609.txt','',3,0);

%Buffer background (with AA)
    buffAAem450=dlmread('buffAAem450.txt','',3,0);
    buffAAem4502=dlmread('buffAAem450_2.txt','',3,0);
    buffAAem350=dlmread('buffAAem350.txt','',3,0);
    buffAAex609=dlmread('buffAAex609.txt','',3,0);

%UV-vis data
cd '/Users/Kristina/Google Drive/Dukovic group/Data/UV-vis/171012'
CdSetolUV= csvread('cdsetol.csv',3,0);
CdSebuffUV = csvread('cdsebuff.csv',3,0);
CdSeAAUV= csvread('CdSeAA.csv',3,0);
cd '/Users/Kristina/Google Drive/Dukovic group/Data/PL/171012'

%% Background correction
%CdSe/CdS in toluene
CdSetolem450bg = CdSetolem450(:,2)- tolem450(:,2);
CdSetolem4502bg = CdSetolem4502(:,2)- tolem4502(:,2);
CdSetolem350bg = CdSetolem350(:,2)- tolem350(:,2);
CdSetolex612bg = CdSetolex612(:,2)- tolex612(:,2);

%CdSe/CdS in buffer
CdSebuffem450bg = CdSebuffem450(:,2) - buffem450(:,2);
CdSebuffem4502bg = CdSebuffem4502(:,2)- buffem4502(:,2);
CdSebuffem350bg = CdSebuffem350(:,2)- buffem350(:,2);
CdSebuffex609bg = CdSebuffex609(:,2)- buffex609(:,2);

%CdSe/CdS in buffer, 100 mM AA
%50 mM buffer subtraction
CdSeAAem450bg = CdSeAAem450(:,2) - buffem450(:,2);
CdSeAAem4502bg = CdSeAAem4502(:,2)- buffem4502(:,2);
CdSeAAem350bg = CdSeAAem350(:,2)- buffem350(:,2);
CdSeAAex609bg = CdSeAAex609(:,2)- buffex609(:,2);

%50 mM buffer with AA subtraction
CdSeAAem450bgsubAA = CdSeAAem450bg - buffAAem450(:,2);

%% Figure 1 Plotting CdSe/CdS in toluene data
hold on
%Emission 1 (450 nm) conditions
plot(CdSetolem450(:,1),CdSetolem450bg,'LineWidth',2)
%Emission 2 (450 nm) conditions
plot(CdSetolem4502(:,1),CdSetolem4502bg,'LineWidth',2)
%Emission 3 (350 nm) conditions
plot(CdSetolem350(:,1),CdSetolem350bg,'LineWidth',2)
%Excitation conditions (612 nm)
plot(CdSetolex612(:,1),CdSetolex612bg,'LineWidth',2)

%Create legend and label axes
legend('CdSe/CdS in toluene, em 450,bg','CdSe/CdS in toluene, em 450-2,bg'...
    ,'CdSe/CdS in toluene, em 350,bg','CdSe/CdS in toluene, ex 612,bg')
xlabel('Wavelength (nm)')
ylabel('CPS')
hold off

%% Figure 2: Plotting CdSe/CdS-S in buffer data
figure(2)
hold on
%Emission 1 (450 nm) conditions
plot(CdSebuffem450(:,1),CdSebuffem450bg,'LineWidth',2)
%Emission 2 (450 nm) conditions
plot(CdSebuffem4502(:,1),CdSebuffem4502bg,'LineWidth',2)
%Emission 2 (350 nm) conditions
plot(CdSebuffem350(:,1),CdSebuffem350bg,'LineWidth',2)
%Excitation conditions (609 nm)
plot(CdSebuffex609(:,1),CdSebuffex609bg,'LineWidth',2)

%Create legend and label axes
legend('CdSe/CdS in buffer, em 450,bg','CdSe/CdS in buffer, em 450-2,bg',...
    'CdSe/CdS in buffer, em 350,bg','CdSe/CdS in buffer, ex 609,bg')
xlabel('Wavelength (nm)')
ylabel('CPS')
hold off


%% Figure 3: Plotting CdSe/CdS-S in buffer with AA data
figure(3)
hold on
%Emission 1 (450 nm) conditions
plot(CdSeAAem450(:,1),CdSeAAem450bg,'LineWidth',2)
%Emission 2 (450 nm) conditions
plot(CdSeAAem4502(:,1),CdSeAAem4502bg,'LineWidth',2)
%Emission 2 (350 nm) conditions
plot(CdSeAAem350(:,1),CdSeAAem350bg,'LineWidth',2)
%Excitation conditions (609 nm)
plot(CdSeAAex609(:,1),CdSeAAex609bg,'LineWidth',2)

% Add AA background data to graphs
%plot all
plot(buffAAem450(:,1),buffAAem450(:,2),'LineWidth',2)
plot(buffAAem4502(:,1),buffAAem4502(:,2),'LineWidth',2)
plot(buffAAem350(:,1),buffAAem350(:,2),'LineWidth',2)
plot(buffAAex609(:,1),buffAAex609(:,2),'LineWidth',2)

%Create legend and label axes
legend('CdSe/CdS in buffer, 100 mM AA, em 450,bg','CdSe/CdS in buffer, 100 mM AA, em 450-2,bg',...
    'CdSe/CdS in buffer, 100 mM AA, em 350,bg','CdSe/CdS in buffer, 100 mM AA, ex 609,bg',...
    'AA buffer, em450','AA buffer,em 450_2','AA buffer,em 350','AA buffer,ex 609')
xlabel('Wavelength (nm)')
ylabel('CPS')
hold off

%% Figure 4 Comparing CdSe/CdS with and without AA
figure(4)
hold on

%Using excitation condition 1 (450 nm) as it's the cleanest data
plot(CdSeAAem450(:,1),CdSeAAem450bg,'LineWidth',2)
plot(CdSebuffem450(:,1),CdSebuffem450bg,'LineWidth',2)
plot(buffAAem450(:,1),buffAAem450(:,2),'LineWidth',2)
plot(CdSebuffem450(:,1),CdSeAAem450bgsubAA,'LineWidth',2)

% %Fit the 2 PL peaks (CdSeAAem450bgsubAA, CdSebuffem450bg)
% CdSebuffem450cut=CdSebuffem450bg(27:127);
% CdSebuffem450cutx=CdSebuffem450((27:127),1);
% CdSebuffem450fit=fit(CdSebuffem450cutx,CdSebuffem450cut,'gauss2')
% CdSeAAem450cut=CdSeAAem450bgsubAA(27:127);
% CdSeAAem450cutx=CdSeAAem450((27:127),1);
% CdSeAAem450fit=fit(CdSeAAem450cutx,CdSeAAem450cut,'gauss2')
% %plot them
% plot(CdSebuffem450fit,CdSebuffem450cutx,CdSebuffem450cut)
% plot(CdSeAAem450fit,CdSeAAem450cutx,CdSeAAem450cut)

legend('CdSe/CdS in buffer, 100 mM AA, em 450,bg',...
    'CdSe/CdS in buffer, em 450,bg','100 mM AA buffer',...
    'CdSe/CdS in buffer minus AA, em 450,bg')%,'Fit-CdSe/CdS in buffer','Fit-CdSe/CdS minus AA')
xlabel('Wavelength (nm)')
ylabel('CPS')
%set axis labels to zoom in
axis([530 750 0 2e4])

hold off
%% Figure 5 UV Vis data
%Plot CdSe/CdS - native in toluene UV-vis data and CdSe/CdS-S in buffer
%UV-vis data
figure(5)
hold on

plot(CdSetolUV(:,1),CdSetolUV(:,2),'LineWidth',2)
plot(CdSebuffUV(:,1),CdSebuffUV(:,2),'LineWidth',2)
plot(CdSeAAUV(:,1),CdSeAAUV(:,2),'LineWidth',2)
legend('CdSe/CdS in toluene','CdSe/CdS-S in buffer',...
    'CdSe/CdS-S in buffer, 100 mM AA')
xlabel('Wavelength (nm)')
ylabel('Absorbance')
axis([250 650 -0.02 0.6])
hold off

%% Figure 6
figure(6)
hold on

%Using excitation condition 1 (450 nm) as it's the cleanest data
plot(CdSebuffem450(:,1),CdSebuffem450bg,'.')
plot(CdSebuffem450(:,1),CdSeAAem450bgsubAA,'.')

%Fit the 2 PL peaks (CdSeAAem450bgsubAA, CdSebuffem450bg)
CdSebuffem450cut=CdSebuffem450bg(27:127);
CdSebuffem450cutx=CdSebuffem450((27:127),1);
CdSebuffem450fit=fit(CdSebuffem450cutx,CdSebuffem450cut,'gauss2')
CdSeAAem450cut=CdSeAAem450bgsubAA(27:127);
CdSeAAem450cutx=CdSeAAem450((27:127),1);
CdSeAAem450fit=fit(CdSeAAem450cutx,CdSeAAem450cut,'gauss2')
%plot them
p1=plot(CdSebuffem450fit);%CdSebuffem450cutx,CdSebuffem450cut);
p2=plot(CdSeAAem450fit);%CdSeAAem450cutx,CdSeAAem450cut);
p1.LineWidth=1;
p2.LineWidth=1;
p1.Color=[0 0 0];
p2.Color=[0 0 0];
legend('CdSe/CdS in buffer, em 450,bg','CdSe/CdS in buffer minus AA, em 450,bg','Fit-CdSe/CdS in buffer','Fit-CdSe/CdS minus AA')
xlabel('Wavelength (nm)')
ylabel('CPS')
%set axis labels to zoom in
axis([530 750 0 1e4])
print(gcf, '-dtiff', '-r500', 'Fig6.tiff')
hold off
