%% TCSPC data import and plot

close all

%Change directory
directory=input('What is the folder? ','s');
path = '/Users/Kristina/Google Drive/Dukovic group/Data/TCSPC/';
fullpath=strcat(path,directory);
cd(fullpath)

%Import data
% Can also use matrix of filenames?
S01=csvread('0_1.csv',1,0);
S02=csvread('0_2.csv',1,0);
S03=csvread('0_3.csv',1,0);
S11=csvread('1_1.csv',1,0);
S12=csvread('1_2.csv',1,0);
S13=csvread('1_3.csv',1,0);

IRF=csvread('IRF.csv',1,0);

filenum=7;

%Graph all assuming common time range?
figure(1)
hold on

p1 = plot(S01(:,1),S01(:,2),S02(:,1),S02(:,2),S03(:,1),S03(:,2),...
    S11(:,1),S11(:,2),S12(:,1),S12(:,2),S13(:,1),S13(:,2),IRF(:,1),IRF(:,2));

for k=1:filenum
    p1(k).LineStyle='none';
    p1(k).Marker='*';
    p1(k).MarkerSize=3;
    k=k+1;
end

set(gca,'YScale','log')
xlabel('Time (ns)')
ylabel('Counts')
legend('S0-1','S0-2','S0-3','S1-1','S1-2','S1-3','IRF')

%set axis range
xmin=input('What is x min?');
xmax=input('What is x max?');
axis([xmin xmax 0 10^5])

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig1.tiff')

hold off