%% TCSPC data import and plot

close all
%Change directory
cd '/Users/Kristina/Google Drive/Dukovic group/Data/TCSPC/171026'

%Import data
%filenames = ['s11b.csv','1_2.csv','1_3.csv'];
numfiles=4;
S11 = csvread('s11b.csv',1,0);
S11b=csvread('1_1.csv',1,0);
S12 = csvread('1_2.csv',1,0);
S13 = csvread('1_3.csv',1,0);
%S21 = csvread(filenames(4,:),1,0);

%adjust time scales if necessary
S11(:,1)=(S11(:,1)/10^3); %converts ns to us

%Plot all data
figure(1)
hold on
pl=plot(S11(:,1),S11(:,2),S11b(:,1),S11b(:,2),S12(:,1),S12(:,2),S13(:,1),S13(:,2));

for k=1:numfiles
    pl(k).Marker='*';
    pl(k).LineStyle='none';
    pl(k).Marker='*';
    pl(k).MarkerSize = 1;
    k=k+1;
end 
set(gca,'YScale','log')

xlabel('Time (us)');
ylabel('Counts');
legend('S1-1','S1-1b','S1-2','S1-3');

%specify export size so it looks ok in Word etc
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3];
print(gcf, '-dtiff', '-r500', 'Fig1.tiff')

hold off
% %%
% figure(2)
% hold on
% pl=plot(S21(:,1),S21(:,2));
% 
% for k=1:1
%     pl(k).Marker='*';
%     pl(k).LineStyle='none';
%     pl(k).Marker='*';
%     pl(k).MarkerSize = 1;
%     k=k+1;
% end 
% set(gca,'YScale','log')
% axis([0 0.8 1 10^5])
% xlabel('Time (us)');
% ylabel('Counts');
% legend('S1-1');
% 
% %specify export size so it looks ok in Word etc
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 5 3];
% print(gcf, '-dtiff', '-r500', 'Fig2.tiff')
% 
% hold off