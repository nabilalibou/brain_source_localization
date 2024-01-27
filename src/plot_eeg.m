function []=plot_eeg(X,n,Fs,names)
%plots EEG data
%
%plot_eeg(X,n,Fs,names,idx)
%
%plots the EEG data given by the matrix X and normalized by the value n.
%
%INPUT: X - EEG data matrix
%       n - normalization factor
%       Fs - sampling rate
%       optional:
%       names - cell array specifying the channel names


%Determine number of electrodes and number of time samples
[N,T]=size(X);

if nargin<4 
    for n1=1:N
        names{n1}=num2str(n1);
    end
end

%construct time vector
t=0:1/Fs:(T-1)/Fs;

%plot EEG data
figure; 
g=plot(t,(X/n+flipud([0:N-1]')*ones(1,T))');
set(g,'Color','b');
axis([0, max(t), -1, N]);
set(gca,'YTick',[0:N-1],'YTickLabel',fliplr(names),'FontSize',14);
xlabel('Time in s','FontSize',16);
ylabel('Channels','FontSize',16);

