% load('data_spiked');
data = importdata('data_spiked.csv');
t = data(:,1); yn = data(:,2);
%
tic;
yfilt = spikeremover(yn,0.03,10);
toc;
plot(t,yn); hold on; plot(t,yfilt); legend('orig','filt'); grid on;
title('Filtering Result')