% load('hyd_1000_FFTL.mat')
% hyd1000 = hyd;
% load('hyd_2000_FFTL.mat')
% load('sei_whale_overlap_test_parm.mat')

hyd1000.detection.call_frequency = hyd1000.detection.call_frequency(1:3530);

start = datenum(hyd.detection.start);
finish = datenum(hyd.detection.end);

one_sec = datenum('01-Jan-2021 12:00:01') - datenum('01-Jan-2021 12:00:00');

figure(10); hold on
subplot(2,1,1)
plot(hyd.detection.call_frequency)
xlabel('Date'); ylabel('Calls per Window')
title('Call Frequency: 2000 FFTL')
axis([0 3530 0 25])
box on; grid on

subplot(2,1,2)
plot(dateTimeVector_num,hyd1000.detection.call_frequency)
xlabel('Date'); ylabel('Calls per Window')
title('Call Frequency: 1000 FFTL')
axis([0 3530 0 25])
box on; grid on
