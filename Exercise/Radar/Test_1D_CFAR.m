close all;
s=randn(1000,1);
s([100, 300, 700])=[30 40 50];
refLength=12;
guardLength=3;
offset=3;
cfarWin=ones((refLength+guardLength)*2+1,1);
cfarWin(refLength+1:refLength+1+2*guardLength)=0;
cfarWin=cfarWin/sum(cfarWin);
noiseLevel=conv(s,cfarWin,'same');
    cfarThreshold=noiseLevel+offset;
    detection = (s > cfarThreshold);
    detection = detection.*s;
figure,plot(s);
hold on;
plot(cfarThreshold,'r--','LineWidth',2);
plot(detection,'g--','LineWidth',2);
legend('Signal','CFAR Threshold','Detection')