clear;
clc;
index=10;
%index=4;
n=2.^(1:index);
w=zeros(n(4));
w(1,1)=0.5;
for i=n
    for j=1:i/2
        w(i,2*j-1)=w(i/2,j)^2;
        w(i,2*j)=2*w(i/2,j)-w(i/2,j)^2;
    end
end
scatter(1:1024,w(1024,1:1024),'.');
axis([0 1024 0 1]);
%scatter(1:16,w(16,1:16),'.');
%axis([0 16 0 1]);
xlabel('channel index');
ylabel('symmetric capacity');
title('Channel polarization N=1024');
set(gcf,'color','white');
grid on