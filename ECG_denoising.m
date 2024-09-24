% This code is for paper"Electrocardiogram denoising using
% Wavelet decomposition and EMD domain filtering" published in the 
% proceedings of 2016 IEEE Region 10 Conference (TENCON). If you are using the
% code for any educational purpose then please cite the paper.
% Author: Ashish Verma
load('100m.mat');
x=val;
Lt=2000;
x1= x(1,1:Lt);
x2= x(2,1:Lt);
y1=((x1-1024)/200);
y2=((x2-1024)/200);
y=y1+y2;

imf1=emd(y);
d1=imf1(1,:)+imf1(2,:)+imf1(3,:);
% Adding White Gaussian Noise to the ECG signal
rng('default');
yn=awgn(y,5,'measured');

% Add baselinewander noise to ECG Signal
% t=0:2*pi/1:22619;
% n=1*sin(2*pi*.000125*t);
% yn2=n+yn1;

% Add Power Line Noise TO the ECG signal
% t=0:2*pi/1:1.2566e+04;
% n1=.5*sin(2*pi*50*t);
% yn=n1+y;

%qrs=dpi_qrs(yn,320,1800,5);%for qrs extraction


 
[threshold,SORH,keepapp]=ddencmp('den','wv',yn);
r=wdencmp('gbl',yn,'db4',2,threshold,SORH,keepapp);

% Using Pan-Tompkin's method for qrs detection
[qrs_amp,qrs,d]=pan_tompkin(r,360,0);

tw=zeros(1,Lt);
w=tukeywin(36,0.5);

for i=1:length(qrs)
  tw(qrs(i)-17:qrs(i)+18)=w;
end

imf=emd(r);
[m,n]=size(imf);
s=0;
d=imf(1,:)+imf(2,:)+imf(3,:);
ws = d.*tw;

for a=1:25

z=0;
IHP=0;

for j= 1:m-1      %  m = number of rows of imf matrix
    
    % Soft-thresholding Method using instantaneous half-period
    MAD(j)   = median(abs(imf(j,:)-median(imf(j,:)))) ;
    sigma(j) = MAD(j)/.6745;
    tau(j)   = sigma(j)*sqrt(2*log10(length(imf(j,:))));
    f=150;
    
    th=a/(2*f); % f= highest freq in signal
    k=1;
    
    for t=1:n-1

        
        if ((imf(j,t)>0 && imf(j,t+1)<=0 ) || (imf(j,t)<=0 && imf(j,t+1)>0))
            
           
                z(j,k) = .002777777778*t;
             
            k=k+1;
      
        end
        z(j,k)=0;
     end
    

        [q,p]=size(z);
        
        for i=1:p
            if z(j,i)==0
                z(j,i)= .002777777778*Lt;
            end
        end
     
        for k=1:p-1
             
            if z(j,k+1)~= .002777777778*Lt;
            IHP(j,k) = z(j,k+1) - z(j,k);
            
            else
                IHP(j,k)=10;
            end
        end
        
             
         for k=1:p-1
             
             if IHP(j,k) < th
                 
                 
                 for t= round(z(j,k)/.002777777778):round(z(j,k+1)/.002777777778)
                     
                         if   imf(j,t) >= tau(j)
                           c(j,t) =  imf(j,t)-tau(j);
            
                         elseif  abs(imf(j,t)) < tau(j)
                            c(j,t) =  0;
            
                         elseif   imf(j,t) <= tau(j)
                            c(j,t) = imf(j,t)+tau(j);
                         end
                 end
             else
                 for t= round(z(j,k)/.002777777778):round(z(j,k+1)/.002777777778)
                 
                         c(j,t)=  imf(j,t);
                         
                 end
             end
         end
        
    end

S=0;
for i=4:m-1               %upto m-1 to remove only WGN
    S=S+c(i,:);           %upto m-3 to remove both WGN and baseline wander noise
end

% X(a,:)= ws+(S+imf(m,:)).*tw.*d+S+imf(m,:);  %add imf(m,:) to RHS to remove WGN only
X(a,:)= ws+S+imf(m,:);
end

for a=1:24
e(a,:) = X(a,:)-X(a+1,:);
cmse(a) = sum(e(a,:).*e(a,:))/length(e(a,:));
end
cmsemin=min(cmse);
for a=1:24
    if cmse(a)==cmsemin
        oth=a/(2*f);
        in=a;
    end
end

% for l=1:length(qrs)
%     X(in,qrs(l)-10:qrs(l)+10)=r(1,qrs(l)-10:qrs(l)+10);
% end
X1(1,:)=X(in,:);
for t=1:Lt
    p1(1,t)=(abs(yn(1,t)-y(1,t)))^2;
    p2(1,t)=(abs(X(1,t)-y(1,t)))^2;
    p3(1,t)=(abs(y(1,t)))^2;
end
sp1 = sum(p1);
sp2 = sum(p2);
sp3 = sum(p3);
SNR = 10*log10(sp1/sp2);

MSE = sp2/Lt;

PRD = 100*sqrt(sp2/sp3);

figure(1)
subplot(4,1,1)
plot(y,'b')
xlabel('Samples')
ylabel('mV')
title('Original ECG signal')

subplot(4,1,2)
plot(yn,'r')
xlabel('samples')
ylabel('mV')
title('Signal with 5 dB WGN')% PLN with 50Hz PLN

subplot(4,1,3)
plot(r,'g')
xlabel('Samples')
ylabel('mV')
title('Wavelet Denoised signal')

subplot(4,1,4)
plot(X(in,:),'g')
xlabel('samples')
ylabel('mV')
title('Final Denoised signal')

%subplot(4,1,4);plot(r,'m');
figure(2)
plot(cmse);
xlabel('Threshold Number');
ylabel('CMSE');

figure(3)

for i=1:m
    subplot(m,2,i);
    plot(imf(i,:))
end

figure(4)

for i=4:m-1
   subplot(m,2,i);
   plot(c(i,:))
   
end

figure(5)
subplot(3,1,1)
plot(d,'r')
xlabel('Samples')
ylabel('mV')
title('d(t) from Wavelet Denoised ECG')

subplot(3,1,2)
plot(tw,'b')
xlabel('samples')
ylabel('mV')
title('Group of Tukey Windows')

subplot(3,1,3)
plot(ws,'g')
xlabel('samples')
ylabel('mV')
title('Windowed Signal')

figure(6)

subplot(4,1,1)
plot(ws,'r')
xlabel('Samples')
ylabel('mV')
title('Windowed Signal(through 1st 3 IMFs)')

subplot(4,1,2)
plot(S,'b')
xlabel('samples')
ylabel('mV')
title('Denoised Signal through Remaining IMFs(IMFs except 1st 3 IMFs) ')

subplot(4,1,3)
plot(imf(m,:),'g')
xlabel('Samples')
ylabel('mV')
title('Residue Signal')

subplot(4,1,4)
plot(X(in,:),'g')
xlabel('Samples')
ylabel('mV')
title('Final Denoised Signal')


