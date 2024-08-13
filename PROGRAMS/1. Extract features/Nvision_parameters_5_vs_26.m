function [x] = Nvision_parameters_5_vs_26 (bdata,rate)
%calculates the desired features and stacks them in a vector need=[ 1 3 4 21 22];
x=zeros(size(bdata,2),26);
%every row is for window

for y=1:size(bdata,2)  
    window_data=bdata(:,y); %only one window.
    % Energy:
    x(y,1) = mean(window_data  .^2);
    % Maximal ampllitude:
    x(y,2) = max(abs(window_data));
    % Amplitude integral:
    c=cumtrapz(window_data);
    x(y,3)=c(end);
    % STD:
    x(y,4)= std(window_data);
    % Skewness:
     x(y,5)=skewness(window_data);
    % Kutrosis:
    x(y,6) = kurtosis(window_data);
    % Curve length:
    x(y,7)= sum(abs(diff(window_data)));
    % Harmonicity:
    h=harmtest (window_data,rate);
    x(y,8) = h;
    % Zero crossing: 
    xaxes = [1:length(window_data)];
    x(y,9) = length(findX(xaxes',window_data,0));

    % fractle dimnsion :
    x(y,10)= calc_FD(window_data'); % calculate the fractle dimnsion of the window using box counting ;
    % Peak:
    x(y,11) = max(abs(window_data))/std(window_data);
    
    % ne features:
    window_data=window_data';
    ne= window_data.^2.- ([window_data(2:length(window_data)) 0]).*([0 window_data(1:(length(window_data)-1))]);
    x(y,12) = mean(ne);
    x(y,13) = std(ne);
    x(y,14) = kurtosis(ne);
    x(y,15) = skewness(ne);
    clear ne;
    
    % fft features:
    L=length(window_data);
    NFFT =2^nextpow2(L);
    f=rate/2*linspace(0,1,NFFT/2+1);
    fftdata = (abs(  fft(window_data,NFFT)/L  ) ).^2;
    allamp=fftdata;
    S_allamp=sum(fftdata);
    
    x(y,16)=std(fftdata);
    x(y,17)=kurtosis(fftdata);
    x(y,18)=skewness(fftdata);
    x(y,19)=median(fftdata); 
    CS=cumsum(fftdata); 
    x(y,20)=CS(end);% sum FFT
    
    % FFT BANDS
    % Delta: 
    st=find(f>0.5,1);
    ed=find(f>3,1);
    x(y,21) = sum(allamp(st:ed-1))/S_allamp;
    %  Theta: 
    st=find(f>3,1);
    ed=find(f>8,1);
    x(y,22) = sum(allamp(st:ed-1))/S_allamp;
    % Alpha: 
    st=find(f>8,1);
    ed=find(f>12,1);
    x(y,23) = sum(allamp(st:ed-1))/S_allamp;
    % Beta
    st=find(f>12,1);
    ed=find(f>20,1);
    x(y,24) = sum(allamp(st:ed-1))/S_allamp;
    % low gamma
    st=find(f>20,1);
    ed=find(f>30,1);
    x(y,25) = sum(allamp(st:ed-1))/S_allamp;
    % high gamma
    st=find(f>30,1);
    ed=find(f>100,1);
    x(y,26) = sum(allamp(st:ed-1))/S_allamp;
end



end 