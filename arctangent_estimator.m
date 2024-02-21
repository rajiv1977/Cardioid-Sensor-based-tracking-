%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Arctangent Estimator                                      %
%                     Copyright @2015_DRDC, version 01_02112015                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv,  and B.Balaji                                      %
%          Defence R&D Canada, 3701 Carling Avenue, Ottawa, ON, K1A 0Z4, Canada.            %
%             rajiv.sithiravel@gmail.com and Bhashyam.Balaji@drdc-rddc.gc.ca                %
%                                                                                           %
%                                   T.Kirubarajan                                           %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                                 kiruba@mcmaster.ca                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ B_hat ] = arctangent_estimator(Meas,T)

Fs = 1/T;
% Arctangent estimator
    m = length(Meas);
    n = pow2(nextpow2(m));
    f = (0:n-1)*(Fs/n);
    
    
    Y1_fft = fft(Meas(1,:),n);
    power1 = Y1_fft.*conj(Y1_fft)/n;
    figure(1);
    plot(f,power1);
    xlabel('Frequency (Hz)');
    ylabel('Power in omni channel');
    
    
    Y2_fft = fft(Meas(2,:),n);
    power2 = Y2_fft.*conj(Y2_fft)/n;
    figure(2);
    plot(f,power2);
    xlabel('Frequency (Hz)');
    ylabel('Power in cosine channel');
    
    Y3_fft = fft(Meas(3,:),n);
    power3 = Y3_fft.*conj(Y3_fft)/n;
    figure(3);
    plot(f,power3);
    xlabel('Frequency (Hz)');
    ylabel('Power in sine channel');
    
    c_hat = 0;
    for mm=1:length(Y1_fft)
       c_hat = c_hat + (Y2_fft(mm) * conj(Y1_fft(mm))) ;
    end
    c_hat = real(c_hat);
    
    s_hat = 0;
    for mm=1:length(Y1_fft)
       s_hat = s_hat + (Y3_fft(mm) * conj(Y1_fft(mm))) ;
    end
    s_hat = real(s_hat);

    B_hat = atan(s_hat/c_hat);
     
end

