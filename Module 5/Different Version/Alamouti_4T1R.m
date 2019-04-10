function ber_ml=Alamoudi_2T1R(SNR,threshold)


Nt=4;                                                          %Number of TX Antennas
Nr=1;                                                          %Number of RX Antennas

no_bit_sym=1;                                                  %Number of bit per symbol
no_it_x_SNR=10000;                                             %Number of iteration per simulation

iter=0;                                                        %Setting up the variables


tot_err_ml = 0;                                                %Number of total errors

while tot_err_ml<threshold                                        %Starting the loop
        
    iter=iter+1;                                                 %Counting the iterations
    
        for i=1:no_it_x_SNR                                                 %Starting the simulation
           
            Data=(2*round(rand(Nt,1))-1)/(sqrt(Nt));                        %Creating random data
            
            %Building the Rayleigh Channel
            
            H = zeros(Nt,Nr);
            r = eye(Nt*Nr);                                               %Correlation matrix. 
            n = randn(Nt*Nr,1)/sqrt(2)+1i*randn(Nt*Nr,1)/sqrt(2);          %Channel coefficients
            H = reshape(r'*n,Nt,Nr);                                      %The complex conjugate matrix

             
            variance = sqrt(0.5/(10^(SNR/10)));                                  %Noise variance 
            
            n = variance * (randn(Nr,Nt) + 1i*randn(Nr,Nt));                    %Noise
         
            %Encoder.We code the data in an Alamouti Matrix

            X = [Data(1) -(Data(2)) -(Data(3)) -(Data(4));
                Data(2) Data(1) Data(4) -(Data(3));
                Data(3) -(Data(4)) Data(1) Data(2);
                Data(4) Data(3) -(Data(2)) Data(1)];              %Coded data
            
            R = H.'*X + n;                                                   %Received matrix
    
            %Combiner
%             s0=(conj(H(1,1))*R(1,1)) + (H(2,1)*conj(R(1,2))) - (H(3,1)*conj(R(1,3))) - (conj(H(4,1))*R(1,4));       
%             s1=(conj(H(2,1))*R(1,1)) - (H(1,1)*conj(R(1,2))) + (conj(H(4,1))*R(1,3)) - (H(3,1)*conj(R(1,4)));
%             s2=(conj(H(3,1))*R(1,1)) + (conj(H(4,1))*R(1,2)) + (H(1,1)*conj(R(1,3))) + (H(2,1)*conj(R(1,4)));
%             
%             S=[s0 s1 s2];

            s0= ;
            s1= conj(R(1,2));
            s2= conj(R(1,3));
            s3= conj(R(1,4));

            S=[s0 s1 s2 s3];
            
       
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %          Decoding              %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                    
            dh = sqrt(2)*[1 -1]/2;
          
            %Computing the distances for the first symbol%
            
            d11=((dh(1)-real(S(1)))^2+(imag(S(1)))^2);
            d12=((dh(2)-real(S(1)))^2+(imag(S(1)))^2);
            d13=((dh(3)-real(S(1)))^2+(imag(S(1)))^2);
              
            
            D1=[d11 d12 d13];                                            %Distances for the received symbol
              
              
            %Building the decisions vector for the first symbol%
                            
              for k=1:3
                  
                  X1_dec(k)=((abs(dh(k)))^2)*sum(sum((abs(H)).^2)-1)+D1(k);    
                  
              end
              
              %Computing the distances for the second symbol%
              
              d21=((dh(1)-real(S(2)))^2+(imag(S(2)))^2);
              d22=((dh(2)-real(S(2)))^2+(imag(S(2)))^2);
              d23=((dh(3)-real(S(2)))^2+(imag(S(2)))^2);
              
              D2=[d21 d22 d23];
              
              %Building the decisions vector for the second symbol%
              
              for x=1:3
                  
                  X2_dec(x)=((abs(dh(k)))^2)*sum(sum((abs(H)).^2)-1)+D2(x);
                  
              end
                  
              %Computing the distances for the third symbol%
              
              d31=((dh(1)-real(S(3)))^2+(imag(S(3)))^2);
              d32=((dh(2)-real(S(3)))^2+(imag(S(3)))^2);
              d33=((dh(3)-real(S(3)))^2+(imag(S(3)))^2);
              
              D3=[d31 d32 d33];
              
              %Building the decisions vector for the third symbol%
              
              for x=1:3
                  
                  X3_dec(x)=((abs(dh(k)))^2)*sum(sum((abs(H)).^2)-1)+D3(x);
                  
              end
              
              %The decisions!! We chose the little one%
              
              [min1, position1]=min(X1_dec);
              [min2, position2]=min(X2_dec);
              [min3, position3]=min(X3_dec);
              
              decoded=[dh(position1) dh(position2) dh(position3)];
              pattern=[min1 min2 min3];
            err_ml = sum(round(Data')~=round(decoded));                 %Computing the errors
            
            tot_err_ml = err_ml + tot_err_ml;                           
            
        end
        
end
        
            ber_ml=tot_err_ml/(no_it_x_SNR*iter*2);                      %Computing the BER
            
        end

