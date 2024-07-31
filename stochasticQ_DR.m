function [quantized,number_of_bits_toSend, R, b]=stochasticQ_DR(quantized,current,prev,bitsToSend, prevRange, prev_b, flag, number_of_bits_toSend)

%flag=0;
%b=bitsToSend;
diff=current - prev;
R=max(abs(diff));
if flag == 0
   b=bitsToSend;
else
    b = log2(1+(2^prev_b-1)*R/prevRange);
    b=ceil(b);
    
    if (b > 32)
        b=32;
    end
end
tau=1/(2^b-1);
number_of_bits_toSend = number_of_bits_toSend+32+5+length(current)*b;% the number of bits to send the value of R.

% Stochastic Quantization
Q=(diff+R)/(2*tau*R);
p=(Q-floor(Q));
for i=1:length(current)
    temp=rand;
    if(temp <=p(i))
        Q(i)=ceil(Q(i));
    else
        Q(i)=floor(Q(i));
    end    
end
quantized=quantized+2*tau*Q*R-R;
%kkk=1;
%b

end    

