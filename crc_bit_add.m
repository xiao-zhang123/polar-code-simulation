function [crc_info_bit]= crc_bit_add(info_bit)
%作用：将冗余序列加入到信息比特中
%返回加入后的比特序列
%output:CRC后的信息比特
    %gCrc24a = [1 0 0 1 1];
    gCrc24a = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    crcLen = length(gCrc24a);
    crc_bit = zeros(1,crcLen-1);
    %state = 0;
    N=length(info_bit);
    %bitMask = (long long int) 1<<crcLen;
    bitMask = bitshift(1,crcLen-1);
    %左移crclen-1位
    % Loop and calculate CRC over each column vector
    state = 0;
    gCrc24a_bit = 0;
    b=zeros(1,crcLen-1);
    y=[info_bit b];
    %最后的state是信息比特的前crclen位
    for n=1:crcLen
      state = bitor(state,bitshift(y(n),crcLen-n));
      gCrc24a_bit = bitor(gCrc24a_bit,bitshift(gCrc24a(n),crcLen-n));
    end

    for n=1:N
       %state = bitshift(state,1);

       if bitand(state,bitMask) >0 %判断state的第25位（从右往左数）是不是位1
          state = bitxor(state,gCrc24a_bit); 
       end
       if(n==N)
           break;
       end
       state =bitshift(state,1);
       f=bitor(state,y(n+crcLen));
       state = bitor(state,y(n+crcLen));
    end
    %将state转化为24位的二进制数
    for n=1:crcLen-1
        crc_bit(crcLen-n) = bitget(state,n);
    end
    % Move CRC-bit into output vector

    crc_info_bit = [info_bit crc_bit];
end