function [crc_info_bit]= crc_bit_add(info_bit)
%���ã����������м��뵽��Ϣ������
%���ؼ����ı�������
%output:CRC�����Ϣ����
    %gCrc24a = [1 0 0 1 1];
    gCrc24a = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    crcLen = length(gCrc24a);
    crc_bit = zeros(1,crcLen-1);
    %state = 0;
    N=length(info_bit);
    %bitMask = (long long int) 1<<crcLen;
    bitMask = bitshift(1,crcLen-1);
    %����crclen-1λ
    % Loop and calculate CRC over each column vector
    state = 0;
    gCrc24a_bit = 0;
    b=zeros(1,crcLen-1);
    y=[info_bit b];
    %����state����Ϣ���ص�ǰcrclenλ
    for n=1:crcLen
      state = bitor(state,bitshift(y(n),crcLen-n));
      gCrc24a_bit = bitor(gCrc24a_bit,bitshift(gCrc24a(n),crcLen-n));
    end

    for n=1:N
       %state = bitshift(state,1);

       if bitand(state,bitMask) >0 %�ж�state�ĵ�25λ���������������ǲ���λ1
          state = bitxor(state,gCrc24a_bit); 
       end
       if(n==N)
           break;
       end
       state =bitshift(state,1);
       f=bitor(state,y(n+crcLen));
       state = bitor(state,y(n+crcLen));
    end
    %��stateת��Ϊ24λ�Ķ�������
    for n=1:crcLen-1
        crc_bit(crcLen-n) = bitget(state,n);
    end
    % Move CRC-bit into output vector

    crc_info_bit = [info_bit crc_bit];
end