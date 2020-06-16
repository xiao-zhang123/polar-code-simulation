function check_result=crc_check_result(output_bit,info_length)
%input:加入冗余校验序列后的序列
%所使用的循环校验多项式
%output:检验结果，1表示正确，0表示错误
    %gCrc24a = [1 0 0 1 1];
    gCrc24a = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    crcLen = length(gCrc24a);
    bitMask = bitshift(1,crcLen-1);
    info_check=0;
    gCrc24a_bit=0;
    for n=1:crcLen
      gCrc24a_bit = bitor(gCrc24a_bit,bitshift(gCrc24a(n),crcLen-n));
    end
    
    for n=1:crcLen
      info_check = bitor(info_check,bitshift(output_bit(n),crcLen-n));%获取前25位
    end
    crc_check=gCrc24a_bit;

    for n=1:info_length
       %state = bitshift(state,1);

       if bitand(info_check,bitMask) >0 %判断state的第25位（从右往左数）是不是位1
          info_check = bitxor(info_check,crc_check); 
       end
       if(n==info_length)
           break;
       end
       info_check =bitshift(info_check,1);
       %f=bitor(info_check,y(n+crcLen));
       info_check = bitor(info_check,output_bit(n+crcLen));
    end
    if(info_check==0)
        check_result=1;
    else
        check_result=0;
    end
end