function check_result=crc_check_result(output_bit,info_length)
%input:��������У�����к������
%��ʹ�õ�ѭ��У�����ʽ
%output:��������1��ʾ��ȷ��0��ʾ����
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
      info_check = bitor(info_check,bitshift(output_bit(n),crcLen-n));%��ȡǰ25λ
    end
    crc_check=gCrc24a_bit;

    for n=1:info_length
       %state = bitshift(state,1);

       if bitand(info_check,bitMask) >0 %�ж�state�ĵ�25λ���������������ǲ���λ1
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