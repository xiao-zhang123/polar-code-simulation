clear
clc
List=[1 2 4 8];
N=256;
K=128;
R=K/N;
crclen=24;
qAry=1;
G = 1;
F = [1 0;1 1];
for j = 1:log2(N)
    G = kron(G,F);
end
EBN0=0:0.5:4;
countCycles=zeros(length(List),length(EBN0)) ;
FER=zeros(length(List),length(EBN0)) ;
BER=zeros(length(List),length(EBN0)) ;
bpskMod = comm.BPSKModulator;
IW=get_BEC_IWi(N);
[info_bit_idx, frozen_bits] = get_info_and_frozen_location(IW, N, K);
Frame=100000;
for Len=1:length(List)
    L=List(Len);
    for time=1:length(EBN0)
        EbNo=EBN0(time);
        error_Frames=0;
        error_Bits=0;
        correct=0;
        EsNo = EbNo + 10*log10(qAry);       
        snrdB = EsNo + 10*log10(R);       % in dB
        noiseVar = 1./(10.^(snrdB/10)); 
        chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
        bpskDemod = comm.BPSKDemodulator('DecisionMethod', ...
    'Approximate log-likelihood ratio','Variance',noiseVar);
        for num=1:Frame
            info_bit = randi([0,1],1,K-crclen);     %获取随机的比特信息
            before_code_bit=zeros(1,N);             %储存未编码之前的总序列           
            [crc_info_bit]= crc_bit_add(info_bit);  %CRC之后的信息序列
            %将信息序列存储到数组before_code_bit中
            before_code_bit(info_bit_idx(:))=crc_info_bit(:);
            S=polar_encode(N,before_code_bit);                   %极化码编码
            %开始调制
            mod = bpskMod(S');
            rSig = chan(mod);
            rxLLR = bpskDemod(rSig); 
            [dec_list] = SCL_CRCdecoder_llr(L,N,rxLLR,noiseVar,info_bit_idx,G);
            %crc校验
            for i=1:L
                after_decode_bit=dec_list(:,i)';
                final_crc_result=crc_check_result(after_decode_bit,length(info_bit));
                if (final_crc_result)
                    break;
                end
            end
        % crc_info_bit
        % %解码之后含CRC冗余的序列
        % after_decode_bit 
            if length(find(crc_info_bit-after_decode_bit))==0 && final_crc_result
                correct=correct+1;
            else
                error_Frames=error_Frames+1;
            end
            countCycles(Len,time)=countCycles(Len,time)+1;
            if error_Frames >=35 && countCycles(Len,time) >= 100
                break;
            end
        end
        FER(Len,time)=error_Frames/countCycles(Len,time);
    end
end
semilogy(EBN0,FER(1,:),'-*b',EBN0,FER(2,:),'-^r',EBN0,FER(3,:),'-sg',EBN0,FER(4,:),'-om');
xlabel('EBN0(dB)')  %x轴坐标描述
ylabel('FER')   %y轴坐标描述
legend('CASCL List=1','CASCL List=2','CASCL List=4','CASCL List=8');   %右上角标注
grid on
title('CA-SCL decode N=256 R=0.5')
set(gcf,'color','white');