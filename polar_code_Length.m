clear
clc
L=1;
code_length=[128 256 512 1024];
R=1/2;
qAry=1;
EBN0=0:0.5:4;
countCycles=zeros(length(code_length),length(EBN0)) ;
FER=zeros(length(code_length),length(EBN0)) ;
BER=zeros(length(code_length),length(EBN0)) ;
bpskMod = comm.BPSKModulator;
Frame=10000;
for Len=1:length(code_length)
    N=code_length(Len);
    K=R*N;
    IW=get_BEC_IWi(N);
    [info_bit_idx, frozen_bits] = get_info_and_frozen_location(IW, N, K);
    G = 1;
    F = [1 0;1 1];
    for j = 1:log2(N)
        G = kron(G,F);
    end
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
            info_bit = randi([0,1],1,K);     %获取随机的比特信息
            %将信息序列存储到数组before_code_bit中
            before_code_bit=zeros(1,N);
            before_code_bit(info_bit_idx(:))=info_bit(:);
            S=polar_encode(N,before_code_bit);                   %极化码编码
            frozen_bits_flag=frozen_bits';         
            %开始调制
            mod = bpskMod(S');
            rSig = chan(mod);
            rxLLR = bpskDemod(rSig); 
            llr=rxLLR;
            [dec_list] = SCL_CRCdecoder_llr(L,N,llr,noiseVar,info_bit_idx,G);
            after_decode_bit=dec_list(:,1)';
            %解码之后的序列
            % after_decode_bit 
            if length(find(info_bit-after_decode_bit))==0
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
ylabel('FER') %y轴坐标描述
legend('N=128','N=256','N=512','N=1024');   %右上角标注
grid on
title('SC decode R=0.5')
set(gcf,'color','white');