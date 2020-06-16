%比较不同码率对polar码的SC解码的影响
clear
clc
L=1;
N=512;
code_rate=[1/8 3/8 5/8 7/8];                             
qAry=1;                                               %qAry=1表示使用BPSK对应的调制阶数
EBN0=-4:0.5:4;                                        %信噪比的范围
countCycles=zeros(length(code_rate),length(EBN0)) ;   %不同码率和信噪比时所运行的帧数
FER=zeros(length(code_rate),length(EBN0)) ;           %误帧率
BER=zeros(length(code_rate),length(EBN0)) ;           %误码率
bpskMod = comm.BPSKModulator;                         %使用BPSK调制
Frame=100000;                                         %总帧数
for Len=1:length(code_rate)
    R=code_rate(Len);
    K=R*N;
    IW=get_BEC_IWi(N);                                %子信道的信道容量
    [info_bit_idx, frozen_bits] = get_info_and_frozen_location(IW, N, K);  %分别表示信息位的下标和冻结位的位置
    G = 1;
    F = [1 0;1 1];
    for j = 1:log2(N)
        G = kron(G,F);
    end                                               %G表示最终的生成矩阵
    for time=1:length(EBN0)
        EbNo=EBN0(time);
        error_Frames=0;
        %error_Bits=0;
        correct=0;
        EsNo = EbNo + 10*log10(qAry);       
        snrdB = EsNo + 10*log10(R);       % in dB
        noiseVar = 1./(10.^(snrdB/10)); 
        chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
        bpskDemod = comm.BPSKDemodulator('DecisionMethod', ...
    'Approximate log-likelihood ratio','Variance',noiseVar);
        for num=1:Frame
            info_bit = randi([0,1],1,K);                          %获取随机的比特信息
            before_code_bit=zeros(1,N);                           %将信息序列存储到数组before_code_bit中
            before_code_bit(info_bit_idx(:))=info_bit(:);
            S=polar_encode(N,before_code_bit);                    %极化码编码
            frozen_bits_flag=frozen_bits';         
            mod = bpskMod(S');                                    %开始调制
            rSig = chan(mod);                                     %过信道
            rxLLR = bpskDemod(rSig);                              %解调
            llr=rxLLR;
            [dec_list] = SCL_CRCdecoder_llr(L,N,llr,noiseVar,info_bit_idx,G);
            after_decode_bit=dec_list(:,1)';                      %选择路径度量值最小的作为译码结果
            if length(find(info_bit-after_decode_bit))==0         %计算误帧率
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
legend('rate=1/8','rate=3/8','rate=5/8','rate=7/8');   %右上角标注
grid on
title('SC decode N=512')
set(gcf,'color','white');