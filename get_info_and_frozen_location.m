function  [info_bits, frozen_bits] = get_info_and_frozen_location(channels, N, K)
%channes:�����ŵ���������
%N:codeword length
%K:info length
[ascend_sort_IWi, ~] = sort(channels, 'ascend');
critical_value = ascend_sort_IWi(N - K);%�ж���Ϣ���غͶ�����ص��ٽ�ֵ
info_bits = zeros(K, 1);
frozen_bits = ones(N, 1);
cnt=1;
for i = 1:N
    if channels(i) >= critical_value && cnt<=K
        info_bits(cnt) = i;
        frozen_bits(i) = 0;
        cnt=cnt+1;
    end
end
end