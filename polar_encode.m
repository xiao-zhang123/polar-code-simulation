function S1=polar_encode(N,X)
%N=input('enter the code length:');
F=[1 0;1 1];
%N=8;
n=log2(N);
%G=kron(F,n);��n��kron�ڻ�
G=F;
for k = 1:1:n-1
    G = kron(G,F);
end
% a=[0:N-1];
% b=bitrevorder(a); %����λ��ת
% BN=zeros(N);
% for i=1:N
%     j=b(1,i)+1;
%     BN(i,j)=1;
% end
% GN=BN*G;
S=X*G;
%��Sת��Ϊ0��1����
% for i=1:N
%     if (mod(S(1,i),2)==1)
%         S(1,i)=1;
%     else
%         S(1,i)=0;
%     end
% end
S1=mod(S,2);