function [RU] = SCL_CRCdecoder_llr(L,N,y,sigma2,A,G)
n = log2(N);
K = length(A);
a = 1;
list = zeros(N,a);
sheet = zeros(N,n+1,a);
PM = zeros(1,a);
sheet(:,n+1,1) = 2*y'/sigma2;
base = zeros(1,n+1);
for q = 1:n+1
    base(q) = 2^(q-1);
end
q = n;
p = [ones(1,n),N+1];

while q<=n
    while q>=1
        for pq = p(q):p(q)+base(q)-1
            if mod(pq-1,base(q))+1==mod(pq-1,base(q+1))+1	% f运算
                com = zeros(2,a);
                com(1,1:a) = abs(sheet(pq,q+1,1:a));
                com(2,1:a) = abs(sheet(pq+base(q),q+1,1:a));
                for l = 1:a
                    minimum = min(com(1,l),com(2,l));
                    sheet(pq,q,l) = sign(sheet(pq,q+1,l))*sign(sheet(pq+base(q),q+1,l))*minimum;
%                     sheet(pq,q,l) = log((exp(sheet(pq,q+1,l)+sheet(pq+base(q),q+1,l))+1)/(exp(sheet(pq,q+1,l))+exp(sheet(pq+base(q),q+1,l))));
                end
            else                                            % g运算
                marks = G(1:base(q),mod(pq-1,base(q))+1);
                Usum = zeros(1,a);
                for v = 1:base(q)
                    if marks(v)==1
                        Usum = mod(Usum+list(v+floor((pq-1)/base(q+1))*base(q+1),:),2);
                    end
                end
                for l = 1:a
                    sheet(pq,q,l) = sheet(pq-base(q),q+1,l)*((-1)^Usum(l))+sheet(pq,q+1,l);
                end
            end
            if q==1
                if sum(ismember(A,pq))==1
                    %----------扩展----------
                    PM0 = zeros(1,2*a);
                    PM0(1:a) = PM;
                    PM0(a+1:2*a) = PM;
                    sheet2(:,:,1:a) = sheet;
                    sheet2(:,:,a+1:2*a) = sheet;
                    ru(1:2*a) = sheet2(pq,1,1:2*a)<0;
                    list_cdd = zeros(N,2*a);
                    list_cdd(1:pq-1,1:a) = list(1:pq-1,:);
                    list_cdd(pq,1:a) = 0;
                    list_cdd(1:pq-1,a+1:2*a) = list(1:pq-1,:);
                    list_cdd(pq,a+1:2*a) = 1;
                    sym = zeros(1,2*a);
                    for l = 1:2*a
                        sym(l) = 1-(ru(l)==list_cdd(pq,l));
                    end
                    PM = PM0;
                    for l = 1:2*a
                        PM(l) = PM0(l)+sym(l)*abs(sheet2(pq,1,l));
                    end
                    %----------竞争----------
                    if 2*a<=L
                        list = list_cdd;
                        sheet = sheet2;
                    else
                        [~,re] = sort(PM);
                        list = list_cdd(:,re(1:L));
                        sheet = sheet2(:,:,re(1:L));
                        PM = PM(re(1:L));
                    end
                    [~,a] = size(list);
                else
                    ru(1:a) = sheet(pq,1,1:a)<0;
                    sym = zeros(1,a);
                    for l = 1:a
                        sym(l) = 1-(ru(l)==0);
                    end
                    for l = 1:a
                        PM(l) = PM(l)+sym(l)*abs(sheet(pq,1,l));
                    end
                end
            end
        end
        p(q) = p(q)+base(q);
        q = q-1;
    end
    while 1
        q = q+1;
        if q>n
            break
        end
        if p(q)~=p(q+1)
            break
        end
    end
end
RU = list(A,:);
% for i = 1:L
%     c = RU(1:K,i);
%     [outdata,res] = detect(det,c);
%     if res==0
%         u = outdata';
%         llr = sheet(:,1,i);
%         break;
%     end
% end
% if res==1
%     c = RU(1:K,1);
%     [outdata,~] = detect(det,c);
%     u = outdata';
%     llr = sheet(:,1,1);
% end