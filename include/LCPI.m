function [ PF_out ] = LCPI( PF_in,W,A,Xt )
Np = size(PF_in,1);
Nf = size(PF_in,3);

for ii = 1:Nf
    PF_out_t = [];
    B=A(ii,:);
    X_hat = mean(PF_in(:,:,ii));
    for jj = 1:Nf
        TF = Xt{ii,jj};
        PF_out_t =  [PF_out_t;PF_in(randsample(Np,ceil(Np.*B(jj)),true,W(:,jj)),:,jj)+ TF];
    end
    PF_out(:,:,ii) = PF_out_t(1:Np,:);
end
end

