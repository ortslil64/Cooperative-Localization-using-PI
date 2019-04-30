function A = GenerateFusionNetwork(Nf,a)
V=rand(Nf);
V=orth(V);
V(:,1)=ones(Nf,1);
V = GramSchmidt(V);
Ei=[1, a.*ones(1,Nf-1)];
Ei = sort(Ei,'descend');
Di = diag(Ei);
A = V*Di*V';
A(A<0)=0;
for ii = 1:Nf
    A(ii,:) =  A(ii,:)./sum(A(ii,:));
end
end

