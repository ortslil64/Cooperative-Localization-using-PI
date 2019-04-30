function [ w12 ] = ParticlesIntersectionK_localization( PF1,w1,PF2,w2,alpha,beta,Xt )
Np = length(w1);
X_hat_f1  = mean(PF1);
X_hat_f2  = mean(PF2);
TF12 = Xt;
particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];


i1 = find(w1>0.00001);
w1 = w1(i1);
mu1 = PF1(i1,1:2);
for kk = 1:size(mu1,1)
    sigma1(:,:,kk)=beta*eye(2);
end
GM1 = gmdistribution(mu1,sigma1,w1);

% PF = [PF1;particles2rel1];
% PF1 = PF(randsample(2*Np,Np),:);


i2 = find(w2>0.00001);
w2 = w2(i2);
mu2 = particles2rel1(i2,1:2);
for kk = 1:size(mu2,1)
    sigma2(:,:,kk)=beta*eye(2);
end
GM2 = gmdistribution(mu2,sigma2,w2);
w12 = (pdf(GM1,PF1(:,1:2)).^(alpha)).*(pdf(GM2,PF1(:,1:2)).^(1-alpha));
end

