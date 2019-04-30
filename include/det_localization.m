function [ alpha ] = det_localization( PF1,w1,PF2,w2,Xt )
Np = length(w1);
X_hat_f1  = mean(PF1);
X_hat_f2  = mean(PF2);
TF12 = Xt;
particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];



i1 = find(w1>0.00001);
mu1 = mean(PF1(i1,1:2));
c1 = cov(PF1(i1,1:2));

i2 = find(w2>0.00001);
mu2 = mean(particles2rel1(i2,1:2));
c2 =  cov(particles2rel1(i2,1:2));


alpha = det(c1)/(det(c1)+det(c2));

end

