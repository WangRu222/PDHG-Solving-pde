function [A]=StiffnessMatrix2D(p,t)
format long;
N = size(p,2);
A = sparse(N,N); 
ve(:,:,1) = (p(:,t(3,:))-p(:,t(2,:)))';
ve(:,:,2) =  (p(:,t(1,:))-p(:,t(3,:)))';
ve(:,:,3) = (p(:,t(2,:))-p(:,t(1,:)))';
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

alpha=1;
for i = 1:3
    for j = 1:3
        Aij = alpha.*((ve(:,1,i).*ve(:,1,j)+ve(:,2,i).*ve(:,2,j)))./(4*area);
        A = A + sparse(t(i,:)',t(j,:)',Aij,N,N);
    end
end
end
