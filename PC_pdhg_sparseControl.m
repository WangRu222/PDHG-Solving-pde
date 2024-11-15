function PC_pdhg_sparseControl(alpha,r,s,n,IterMax)
%prediction-correction scheme
%alpha=1e-5,r=4.2e3,s=1e-1 for sparse control
tic
beta=1e-3;
[point,edge,t]=getmesh(n);
n1=size(point,2);
n2=size(edge,2);
M=MassMatrix2D(point',t'); %%Mass matrix
A=StiffnessMatrix2D(point,t);
W=Lump_massMatrix2D(n);
bdNode = unique([edge(1,:) edge(2,:)]);
freenode=setdiff(1:n1,bdNode);
M=M(freenode,freenode);
A=A(freenode,freenode);
utt=1/sqrt(n1-n2)*ones(n1-n2,1);
norm(utt)
ytt=A\(M*utt);
ptt=A\(M*ytt);
norm(ptt)
W=W(freenode,freenode);
u=1e-4*ones(n1-n2,1);
p=1e-4*ones(n1-n2,1);
a=-30*ones(n1-n2,1);
b=30*ones(n1-n2,1);
yd1=sin(2*pi*point(1,:)).*sin(2*pi*point(2,:)).*exp(2*point(1,:))*(1/6);
%  yd1=sin(2*pi*point(1,:)).*sin(2*pi*point(2,:));
yd=yd1(freenode);
err=1;
tol=1e-5;
iter=0;
[LA,DA,PA]=ldl(A);
while(err>tol && iter <=IterMax)
   iter=iter+1;
   q=PA*(LA'\(DA\(LA\(PA'*(M*p)))));
    u_new=wthresh((u-r*q)/(r*alpha+1),'s',r*beta/(r*alpha+1));
   u_new=max(a,min(b,u_new));
   u_bar=u_new+ 1*(u_new-u);
%    (sqrt(5)-1)/2
   y=PA*(LA'\(DA\(LA\(PA'*(M*u_bar)))));
   p_new=(s*y+1*p-s*yd')/(s+1);
   u_new1=u+1.8*(u_new-u);
   p_new1=p+1.8*(p_new-p);
%    (2/(sqrt(5)-1))
   dual_residual(iter)=sqrt((p_new-p)'*M*(p_new-p))/sqrt(p'*M*p);
   primal_residual(iter)=sqrt((u_new-u)'*M*(u_new-u))/sqrt(u'*M*u);
   err=max(primal_residual(iter),dual_residual(iter));
    y_new=y;
    objective_value(iter)=0.5*(y_new-yd')'*M*(y_new-yd')+0.5*alpha*u_new'*M*u_new+beta*norm(W*u_new,1);
    u=u_new1;
    p=p_new1;
end
toc
distance_PC=sqrt((y_new-yd')'*M*(y_new-yd'))
iter
Obj_PC=objective_value(iter)
yy(freenode)=y_new;
yy(bdNode)=0;
uu(freenode)=u_new;
uu(bdNode)=0;
figure
loglog(1:iter,dual_residual,'-*')
xlim([0 iter])
hold on
loglog(1:iter,primal_residual,'r-square')
hold on
loglog(1:iter,objective_value,'g-diamond')
xlim([0 iter])
% xlabel('iteration')
% ylabel('Objective value')
hold off
legend('Dual residual', 'Primal residual','Objective value')
% t=[t;ones(1,size(t,2))];
figure ,colormap jet;trimesh(t',point(1,:),point(2,:),yy-yd1,'LineWidth',1);
figure ,colormap jet;trimesh(t',point(1,:),point(2,:),yy,'LineWidth',1);
figure ,colormap jet;trimesh(t',point(1,:),point(2,:),uu,'LineWidth',1);
% figure
% view(-38,40)
% pdemesh(point,edge,t,yy-yd1)
% figure
% view(-38,40)
% pdesurf(point,t,yy')
% figure
% view(-38,40)
% pdesurf(point,t,uu')
% figure
% view(-38,40)
% pdesurf(point,t,yy(:,k)-True_Solution(:,k))
% figure
% view(-38,40)
% pdesurf(point,t,yd1(:,N+1))
end
