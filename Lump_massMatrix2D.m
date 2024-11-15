function [W]=Lump_massMatrix2D(n)
[p,e,t]=getmesh(n);
p=p';
t=t';
np=size(p,1);
% 1. triangles area
x21=p(t(:,2),1)-p(t(:,1),1); y12=p(t(:,1),2)-p(t(:,2),2);
x32=p(t(:,3),1)-p(t(:,2),1); y23=p(t(:,2),2)-p(t(:,3),2);
x13=p(t(:,1),1)-p(t(:,3),1); y31=p(t(:,3),2)-p(t(:,1),2);
tarea=(x21.*y31-x13.*y12)/2;
W=sparse(np,np);
for i=1:3
        W=W+sparse(t(:,i),t(:,i),tarea/3,np,np);
end
end