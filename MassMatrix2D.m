function [M]=MassMatrix2D(p,t)
np=size(p,1);
% 1. triangles area
x21=p(t(:,2),1)-p(t(:,1),1); y12=p(t(:,1),2)-p(t(:,2),2);
x32=p(t(:,3),1)-p(t(:,2),1); y23=p(t(:,2),2)-p(t(:,3),2);
x13=p(t(:,1),1)-p(t(:,3),1); y31=p(t(:,3),2)-p(t(:,1),2);
tarea=(x21.*y31-x13.*y12)/2;
Ah=sparse(np,np);
for i=1:3
    for j=1:i
        Ah=Ah+sparse(t(:,i),t(:,j),tarea/12,np,np);
    end
end
Ah=Ah+Ah.';
M=Ah;
end