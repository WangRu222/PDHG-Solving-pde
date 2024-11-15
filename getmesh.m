function [p,e,t]=getmesh(n)
m=n;
left=0;
right=1;
h=(right-left)/n;
top=1;
down=0;
A=zeros(m+1,n+1);
p=zeros(3,(m+1)*(n+1));
t=zeros(3,2*m*n);
e=zeros(2,2*(n+m));
N=1;
for k=(m+n+2):-1:2
    if mod(k,2)==0
        I=max(k-(n+1),1):min(k-1,m+1);
    else
        I=min(k-1,m+1):-1:max(k-(n+1),1);
    end
    for i=I
        j=k-i;
        p(1,N)=left+(j-1)*h;
        p(2,N)=top-(i-1)*h;
        if (i-1)*(i-(m+1))*(j-1)*(j-(n+1))
            p(3,N)=0;
        else
            p(3,N)=1;
        end
        A(i,j)=N;
        N=N+1;
    end
end
N=1;
for k=(m+n+1):-1:2
    if mod(k,2)==1
        I=max(k-(n+1),1):min(k-1,m+1);
        for i=I
            j=k-i;
            if (i-1)*(j-(n+1))~=0
                t(1,N)=A(i,j);
                t(2,N)=A(i,j+1);
                t(3,N)=A(i-1,j+1);
                N=N+1;
            end
            if (i-(m+1))*(j-(n+1))~=0
                t(1,N)=A(i,j);
                t(2,N)=A(i+1,j);
                t(3,N)=A(i,j+1);
                N=N+1;
            end
        end
    else
        I=min(k-1,m+1):-1:max(k-(n+1),1);
        for i=I
            j=k-i;             
            if (i-(m+1))*(j-(n+1))~=0
                t(1,N)=A(i,j);
                t(2,N)=A(i+1,j);
                t(3,N)=A(i,j+1);
                N=N+1;
            end
            if (i-1)*(j-(n+1))~=0
                t(1,N)=A(i,j);
                t(2,N)=A(i,j+1);
                t(3,N)=A(i-1,j+1);
                N=N+1;
            end
        end
    end
        
end
N=1;
for k=1:n
    e(1,N)=A(1,k);
    e(2,N)=A(1,k+1);
    N=N+1;
end
for k=1:m
    e(1,N)=A(k,n+1);
    e(2,N)=A(k+1,n+1);
    N=N+1;
end
for k=1:n
    e(1,N)=A(m+1,k);
    e(2,N)=A(m+1,k+1);
    N=N+1;
end
for k=1:m
    e(1,N)=A(k,1);
    e(2,N)=A(k+1,1);
    N=N+1;
end
%[s t]=sort_edge(e,t);
% % save 
% % hold on
% % for i=1:size(p,2)
% %     temp=num2str(i);
% %     text(p(1,i)+0.1*h,p(2,i)-0.1*h,temp,'Color','blue');
% %     plot(p(1,i),p(2,i),'r*');
% % end
% % for i=1:size(s,2)
% %     temp=num2str(i);
% %     text(sum(p(1,s(1:2,i)))/2,sum(p(2,s(1:2,i)))/2,temp,'Color','green');
% % end
% % for i=1:size(t,2)
% %     temp=num2str(i);
% %     text(sum(p(1,t(1:3,i)))/3,sum(p(2,t(1:3,i)))/3,temp,'Color','red');
% %     plot(p(1,t([1 2],i)),p(2,t([1 2],i)));
% %     plot(p(1,t([2 3],i)),p(2,t([2 3],i)));
% %     plot(p(1,t([1 3],i)),p(2,t([1 3],i)));
% % end
% % hold off

% % p1=p(1:2,:);
% % s1=s;
% % t1=t;

p=[size(p,2) 0 0 0;(1:size(p,2))' p'];
t=[size(t,2) 0 0 0;(1:size(t,2))' t(1:3,:)'];

p=p';
p=p(2:3,2:end);
t=t';
t=t(2:4,2:end);
% pdemesh(p,e,t)
end