%------------------------------------------------------------------
% Tri_Quad.m
%
% The main program for calculating plane problems, including triangular element, 
% four node isoparametric element and hybrid element
% 
% This program is only used for learning and communication. 
% Please do not use it for business without permission. 
% Please indicate the source for reprint
%
% Written by Jia Xu & Ruoxi Liu - 04/13/2022
% Contact: liuroci@163.com
%------------------------------------------------------------------

%------------------------------- Mix -------------------------------
%-------------------------- Main Program ---------------------------
function [displacement,stress,stress1]=Tri_Quad(n_node,n_element,type,element_type,x,y,element_material,element_nodes,material,constraint,load)
% Integrate global stiff matrix
stiffGlobal=zeros(2*n_node,2*n_node);
for iE=1:n_element
    if element_type(iE)==3
       n_elementnode=3;
       stiffLocal=Tri_Local_stiff(iE,x,y,element_material,element_nodes,material,type);  % Stiff matrix of each element
       element_node=element_nodes(:,1:n_elementnode);
       location=zeros(1,6);
       location(1:2:2*n_elementnode-1)=2*element_node(iE,:)-1;
       location(2:2:2*n_elementnode)=2*element_node(iE,:);
       stiffGlobal(location,location)=stiffGlobal(location,location)+stiffLocal;
    else
       n_gauss=4;
       n_elementnode=4;
       stiffLocal=Quad_Local_stiff(iE,x,y,element_material,element_nodes,material,type,n_gauss);
       element_node=element_nodes(:,1:n_elementnode);
       location=zeros(1,8);
       location(1:2:2*n_elementnode-1)=2*element_node(iE,:)-1;
       location(2:2:2*n_elementnode)=2*element_node(iE,:);
       stiffGlobal(location,location)=stiffGlobal(location,location)+stiffLocal;
    end
end
% Integrate global load matrix
loadGlobal=zeros(2*n_node,1);
for i=1:size(load,1)
    loadGlobal(load(i,1))=load(i,2);
end
% Solve displacement
[stiffGlobal2,loadGlobal2]=Replace01(stiffGlobal,loadGlobal,constraint'); % Process boundary conditions
displacement=stiffGlobal2\loadGlobal2;  
%Solve Stress
stress=zeros(n_element,3);
stress1=zeros(n_element,3);
for iE=1:n_element
    if element_type(iE)==3
        [stress(iE,:),stress1(iE,:)]=Tri_stress_solve(iE,displacement,x,y,element_material,element_nodes,material,type); 
                                     
    else
        [stress(iE,:),stress1(iE,:)]=Quad_stress_solve(iE,displacement,x,y,element_nodes);
    end                            
end
end

%----------------------------Triangular-----------------------------
%-------------------- Local stiffness matrix  ----------------------
function stiffLocal=Tri_Local_stiff(iE,x,y,element_material,element_nodes,material,type)
    n1=element_nodes(iE,1);
    n2=element_nodes(iE,2);
    n3=element_nodes(iE,3);
    n=element_material(iE);
    if type>0
        material(:,1)=material(:,1)/(1-material(:,2)^2);
        material(:,2)=material(:,2)/(1-material(:,2));
    end
    E=material(n,1);um=material(n,2);t=material(n,3);
    c(3)=x(n2)-x(n1);c(2)=x(n1)-x(n3);c(1)=-c(2)-c(3);
    b(3)=y(n1)-y(n2);b(2)=y(n3)-y(n1);b(1)=-b(2)-b(3);
    A=(b(2)*c(3)-b(3)*c(2))/2;
    a1=E*t/4/(A-um^2*A);u1=(1-um)/2;
    st=zeros(6);
    for ii=1:3
        for jj=1:3
            bb=b(ii)*b(jj);cc=c(ii)*c(jj);
            cb=c(ii)*b(jj);bc=b(ii)*c(jj);
            i2=2*ii;j2=2*jj;
            i1=i2-1;j1=j2-1;
            st(i1,j1)=a1*(bb+u1*cc);
            st(i1,j2)=a1*(um*bc+u1*cb);
            st(i2,j1)=a1*(um*cb+u1*bc);
            st(i2,j2)=a1*(cc+u1*bb);
        end
    end
    stiffLocal=st;
end

%--------------------  Solve element stress  -----------------------
function [stress,stress1]=Tri_stress_solve(iE,displacement,x,y,element_material,element_nodes,material,type)
if type>0
    material(:,1)=material(:,1)/(1-material(:,2)^2);
    material(:,2)=material(:,2)/(1-material(:,2));
end
    n1=element_nodes(iE,1);
    n2=element_nodes(iE,2);
    n3=element_nodes(iE,3);
    n=element_material(iE);
    E=material(n,1);um=material(n,2);t=material(n,3);
    c(3)=x(n2)-x(n1);
    b(3)=y(n1)-y(n2);
    c(2)=x(n1)-x(n3);
    b(2)=y(n3)-y(n1);
    c(1)=-c(2)-c(3);
    b(1)=-b(2)-b(3);
    A=(b(2)*c(3)-b(3)*c(2))/2;
    aa=E/(2*A*(1-um^2));
    location(1:2:5)=2*element_nodes(iE,1:3)-1;
    location(2:2:6)=2*element_nodes(iE,1:3);
    w=displacement(location);
    s1=zeros(3,6);
    for ii=1:3
        j1=2*ii-1;j2=j1+1;
        s1(1,j1)=b(ii)*aa;
        s1(2,j1)=s1(1,j1)*um;
        s1(3,j2)=s1(1,j1)*(1-um)/2.0;
        s1(2,j2)=c(ii)*aa;
        s1(1,j2)=s1(2,j2)*um;
        s1(3,j1)=s1(2,j2)*(1-um)/2.0;
    end
    st=zeros(1,3);
    for ii=1:3
        for jj=1:6
        st(ii)=st(ii)+s1(ii,jj)*w(jj);
        end
    end
    stress(1,:)=st;
    stress1(1,1)=(st(1)+st(2))/2+sqrt(((st(1)-st(2))/2)^2+st(3)^2);
    stress1(1,2)=(st(1)+st(2))/2-sqrt(((st(1)-st(2))/2)^2+st(3)^2); 
    stress1(1,3)=90-1/2/pi*180*atan(2*st(3)/(st(1)-st(2)));         
end


%-------------------------Isoparametric_Fournode--------------------------
%--------------------- Gauss point and weight  ---------------------
function [xg,w]=gauss(n_gauss)
    x_gauss(2,1:2)=[-sqrt(3)/3 sqrt(3)/3];
    weight(2,1:2)=[1 1];
    x_gauss(3,1:3)=[-sqrt(15)/5 0 sqrt(15)/5];
    weight(3,1:3)=[5/9 8/9 5/9];
    x_gauss(4,1:4)=[-0.8611363115940520 -0.3399810435848560 0.3399810435848560 0.8611363115940520];
    weight(4,1:4)=[0.3478548451374530 0.6521451548625460 0.6521451548625460 0.3478548451374530];
    xg=x_gauss(n_gauss,1:n_gauss);
    w=weight(n_gauss,1:n_gauss);
end

%------------------ Basis function & Derivative  -------------------
% function N=Basis(ks,et)  % ks-¦Î,et-¦Ç
% N(1)=(1-ks)*(1-et)/4; 
% N(2)=(1+ks)*(1-et)/4; 
% N(3)=(1+ks)*(1+et)/4; 
% N(4)=(1-ks)*(1+et)/4; 
% end

function [N_ks,N_et]=Basis_Derivative(ks,et)
N_ks(1)=-1/4*(1-et); 
N_ks(2)=1/4*(1-et); 
N_ks(3)=1/4*(1+et); 
N_ks(4)=-1/4*(1+et); 
N_et(1)=-1/4*(1-ks); 
N_et(2)=-1/4*(1+ks); 
N_et(3)=1/4*(1+ks); 
N_et(4)=1/4*(1-ks); 
end

%--------------------------- Jacobi  -------------------------------
function J=Jacobi(iE,ks,et,x,y,element_nodes)
nodes=element_nodes(iE,:);
x_element=x(nodes);y_element=y(nodes);
[N_ks,N_et]=Basis_Derivative(ks,et);
J=[N_ks;N_et]*[x_element y_element];
end

%------------------------- Strain matrix  --------------------------
function B=strain(iE,ks,et,x,y,element_nodes)
[N_ks,N_et]=Basis_Derivative(ks,et);
J=Jacobi(iE,ks,et,x,y,element_nodes);
for i=1:4
    B1=J(2,2)*N_ks(i)-J(1,2)*N_et(i); 
    B2=-J(2,1)*N_ks(i)+J(1,1)*N_et(i); 
    B(1:3,2*i-1:2*i)=[B1 0;0 B2;B2 B1]/det(J); 
end
end

%-------------------- Local stiffness matrix  ----------------------
function stiffLocal=Quad_Local_stiff(iE,x,y,element_material,element_nodes,material,type,n_gauss)
global D
    n=element_material(iE);
    if type>0
        material(:,1)=material(:,1)/(1-material(:,2)^2);
        material(:,1)=material(:,2)/(1-material(:,2));
    end
    E=material(n,1);um=material(n,2);t=material(n,3);
    D=(E/(1-um*um))*[1 um 0;um 1 0;0 0 (1-um)/2]; % Elastic matrix
    [xg,w]=gauss(n_gauss);
    st=zeros(8,8);
    for i=1:n_gauss
        for j=1:n_gauss
            B=strain(iE,xg(i),xg(j),x,y,element_nodes);
            J=Jacobi(iE,xg(i),xg(j),x,y,element_nodes);
            st=st+w(i)*w(j)*B'*D*B*t*det(J);
        end
    end
    stiffLocal=st;
end

%----------------------  Solve element stress  ---------------------
function [stress,stress1]=Quad_stress_solve(iE,displacement,x,y,element_nodes)
global D
    location(1:2:7)=2*element_nodes(iE,:)-1;
    location(2:2:8)=2*element_nodes(iE,:);
    w=displacement(location);
    J=Jacobi(iE,0,0,x,y,element_nodes);
    [N_ks,N_et]=Basis_Derivative(0,0);
    for i=1:4
        B1=J(2,2)*N_ks(i)-J(1,2)*N_et(i); 
        B2=-J(2,1)*N_ks(i)+J(1,1)*N_et(i); 
        B(1:3,2*i-1:2*i)=[B1 0;0 B2;B2 B1]/det(J); 
    end
    stress(1,:)=(D*B*w)';
    st=stress;
    stress1(1,1)=(st(1)+st(2))/2+sqrt(((st(1)-st(2))/2)^2+st(3)^2);
    stress1(1,2)=(st(1)+st(2))/2-sqrt(((st(1)-st(2))/2)^2+st(3)^2); 
    stress1(1,3)=90-1/2/pi*180*atan(2*st(3)/(st(1)-st(2)));   
end


%-------------------- 0-1 replacement method  ----------------------
function [stiffGlobal2,loadGlobal2]=Replace01(stiffGlobal,loadGlobal,constraint)
for i=1:size(constraint,2)
    stiffGlobal(constraint(i),:)=zeros();
    stiffGlobal(:,constraint(i))=zeros();
    stiffGlobal(constraint(i),constraint(i))=1;
    loadGlobal(constraint(i))=0;
    stiffGlobal2=stiffGlobal;
    loadGlobal2=loadGlobal;
end
end
