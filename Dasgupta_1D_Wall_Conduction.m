clc
close all
%Steady one dimensional heat conduction program
clear all

%Input data
L=1.5; %[m]
%N=8; %number of grid points
N_list=4:4:16; 
N=N_list(end); %Number of grid points
delta_x=L/(N_list(end)-2); %Thickness for each grid point
alpha=20.0; %thermal diffusivity [m^2/s]
T_left=350.0; %Temperature left wall [K]
T_right=300.0; %Temperature right wall [K]
Qt=50; %Heat source proportional to temperature, i.e.Qt*T, [1/s]
maxIter=1e6; %Maximum number of iterations for the solver
maxRes=1e-6; %Solver residual |Ax-b|<maxRes

Error=zeros(1,length(N_list));
cas=1;
for N=N_list
%Mesh
%posX=0:L/(N-1):L;
x=linspace(0,L,N);

%Preallocate coefficients & mesh distances
ap=ones(1,N);
ae=zeros(1,N);
aw=zeros(1,N);
b=zeros(1,N);
Dxe=zeros(1,N);
Dxw=zeros(1,N);
Dx=zeros(1,N);

%Preallocate vector for local residual of solver
loc_res=zeros(1,N);

%Inner coefficients
for iX=2:N-1
  Dxe(iX)=(x(iX+1)-x(iX));
  Dxw(iX)=(x(iX)-x(iX-1));
  Dx(iX)=(Dxe(iX)+Dxw(iX))/2;
  ae(iX)=alpha/Dxe(iX);
  aw(iX)=alpha/Dxw(iX);
  ap(iX)=-(ae(iX)+aw(iX)+Qt*Dx(iX));
  b(iX)=0;
end


%Boundary coefficients
b(1)=T_left;
b(end)=T_right;

%Initialize temperature vector
%T=ones(size(x))*((Tw+Te)*0.5);
T=zeros(size(x));
T_reset=zeros(size(x)); %previous iteration array for Jacobi method
T(1)=T_left;
T(end)=T_right;

%Solver
res=maxRes+1;
ite=0;
tic;
while res>maxRes && ite<maxIter
  %Jacobi iteration
  T_reset=T;
  T(1)=(-ae(1)*T_reset(2)+b(1))/ap(1);
  
  for iX=2:numel(T)-1
    T(iX)=(-aw(iX)*T_reset(iX-1)-ae(iX)*T_reset(iX+1)+b(iX))/ap(iX);
  end
  T(end)=(-aw(end)*T_reset(end-1)+b(end))/ap(end);
  
  %Calculation of the solver residual res=|Ax-b|
  for iX=2:numel(T)-1
    loc_res(iX)=aw(iX)*T(iX-1)+ap(iX)*T(iX)+ae(iX)*T(iX+1)-b(iX);
  end
  res=max(abs(loc_res));

  ite=ite+1;
  if mod(ite,10000)==0
    fprintf('ite: %d solver residual: %e\n',ite,res);
  end
end
if ite==maxIter
  warning(['Maximum number of iterations reached (lastRes=',num2str(res),').']);
end

%Matlab solver for systems of linear equations
%A is the full matrix. Try full(A) to see the whole matrix
A=spdiags([[aw(2:end)'; 0] ap(:) [0; ae(1:end-1)']],[-1 0 1],N,N);
T2=A\b(:);
toc;

%Tana=@(x) (-0.5*Q/lambda)*x.^2+((T_right-T_left+0.5*Q/lambda*L^2)/L)*x+T_left; %analytic solution
k=sqrt(Qt/alpha);
C2 = (T_right-T_left*exp(k*L))/(exp(-k*L)-exp(k*L));
C1 = T_left-C2;
Tana=@(x) C1*exp(k*x)+C2*exp(-k*x); %analytic solution

%plot(x,T);
%hold on
%plot(x,T2);
[AX,H1,H2] = plotyy(x,T,x,feval(Tana,x)-T);
legend('T numerical','Error')
ylabel(AX(1),'Temperature') % left y-axis 
ylabel(AX(2),'Error') % right y-axis
title('Temperature distribution (steady case)');
drawnow

Error(cas)=max(abs(feval(Tana,x)-T));

cas=cas+1;
end
%% Plotting Steady
figure(2);
%label={'Grids=4','Grids=8','Grids=12','Grids=16'};
loglog(L./N_list,Error,'-o','LineWidth',1);
%text(L./N_list,Error,label)
ylabel('Error');
xlabel('Dx')
title('Truncation error analysis (steady case)');

%print('tempDist','-dpng');

%% Transient code for the same case 
%% Unsteady Numerical
trun_trans=ones(1,length(N_list)); %zeros(1,length(N_list))
tas=1;

%N=8;
%for N=N_list
Time_start=0; 
Time_end=5;
Time_step=0.20*(delta_x^2/alpha); %requirement for convergence---> t<=0.25*(delta_x^2/alpha)
x=linspace(0,L,N);
A=50;
beta=12.56637062;
T_prev=C1*exp(k*x)+C2*exp(-k*x)+A*sin(beta*x);%ones(1,N);

Unsteady_Residue=10; %initial residue to start the process
j=1; %indexing for time in rows
t=Time_start;
while Unsteady_Residue>maxRes && t<Time_end %specifying accuracy and end time
    Time(j,1)=t; %time step array
    T_unsteady(j,1)=b(1);
    for i=2:1:N-1 %indexing for grid points in columns
        T_unsteady(j,i)=T_prev(i)+(Time_step/delta_x)*(ap(i)*T_prev(i)+aw(i)*T_prev(i-1)+ae(i)*T_prev(i+1));
    end
    T_unsteady(j,N)=b(end);

    Unsteady_Residue=max(abs(T_unsteady(j,:)-T_prev));

    T_prev=T_unsteady(j,:);
    t=t+Time_step;
    j=j+1;
end

%% Unsteady Analytical

Tana_unst=@(a,b) C1*exp(k*a)+C2*exp(-k*a)+A*exp(-(alpha*beta^2+Qt)*b)*sin(beta*a); %analytic solution function
TUnst_ana=[];
Error_transient=zeros(1,N);
z=1;
for l=Time' %t=0 to t=end 
    w=1;
    L_step=x;
    for k=L_step %L=0 to L=1.5
        TUnst_ana(z,w)=Tana_unst(k,l); %time in rows and length in column 
        w=w+1;
    end
    Error_transient(z,:)=(TUnst_ana(z,:))-(T_unsteady(z,:)); % Analytical - Numerical in transient
    
    z=z+1;
end 
% trun_trans(tas,:)=max(abs(T_unsteady-TUnst_ana),[],'all');
% tas=tas+1;
% end 
%% Plotting Unsteady Numerical 
%Plotting for Transient (x,0.0002)
for plt=1:10:length(T_unsteady)
   figure(3)
   plot_T=T_unsteady(plt,:); %Numerical value of each time step
   plot_E=Error_transient(plt,:);
   [AX,H1,H2] = plotyy(x,plot_T,x,plot_E);
   legend('T numerical','Error')
   ylabel(AX(1),'Temperature') % left y-axis 
   ylabel(AX(2),'Error') % right y-axis
   title('Temperature distribution (transient case)');
   drawnow
end 
%% Plotting Unsteady Aanlytical 
figure(4) %Temperature at t=0.0000s
plot (L_step,TUnst_ana(1,:),'-o')
hold on
plot (x,T_unsteady(1,:))
ylabel('Temperature (K) at 0 sec')
xlabel('Distance (m)')
title ('Instantaneous Transient temperature at 0 sec')
legend('T Analytical','T numerical')

figure(5) %Temperature at t=0.0002s
plot (L_step,TUnst_ana(3,:),'-o')
hold on
plot (x,T_unsteady(3,:))
hold off
ylabel('Temperature (K) at 0.0002 sec')
xlabel('Distance (m)')
title ('Instantaneous Transient temperature at 0.0002 sec')
legend('T Analytical','T numerical')

figure(6) %Temperature at t=0.0004s
plot (L_step,TUnst_ana(5,:),'-o')
hold on
plot (x,T_unsteady(5,:))
hold off
ylabel('Temperature (K) at 0.0004 sec')
xlabel('Distance (m)')
title ('Instantaneous Transient temperature at 0.0004 sec')
legend('T Analytical','T numerical')

figure(7) %Temperature at t=0.0006s
plot (L_step,TUnst_ana(7,:),'-o')
hold on
plot (x,T_unsteady(7,:))
hold off
ylabel('Temperature (K) at 0.0006 sec')
xlabel('Distance (m)')
title ('Instantaneous Transient temperature at 0.0006 sec')
legend('T Analytical','T numerical')


figure(8) %Temperature after reaching steady condition
p=1;
for k=x; %L=0 to L=1.5
        TUnstana_plot(:,p)=Tana_unst(k,Time(end)); %time in rows and length in column
        p=p+1;
end 

plot (x,TUnstana_plot,'-o') %plotting transient analytical solution
hold on 
plot (x,T_unsteady(end,:),'LineWidth',2) %plotting transient numerical solution 
hold off
ylabel('Temperature')
xlabel('Distance (m)')
title ('Transient Analytical vs Numerical Solution')
legend('T Analytical','T numerical')

% %% Plotting truncation error (transient) 
% figure(9)
% label={'Grids=4','Grids=8','Grids=12','Grids=16'};
% loglog(L./N_list,trun_trans,'-o','LineWidth',1);
% text(L./N_list,trun_trans,label)
% ylabel('Error');
% xlabel('Dx')
% title('Truncation error analysis (Transient case)');
