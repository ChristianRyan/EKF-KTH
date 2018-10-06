clear all
load('GPSdata.mat')
est = NonLinearLeastSquares(gps_data,ref_data_struct);
%plot3(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:), ref_data_struct.traj_ned(3,:))
%hold on
%figure
%plot(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:))
% plot3(est.x_h(1,:), est.x_h(2,:), est.x_h(3,:))
% hold on
% plot3(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:), ref_data_struct.traj_ned(3,:))
%plot(est.x_h(1,:), est.x_h(2,:))

sat_with_val = [1, 8, 10, 19, 20, 25, 26];
conf=[1,3,5,6,7];

N=length(gps_data(1).PseudoRange);
M = length(gps_data);
y = [];
satpos = [];
for i = 1:length(sat_with_val)
    y = [y; gps_data(sat_with_val(i)).PseudoRange];
    satpos = [satpos; gps_data(sat_with_val(i)).Satellite_Position_NED];
end

T = ref_data_struct.Ts;
c = ref_data_struct.c;
R_k = ref_data_struct.s2r;

%est_x =[-0.21,2,0.0050,2,-0.001,0.00000667,0.000001]'%zeros(7,1);
%est_x=[1.9873,-0.01,-1.3992,-0.01,1.1003,0.000000000001,0]'
est_x = [ref_data_struct.traj_ned(1,1), 0, ref_data_struct.traj_ned(1,2),0 , ref_data_struct.traj_ned(1,3), ref_data_struct.x_clk(1,1), ref_data_struct.x_clk(1,2)]'%ones(7,1)-[0,1,0,1,0,0.99999,0.99999]'%zeros(7,1);
%est_x=zeros(7,1)
%P_k = blkdiag(1,1,1,1,1,0*eye(2))%zeros(7,7)ones(7,7);0.01*eye(7);
F_2 = [1 T; 0 1];
Q_2 = [T^4/3 T^3/2; T^3/2 T^2];
Q_clk = [ref_data_struct.PSD_clk(1)*T + ref_data_struct.PSD_clk(2)*T^3/3, T^2*ref_data_struct.PSD_clk(2); T^2*ref_data_struct.PSD_clk(2), ref_data_struct.PSD_clk(2)*T];
S_x = 3000*10^-5;
S_y = S_x;
S_z = 0.01*10^-8;
Q_k = blkdiag((S_x/T).*Q_2,(S_y/T).*Q_2,S_z/T,c^2* Q_clk);
G_k =  eye(7);
F_k = blkdiag(F_2, F_2, 1, F_2);
%H_k = 2/(abs(satpos(1)-est_x))*[est_x(1)-satpos(1,1), 0,est_x(2)-satpos(1,2), 0,est_x(3)-satpos(1,3), c, 0]; 
S_k = zeros(7,1);
P_k=Q_k;
for n = 1:N
    H_k = [];
    e_k = [];
    for a = 1:length(sat_with_val)
        H_k_1 = (1/(norm(gps_data(sat_with_val(a)).Satellite_Position_NED(:,n)- est_x([1,3,5],n))))*[est_x(1,n)-gps_data(sat_with_val(a)).Satellite_Position_NED(1,n), 0,est_x(3,n)-gps_data(sat_with_val(a)).Satellite_Position_NED(2,n), 0,est_x(5,n)-gps_data(sat_with_val(a)).Satellite_Position_NED(3,n), 0, 0]+[0,0,0,0,0,1,0]; 
        H_k = [H_k; H_k_1];
        e_k = [e_k; y(a,n)-(norm(gps_data(sat_with_val(a)).Satellite_Position_NED(:,n)- est_x([1,3,5],n)))-est_x(6,n)];
    end
    %e_k = y(:,n) -H_k*est_x(:,n);
    
        R_ek = H_k*P_k*H_k' + R_k*eye(7);
        K_k = (F_k*P_k*H_k')/(R_ek);
    

    oo(:,n)=[ref_data_struct.traj_ned(:,n)- est_x([1,3,5],n);c*ref_data_struct.x_clk(:,n)-est_x([6,7],n)];
    for kkk=1:5
        if kkk<=3
        oop(kkk,n)=oo(kkk,n)+norm(3*sqrt(P_k(conf(kkk),conf(kkk))));
        ool(kkk,n)=oo(kkk,n)-norm(3*sqrt(P_k(conf(kkk),conf(kkk))));
        else 
        oop(kkk,n)=(oo(kkk,n)+norm(3*sqrt(P_k(conf(kkk),conf(kkk)))));
        ool(kkk,n)=(oo(kkk,n)-norm(3*sqrt(P_k(conf(kkk),conf(kkk)))));
        end
    end
    est_x(:,n+1) = F_k*est_x(:,n) + K_k*e_k; 

    
    P_k = F_k*P_k*F_k' + G_k*Q_k*G_k' - K_k*R_ek*K_k';
    
    b(:,n)=est_x(:,n+1)+P_k*H_k'/R_ek*e_k;
end
%%
plot(b(1,500:end), b(3,500:end));
hold on
plot(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:))
plot(est.x_h(1,:),est.x_h(2,:))
hold off
%plot3(est_x(1,:), est_x(3,:),est_x(5,:)) 
%plot3(b(1,:), b(3,:),b(5,:)) 
%plot3(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:), ref_data_struct.traj_ned(3,:))
%figure 
%plot(est_x(1,:), est_x(3,:))


%% x
figure(2)




plot(ref_data_struct.traj_ned(1,:),ref_data_struct.traj_ned(2,:))
hold on
plot(est.x_h(1,:),est.x_h(2,:))
plot(b(1,500:end), b(3,500:end));
legend('True','Non-linear LS','EKF estimation')
title('The estimation of the vehicle position for large covariance')
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off



figure (3)
errorx=ref_data_struct.traj_ned(1,:)-est_x(1,1:2015);
sqerrorx=ref_data_struct.traj_ned(1,:)-est.x_h(1,:);
plot(errorx)
hold on
plot(sqerrorx)
title('Errors of X axis')
xlabel('index')
ylabel('error')
plot(ool(1,:))
plot(oop(1,:))
legend('EKF estimation','non-linear LS','3\sigma interval','3\sigma interval')

hold off

figure (4)
errory=ref_data_struct.traj_ned(2,:)-est_x(3,1:2015);
sqerrory=ref_data_struct.traj_ned(2,:)-est.x_h(2,:);
plot(errory)
hold on
plot(sqerrory)
title('Errors of Y axis')
xlabel('index')
ylabel('error')
plot(ool(2,:))
plot(oop(2,:))
legend('EKF estimation','non-linear LS','3\sigma interval','3\sigma interval')
hold off


figure (5)
errorz=ref_data_struct.traj_ned(3,:)-est_x(5,1:2015);
sqerrorz=ref_data_struct.traj_ned(3,:)-est.x_h(3,:);
plot(errorz)
hold on
plot(sqerrorz)
title('Errors of Z axis')
xlabel('index')
ylabel('error')
plot(ool(3,:))
plot(oop(3,:))
legend('EKF estimation','non-linear LS','3\sigma interval','3\sigma interval')

hold off

figure (6)
errort=ref_data_struct.x_clk(1,:).*c-est_x(6,1:2015);
sqerrort=ref_data_struct.x_clk(1,:).*c-est.x_h(4,:);
plot(errort)
title('Errors of clock offset (meter)')
xlabel('index')
ylabel('error')
hold on
plot(sqerrort)
plot(ool(4,:))
plot(oop(4,:))
legend('EKF estimation','non-linear LS','3\sigma interval','3\sigma interval')
hold off

figure (7)
errort2=ref_data_struct.x_clk(2,:).*c-est_x(7,1:2015);
plot(errort2)
title('Errors of clock drift(meter)')
xlabel('index')
ylabel('error')
hold on
plot(ool(5,:))
plot(oop(5,:))
legend('EKF estimation','3\sigma interval','3\sigma interval')
hold off