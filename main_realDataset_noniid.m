clc; clear; close all;
%% initial workspace
rng('default')
rng('shuffle')
linewidth = 2.5;
fontsize = 14;
MarkerSize = 10;

SNR=10;
P=1E-3;
bandwidth = 1e6;
rate = bandwidth*log2(1+SNR);

%load('Datasets/a1a.mat'); %rho = 500
%load('Datasets/w7a.mat'); %rho = 500
%load('Datasets/w8a.mat'); %rho = 500
%load('Dataset_noniid/phishing.mat');

load('Dataset_noniid/a9a.mat'); %rho = 500

%load('datasets_clean/a1a_clean.mat'); %rho = 500
%load('datasets_clean/w7a_clean.mat'); %rho = 500
%load('datasets_clean/w8a_clean.mat'); %rho = 500
%load('datasets_clean/phishing_clean.mat');
%load('datasets_clean/a9a_clean.mat'); %rho = 500

%load('Datasets/a1a.mat'); %rho = 500
%load('Datasets/w7a.mat'); %rho = 500
%load('Datasets/w8a.mat'); %rho = 500
%load('Datasets/phishing.mat');

%load('Datasets/a9a.mat'); %rho = 500

% for i=1:size(X,1)
%     for j=2:size(X,2)
%         X(i,j)=X(i,j)-0.5+rand();
%     end
% end

X = reshape(X,size(X,1)*size(X,2),size(X,3));
for j=2:size(X,2)
    temp1=abs(X(:,j));
    temp=max(temp1);
        for i=1:size(X,1)
            %XXXX(i,j)=(XXX(i,j)-mean(XXX(:,j)))/(max(XXX(:,j))-min(XXX(:,j)));
            X(i,j)=X(i,j)/temp;
        end
end

num_feature=size(X,2);
total_sample=size(y',1);

%% setup hyperparameters
para.dimx=num_feature+1;%size(X,2);%20;
para.node=100;
para.bs=ones(para.node,1)*300;%400;%round(randi([10,500],para.node,1));
para.iter=1000;
para.repeat=1;

%K=sum(para.bs);
bitsToSend = 2;
%% setup data
%temp=X';
data.features = [X';ones(1,size(X',2))];%[randn(para.dimx-1,K);ones(1,K)];
data.labels=y;%zeros(K,1);



%data.labels=data.labels';



%% setup parameters (stepsizes)


%% setup logs
log_SGD = zeros(para.iter,2,para.repeat);
log_QSGD = zeros(para.iter,2,para.repeat);
log_GD = zeros(para.iter,2,para.repeat);
log_GD_1 = zeros(para.iter,2,para.repeat);
log_VR = zeros(para.iter,2,para.repeat);
log_QVR = zeros(para.iter,2,para.repeat);
log_LGD = zeros(para.iter,2,para.repeat);
log_LSGD = zeros(para.iter,2,para.repeat);
log_Prox = zeros(para.iter,2,para.repeat);
log_QGD = zeros(para.iter,2,para.repeat);
log_QGD2 = zeros(para.iter,2,para.repeat);
log_QGD_1 = zeros(para.iter,2,para.repeat);

log_GD_2 = zeros(para.iter,2,para.repeat);
log_QGD_2 = zeros(para.iter,2,para.repeat);
log_QGD_2_2 = zeros(para.iter,2,para.repeat);

g=@gradx;
f=@fx;
%% main iteration
for loop_i=1:para.repeat
    init_x=zeros(para.node,para.dimx);
    %Deterministic
    para.W=ones(para.node)/para.node;
    para.stepsize = 4;
    
    para.I = 1;
    para.R = 1;
    para.VR=0;
    para.mini_bs = 0.0025;
    para.localiter = para.iter;
    para.stepsize = 4*2;
    [log_SGD(:,:,loop_i),itrSGD, number_of_bits_toSend_nq0, Energy_nq0] = Fed_PD(para,data,g,f,init_x, P, rate);
    
    [log_QSGD(:,:,loop_i),itrQSGD, number_of_bits_toSend_q0, Energy_q0] = QFed_PD(para,data,g,f,init_x, bitsToSend, P, rate);
    
    para.VR=0;
    para.mini_bs = 1;
    para.localiter = 8;
    para.R = 2;%1.5;%2;
    [log_GD(:,:,loop_i), itr1, number_of_bits_toSend_nq1, Energy_nq1] = Fed_PD(para,data,g,f,init_x, P, rate);
    
    [log_QGD(:,:,loop_i), itr1_q, number_of_bits_toSend_q1, Energy_q1] = QFed_PD(para,data,g,f,init_x, bitsToSend, P, rate);
    
% %     bitsToSend = 4;
% %     
% %     [log_QGD2(:,:,loop_i), itr1_q2] = QFed_PD(para,data,g,f,init_x, bitsToSend);
% %     
    para.VR=0;
    para.mini_bs = 1;
    %para.localiter = 8;
    para.R = 1;
    [log_GD_1(:,:,loop_i),itr2, number_of_bits_toSend_nq2, Energy_nq2] = Fed_PD(para,data,g,f,init_x, P, rate);
    
    [log_QGD_1(:,:,loop_i),itr2_q, number_of_bits_toSend_q2, Energy_q2] = QFed_PD(para,data,g,f,init_x, bitsToSend, P, rate);
% %     
    para.I = 100;
    para.R =1;
    para.VR=1;
    para.mini_bs = 0.0025;
    para.localiter = 2;
    [log_VR(:,:,loop_i), itrVR, number_of_bits_toSend_nq_vr, Energy_nq_vr] = Fed_PD(para,data,g,f,init_x, P, rate);
    [log_QVR(:,:,loop_i), itrQVR, number_of_bits_toSend_q_vr, Energy_q_vr] = QFed_PD(para,data,g,f,init_x, bitsToSend, P, rate);
%     
    para.mini_bs = 1;
    para.localiter = 8;
    para.stepsize = 4/para.localiter;
    log_LGD(:,:,loop_i) = mainiter_SGD(para,data,g,f,init_x);
    
    para.mini_bs = 0.0025;
    para.localiter = para.iter;
    para.stepsize = 4;
    log_LSGD(:,:,loop_i) = mainiter_SGD(para,data,g,f,init_x);


%% draw figures
para.VR=1;
para.mini_bs = 0.0025;
para.localiter = 8;
para.R = 1;
para.mu = 0.1;
log_Prox(:,:,1) = Fed_Prox(para,data,g,f,init_x);



para.VR=0;
para.mini_bs = 1;
para.localiter = 8;
para.R = -1;
[log_GD_2(:,:,loop_i),itr3, number_of_bits_toSend_nq3, Energy_nq3] = Fed_PD(para,data,g,f,init_x, P, rate);

%bitsToSend = 2;


para.VR=0;
para.mini_bs = 1;
%para.localiter = 8;
para.R = -1;
[log_QGD_2(:,:,loop_i),itr3_q, number_of_bits_toSend_q3, Energy_q3] = QFed_PD(para,data,g,f,init_x, bitsToSend, P, rate);

end
number_of_bits_toSend_nq=zeros(para.iter,1);
number_of_bits_toSend_nq(1)=32*para.dimx*para.node;

for i=2:para.iter
    number_of_bits_toSend_nq(i)=number_of_bits_toSend_nq(i-1)+32*para.dimx*para.node;
end

% bitsToSend = 4;
% 
% 
% %para.VR=0;
% %para.mini_bs = 1;
% %para.localiter = 8;
% %para.R = -1;
% [log_QGD_2_2(:,:,loop_i),itr3_q_2] = QFed_PD(para,data,g,f,init_x, bitsToSend);

%print('plotting figures');

%% loss



figure(1);
% idx_xfilter=cumsum(log_xFilter(:,1,R=1))*400;
%idx = 0:para.iter-1;

semilogy(cumsum(itr1)-1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
hold on
semilogy(cumsum(itr1_q)-1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
%semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');

semilogy(cumsum(itr2)-1, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
semilogy(cumsum(itr2_q)-1, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');

semilogy(cumsum(itr3)-1, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');

semilogy(cumsum(itr3_q)-1, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');
hold off
xl = xlabel('Communication r','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T');%'EXTRA','xFilter','DSGD','GNSD' );
xlim([0 500])







figure(2);

semilogy(number_of_bits_toSend_nq1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
hold on
semilogy(number_of_bits_toSend_q1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
%semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');

semilogy(number_of_bits_toSend_nq2, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
semilogy(number_of_bits_toSend_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');

semilogy(number_of_bits_toSend_nq3, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');

semilogy(number_of_bits_toSend_q3, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');

hold off
xl = xlabel('number of transmitted bits','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T');%'EXTRA','xFilter','DSGD','GNSD' );
%xlim([0 500])







figure(3);

semilogy(Energy_nq1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
hold on
semilogy(Energy_q1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
%semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');

semilogy(Energy_nq2, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
semilogy(Energy_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');

semilogy(Energy_nq3, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');

semilogy(Energy_q3, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');

hold off
xl = xlabel('Total energy bill','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T');%'EXTRA','xFilter','DSGD','GNSD' );
%xlim([0 500])

iter_g = 0:para.iter-1;

% save non_iid_fixedQR.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
%     number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
%     number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
%     itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
%     number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq

save non_iid_v2.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
    number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
    number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
    itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
    number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq

figure(4);
% idx_xfilter=cumsum(log_xFilter(:,1,R=1))*400;
%idx = 0:para.iter-1;
semilogy(cumsum(itrSGD)-1, mean(log_SGD(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'b');
hold on
semilogy(cumsum(itrQSGD)-1, mean(log_QSGD(:,2,:),3),'linestyle', '--','linewidth',linewidth,'color', 'b');

semilogy(cumsum(itr1)-1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
hold on
semilogy(cumsum(itr1_q)-1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
%semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');

semilogy(cumsum(itr2)-1, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
semilogy(cumsum(itr2_q)-1, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');

semilogy(cumsum(itr3)-1, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');

semilogy(cumsum(itr3_q)-1, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');

%plot(cumsum(itr3_q_2)-1, mean(log_QGD_2_2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'k');
% hold on
% plot(cumsum(itr3)-1, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'g');
% hold on
semilogy(cumsum(itrVR)-1, mean(log_VR(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'k');
semilogy(cumsum(itrQVR)-1, mean(log_QVR(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'k');
% hold on
semilogy((iter_g), mean(log_LGD(:,2,:),3),'linestyle', '--','linewidth',linewidth,'color', 'm');
% hold on
semilogy((iter_g), mean(log_LSGD(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'm');
% hold on
semilogy((iter_g), mean(log_Prox(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');
hold off
xl = xlabel('Communication r','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-SGD','FedPD-QSGD','FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T','FedPD-VR','FedPD-QVR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
xlim([0 500])







figure(5);
semilogy(number_of_bits_toSend_nq0, mean(log_SGD(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'b');
hold on
semilogy(number_of_bits_toSend_q0, mean(log_QSGD(:,2,:),3),'linestyle', '--','linewidth',linewidth,'color', 'b');

semilogy(number_of_bits_toSend_nq1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
hold on
semilogy(number_of_bits_toSend_q1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
%semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');

semilogy(number_of_bits_toSend_nq2, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
semilogy(number_of_bits_toSend_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');

semilogy(number_of_bits_toSend_nq3, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');

semilogy(number_of_bits_toSend_q3, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');

%plot(cumsum(itr3_q_2)-1, mean(log_QGD_2_2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'k');
% hold on
% plot(cumsum(itr3)-1, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'g');
% hold on
semilogy(number_of_bits_toSend_nq_vr, mean(log_VR(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'k');
semilogy(number_of_bits_toSend_q_vr, mean(log_QVR(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'k');
% hold on
semilogy(number_of_bits_toSend_nq, mean(log_LGD(:,2,:),3),'linestyle', '--','linewidth',linewidth,'color', 'm');
% hold on
semilogy(number_of_bits_toSend_nq, mean(log_LSGD(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'm');
% hold on
semilogy(number_of_bits_toSend_nq, mean(log_Prox(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');
hold off
xl = xlabel('number of transmitted bits','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-SGD','FedPD-QSGD','FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T','FedPD-VR','FedPD-QVR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
%xlim([0 500])




% load nono_iid.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
%     number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
%     number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
%     itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
%     number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq
% 



