clc; clear; close all;
%% initial workspace
rng('default')
rng('shuffle')
linewidth = 2.5;
fontsize = 14;
MarkerSize = 10;

%% setup hyperparameters
para.dimx=20;
para.node=100;
para.bs=ones(para.node,1)*400;%*400;%round(randi([10,500],para.node,1));
para.iter=600;
para.repeat=1;

K=sum(para.bs);
bitsToSend = 2;
%% setup data

% data.features=[randn(para.dimx-1,K);ones(1,K)];
% data.labels=randi([0,1],K,1);
% data.labels(data.labels==0) = -1;
% % figure(1);
% % subplot(1,2,1)
% % plot(data.labels)


data.features = [randn(para.dimx-1,K);ones(1,K)];
data.labels=zeros(K,1);

past=1;
for i=1:para.node
    %x=rand(para.dimx,1)*40+randi([-10,10],para.dimx,1);
    pp=rand();
    pp=round(pp);
    if(pp==0)
        pp=-1;
    end
    current=para.bs(i);
    data.labels(past:past+current-1) = pp;%data.features(:,past:past+current-1)'*x;
    past=past+current;
end


% past=1;
% for i=1:para.node
%     x=rand(para.dimx,1)*40+randi([-10,10],para.dimx,1);
%     current=para.bs(i);
%     data.labels(past:past+current-1) = data.features(:,past:past+current-1)'*x;
%     past=past+current;
% end
% 
% threshold = randn(K,1)*max(data.labels);
% data.labels((data.labels>0)&(data.labels>threshold)) = 1;
% data.labels(data.labels>0&data.labels<threshold) = 2;
% data.labels(data.labels<0&(-data.labels)>threshold) = 2;
% data.labels(data.labels<0&(-data.labels)<threshold) = 1;
% data.labels(data.labels==2) = -1;

% % subplot(1,2,2)
% % plot(data.labels)

%data.labels

data.labels=data.labels';



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
    [log_SGD(:,:,loop_i),itrSGD, number_of_bits_toSend_nq0] = Fed_PD(para,data,g,f,init_x);
    
    [log_QSGD(:,:,loop_i),itrQSGD, number_of_bits_toSend_q0] = QFed_PD(para,data,g,f,init_x, bitsToSend);
    
    para.VR=0;
    para.mini_bs = 1;
    para.localiter = 8;
    para.R = 1.5;%2;
    [log_GD(:,:,loop_i), itr1, number_of_bits_toSend_nq1] = Fed_PD(para,data,g,f,init_x);
    
    [log_QGD(:,:,loop_i), itr1_q, number_of_bits_toSend_q1] = QFed_PD(para,data,g,f,init_x, bitsToSend);
    
% %     bitsToSend = 4;
% %     
% %     [log_QGD2(:,:,loop_i), itr1_q2] = QFed_PD(para,data,g,f,init_x, bitsToSend);
% %     
    para.VR=0;
    para.mini_bs = 1;
    %para.localiter = 8;
    para.R = 1;
    [log_GD_1(:,:,loop_i),itr2, number_of_bits_toSend_nq2] = Fed_PD(para,data,g,f,init_x);
    
    [log_QGD_1(:,:,loop_i),itr2_q, number_of_bits_toSend_q2] = QFed_PD(para,data,g,f,init_x, bitsToSend);
% %     
    para.I = 100;
    para.R =1;
    para.VR=1;
    para.mini_bs = 0.0025;
    para.localiter = 2;
    [log_VR(:,:,loop_i), itrVR, number_of_bits_toSend_nq_vr] = Fed_PD(para,data,g,f,init_x);
    [log_QVR(:,:,loop_i), itrQVR, number_of_bits_toSend_q_vr] = QFed_PD(para,data,g,f,init_x, bitsToSend);
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
[log_GD_2(:,:,loop_i),itr3, number_of_bits_toSend_nq3] = Fed_PD(para,data,g,f,init_x);

%bitsToSend = 2;


para.VR=0;
para.mini_bs = 1;
%para.localiter = 8;
para.R = -1;
[log_QGD_2(:,:,loop_i),itr3_q, number_of_bits_toSend_q3] = QFed_PD(para,data,g,f,init_x, bitsToSend);

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

% figure(1);
% % idx_xfilter=cumsum(log_xFilter(:,1,R=1))*400;
% idx = 0:para.iter-1;
% 
% semilogy(cumsum(itr1)-1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
% hold on
% semilogy(cumsum(itr1_q)-1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
% %semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');
% 
% semilogy(cumsum(itr2)-1, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
% semilogy(cumsum(itr2_q)-1, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');
% 
% semilogy(cumsum(itr3)-1, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');
% 
% semilogy(cumsum(itr3_q)-1, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');
% hold off
% xl = xlabel('Communication r','FontSize',fontsize,'interpreter','latex');
% yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
% %le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
% le = legend('FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T');%'EXTRA','xFilter','DSGD','GNSD' );
% xlim([0 500])
% 
% 
% 
% 
% 
% 
% 
% figure(2);
% 
% semilogy(number_of_bits_toSend_nq1, mean(log_GD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'b');
% hold on
% semilogy(number_of_bits_toSend_q1, mean(log_QGD(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
% %semilogy(cumsum(itr1_q2)-1, mean(log_QGD2(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');
% 
% semilogy(number_of_bits_toSend_nq2, mean(log_GD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
% semilogy(number_of_bits_toSend_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'g');
% 
% semilogy(number_of_bits_toSend_nq3, mean(log_GD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'c');
% 
% semilogy(number_of_bits_toSend_q3, mean(log_QGD_2(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'm');
% 
% hold off
% xl = xlabel('number of transmitted bits','FontSize',fontsize,'interpreter','latex');
% yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
% %le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
% le = legend('FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T');%'EXTRA','xFilter','DSGD','GNSD' );
% %xlim([0 500])
% 
% 






figure(1);
% idx_xfilter=cumsum(log_xFilter(:,1,R=1))*400;
idx = 0:para.iter-1;
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
semilogy((0:para.iter-1), mean(log_LGD(:,2,:),3),'linestyle', '--','linewidth',linewidth,'color', 'm');
% hold on
semilogy((0:para.iter-1), mean(log_LSGD(:,2,:),3),'linestyle', ':','linewidth',linewidth,'color', 'm');
% hold on
semilogy((0:para.iter-1), mean(log_Prox(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'm');
hold off
xl = xlabel('Communication r','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-SGD','FedPD-QSGD','FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T','FedPD-VR','FedPD-QVR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
xlim([0 500])







figure(2);
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


