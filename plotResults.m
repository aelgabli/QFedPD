clc; clear; close all;
%% initial workspace
rng('default')
rng('shuffle')
linewidth = 2.5;
fontsize = 14;
MarkerSize = 10;

load non_iid_v2.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
    number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
    number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
    itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
    number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq

% load iid_v5.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
%     number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
%     number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
%     itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
%     number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq




figure(1);
subplot(1,3,1)
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
xl = xlabel({'Communication r';'(a)'},'FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('FedPD-SGD','FedPD-QSGD','FedPD-GD p=0.5','QFedPD-GD p=0.5','FedPD-GD p=0','QFedPD-GD p=0','FedPD-GD p=(T-r)/T','QFedPD-GD p=(T-r)/T','FedPD-VR','FedPD-QVR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
xlim([0 1000])







subplot(1,3,2)
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



subplot(1,3,3)

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


%----------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------

figure(2);
subplot(1,3,1)
semilogy(cumsum(itr2_q)-1, mean(log_QGD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
hold on



subplot(1,3,2)
semilogy(number_of_bits_toSend_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
hold on



subplot(1,3,3)
semilogy(Energy_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-','linewidth',linewidth,'color', 'k');
hold on


load non_iid_fixedQR.mat itr1 log_GD itr1_q log_QGD itr2 log_GD_1 itr2_q log_QGD_1 itr3 log_GD_2 itr3_q log_QGD_2 ...
    number_of_bits_toSend_nq1 number_of_bits_toSend_q1 number_of_bits_toSend_nq2 number_of_bits_toSend_q2 number_of_bits_toSend_nq3 ...
    number_of_bits_toSend_q3 Energy_nq1 Energy_q1 Energy_nq2 Energy_q2 Energy_nq3 Energy_q3 ...
    itrSGD log_SGD itrQSGD log_QSGD itrVR log_VR itrQVR log_QVR iter_g log_LGD log_LSGD log_Prox ...
    number_of_bits_toSend_nq0 number_of_bits_toSend_q0 number_of_bits_toSend_nq_vr number_of_bits_toSend_q_vr number_of_bits_toSend_nq

subplot(1,3,1)
semilogy(cumsum(itr2_q)-1, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
xl = xlabel({'Communication r';'(a)'},'FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('QFedPD-GD-DR p=0','QFedPD-GD-2bits p=0');%'EXTRA','xFilter','DSGD','GNSD' );
xlim([0 1000])


subplot(1,3,2)
semilogy(number_of_bits_toSend_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
xl = xlabel('number of transmitted bits','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
%le = legend('FedPD-SGD','FedPD-GD p=0.5','FedPD-GD p=0','FedPD-GD p=(T-r)/T','FedPD-VR','FedAvg-GD','FedAvg-SGD','FedProx-VR');%'EXTRA','xFilter','DSGD','GNSD' );
le = legend('QFedPD-GD-DR p=0','QFedPD-GD-2bits p=0');

subplot(1,3,3)
semilogy(Energy_q2, mean(log_QGD_1(:,2,:),3),'linestyle', '-.','linewidth',linewidth,'color', 'r');
xl = xlabel('Total energy bill','FontSize',fontsize,'interpreter','latex');
yl = ylabel('$\Vert\nabla f(x^r_0) \Vert^2$','FontSize',fontsize,'interpreter','latex');
le = legend('QFedPD-GD-DR p=0','QFedPD-GD-2bits p=0');


