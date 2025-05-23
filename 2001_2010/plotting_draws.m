load('param2001.mat');
subplot(2,5,1);plot(prop_educ(:,1),'LineWidth',2.5);title('prop category 1 (education)');
subplot(2,5,2);plot(prop_educ(:,2),'LineWidth',2.5);title('prop category 2 (education)');
subplot(2,5,3);plot(prop_educ(:,3),'LineWidth',2.5);title('prop category 3 (education)');
subplot(2,5,4);plot(prop_educ(:,4),'LineWidth',2.5);title('prop category 4 (education)');
subplot(2,5,5);plot(prop_educ(:,5),'LineWidth',2.5);title('prop category 5 (education)');
subplot(2,5,6);plot(prop_happiness(:,1),'LineWidth',2.5);title('prop category 1 (happiness)');
subplot(2,5,7);plot(prop_happiness(:,2),'LineWidth',2.5);title('prop category 2 (happiness)');
subplot(2,5,8);plot(prop_happiness(:,3),'LineWidth',2.5);title('prop category 3 (happiness)');
subplot(2,5,9);plot(prop_happiness(:,4),'LineWidth',2.5);title('prop category 4 (happiness)');
subplot(2,5,10);plot(prop_happiness(:,5),'LineWidth',2.5);title('prop category 5 (happiness)');

subplot(2,3,1);plot(theta_gauss(:,1),'LineWidth',2);title('r12');
subplot(2,3,2);plot(theta_gauss(:,2),'LineWidth',2);title('r13');
subplot(2,3,3);plot(theta_gauss(:,3),'LineWidth',2);title('r14');
subplot(2,3,4);plot(theta_gauss(:,4),'LineWidth',2);title('r23');
subplot(2,3,5);plot(theta_gauss(:,5),'LineWidth',2);title('r24');
subplot(2,3,6);plot(theta_gauss(:,6),'LineWidth',2);title('r34');

subplot(3,3,1);plot(w_gamma(:,1),'LineWidth',2);title('\xi_{1G}');
subplot(3,3,2);plot(w_gamma(:,2),'LineWidth',2);title('\xi_{2G}');
subplot(3,3,3);plot(w_gamma(:,3),'LineWidth',2);title('\xi_{3G}');
subplot(3,3,4);plot(m_gamma(:,1),'LineWidth',2);title('\mu_{1}');
subplot(3,3,5);plot(m_gamma(:,2),'LineWidth',2);title('\mu_{2}');
subplot(3,3,6);plot(m_gamma(:,3),'LineWidth',2);title('\mu_{3}');
subplot(3,3,7);plot(v_gamma(:,1),'LineWidth',2);title('v_{1}');
subplot(3,3,8);plot(v_gamma(:,2),'LineWidth',2);title('v_{2}');
subplot(3,3,9);plot(v_gamma(:,3),'LineWidth',2);title('v_{3}');

subplot(3,3,1);plot(w_beta(:,1),'LineWidth',2);title('\xi_{1B}');
subplot(3,3,2);plot(w_beta(:,2),'LineWidth',2);title('\xi_{2B}');
subplot(3,3,3);plot(w_beta(:,3),'LineWidth',2);title('\xi_{3B}');
subplot(3,3,4);plot(m_beta(:,1),'LineWidth',2);title('m_{1}');
subplot(3,3,5);plot(m_beta(:,2),'LineWidth',2);title('m_{2}');
subplot(3,3,6);plot(m_beta(:,3),'LineWidth',2);title('m_{3}');
subplot(3,3,7);plot(s_beta(:,1),'LineWidth',2);title('s_{1}');
subplot(3,3,8);plot(s_beta(:,2),'LineWidth',2);title('s_{2}');
subplot(3,3,9);plot(s_beta(:,3),'LineWidth',2);title('s_{3}');

load('param2010.mat');
subplot(2,5,1);plot(prop_educ(:,1),'LineWidth',2.5);title('prop category 1 (education)');
subplot(2,5,2);plot(prop_educ(:,2),'LineWidth',2.5);title('prop category 2 (education)');
subplot(2,5,3);plot(prop_educ(:,3),'LineWidth',2.5);title('prop category 3 (education)');
subplot(2,5,4);plot(prop_educ(:,4),'LineWidth',2.5);title('prop category 4 (education)');
subplot(2,5,5);plot(prop_educ(:,5),'LineWidth',2.5);title('prop category 5 (education)');
subplot(2,5,6);plot(prop_happiness(:,1),'LineWidth',2.5);title('prop category 1 (happiness)');
subplot(2,5,7);plot(prop_happiness(:,2),'LineWidth',2.5);title('prop category 2 (happiness)');
subplot(2,5,8);plot(prop_happiness(:,3),'LineWidth',2.5);title('prop category 3 (happiness)');
subplot(2,5,9);plot(prop_happiness(:,4),'LineWidth',2.5);title('prop category 4 (happiness)');
subplot(2,5,10);plot(prop_happiness(:,5),'LineWidth',2.5);title('prop category 5 (happiness)');

subplot(2,3,1);plot(theta_gauss(:,1),'LineWidth',2);title('r12');
subplot(2,3,2);plot(theta_gauss(:,2),'LineWidth',2);title('r13');
subplot(2,3,3);plot(theta_gauss(:,3),'LineWidth',2);title('r14');
subplot(2,3,4);plot(theta_gauss(:,4),'LineWidth',2);title('r23');
subplot(2,3,5);plot(theta_gauss(:,5),'LineWidth',2);title('r24');
subplot(2,3,6);plot(theta_gauss(:,6),'LineWidth',2);title('r34');

subplot(2,3,1);plot(w_gamma(:,1),'LineWidth',2);title('\xi_{1G}');
subplot(2,3,2);plot(w_gamma(:,2),'LineWidth',2);title('\xi_{2G}');
subplot(2,3,3);plot(m_gamma(:,1),'LineWidth',2);title('\mu_{1}');
subplot(2,3,4);plot(m_gamma(:,2),'LineWidth',2);title('\mu_{2}');
subplot(2,3,5);plot(v_gamma(:,1),'LineWidth',2);title('v_{1}');
subplot(2,3,6);plot(v_gamma(:,2),'LineWidth',2);title('v_{2}');

subplot(3,3,1);plot(w_beta(:,1),'LineWidth',2);title('\xi_{1B}');
subplot(3,3,2);plot(w_beta(:,2),'LineWidth',2);title('\xi_{2B}');
subplot(3,3,3);plot(w_beta(:,3),'LineWidth',2);title('\xi_{3B}');
subplot(3,3,4);plot(m_beta(:,1),'LineWidth',2);title('m_{1}');
subplot(3,3,5);plot(m_beta(:,2),'LineWidth',2);title('m_{2}');
subplot(3,3,6);plot(m_beta(:,3),'LineWidth',2);title('m_{3}');
subplot(3,3,7);plot(s_beta(:,1),'LineWidth',2);title('s_{1}');
subplot(3,3,8);plot(s_beta(:,2),'LineWidth',2);title('s_{2}');
subplot(3,3,9);plot(s_beta(:,3),'LineWidth',2);title('s_{3}');