function y =plot_graphics(dlambda,cas,VAR_HAM,COMPW,COMPTh,itermax,h)

%     figure(1)
%     title('Hamiltonian variation','fontsize',12)
%     hold on
%     xlabel('Time','fontsize',25)
%     ylabel('abs(H(W)-H(W0))','fontsize',12)
%     plot(h*(1:itermax),VAR_HAM,'.-');
%     set(gca, 'FontSize',12)
%     W12=zeros(itermax,1);
%     W23=zeros(itermax,1);
%     for i=1:1:itermax
%         W12(i)=COMP(i,3,2);
%         W23(i)=COMP(i,3,1);
%     end
%     figure(2)
%     title('Phase portrait','fontsize',12)
%     hold on
%     xlabel('W12','fontsize',12)
%     ylabel('W23','fontsize',12)
%     plot(W12,W23);
%     y = W12;
    figure(1)
    title('Eigenvalues of $\Theta$ variation','fontsize',12,'Interpreter','latex')
    hold on
    xlabel('t','fontsize',14,'Interpreter','latex')
   % ylabel('abs(lambda(W)-lambda(W0))','fontsize',12)
    p1=plot(h*(1:itermax),dlambda,'.-');
    saveas(p1,'spec_Kir_LSK','epsc');
    set(gca, 'FontSize',12)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    title('Cross-helicity variation','fontsize',12,'Interpreter','latex')
    hold on
    xlabel('t','fontsize',14,'Interpreter','latex')
   % ylabel('abs(lambda(W)-lambda(W0))','fontsize',12)
    p1=plot(h*(1:itermax),cas,'.-');
    saveas(p1,'cas_Kir_LSK','epsc');
    set(gca, 'FontSize',12)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure(3)
    title('Hamiltonian variation','fontsize',12,'Interpreter','latex')
    hold on
    xlabel('t','fontsize',14,'Interpreter','latex')
   % ylabel('abs(lambda(W)-lambda(W0))','fontsize',12)
    p1=plot(h*(1:itermax),VAR_HAM,'.-');
    saveas(p1,'ham_Kir_LSK','epsc');
    set(gca, 'FontSize',12)
    W12=zeros(itermax,1);
    W23=zeros(itermax,1);
    for i=1:1:itermax
        W12(i)=COMPW(i,1,3);
        W23(i)=COMPW(i,3,2);
    end
    
    
    figure(4)
    title('Phase portrait for $W$','fontsize',12,'Interpreter','latex')
    hold on
    xlabel('$W_{13}$','fontsize',12,'Interpreter','latex')
    ylabel('$W_{32}$','fontsize',12,'Interpreter','latex',"Rotation",0)
    p1=plot(W12,W23);
    saveas(p1,'phaseW_Kir_LSK','epsc');
    for i=1:1:itermax
        W12(i)=COMPTh(i,1,3);
        W23(i)=COMPTh(i,3,2);
    end
    
    
    figure(5)
    title('Phase portrait for $\Theta$','fontsize',12,'Interpreter','latex')
    hold on
    xlabel('$\Theta_{13}$','fontsize',12,'Interpreter','latex')
    ylabel('$\Theta_{32}$','fontsize',12,'Interpreter','latex',"Rotation",0)
    p1=plot(W12,W23);
    saveas(p1,'phaseTh_Kir_LSK','epsc');
    