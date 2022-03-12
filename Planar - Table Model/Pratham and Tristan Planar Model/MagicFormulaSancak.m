%% Magic Formula Tire Model by Faruk Sancak
%Slip angle takes values between -20 and 20 degrees.
%Slip ratio takes values between -100% and 100%.
ax1=-22.3; ax2=1144; ax3=49.6; ax4=226; ax5=0.069; ax6=-0.006; ax7=0.056;
ax8=0.486;
ay1=-22.1; ay2=1011; ay3=1078; ay4=1.82; ay5=0.208; ay6=0; ay7=-0.354;
ay8=0.707;
Cx=1.65;
Cy=1.3;
%%
j=0;
for Fz=2:2:8
    i=0;
    j=j+1;
    for alpha=-20:0.01:20
        i=i+1;
        Dx=(ax1*(Fz^2))+(ax2*Fz);
        BCDx=((ax3*(Fz^2))+(ax4*(Fz^2)))/(exp(ax5*Fz));
        Bx=BCDx/(Cx*Dx);
        Ex=(ax6*(Fz^2))+ax7*Fz+ax8;
        %%
        Dy=(ay1*(Fz^2))+(ay2*Fz);
        BCDy=ay3*sind(ay4*atand(ay5*Fz));
        By=BCDy/(Cy*Dy);
        Ey=(ay6*(Fz^2))+ay7*Fz+ay8;
        %%
        Fx(i,j)=Dx*sind(Cx*atand(Bx*(alpha)));
        Fy(i,j)=Dy*sind(Cy*atand(By*(alpha)));
    end
end
%%
ratio=-100:0.05:100;
figure
plot(ratio,Fx(:,1),ratio,Fx(:,2),ratio,Fx(:,3),ratio,Fx(:,4),'linewidth',2)
grid on
title('Magic Formula Tire Model','fontsize',12,'fontweight','b')
legend('F_z = 2 kN','F_z = 4 kN','F_z = 6 kN','F_z = 8kN','Location','NorthWest')
xlabel('Slip %','fontsize',12,'fontweight','b')
ylabel('F_x [N]','fontsize',12,'fontweight','b')

alpha=-20:0.01:20;
figure
plot(alpha,Fy(:,1),alpha,Fy(:,2),alpha,Fy(:,3),alpha,Fy(:,4),'linewidth',2)
grid on
title('Magic Formula Tire Model','fontsize',12,'fontweight','b')
legend('F_z = 2 kN','F_z = 4 kN','F_z = 6 kN','F_z = 8kN','Location','NorthWest')
xlabel('Slip angle [deg]','fontsize',12,'fontweight','b')
ylabel('F_y [N]','fontsize',12,'fontweight','b')