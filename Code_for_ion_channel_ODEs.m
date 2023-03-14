%% This code complements the article "Modelling the flow through ion 
% channels at the cell membrane" by Pamela Burrage and Adrianne Jenner. 

%% In this section, we plot the steady state values Cstar and Pstar for 
% different k1, k2 values

%setting up the k1 and k2 values to calculate Cstar and Pstar
k1 = linspace(1e-6,100,10);
k2 = k1;
[k1m,k2m] = meshgrid(k1,k2);

N = 50; %fixing N = 50

%calculating Pstar and Cstar
Pstar = k1m./(k1m+k2m).*N;
Cstar = k2m./(k1m+k2m).*N;

% setting the colouring
colmap = cbrewer2('YlGnBu');

% plotting the value of Pstar and Cstar for the different k1 and k2 values
surf(k1m,k2m,Pstar,'EdgeColor','none')
colormap(colmap)
colorbar
title('P^*')
set(gca,'FontSize',18)
xlabel('k_1')
ylabel('k_2')

figure
surf(k1m,k2m,Cstar,'EdgeColor','none')
colormap(colmap)
colorbar
title('C^*')
set(gca,'FontSize',18)
xlabel('k_1')
ylabel('k_2')

%% In this section, we solve the ODE model and plot the Cstar and Pstar

% setting the model parameters
k1 = 1;
k2 = 2;
N = 50;
C0 = 5;

t0 = 0;
T = 10;

P0 = N - C0;

%calculating the steady state values
Pstar = k1./(k1+k2).*N;
Cstar = k2./(k1+k2).*N;

%solving the ODE
[t,y] = ode45(@(t,y)ionChannelFn(t,y,k1,k2),[t0 T],[C0 P0]);

% plotting the ODE
figure
hold on
plot([0 T],[Cstar, Cstar],'Color',[242,184,0]/255,'LineWidth',2)
plot([0 T],[Pstar, Pstar],'Color',[100 36 228]/255,'LineWidth',2)
plot(t,y(:,1),'*--','Color',[242,184,0]/255) % y1, y2 vs time
plot(t,y(:,2),'*--','Color',[100 36 228]/255) % y1, y2 vs time
legend('Closed C(t)','Open P(t)')

xlabel('Time (t)')
ylabel('Channels')
set(gca,'FontSize',22)
ylim([0 50])

%% Simulating the SSA

% initialising the index
i = 1;

%initialising the vectors for storing the value of C and P
C_vec(i) = C0;
P_vec(i) = P0;

t_vec(i) = 0;

while max(t_vec)<T
    i = i+1;

    q1 = k1*C_vec(i-1);
    q2 = k2*P_vec(i-1);

    r = q1+q2;

    u = rand;

    delta_t = -1/r*log(u);

    t_vec(i) = t_vec(i-1)+delta_t;

u = rand;
if u<q1/r
C_vec(i) = C_vec(i-1)-1;
P_vec(i) = P_vec(i-1)+1;
else
C_vec(i) = C_vec(i-1)+1;
P_vec(i) = P_vec(i-1)-1;

end

end

figure
hold on 
plot(t_vec,C_vec,'Color',[242,184,0]/255,'LineWidth',0.5)
plot(t_vec,P_vec,'Color',[100 36 228]/255,'LineWidth',0.5)
plot(t,y(:,1),':','Color',[242,184,0]/255,'LineWidth',3) % y1, y2 vs time
plot(t,y(:,2),':','Color',[100 36 228]/255,'LineWidth',3) % y1, y2 vs time
legend('C(t) - SSA','Open P(t)- SSA','Closed C(t) - ODE','Open P(t) - ODE')

xlabel('Time (t)')
ylabel('Channels')
set(gca,'FontSize',18)
ylim([0 50])

%%

function dydt = ionChannelFn(t,y,k1,k2)

    C = y(1);
    P = y(2);

    dCdt = -k1*C+k2*P;
    dPdt = k1*C-k2*P;

    dydt = [dCdt;dPdt];

end
